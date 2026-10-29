"""Tests for the qubit chemistry mapping.

The mapping must obey fermionic law on the qubits, conserve particles,
and reproduce the chemistry exactly (FCI agreement). The physics checks
are self-consistent: the code's own H atom sets the dissociation limit.
© Crimson Rose LLC, Curtis Ray Dyess 2026-27
"""

import numpy as np
import qubit_chemistry as qc


def check(name, cond):
    assert cond, f"FAILED: {name}"
    print(f"pass: {name}")


nq = 4
cp, ap = qc.build_operators(nq)
I = np.eye(1 << nq)

# --- Fermionic law survives the mapping ---------------------------------------
for p in range(nq):
    for q in range(nq):
        anti = ap[p] @ cp[q] + cp[q] @ ap[p]
        want = I if p == q else np.zeros_like(I)
        check(f"anticommutation {{a_{p}, a^dag_{q}}} = delta",
              np.allclose(anti, want, atol=1e-12))
        check(f"Pauli exclusion a^dag_{p} a^dag_{p} = 0",
              np.allclose(cp[p] @ cp[p], 0, atol=1e-12) if p == q else True)

# --- The Hamiltonian is honest --------------------------------------------------
h1, eri, e_nuc, mol, mf = qc.h2_integrals(0.74)
h, g = qc.spin_orbital_integrals(h1, eri)
H = qc.qubit_hamiltonian(h, g, e_nuc)
check("4 spin-orbitals -> 4 qubits -> 16x16", H.shape == (16, 16))
check("Hamiltonian is Hermitian", np.allclose(H, H.T, atol=1e-12))
N = qc.number_operator(nq)
check("particle number conserved [H, N] = 0",
      np.allclose(H @ N - N @ H, 0, atol=1e-10))

# --- Qubit answer == FCI answer (same physics, two representations) --------------
for R in (0.5, 0.74, 1.0, 2.0):
    e_qubit, _, (m, mf_) = qc.ground_energy(R)
    e_fci = qc.fci_reference(m, mf_)
    check(f"qubit == FCI at R={R} A",
          abs(e_qubit - e_fci) < 1e-8)
    print(f"      R={R}: qubit {e_qubit:.8f}  FCI {e_fci:.8f}")

# --- The chemistry plays out correctly -------------------------------------------
e_mid, _, _ = qc.ground_energy(0.74)
e_short, _, _ = qc.ground_energy(0.4)
e_long, _, _ = qc.ground_energy(1.5)
check("bond minimum near 0.74 A", e_mid < e_short and e_mid < e_long)

e_diss, _, _ = qc.ground_energy(4.0)
e_h = qc.h_atom_energy()
check("dissociation -> two isolated H atoms",
      abs(e_diss - 2 * e_h) < 0.02)
print(f"      E(4.0 A) = {e_diss:.6f}, 2*E(H) = {2*e_h:.6f}")

# --- Interference really flows ------------------------------------------------------
ts, flow = qc.play_out_interference(0.74, t_max=30.0, nsteps=60)
check("correlation flows out of the single determinant", flow.max() > 0.01)
check("flow stays a probability", np.all((flow >= 0) & (flow <= 1)))

print("\nAll mapping tests pass: fermions obey the law on qubits,")
print("the qubit chemistry matches FCI exactly, and the bond behaves like H2.")
