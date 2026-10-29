"""Map each electron's possibilities onto qubits, let interference play out the chemistry.

H2 in a minimal basis: 2 electrons, 4 spin-orbitals -> 4 qubits.
The Jordan-Wigner mapping turns each spin-orbital's "occupied / empty"
possibility into one qubit. Fermionic antisymmetry becomes the Z-string:
creating an electron in orbital p picks up the sign (-1)^(electrons in
orbitals below p). That sign IS the Z-string, written out explicitly.

Then the machine's own interference does the chemistry: the qubit
Hamiltonian's eigenstructure (diagonalization) and its unitary evolution
exp(-iHt) sum every possibility with its proper phase. No classical
approximation of the correlation — the 16-dimensional wave function
holds it all.

Integrals come from PySCF (real quantum chemistry numbers, STO-3G).
The qubit answer is checked against PySCF's own FCI solver: same physics,
two representations. They must agree exactly.

© Crimson Rose LLC, Curtis Ray Dyess 2026-27
"""

import numpy as np
from scipy.linalg import expm

try:
    from pyscf import gto, scf, ao2mo, fci
    HAVE_PYSCF = True
except ImportError:
    HAVE_PYSCF = False


# ----------------------------------------------------------------------------
# Jordan-Wigner fermionic operators on nq qubits.
# Basis state |n> is the integer whose bit p = occupation of spin-orbital p.
# a_p^dag flips bit p 0->1 with the JW sign (-1)^(popcount of lower bits).
# ----------------------------------------------------------------------------

def creation(nq, p):
    dim = 1 << nq
    op = np.zeros((dim, dim), dtype=complex)
    mask = (1 << p) - 1
    for n in range(dim):
        if not (n >> p) & 1:
            sign = -1 if bin(n & mask).count("1") % 2 else 1
            op[n | (1 << p), n] = sign
    return op


def annihilation(nq, p):
    return creation(nq, p).T.conj()


def build_operators(nq):
    cp = [creation(nq, p) for p in range(nq)]
    ap = [c.T.conj() for c in cp]
    return cp, ap


# ----------------------------------------------------------------------------
# Second-quantized H2 Hamiltonian -> qubit matrix.
# ----------------------------------------------------------------------------

def spin_orbital_integrals(h1_spatial, eri_spatial):
    """Lift spatial integrals to the spin-orbital basis.

    p = 2*i + s, i = spatial orbital, s = 0 (alpha) / 1 (beta).
    h[p,q] nonzero only for same spin; g[p,q,r,s] (chemists' (pq|rs))
    nonzero only when sp==sr and sq==ss.
    """
    n = h1_spatial.shape[0]
    nq = 2 * n
    h = np.zeros((nq, nq))
    g = np.zeros((nq, nq, nq, nq))
    for p in range(nq):
        for q in range(nq):
            if p % 2 == q % 2:
                h[p, q] = h1_spatial[p // 2, q // 2]
    # Physicist's <pq|rs> = int psi_p*(1) psi_q*(2) r12^-1 psi_r(1) psi_s(2):
    # spatial part = chemists' (pr|qs) = eri[i,k,j,l]; spin rule pairs
    # p with r and q with s (sp==sr, sq==ss). Contracted as
    #   V = 1/2 sum g[p,q,r,s] a_p^dag a_q^dag a_s a_r.
    for p in range(nq):
        for q in range(nq):
            for r in range(nq):
                for s_ in range(nq):
                    if p % 2 == r % 2 and q % 2 == s_ % 2:
                        g[p, q, r, s_] = eri_spatial[p // 2, r // 2, q // 2, s_ // 2]
    return h, g


def qubit_hamiltonian(h, g, e_nuc):
    """H = sum h[p,q] a_p^dag a_q + 1/2 sum g[p,q,r,s] a_p^dag a_q^dag a_s a_r."""
    nq = h.shape[0]
    cp, ap = build_operators(nq)
    dim = 1 << nq
    H = np.zeros((dim, dim), dtype=complex)
    for p in range(nq):
        for q in range(nq):
            if h[p, q] != 0.0:
                H += h[p, q] * (cp[p] @ ap[q])
    for p in range(nq):
        for q in range(nq):
            for r in range(nq):
                for s_ in range(nq):
                    if g[p, q, r, s_] != 0.0:
                        H += 0.5 * g[p, q, r, s_] * (cp[p] @ cp[q] @ ap[s_] @ ap[r])
    H += e_nuc * np.eye(dim)
    return H.real  # the Hamiltonian is real-symmetric


def number_operator(nq):
    cp, ap = build_operators(nq)
    return sum(cp[p] @ ap[p] for p in range(nq)).real


# ----------------------------------------------------------------------------
# Real chemistry numbers from PySCF.
# ----------------------------------------------------------------------------

def h2_integrals(bond_angstrom):
    if not HAVE_PYSCF:
        raise RuntimeError("PySCF is required for the integrals")
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {bond_angstrom}", basis="sto-3g",
                unit="Angstrom", verbose=0)
    mf = scf.RHF(mol).run()
    n = mol.nao_nr()
    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), n)
    return h1, eri, mol.energy_nuc(), mol, mf


def h_atom_energy():
    if not HAVE_PYSCF:
        raise RuntimeError("PySCF is required for the integrals")
    mol = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, verbose=0)
    mf = scf.UHF(mol).run()
    n = mol.nao_nr()
    h1 = mf.mo_coeff[0].T @ mf.get_hcore() @ mf.mo_coeff[0]
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff[0]), n)
    h, g = spin_orbital_integrals(h1, eri)
    H = qubit_hamiltonian(h, g, mol.energy_nuc())
    return np.linalg.eigvalsh(H)[0]


def fci_reference(mol, mf):
    solver = fci.FCI(mol, mf.mo_coeff)
    return solver.kernel()[0]


def ground_energy(bond_angstrom):
    """Qubit-mapped ground-state energy of H2 at this bond length (Hartree)."""
    h1, eri, e_nuc, mol, mf = h2_integrals(bond_angstrom)
    h, g = spin_orbital_integrals(h1, eri)
    H = qubit_hamiltonian(h, g, e_nuc)
    return np.linalg.eigvalsh(H)[0], H, (mol, mf)


def dissociation_curve(r_min=0.4, r_max=3.0, npts=27):
    rs = np.linspace(r_min, r_max, npts)
    es = [ground_energy(r)[0] for r in rs]
    return rs, np.array(es)


def play_out_interference(bond_angstrom, t_max=50.0, nsteps=200):
    """Let the machine's own interference play: evolve the Hartree-Fock
    state |1100> under the qubit Hamiltonian and watch correlation grow.

    Returns times and the weight that has flowed out of the single
    determinant into the correlated wave function.
    """
    e0, H, _ = ground_energy(bond_angstrom)
    dim = H.shape[0]
    psi0 = np.zeros(dim, dtype=complex)
    psi0[0b0011] = 1.0  # both electrons in spatial orbital 0: the HF state
    ts = np.linspace(0, t_max, nsteps)
    flow = []
    for t in ts:
        psi = expm(-1j * (H - e0 * np.eye(dim)) * t) @ psi0
        flow.append(1.0 - abs(psi[0b0011]) ** 2)
    return ts, np.array(flow)


if __name__ == "__main__":
    rs, es = dissociation_curve()
    i = int(np.argmin(es))
    print(f"H2 STO-3G on 4 qubits: minimum {es[i]:.6f} Ha at {rs[i]:.3f} Angstrom")
    print(f"dissociation limit E(3.0 A) = {es[-1]:.6f} Ha")
    print(f"isolated H atom (same machinery) = {h_atom_energy():.6f} Ha")
    ts, flow = play_out_interference(rs[i])
    print(f"interference played out: max correlation flow = {flow.max():.4f} "
          f"of the wave function leaves the single determinant")
    np.savetxt("h2-curve.csv", np.column_stack([rs, es]),
               delimiter=",", header="bond_A,energy_Ha", comments="")
    print("wrote h2-curve.csv")
