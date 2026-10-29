# Quantum Chemistry on Qubits

H₂ mapped onto 4 qubits via the Jordan–Wigner transform — Curtis's order 2026-10-29: *"Map each electrons possibities onto qubits and let the machines own interference play out the chemistry."*

**Status: DRAFT / GATHERED.** The code is written and the chemistry is real; nothing here stands until Curtis says it stands.

## The idea, plainly

A hydrogen molecule has 2 electrons and, in the minimal picture, 2 spatial orbitals. Each electron can sit in either orbital with either spin: 4 spin-orbitals, 4 possibilities per electron.

**The map:** one qubit per spin-orbital. |1100⟩ means both electrons piled in orbital 0 (the simple Hartree–Fock picture). |0011⟩ means both jumped to orbital 1. The real molecule lives in ALL of these at once, weighted — and the machine's own interference (diagonalizing the 16×16 qubit Hamiltonian) plays out the chemistry.

**The fermionic law survives the trip.** Pauli exclusion and exchange antisymmetry become the Jordan–Wigner Z-string, written out explicitly in the code. Tested: {a_p, a_q†} = δ_pq exactly on the qubits, and a†a† = 0 (no two electrons in one state, ever).

## What the build proved

- 4 spin-orbitals → 4 qubits → 16-dimensional wave function. Exact.
- Qubit ground energy matches PySCF's FCI solver to 8 decimals at R = 0.5, 0.74, 1.0, 2.0 Å (agreement < 1e-8).
- Bond minimum near 0.74 Å; dissociation → two isolated H atoms, exact to the model's precision.
- Particle number conserved: [H, N] = 0.
- Correlation is real: up to 4.3% of the wave function leaves the single determinant during time evolution — the part Hartree–Fock can never see.
- Recorded curve: `h2-curve.csv` (energies), `h2-curve.html` (plot), `h2-qubit-schematic.html` (interactive schematic).

## Honest limits

- STO-3G is the smallest basis in the business: qualitatively right, not chemically accurate (real H₂ binds at −1.174 Ha, 0.741 Å; we get −1.137 at ~0.70–0.74 Å).
- 4 qubits fit on a laptop. A useful molecule (say, caffeine) needs hundreds of qubits and error correction.
- This ran on a classical machine simulating the qubits exactly. A real quantum computer would prepare the state and let ITS interference do it. The map is identical; only the stage changes.

## Run it

```bash
pip install pyscf numpy scipy   # PySCF supplies the one- and two-electron integrals
python3 test_qubit_chemistry.py # operator identities run without PySCF; energy tests need it
```

`qubit_chemistry.py` raises a clear error if PySCF is missing — the operator-level tests (anticommutation, Pauli exclusion) run anywhere.

## Files

- `qubit_chemistry.py` — the map: integrals, Jordan–Wigner encoding, 16×16 Hamiltonian, time evolution
- `test_qubit_chemistry.py` — tests: fermionic identities, Pauli exclusion, energy vs FCI
- `MAP.md` — the full map document, plain language
- `h2-curve.csv` / `h2-curve.html` — dissociation curve data and plot
- `h2-qubit-schematic.html` — interactive schematic of the H₂ qubit map

---

© Crimson Rose LLC, Curtis Ray Dyess 2026-27
