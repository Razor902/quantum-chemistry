# The Map: Each Electron's Possibilities Onto Qubits

Curtis's order 2026-10-29: "Map each electrons possibities onto qubits and let the machines own interference play out the chemistry." Built the same night.

STATUS: DRAFT / GATHERED. The code is tested; the chemistry is real (PySCF integrals, FCI agreement). Nothing here stands until Curtis says it stands.
Watermark: © Crimson Rose LLC, Curtis Ray Dyess 2026-27

## The map, in plain language

A hydrogen molecule has 2 electrons and, in the minimal picture, 2 spatial orbitals (one on each atom). Each electron can sit in either orbital with either spin: 4 spin-orbitals, 4 possibilities per electron.

**The map:** one qubit per spin-orbital. Qubit 0 = "is there an electron in orbital 0, spin up?" Qubit 1 = orbital 0 spin down. Qubits 2, 3 = orbital 1 up/down. |1100⟩ means both electrons piled in orbital 0 — the simple Hartree-Fock picture. |0011⟩ means both jumped to orbital 1. The real molecule lives in ALL of these at once, weighted.

**The fermionic law survives the trip.** Electrons refuse to share a state (Pauli) and pick up a minus sign when swapped (antisymmetry). On qubits this becomes the Jordan-Wigner Z-string: creating an electron in orbital p flips the sign for every electron sitting in orbitals below p. The code writes this sign out explicitly — tested: {a_p, a_q†} = δ_pq holds exactly on the qubits, and a†a† = 0 (no two electrons in one state, ever).

**The machine's own interference plays out the chemistry.** The qubit Hamiltonian is a 16×16 matrix holding every possibility and every phase. Diagonalizing it sums all paths at once — the ground state comes out -1.13728383 Ha at 0.74 Å, matching PySCF's FCI solver to 8 decimals. Same physics, two representations. And evolving |1100⟩ under exp(-iHt) shows the wave function flowing out of the single determinant into correlated configurations — interference doing chemistry in real time.

## What the build proved

- 4 spin-orbitals → 4 qubits → 16-dimensional wave function. Exact.
- Qubit energy == FCI energy at R = 0.5, 0.74, 1.0, 2.0 Å (agreement < 1e-8).
- Bond minimum near 0.74 Å; dissociation → two isolated H atoms (-0.933171 vs 2×(-0.466582), exact to the model's precision).
- Particle number conserved: [H, N] = 0.
- Correlation is real: up to 4.3% of the wave function leaves the single determinant during evolution — the part Hartree-Fock can never see.

## Honest limits

- STO-3G is the smallest basis in the business: the numbers are qualitatively right, not chemically accurate (real H2 binds at -1.174 Ha, 0.741 Å; we get -1.137 at ~0.70-0.74).
- 4 qubits fit on a laptop. A useful molecule (say, caffeine) needs hundreds of qubits and error correction — Phase 4/5 territory.
- This ran on a classical machine simulating the qubits exactly. A real quantum computer would prepare the state and let ITS interference do it. The map is identical; only the stage changes.

## Files

- `qubit_chemistry.py` — the map, the Hamiltonian, the evolution, the curve.
- `test_qubit_chemistry.py` — fermionic law, Hermiticity, number conservation, FCI agreement, bond physics, interference flow. All pass.
- `h2-curve.csv` — the dissociation curve, 0.4 to 3.0 Å.
- `h2-curve.html` — the curve plotted, opens on his phone.

## Next rungs (his call)

1. Bigger basis (cc-pVDZ): 20 spin-orbitals, 20 qubits — still exactly simulable, closer to the real bond.
2. A real VQE loop: parameterized circuit on the 4 qubits, classical optimizer, watch it find -1.137.
3. LiH (4 electrons, 12 spin-orbitals): the next molecule up the staircase.
