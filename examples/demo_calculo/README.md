[**English**](README.md) · [Español](README.es.md)

# Silicon: convergence, equation of state and elastic constants

The three calculation modules (`converge`, `eos`, `elastic`) run on silicon
with LDA and a norm-conserving pseudopotential.

    olla-dft converge Si.cif --kind ecutwfc --values 20,30,40,50,60,70 --run
    olla-dft eos Si.cif --ecutwfc 60 --run
    olla-dft elastic Si_eq.cif --ecutwfc 60 --run

Results: the cutoff converges at 50 Ry (1 meV/atom); a0 = 5.402 Å and
B0 = 94.2 GPa; C11/C12/C44 = 159.9/61.7/76.6 GPa, with B = 94.45 GPa
computed from the Cij (0.25 % off the EOS value). Experiment:
165.8/63.9/79.6 GPa and B ≈ 98 GPa.

### Files

| File | What it is |
|---|---|
| `CONVERGENCIA.dat`, `CONVERGENCIA.txt` | E(ecutwfc) table and the convergence report with the recommended cutoff |
| `convergencia.png` | convergence figure |
| `EOS.dat`, `EOS.txt` | E(V) points and the Birch–Murnaghan fit (V0, B0, B0') |
| `eos.png` | equation-of-state figure |
| `ELASTIC_C.dat`, `ELASTIC.txt` | elastic matrix Cij and the report (Voigt–Reuss–Hill moduli, stability) |
| `elastic.png` | stress–strain figure of the applied deformations |
