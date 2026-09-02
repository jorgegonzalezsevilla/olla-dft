# bcc iron: spin-polarised DOS

Ferromagnetic bcc iron: input generation with an initial magnetic moment
and a spin-resolved density of states.

    olla-dft gen Fe.cif -p dos -o . --mag Fe=0.7
    olla-dft dos . -o . --journal aps --emin -10 --emax 6

pw.x converges to 2.28 μB/cell; integrating the spin-resolved DOS that
Olla-DFT exports recovers 2.27 μB independently.

### Files

| File | What it is |
|---|---|
| `scf.in`, `nscf.in` | pw.x inputs written by `olla-dft gen` (scf with `starting_magnetization` and a dense nscf for the DOS) |
| `Fe_dos.pdf`, `Fe_dos.png` | spin-resolved DOS (up/down) in APS style |
