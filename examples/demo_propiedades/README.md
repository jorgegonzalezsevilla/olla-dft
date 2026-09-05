[**English**](README.md) · [Español](README.es.md)

# Silicon: optical functions, phonons and effective mass

The `optics`, `phonons` and `effmass` modules run on relaxed silicon (LDA,
norm-conserving, QE 6.6). The structure is `Si_relajado.cif` (a0 = 5.402 Å,
output of the EOS in `demo_calculo`).

    olla-dft optics  Si_relajado.cif --ecutwfc 60 --run
    olla-dft optics  Si_relajado.cif --collect --scissor 0.65
    olla-dft phonons Si_relajado.cif --qgrid 2x2x2 --run
    olla-dft phonons Si_relajado.cif --collect

### Optics (`opticas_Si.png`, `OPTICS.dat`)

14×14×14 nscf without symmetry, 27 bands, 0.1 eV broadening.

Without a scissor the calculation gives ε1(0) = 16.13, well above the
experimental 11.7: this is the direct consequence of the LDA gap being
0.52 eV instead of 1.17 eV. The spectrum itself is correct and was checked
twice:

- f-sum rule: ∫ E·ε2(E) dE = 451.5 eV² against (π/2)(ħωp)² = 451.5 eV² with
  the plasmon frequency reported by epsilon.x (16.95 eV; experimental
  16.7). Factor 1.000.
- Kramers–Kronig of ε2 reproduces the ε1 of epsilon.x: 16.09 vs 16.13
  (0.3 %).

With `--scissor 0.65` (1.17 − 0.52) the ε2 peak lands at 4.30 eV, exactly
the E2 critical point of silicon, and ε1(0) drops to 10.44.

### Phonons (`fonones_Si.png`, `FONONES_DOS.dat`, `FONONES_TERMO.dat`)

2×2×2 q mesh (8 points), DOS interpolated on 12×12×12.

Frequencies (cm⁻¹) against inelastic neutron scattering:

    point        Olla-DFT   experimental
    Γ     TO/LO  508.9      517
    X     TA     140.8      150
    X     LA/LO  406.5      410
    X     TO     455.9      463
    L     TA     107.7      114
    L     LA     372.3      378
    L     LO     408.5      417
    L     TO     484.8      490

All 1 to 6 % low, which is what LDA with such a small q mesh is expected to
give, and with the correct degeneracies. No imaginary frequencies: the
structure is well relaxed.

Harmonic thermodynamics per cell (2 atoms):
- zero-point energy 122.25 meV
- at 300 K: C_v = 0.411 meV/K, S = 0.422 meV/K, F = 65.0 meV

The experimental C_v of Si at 300 K is 20 J/(mol·K), which per two-atom
cell is 0.415 meV/K: a 1 % difference.

### Effective mass (`MASA_EFECTIVA.dat`)

    olla-dft gen Si_relajado.cif --preset bands -o bandas    # then run it
    olla-dft effmass Si_relajado.cif --bands-dir bandas -o masa --run

Fine path of 6 lines, ±0.06 Å⁻¹, 21 points each.

    mass                 Olla-DFT   reference
    electron long.       0.949      0.916  (experimental)
    electron transv.     0.193      0.190  (experimental)
    heavy hole [100]     0.269      0.277  (Luttinger)
    heavy hole [111]     0.670      0.718  (Luttinger)

The two transverse masses come out identical, as valley symmetry demands.
The quick fit on the ordinary band path is NOT publishable: its window is
0.35 Å⁻¹, outside the parabolic regime, and the report flags it as
unreliable.

### Files

| File | What it is |
|---|---|
| `Si_relajado.cif` | starting structure (relaxed silicon, a0 = 5.402 Å) |
| `OPTICS.dat` | ε1, ε2, n, k, α and R versus energy (isotropic average) |
| `opticas_Si.png` | optical-functions figure |
| `FONONES_DOS.dat` | phonon density of states |
| `FONONES_TERMO.dat` | F, U, S and Cv versus T in the harmonic approximation |
| `fonones_Si.png` | phonon dispersion and DOS |
| `MASA_EFECTIVA.dat` | fitted effective masses with their window and fit quality |
