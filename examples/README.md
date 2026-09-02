# Olla-DFT examples

Each subfolder carries REAL Quantum ESPRESSO output, not mock-ups: the
data, the figures and a `README.md` with the exact `olla-dft` commands
that produced them and the comparison with experiment. The commands in
every README are validated against the CLI in `tests/test_examples.py`,
so they cannot go stale without a test failing.

| Folder | What it demonstrates |
|---|---|
| [`demo_Si/`](demo_Si/) | silicon bands, DOS and PDOS: `gen`, `run.sh`, `plot` |
| [`demo_Fe/`](demo_Fe/) | spin-polarized bcc iron: `gen --mag`, `dos` |
| [`demo_calculo/`](demo_calculo/) | convergence, equation of state and elastic constants: `converge`, `eos`, `elastic` |
| [`demo_propiedades/`](demo_propiedades/) | optical functions, phonons and effective mass of silicon: `optics`, `phonons`, `effmass` |
| [`demo_derivadas/`](demo_derivadas/) | Debye/Slack from the Cij, quasi-harmonic expansion and material datasheet: `derived`, `qha`, `datasheet` |
| [`demo_laminar/`](demo_laminar/) | layers, diffractogram and exfoliation of graphite: `layers`, `xrd`, `exfoliate` |
| [`demo_espectros_avanzados/`](demo_espectros_avanzados/) | XANES, Hubbard U, electron-phonon, band unfolding and VDOS: `xanes`, `hubbard`, `elph`, `unfold`, `md` |
| [`demo_tddft_balistico/`](demo_tddft_balistico/) | TDDFPT of ethylene and ballistic conductance of an Al wire: `tddft`, `ballistic` |
| [`plantillas/`](plantillas/) | the same figure in every visual template: `templates`, `plot -t` |

### Loose structure files

Structure files to try any command without hunting for one of your own:

| File | What it is |
|---|---|
| `grafito.cif` | hexagonal graphite (4 atoms, 2 layers); used by `demo_laminar/` |
| `hbn.cif` | hexagonal boron nitride, another layered material |
| `ZnO.cif` | wurtzite zinc oxide (4 atoms), the example used in the main README |
| `POSCAR_NaCl` | sodium chloride in POSCAR format (8-atom conventional cell): Olla-DFT reads VASP as well as CIF |

Some demo commands name the starting structure they were run with
(`Si.cif`, `NiO.cif`, `Al.cif`, `c2h4.cif`...), which is not included: any
equivalent CIF works (Materials Project, COD...) or one of the files in
this folder.

For instance:

    olla-dft info POSCAR_NaCl
    olla-dft prim POSCAR_NaCl -o NaCl_prim.cif
    olla-dft layers hbn.cif
    olla-dft kpath ZnO.cif
    olla-dft gen ZnO.cif -p all -o ZnO_run --insulator

For guided start-to-finish sessions, see `olla-dft recetas`.
