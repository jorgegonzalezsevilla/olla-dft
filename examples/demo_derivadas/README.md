[**English**](README.md) · [Español](README.es.md)

# Silicon: Debye and Slack from the Cij, quasi-harmonic and material datasheet

Derived modules (`derived`, `qha`, `datasheet`): properties obtained from
results you already have, with no new Quantum ESPRESSO calculation.

Everything in this folder is POST-PROCESSING: it costs no new Quantum
ESPRESSO run, it comes out of results you already had.

### What to run

Debye temperature, sound velocities and Slack thermal conductivity from the
elastic constants:

    olla-dft derived ../demo_propiedades/Si_relajado.cif --cij ELASTIC_C.dat

Quasi-harmonic thermal expansion: V(T), α(T), a(T), Cv and Cp:

    olla-dft qha QHA_entrada_Si.dat --cells 8 --natoms 2 --cubic -o salida_qha

Material datasheet (gathers everything found in a project folder into one
Markdown file):

    olla-dft datasheet carpeta_del_proyecto/ -o .

### What to look at

In `derived`, note that the elastic Debye temperature (~635 K for silicon,
645 K experimental) is NOT the one that comes out of the phonon DOS: they
are two different definitions, and the report says so.

In `qha`, note that α(T) comes out NEGATIVE below ~165 K. It is not a bug:
real silicon has negative thermal expansion below ~120 K. That the model
reproduces it — with the zero crossing shifted, because the frequencies
come from a machine-learned potential and not from DFPT — is precisely the
sign that the quasi-harmonic approximation is doing its job.

### Files

| File | What it is |
|---|---|
| `ELASTIC_C.dat` | elastic matrix of silicon (LDA, primitive cell), exactly as `olla-dft elastic --collect` writes it |
| `QHA_entrada_Si.dat` | input table for `olla-dft qha`: one line per volume, with V (Å³), E (eV) and the 48 frequencies (cm⁻¹) of a 2×2×2 supercell (8 primitive cells). Generated with MACE-MP-0 in seconds, not with DFPT |
| `qha_Si.png` | the figure produced by the `qha` command above |
| `ficha_Si.md` | sample output of `olla-dft datasheet` |
