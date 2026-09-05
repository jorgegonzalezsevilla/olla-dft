[**English**](README.md) · [Español](README.es.md)

# Silicon: bands, DOS and PDOS from start to finish

Complete silicon example: generate the Quantum ESPRESSO inputs, run them,
and produce journal-ready band-structure and density-of-states figures.

**Generate the inputs:**

    olla-dft gen Si.cif -p all -o . --insulator

**Run** (pw.x, dos.x, projwfc.x and bands.x in order):

    ./run.sh

**Figures:**

    olla-dft plot . -o . --gap-label            # -> Si_bandas_dos
    olla-dft plot . -o . --gap-label --mono     # -> monochrome version

Result: indirect gap of 0.524 eV, VBM at Γ and CBM along Γ–X. The PDFs are
exactly the requested column width and carry the fonts embedded as
TrueType, ready to send to a journal.

### Files

| File | What it is |
|---|---|
| `scf.in`, `nscf.in`, `bands.in` | pw.x inputs written by `olla-dft gen` (scf, dense nscf and band path) |
| `KPATH.txt` | the high-symmetry path used in `bands.in`, with labels |
| `run.sh` | generated run script; it also expects `dos.in`, `projwfc.in` and `bands_pp.in`, which `gen -p all` writes and are not included here |
| `Si_bandas_dos.pdf`, `.png` | bands + DOS in a single figure, with the gap labelled |
| `Si_bandas_dos_mono.pdf` | the same figure in monochrome (`--mono`) |
| `Si_bandas.pdf`, `Si_dos.pdf` | bands and DOS as separate figures |
