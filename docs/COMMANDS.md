# Olla-DFT command reference

The 80 `olla-dft` subcommands, grouped by area, with their options. Generated from the code itself with `python tools/build_docs.py`; the same information is printed by `olla-dft COMMAND --help --language en` and, as a browsable page, by `olla-dft docs --language en`.

## Index

- **Getting started**: [`start`](#start), [`wizard`](#wizard), [`recetas`](#recetas), [`teoria`](#teoria), [`docs`](#docs), [`sistema`](#sistema), [`selftest`](#selftest), [`update`](#update)
- **Structures and inputs**: [`gen`](#gen), [`info`](#info), [`kpath`](#kpath), [`prim`](#prim), [`conv`](#conv), [`supercell`](#supercell), [`convert`](#convert)
- **Electronic structure**: [`bands`](#bands), [`dos`](#dos), [`plot`](#plot), [`gap`](#gap), [`fermi`](#fermi), [`effmass`](#effmass), [`wannier`](#wannier), [`unfold`](#unfold), [`topology`](#topology), [`hubbard`](#hubbard)
- **Spectra and response**: [`optics`](#optics), [`tddft`](#tddft), [`xanes`](#xanes), [`xps`](#xps), [`corehole`](#corehole), [`charge`](#charge), [`charges`](#charges), [`wf`](#wf), [`berry`](#berry)
- **Phonons, transport and thermal**: [`phonons`](#phonons), [`elph`](#elph), [`transport`](#transport), [`ballistic`](#ballistic), [`kappa`](#kappa), [`qha`](#qha), [`thermochem`](#thermochem), [`md`](#md), [`derived`](#derived)
- **Mechanics and stability**: [`converge`](#converge), [`eos`](#eos), [`elastic`](#elastic), [`strain`](#strain), [`layers`](#layers), [`xrd`](#xrd), [`exfoliate`](#exfoliate), [`gamma`](#gamma)
- **Surfaces, defects and chemistry**: [`surface`](#surface), [`defect`](#defect), [`interface`](#interface), [`adsorb`](#adsorb), [`eform`](#eform), [`align`](#align), [`esm`](#esm), [`echem`](#echem), [`neb`](#neb), [`amorphous`](#amorphous)
- **Automation and quality**: [`doctor`](#doctor), [`audit`](#audit), [`crosscheck`](#crosscheck), [`cost`](#cost), [`db`](#db), [`hull`](#hull), [`mlip`](#mlip), [`suggest`](#suggest), [`datasheet`](#datasheet), [`report`](#report), [`compare`](#compare), [`tune`](#tune), [`results`](#results), [`campaign`](#campaign), [`pseudos`](#pseudos)
- **Project**: [`project`](#project), [`resilient`](#resilient)
- **Appearance and configuration**: [`templates`](#templates), [`config`](#config)

## Getting started

### `start`

guided start to create a project without knowing the CLI

**Usage:** `olla-dft start [-h] [--project PROJECT] [--structure STRUCTURE] [--goal GOAL] [--name NAME] [--non-interactive] [--no-validate] [--language {es,en}]`

**Options:**

| Option | Description |
|---|---|
| `--project` | project folder (default: `.`) |
| `--structure` | CIF, POSCAR or pw.x input |
| `--goal` | relax, gap, dos, phonons, optics or scf |
| `--name` | display name of the project |
| `--non-interactive` | do not ask; requires --structure in a new project |
| `--no-validate` | do not run the initial validation |
| `--language {es,en}` | language of the guided start (default: es) |

### `wizard`

assistant: tell me WHAT you want to know and I tell you what to run, in order and with the commands

**Usage:** `olla-dft wizard [-h] [--goal GOAL] [--ask TEXTO] [--list] [--term TERM] [--no-glossary] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Arguments:**

- `file` — your structure (optional)

**Options:**

| Option | Description |
|---|---|
| `--goal` | goal key; they are listed with --list |
| `--ask TEXTO` | describe it in your own words (in Spanish, which is what the matcher understands), e.g. 'quiero saber si absorbe luz' |
| `--list` | list everything the assistant knows how to do |
| `--term` | what a term means |
| `--no-glossary` | do not explain the technical terms at the end of the answer |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |

### `recetas`

complete sessions from start to finish: which command comes after which and which file they pass to each other

**Usage:** `olla-dft recetas [-h] [--buscar TEXTO] [--script [ARCHIVO]] [receta]`

**Arguments:**

- `receta` — recipe key; with nothing, lists them all

**Options:**

| Option | Description |
|---|---|
| `--buscar TEXTO` | search for it in your own words, without knowing the key |
| `--script ARCHIVO` | write the recipe as a commented shell script, ready to edit |

### `teoria`

the physics behind a command: what it answers, the formulas it implements, which module they come from and where each number comes from

**Usage:** `olla-dft teoria [-h] [--all] [-o ARCHIVO.md] [comando]`

**Arguments:**

- `comando` — command to explain; with nothing, the index

**Options:**

| Option | Description |
|---|---|
| `--all` | the full document (all areas) |
| `-o, --output ARCHIVO.md` | save it as Markdown instead of printing it |

### `docs`

browsable reference of all subcommands, generated from the code itself

**Usage:** `olla-dft docs [-h] [-o OUTPUT] [--open] [--language {es,en}] [--both]`

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output HTML file (default: `olla-dft-docs.html`) |
| `--open` | open it in the browser when done |
| `--language {es,en}` | language of the reference interface (default: es) |
| `--both` | generate separate Spanish and English references |

### `sistema`

what Olla-DFT sees on this machine: encoding, where it stores the configuration, which QE binaries it finds and how to launch calculations here

**Usage:** `olla-dft sistema [-h]`

### `selftest`

check Olla-DFT against published values, not against itself

**Usage:** `olla-dft selftest [-h] [--full] [--mlip] [--only ONLY] [--list] [--pseudo-dir PSEUDO_DIR] [--pw-cmd PW_CMD] [--nproc NPROC] [-j JOBS] [--keep CARPETA]`

**Options:**

| Option | Description |
|---|---|
| `--full` | include the tests that actually run pw.x (about ten minutes) |
| `--mlip` | include separately the machine-learned potential test (requires MACE) |
| `--only` | only these tests, comma-separated |
| `--list` | list the tests and their references, without running anything |
| `--pseudo-dir` | pseudopotentials for the --full tests |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `-j, --jobs JOBS` | simultaneous tests (default: 1) |
| `--keep CARPETA` | keep the calculations here instead of deleting them |

**Physics:** [`olla-dft teoria selftest`](THEORY.md)

### `update`

check whether a newer Olla-DFT is published and, if so, install it after a confirmation; it never runs on its own

**Usage:** `olla-dft update [-h] [--check] [--yes] [--version TAG]`

**Options:**

| Option | Description |
|---|---|
| `--check` | only check and report, install nothing |
| `--yes` | do not ask; install directly if a newer version exists |
| `--version TAG` | install a specific version (e.g. v1.0.1) instead of the latest |

## Structures and inputs

### `gen`

generate pw.x inputs and post-processing

**Usage:** `olla-dft gen [-h] [-p {scf,relax,vc-relax,nscf,bands,dos,all,md}] [-o OUTDIR] [-k {coarse,fine,gamma,medium,very-fine}] [--kspacing KSPACING] [--kgrid N N N] [--band-points BAND_POINTS] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--insulator] [--primitive] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--prefix PREFIX] [--nspin {1,2}] [--mag MAG] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--soc] [--hubbard EL=U] [--hubbard-style {legacy,card}] [--charge Q] [--dipole [EJE]] [--nosym] [--functional {b3lyp,gaupbe,hse,pbe0}] [--exx-grid NxNxN] [--exx-fraction EXX_FRACTION] [--dt FS] [--nstep NSTEP] [--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}] [-T TEMPERATURE] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input, ...)

**Options:**

| Option | Description |
|---|---|
| `-p, --preset {scf,relax,vc-relax,nscf,bands,dos,all,md}` | calculation type (default: scf) |
| `-o, --outdir` | output folder (default: `.`) |
| `-k, --klevel {coarse,fine,gamma,medium,very-fine}` | k-mesh density (gamma/coarse/medium/fine/very-fine) |
| `--kspacing KSPACING` | k spacing in Å^-1 (overrides --klevel) |
| `--kgrid N` | explicit k-grid for scf/relax (three integers; overrides --kspacing and --klevel) |
| `--band-points BAND_POINTS` | points per segment of the k-path |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--insulator` | occupations='fixed' (insulators; default: smearing) |
| `--primitive` | reduce to the standardized primitive cell before generating |
| `--pseudo-dir` | pseudopotential folder (overrides config) |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--prefix` | calculation prefix (default: formula) |
| `--nspin {1,2}` | 2 turns on spin polarization (default: `1`) |
| `--mag` | starting magnetization: a number (0.5) or per element (Fe=0.7,O=0). Implies --nspin 2 |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion (van der Waals) correction |
| `--soc` | spin-orbit coupling: noncollinear calculation with lspinorb (requires fully relativistic pseudopotentials) |
| `--hubbard EL=U` | Hubbard U in eV per element, e.g. Ni=4.1. Can be repeated. To compute it instead of guessing it:  olla-dft hubbard --cycle |
| `--hubbard-style {legacy,card}` | legacy = lda_plus_u (QE <= 7.0), card = HUBBARD card (QE >= 7.1) (default: `legacy`) |
| `--charge Q` | total charge of the cell (tot_charge): +1 removes one electron, -1 adds one |
| `--dipole EJE` | dipole correction for polar slabs; with no value it uses the c axis. Places the sawtooth inside the vacuum |
| `--nosym` | turn off symmetry (nosym and noinv) |
| `--functional {b3lyp,gaupbe,hse,pbe0}` | hybrid functional: hse, pbe0, b3lyp or gaupbe. Costs one to two orders of magnitude more than PBE, and the report says so with numbers |
| `--exx-grid NxNxN` | q mesh for exact exchange (default 1x1x1). Must divide the k mesh |
| `--exx-fraction EXX_FRACTION` | fraction of exact exchange, if you want to change the functional's own |
| `--dt FS` | MD time step in fs (default: 1.0) |
| `--nstep NSTEP` | MD steps (default: 1000) |
| `--thermostat {none,rescaling,berendsen,andersen,initial,reduce-history}` | MD thermostat; none = NVE (default) |
| `-T, --temperature TEMPERATURE` | MD target temperature in K (default: 300) |

**Physics:** [`olla-dft teoria gen`](THEORY.md)

### `info`

structure and symmetry information

**Usage:** `olla-dft info [-h] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Physics:** [`olla-dft teoria info`](THEORY.md)

### `kpath`

high-symmetry path (seekpath)

**Usage:** `olla-dft kpath [-h] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Physics:** [`olla-dft teoria kpath`](THEORY.md)

### `prim`

standardized primitive cell

**Usage:** `olla-dft prim [-h] [-o OUTPUT] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output structure file (by default, a .cif named after the command) |

**Physics:** [`olla-dft teoria prim`](THEORY.md)

### `conv`

standardized conventional cell

**Usage:** `olla-dft conv [-h] [-o OUTPUT] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output structure file (by default, a .cif named after the command) |

**Physics:** [`olla-dft teoria conv`](THEORY.md)

### `supercell`

build a supercell

**Usage:** `olla-dft supercell [-h] [-o OUTPUT] file nx ny nz`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)
- `nx` — cell repetitions along a
- `ny` — cell repetitions along b
- `nz` — cell repetitions along c

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output structure file (by default, a .cif named after the command) |

**Physics:** [`olla-dft teoria supercell`](THEORY.md)

### `convert`

convert format (CIF/POSCAR/XYZ)

**Usage:** `olla-dft convert [-h] [-o OUTPUT_FLAG] file [output]`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)
- `output` — destination file; the format is inferred from the extension (.cif, .vasp, .xyz...)

**Options:**

| Option | Description |
|---|---|
| `-o, --output-flag` | output file (alternative to giving it positionally) |

## Electronic structure

### `bands`

analyze and plot the band structure

**Usage:** `olla-dft bands [-h] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--fat SELECTOR] [--fat-scale FAT_SCALE] [--projwfc ARCHIVO] [path]`

**Arguments:**

- `path` — calculation folder (or path to the .xml)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--prefix` | calculation prefix (detected automatically) |
| `--ref {auto,fermi,vbm,none}` | energy origin (default: auto) |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-6.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `6.0`) |
| `--no-plot` | only export the data, without generating the plot |
| `--dpi DPI` | resolution of the bitmap formats (default: `600`) |
| `--format` | comma-separated formats: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | visual template: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (or the path to your own JSON) |
| `--size {paper,poster,presentation}` | type scale: paper / presentation / poster |
| `--font {sans,serif,latex}` | font family (latex = Computer Modern) |
| `--usetex` | render the text with real LaTeX |
| `--palette` | palette: grayscale, okabe-ito, okabe-ito-dark, or comma-separated hex colors |
| `--background` | background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | publisher's column widths (default: `generic`) |
| `--width` | width: single / onehalf / double, or a number in mm |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | monochrome: black ink and line patterns (for journals that charge for color) |
| `--dashes {auto,always,never}` | line patterns as a secondary encoding (default: `auto`) |
| `--title` | title inside the figure (none by default: in a paper the text goes in the caption) |
| `--gap-label` | annotate the gap value inside the plot |
| `--panel` | panel label, e.g. '(a)' |
| `--fat SELECTOR` | fatbands: weight of an orbital on each band. For example Ni-d, Si-p, O, d or atomo:3. Needs the projwfc.x output of the SAME bands calculation |
| `--fat-scale FAT_SCALE` | size of the fatband markers (default: `55.0`) |
| `--projwfc ARCHIVO` | projwfc.x output (by default projwfc.out in the same folder) |

**Physics:** [`olla-dft teoria bands`](THEORY.md)

### `dos`

analyze and plot DOS and PDOS

**Usage:** `olla-dft dos [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [--dband EL[-ORB]] [--dband-emax eV] [path]`

**Arguments:**

- `path` — calculation folder (or path to the .xml)

**Options:**

| Option | Description |
|---|---|
| `--mode {orbital,element,total}` | how to decompose the PDOS (default: `orbital`) |
| `-o, --outdir` | output folder (default: `.`) |
| `--prefix` | calculation prefix (detected automatically) |
| `--ref {auto,fermi,vbm,none}` | energy origin (default: auto) |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-6.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `6.0`) |
| `--no-plot` | only export the data, without generating the plot |
| `--dpi DPI` | resolution of the bitmap formats (default: `600`) |
| `--format` | comma-separated formats: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | visual template: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (or the path to your own JSON) |
| `--size {paper,poster,presentation}` | type scale: paper / presentation / poster |
| `--font {sans,serif,latex}` | font family (latex = Computer Modern) |
| `--usetex` | render the text with real LaTeX |
| `--palette` | palette: grayscale, okabe-ito, okabe-ito-dark, or comma-separated hex colors |
| `--background` | background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | publisher's column widths (default: `generic`) |
| `--width` | width: single / onehalf / double, or a number in mm |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | monochrome: black ink and line patterns (for journals that charge for color) |
| `--dashes {auto,always,never}` | line patterns as a secondary encoding (default: `auto`) |
| `--title` | title inside the figure (none by default: in a paper the text goes in the caption) |
| `--gap-label` | annotate the gap value inside the plot |
| `--panel` | panel label, e.g. '(a)' |
| `--dband EL[-ORB]` | center, width and filling of a projected band, e.g. Pt (uses d) or Ni-p. This is the descriptor that correlates with the adsorption energy |
| `--dband-emax eV` | upper cutoff of the integral, relative to the Fermi level |

**Physics:** [`olla-dft teoria dos`](THEORY.md)

### `plot`

combined bands + DOS plot

**Usage:** `olla-dft plot [-h] [--mode {orbital,element,total}] [-o OUTDIR] [--prefix PREFIX] [--ref {auto,fermi,vbm,none}] [--emin EMIN] [--emax EMAX] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--dashes {auto,always,never}] [--title TITLE] [--gap-label] [--panel PANEL] [path]`

**Arguments:**

- `path` — calculation folder (or path to the .xml)

**Options:**

| Option | Description |
|---|---|
| `--mode {orbital,element,total}` | how to decompose the PDOS (default: `orbital`) |
| `-o, --outdir` | output folder (default: `.`) |
| `--prefix` | calculation prefix (detected automatically) |
| `--ref {auto,fermi,vbm,none}` | energy origin (default: auto) |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-6.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `6.0`) |
| `--no-plot` | only export the data, without generating the plot |
| `--dpi DPI` | resolution of the bitmap formats (default: `600`) |
| `--format` | comma-separated formats: pdf,png,svg,eps,tif (default: `pdf,png`) |
| `-t, --template` | visual template: dark, journal, latex, latex-true, minimal, mono, mono-latex, poster, slides (or the path to your own JSON) |
| `--size {paper,poster,presentation}` | type scale: paper / presentation / poster |
| `--font {sans,serif,latex}` | font family (latex = Computer Modern) |
| `--usetex` | render the text with real LaTeX |
| `--palette` | palette: grayscale, okabe-ito, okabe-ito-dark, or comma-separated hex colors |
| `--background` | background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | publisher's column widths (default: `generic`) |
| `--width` | width: single / onehalf / double, or a number in mm |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | monochrome: black ink and line patterns (for journals that charge for color) |
| `--dashes {auto,always,never}` | line patterns as a secondary encoding (default: `auto`) |
| `--title` | title inside the figure (none by default: in a paper the text goes in the caption) |
| `--gap-label` | annotate the gap value inside the plot |
| `--panel` | panel label, e.g. '(a)' |

**Physics:** [`olla-dft teoria plot`](THEORY.md)

### `gap`

band gap report only (fast)

**Usage:** `olla-dft gap [-h] [--prefix PREFIX] [path]`

**Arguments:**

- `path` — calculation folder (or path to the .xml)

**Options:**

| Option | Description |
|---|---|
| `--prefix` | calculation prefix (detected automatically) |

**Physics:** [`olla-dft teoria gap`](THEORY.md)

### `fermi`

export the Fermi surface as BXSF

**Usage:** `olla-dft fermi [-h] [-o OUTDIR]`

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `transporte`) |

**Physics:** [`olla-dft teoria fermi`](THEORY.md)

### `effmass`

effective mass by parabolic fit of the bands

**Usage:** `olla-dft effmass [-h] [-o OUTDIR] [--bands-dir BANDS_DIR] [--collect] [--run] [--half-width HALF_WIDTH] [--points POINTS] [--window WINDOW] [--min-points MIN_POINTS] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `masa_efectiva`) |
| `--bands-dir` | folder with a bands calculation already done (VBM and CBM come from there) |
| `--collect` | read the fine calculation that has already run |
| `--run` | run the fine calculation as soon as it is prepared |
| `--half-width HALF_WIDTH` | half-width of each line in Å⁻¹ (default: `0.06`) |
| `--points POINTS` | k-points per line (odd) (default: `21`) |
| `--window WINDOW` | half-width of the quick fit along the path, in Å⁻¹ on each side of the extremum (by default, half the parabolic limit: ±0.06) |
| `--min-points MIN_POINTS` | minimum points for the quick fit (default: `7`) |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--timeout TIMEOUT` | time limit in seconds for each pw.x run |

**Physics:** [`olla-dft teoria effmass`](THEORY.md)

### `wannier`

Wannier functions: interpolate bands, centers and spread, without needing wannier90

**Usage:** `olla-dft wannier [-h] [-o OUTDIR] [-g NxNxN] [-p SITIO:ORBITAL] [--bands BANDS] [--exclude 5-8] [--window MIN:MAX] [--frozen MIN:MAX] [--no-minimize] [--iterations ITERATIONS] [--points POINTS] [--dft-bands DIR] [--no-dft-bands] [--dos N] [--sigma SIGMA] [--run] [--collect] [--pw-cmd PW_CMD] [--pw2wan-cmd PW2WAN_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `wannier`) |
| `-g, --grid NxNxN` | FULL k-point mesh (default 4x4x4). It sets the quality of the interpolation |
| `-p, --projections SITIO:ORBITAL` | trial orbitals: 'Si:sp3', 'O:p;Ti:d', 'f=0.125,0.125,0.125:s'. Several separated by ';'. With 'auto', s and p are placed on every atom (default: `auto`) |
| `--bands BANDS` | nscf bands (default: as many as needed) |
| `--exclude 5-8` | bands that do NOT enter the wannierization |
| `--window MIN:MAX` | outer disentanglement window in eV: which bands the subspace can be chosen from. Needed when the bands are entangled with others (conduction, metals) |
| `--frozen MIN:MAX` | frozen window in eV: the bands inside are reproduced EXACTLY. Usually the valence plus the piece of conduction you care about |
| `--no-minimize` | stay in the projection gauge, without minimizing the spread |
| `--iterations ITERATIONS` | minimization steps (default 500) |
| `--points POINTS` | points per segment of the interpolated path (default: `30`) |
| `--dft-bands DIR` | folder with the DFT bands calculation to compare against; without it there is no real validation |
| `--no-dft-bands` | with --run, skip step 4 (bands) |
| `--dos N` | also, interpolated DOS on an NxNxN mesh |
| `--sigma SIGMA` | broadening of the interpolated DOS (eV) (default: `0.05`) |
| `--run` | launch the four steps in order |
| `--collect` | analyze what has already run |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--pw2wan-cmd` | pw2wannier90.x executable (default: next to pw.x) |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--timeout TIMEOUT` | time limit in seconds for each pw.x run |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kgrid NxNxN` | mesh of the initial scf |
| `--insulator` | occupations='fixed' (insulators; default: smearing) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria wannier`](THEORY.md)

### `unfold`

unfold the bands of a supercell onto the primitive Brillouin zone

**Usage:** `olla-dft unfold [-h] [-o OUTDIR] [--prefix PREFIX] [--bands BANDS] [--spin {up,dw}] [--emin EMIN] [--emax EMAX] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path primitive`

**Arguments:**

- `path` — folder of the supercell's bands calculation
- `primitive` — structure of the PRIMITIVE cell

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--prefix` | calculation prefix (detected automatically) |
| `--bands BANDS` | how many bands to unfold (from the lowest) |
| `--spin {up,dw}` | spin channel to unfold if the calculation is lsda (ONE channel is unfolded per run; by default, up) |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-6.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `6.0`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria unfold`](THEORY.md)

### `topology`

Chern number and Wilson loops of a Wannier model

**Usage:** `olla-dft topology [-h] (--occupied N | --fermi EV) [-g NxN] [--plane {xy,xz,yz}] [--fixed K] [--gap-tol EV] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] MODELO`

**Arguments:**

- `MODELO` — *_hr.dat file or folder containing WANNIER_hr.dat

**Options:**

| Option | Description |
|---|---|
| `--occupied N` | number of occupied bands of the isolated subspace |
| `--fermi EV` | Fermi level; rejected if it crosses a band |
| `-g, --grid NxN` | periodic mesh of the 2D section (default: 40x40) |
| `--plane {xy,xz,yz}` | oriented plane of the BZ section (default: xy) |
| `--fixed K` | perpendicular fractional coordinate (default: 0) |
| `--gap-tol EV` | minimum direct gap to accept the invariant (default: 1e-8) |
| `-o, --outdir` | output folder (default: `topology`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria topology`](THEORY.md)

### `hubbard`

Hubbard U by linear response (hp.x), instead of copying it from a paper

**Usage:** `olla-dft hubbard [-h] [-o OUTDIR] [--species SPECIES] [--qgrid QGRID] [--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}] [--hubbard-style {legacy,card}] [--cycle] [--max-iter MAX_ITER] [--tol TOL] [--mixing MIXING] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--intersite] [--v-threshold eV] file`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `hubbard`) |
| `--species` | species to perturb, comma-separated. By default, the transition metals and rare earths in the structure |
| `--qgrid` | q mesh of the linear response; equivalent to a supercell of nq1*nq2*nq3 cells (default: `2x2x2`) |
| `--projection {atomic,ortho-atomic,norm-atomic,wannier,pseudo}` | projection scheme. The U is ONLY valid with the same scheme it was computed with (default: `ortho-atomic`) |
| `--hubbard-style {legacy,card}` | DFT+U syntax of the scf: legacy = lda_plus_u (QE <= 7.0), card = HUBBARD card (QE >= 7.1, where the old syntax is an error) (default: `legacy`) |
| `--cycle` | full self-consistency cycle: scf -> hp.x -> scf with the new U, until it stops changing |
| `--max-iter MAX_ITER` | maximum iterations of the scf -> hp.x -> scf cycle with --cycle (default: 6) |
| `--tol TOL` | change in eV below which it is considered converged (default: `0.05`) |
| `--mixing MIXING` | step damping; lower it to 0.5 if it oscillates (default: `1.0`) |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--metal` | metallic system: smeared occupations instead of fixed ones |
| `--nspin {1,2}` | 2 turns on spin polarization (default: 1) |
| `--mag` | starting magnetization: a number (0.5) or per element (Fe=0.7,O=0). Implies --nspin 2 |
| `--intersite` | besides the U values, read the intersite V that hp.x already writes and generate the HUBBARD card for QE >= 7.1 |
| `--v-threshold eV` | V below this is neither listed nor written (default: `0.01`) |

**Physics:** [`olla-dft teoria hubbard`](THEORY.md)

## Spectra and response

### `optics`

ε(ω), absorption and Tauc with epsilon.x (NC pseudopotentials)

**Usage:** `olla-dft optics [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--wmax WMAX] [--smear SMEAR] [--metal] [--suite] [--tauc {direct,indirect}] [--scissor SCISSOR] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `opticas`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--wmax WMAX` | maximum energy of the spectrum (eV) (default: `20.0`) |
| `--smear SMEAR` | interband broadening (eV) (default: `0.1`) |
| `--metal` | metallic system (smeared occupations) |
| `--suite` | also export an interchange JSON for the other apps in the suite |
| `--tauc {direct,indirect}` | transition type for the Tauc plot (default: `direct`) |
| `--scissor SCISSOR` | rigid gap shift in eV (experimental or GW gap minus the computed gap); shifts ε2 and rebuilds ε1 via Kramers-Kronig |

**Physics:** [`olla-dft teoria optics`](THEORY.md)

### `tddft`

optical absorption with TDDFPT: lets the excited electron and its hole see each other

**Usage:** `olla-dft tddft [-h] [-o OUTDIR] [--method {lanczos,davidson}] [--iter ITER] [--pol {1,2,3,4}] [--states STATES] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--scissor SCISSOR] [--extrapolation {no,constant,osc}] [--tamm-dancoff] [--rpa] [--gamma] [--gap GAP] [--compare OPTICS.dat] [--nbnd NBND] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `tddft`) |
| `--method {lanczos,davidson}` | lanczos gives the whole spectrum; davidson gives the first excitations one by one (default: `lanczos`) |
| `--iter ITER` | Lanczos iterations: they set the resolution (default: `500`) |
| `--pol {1,2,3,4}` | 1/2/3 = xx/yy/zz, 4 = full tensor (default: `4`) |
| `--states STATES` | excitations to look for (davidson) (default: `10`) |
| `--emin EMIN` | lower limit of the energy axis (eV) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `15.0`) |
| `--broadening BROADENING` | broadening in eV (default 0.05). With --collect it sets the exciton detection threshold; if omitted it is read from spectrum.in |
| `--scissor SCISSOR` | rigid shift of the empty bands in eV (lanczos only): compensates the underestimated gap |
| `--extrapolation {no,constant,osc}` | extrapolation of the Lanczos chain in the spectrum: no, constant or osc (default: osc) |
| `--tamm-dancoff` | Tamm-Dancoff approximation: cheaper, not exact |
| `--rpa` | turn off the xc kernel, to see how much it contributes |
| `--gamma` | force K_POINTS gamma. Detected automatically when the structure is a molecule |
| `--gap GAP` | independent-particle gap in eV, to detect whether there is a bound exciton |
| `--compare OPTICS.dat` | overlay the spectrum from 'olla-dft optics' |
| `--nbnd NBND` | number of bands in the scf (by default, whatever pw.x decides) |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF; can be repeated |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--metal` | metallic system: smeared occupations instead of fixed ones |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria tddft`](THEORY.md)

### `xanes`

XANES/NEXAFS: near-edge X-ray absorption (xspectra.x)

**Usage:** `olla-dft xanes [-h] [-o OUTDIR] [--element ELEMENT] [--site SITE] [--edge EDGE] [--core-hole UPF] [--polarization POLARIZATION] [--average] [--emin EMIN] [--emax EMAX] [--broadening BROADENING] [--r-paw R_PAW] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `xanes`) |
| `--element` | absorbing element |
| `--site SITE` | which atom of that element (from 0) |
| `--edge` | edge: K, L1, L2, L3 or L23 (the ones xspectra.x computes; M edges are not) (default: `K`) |
| `--core-hole UPF` | core-hole pseudopotential (olla-dft corehole) |
| `--polarization` | electric field direction, e.g. '0 0 1' (default: `1 0 0`) |
| `--average` | three orthogonal directions and their average: what corresponds to a powder sample |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-10.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `30.0`) |
| `--broadening BROADENING` | broadening in eV (xgamma) (default: `0.8`) |
| `--r-paw R_PAW` | PAW sphere radius of the absorber for xspectra.x, in bohr (default: 3.0) |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--metal` | metallic system: smeared occupations instead of fixed ones |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria xanes`](THEORY.md)

### `xps`

core-level shifts (initial state)

**Usage:** `olla-dft xps [-h] [-o OUTDIR] [--core-hole EL=UPF] [--collect] [--suite] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `xps`) |
| `--core-hole EL=UPF` | core-hole pseudopotential, e.g. Si=Si.star1s.UPF. Can be repeated. Without it, initial_state.x returns a table of zeros |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--suite` | also export an interchange JSON for the other apps in the suite |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--metal` | metallic system: smeared occupations instead of fixed ones |

**Physics:** [`olla-dft teoria xps`](THEORY.md)

### `corehole`

generate the normal + core-hole pseudopotential pair (ld1.x), for XPS and XANES

**Usage:** `olla-dft corehole [-h] [-o OUTDIR] [--edge EDGE] [--functional FUNCTIONAL] [--rcut RCUT] [--rel {0,1,2}] [--semicore] [--pseudotype {1,2,3}] [--plain] [--only-inputs] [--projectors {1,2}] [--ld1-cmd LD1_CMD] [--core-wfc UPF] [--orbital ORBITAL] [--output OUTPUT] [element]`

**Arguments:**

- `element` — element symbol, e.g. Si

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `pseudos`) |
| `--edge` | edge/level of the hole: K (1s), L1 (2s), L23 (2p), M1, M23, M45 (default: `K`) |
| `--functional` | functional of the pseudopotential; it must be the same one you will run pw.x with (default: `PBE`) |
| `--rcut RCUT` | cutoff radius in bohr (by default, one per row of the periodic table) |
| `--rel {0,1,2}` | 0 nonrelativistic, 1 scalar, 2 full |
| `--semicore` | put the (n-1)s(n-1)p shell in the valence |
| `--pseudotype {1,2,3}` | 1 and 2 are norm-conserving, 3 is ultrasoft (default: `2`) |
| `--plain` | generate ONLY the normal pseudopotential, without the core-hole one. Useful to get a consistent pseudo for an element that does not support it |
| `--only-inputs` | write the ld1.x inputs without running them |
| `--projectors {1,2}` | GIPAW projectors per channel. XSpectra recommends 2, but with 2 the pseudo comes out ultrasoft and --rcut usually has to be tuned by hand (default: `1`) |
| `--ld1-cmd` | path to ld1.x |
| `--core-wfc UPF` | instead of generating: extract the core wavefunction from a UPF in the format read by xspectra.x |
| `--orbital` | orbital to check, e.g. 1S |
| `--output` | output file for --core-wfc |

**Physics:** [`olla-dft teoria corehole`](THEORY.md)

### `charge`

charge density / ELF / spin with pp.x

**Usage:** `olla-dft charge [-h] [-o OUTDIR] [--field {density,elf,spin,potential,vtotal}] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Arguments:**

- `path` — calculation folder

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--field {density,elf,spin,potential,vtotal}` | field to compute with pp.x: density, elf, spin, potential or vtotal (default: density) |
| `--axis` | axis of the planar profile (default: `c`) |
| `--rerun` | rerun pp.x even if the cube file already exists |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria charge`](THEORY.md)

### `charges`

Löwdin/Bader charges and density difference

**Usage:** `olla-dft charges [-h] [--lowdin LOWDIN] [--bader BADER] [--difference CUBE [CUBE ...]] [--pseudo-dir PSEUDO_DIR] [--axis {0,1,2}] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Arguments:**

- `file` — structure (for Bader)

**Options:**

| Option | Description |
|---|---|
| `--lowdin` | projwfc.x output |
| `--bader` | density cube (plot_num=0) |
| `--difference CUBE` | total.cube part1.cube part2.cube ... |
| `--pseudo-dir` | folder with the calculation's UPF files: Z_valence for the 'neta' (net) column comes from there (overrides config) |
| `--axis {0,1,2}` | axis of the planar profile of the density difference: 0, 1 or 2 (default: 2) |
| `-o, --outdir` | output folder (default: `.`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria charges`](THEORY.md)

### `wf`

work function from a calculation with vacuum

**Usage:** `olla-dft wf [-h] [-o OUTDIR] [--axis AXIS] [--rerun] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Arguments:**

- `path` — calculation folder

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--axis` | vacuum axis: a/b/c (default c) |
| `--rerun` | rerun pp.x even if the cube already exists |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria wf`](THEORY.md)

### `berry`

Berry-phase polarization: spontaneous ΔP, Born charges and ferroelectricity

**Usage:** `olla-dft berry [-h] [-o OUTDIR] [--gdir {1,2,3}] [--nppstr NPPSTR] [--kperp NxN] [-r ARCHIVO] [--displace ATOMO:dx,dy,dz] [--nlambda NLAMBDA] [--run] [--collect] [--redo] [--pw-cmd PW_CMD] [--nproc NPROC] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kgrid NxNxN] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Arguments:**

- `file` — structure (the polar one, if there is a path)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `berry`) |
| `--gdir {1,2,3}` | direction: reciprocal-lattice vector (default 3) |
| `--nppstr NPPSTR` | points per k string (default 9); raise it until the phase stops changing |
| `--kperp NxN` | mesh perpendicular to the string (default 6x6) |
| `-r, --reference ARCHIVO` | reference structure, normally the centrosymmetric one: an adiabatic path to the polar one is interpolated and ΔP is the spontaneous polarization |
| `--displace ATOMO:dx,dy,dz` | displacement path of one atom, in Å; the slope of P gives the Born effective charge |
| `--nlambda NLAMBDA` | points along the path (default 5) |
| `--run` | run the calculations as soon as the inputs are prepared |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--redo` | also redo the calculations that were already finished |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--timeout TIMEOUT` | time limit in seconds for each pw.x run |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kgrid NxNxN` | k mesh of the scf, e.g. 6x6x6 (by default, from the k spacing) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria berry`](THEORY.md)

## Phonons, transport and thermal

### `phonons`

DFPT phonons: dispersion, DOS, thermodynamics, IR

**Usage:** `olla-dft phonons [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--qgrid QGRID] [--gamma] [--raman] [--laser LASER] [--suite] [--tscan T1,T2,...] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `fonones`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--qgrid` | q mesh, e.g. 2x2x2 |
| `--gamma` | Γ only with dynmat.x: frequencies and IR activities |
| `--raman` | also Raman tensors and intensities at Γ (lraman; norm-conserving pseudopotentials only, and quite a bit more expensive) |
| `--laser LASER` | laser wavelength in nm to simulate the Raman spectrum (default: `532.0`) |
| `--suite` | also export an interchange JSON (only with --gamma) for the FTIR and Raman apps |
| `--tscan T1,T2,...` | ELECTRONIC temperature sweep in K: repeats the phonons with fermi-dirac smearing at each one and checks whether an imaginary mode stabilizes on heating (charge density waves, structural transitions) |

**Physics:** [`olla-dft teoria phonons`](THEORY.md)

### `elph`

electron-phonon coupling: lambda, Tc and a real tau for transport

**Usage:** `olla-dft elph [-h] [-o OUTDIR] [--qgrid QGRID] [--kgrid KGRID] [--kgrid-nscf KGRID_NSCF] [--nsigma NSIGMA] [--sigma SIGMA] [--degauss DEGAUSS] [--debye DEBYE] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `elph`) |
| `--qgrid` | DFPT q mesh, e.g. 2x2x2 (default: 2x2x2) |
| `--kgrid` | k mesh of the scf |
| `--kgrid-nscf` | k mesh of the dense nscf; by default, twice the scf one rounded to a multiple of the q mesh |
| `--nsigma NSIGMA` | how many broadenings ph.x sweeps for lambda (el_ph_nsigma; default: 10) |
| `--sigma SIGMA` | step of the broadening sweep, in Ry (default: `0.005`) |
| `--degauss DEGAUSS` | scf smearing in Ry (default: 0.02) |
| `--debye DEBYE` | Debye temperature in K, to mark the regime where the tau formula is valid |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria elph`](THEORY.md)

### `transport`

Seebeck, sigma/tau and power factor (CRTA)

**Usage:** `olla-dft transport [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--grid GRID] [--temperatures TEMPERATURES] [--mu-span MU_SPAN] [--metal] [--nspin {1,2}] [--mag MAG] [--spin-resolved] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `transporte`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--grid` | nscf mesh, e.g. 16x16x16 |
| `--temperatures` | comma-separated temperatures in K (default: `300`) |
| `--mu-span MU_SPAN` | chemical-potential range around E_F (eV) (default: `1.0`) |
| `--metal` | metallic system: smeared occupations instead of fixed ones |
| `--nspin {1,2}` | 2 turns on spin polarization in scf and nscf (needed for --spin-resolved) (default: `1`) |
| `--mag` | starting magnetization, e.g. Fe=0.7 (implies --nspin 2) |
| `--spin-resolved` | separate the two spin channels (two-current model) and give the conductivity polarization and the spin thermopower |

**Physics:** [`olla-dft teoria transport`](THEORY.md)

### `ballistic`

Landauer ballistic conductance (pwcond.x), for nanocontacts and molecules between electrodes

**Usage:** `olla-dft ballistic [-h] [--scatterer SCATTERER] [-o OUTDIR] [--ikind {0,1}] [--emin EMIN] [--emax EMAX] [--points POINTS] [--nz1 NZ1] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [file]`

**Arguments:**

- `file` — electrode: the cell periodic in z

**Options:**

| Option | Description |
|---|---|
| `--scatterer` | scattering region (the molecule or the defect). Without it only the complex bands come out |
| `-o, --outdir` | output folder (default: `balistico`) |
| `--ikind {0,1}` | 0 = complex bands only, 1 = conductance with the same electrode on both sides (default: 1 if there is --scatterer, 0 otherwise). Different electrodes (ikind=2 of pwcond.x) are not supported |
| `--emin EMIN` | lower limit of the energy axis (eV) (default: `-3.0`) |
| `--emax EMAX` | upper limit of the energy axis (eV) (default: `3.0`) |
| `--points POINTS` | number of energies in the conductance sweep (default: 61) |
| `--nz1 NZ1` | z subdivisions of each pwcond.x slice (nz1; default: 3) |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF; can be repeated |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria ballistic`](THEORY.md)

### `kappa`

lattice thermal conductivity: fc3, phonon Boltzmann equation and mean free path

**Usage:** `olla-dft kappa [-h] [-o OUTDIR] [--dim NxNxN] [--dim-fc2 NxNxN] [--distance DISTANCE] [--mesh MESH] [--temps TEMPS] [--isotopes] [--grain UM] [--model MODEL] [--collect] [--force] [--metal] [--pseudo-dir PSEUDO_DIR] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file`

**Arguments:**

- `file` — structure (primitive cell)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `kappa`) |
| `--dim NxNxN` | fc3 supercell (default 2x2x2). This is the expensive part: the number of configurations grows quickly |
| `--dim-fc2 NxNxN` | LARGER supercell for the harmonic part only, which is cheap and needs longer range |
| `--distance DISTANCE` | finite displacement in Å (default 0.03) |
| `--mesh MESH` | q mesh for the Boltzmann equation (default 13) |
| `--temps` | temperatures: 100:800:8 or 300,500,700 (default: `100:800:8`) |
| `--isotopes` | add isotope scattering with natural abundances (in Si it is ~10 %%) |
| `--grain UM` | grain size in µm: adds boundary scattering |
| `--model` | compute the forces with a machine-learned potential (mace, chgnet, m3gnet) instead of pw.x: seconds instead of hours, but the absolute value may be far off |
| `--collect` | read the forces already computed and solve |
| `--force` | write the inputs even if there are a huge number of them |
| `--metal` | metallic system (smeared occupations in the fc2/fc3 scf runs). Without it occupations='fixed' is used, which is right for insulators |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 (default: `0.35`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria kappa`](THEORY.md)

### `qha`

quasi-harmonic: thermal expansion and a(T)

**Usage:** `olla-dft qha [-h] [-o OUTDIR] [--natoms NATOMS] [--cells CELLS] [--cubic] [--structure CIF] [--tmax TMAX] [--dt DT] [--temp TEMP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] data`

**Arguments:**

- `data` — table: V(A^3) E(eV) w1 w2 ... per volume

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--natoms NATOMS` | atoms per primitive cell, for the per-atom quantities (default: 1) |
| `--cells CELLS` | primitive cells per supercell of the modes (default: `1`) |
| `--cubic` | also a(T). Without --structure it is only V_prim^(1/3) |
| `--structure CIF` | structure of the material: with it a(T) is converted to the CONVENTIONAL lattice parameter (factor 4 in fcc/diamond, 2 in bcc) and cubic symmetry is detected |
| `--tmax TMAX` | maximum temperature of the T grid in K (default: 1000) |
| `--dt DT` | step: integration step in fs (amorphous) or temperature step in K (qha) (default: `5.0`) |
| `--temp TEMP` | working temperature in K (default: `300.0`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria qha`](THEORY.md)

### `thermochem`

ZPE, entropy and free energy: from a DFT energy to one comparable with experiment

**Usage:** `olla-dft thermochem [-h] [--phase {solido,adsorbato,gas,transicion}] [--structure STRUCTURE] [--temp TEMP] [--pressure PRESSURE] [--symmetry SYMMETRY] [--multiplicity MULTIPLICITY] [--floor FLOOR] [--energy ENERGY] [-o OUTDIR] freqs`

**Arguments:**

- `freqs` — file of frequencies in cm-1, or the comma-separated list

**Options:**

| Option | Description |
|---|---|
| `--phase {solido,adsorbato,gas,transicion}` | gas adds translations and rotations; transicion requires exactly one imaginary frequency (default: `solido`) |
| `--structure` | structure (needed for the gas phase) |
| `--temp TEMP` | working temperature in K (default: `298.15`) |
| `--pressure PRESSURE` | in bar (default: `1.0`) |
| `--symmetry SYMMETRY` | symmetry number of the point group: 2 for H2O and O2, 3 for NH3, 12 for CH4 (default: `1`) |
| `--multiplicity MULTIPLICITY` | spin multiplicity of the ground state (default: `1`) |
| `--floor FLOOR` | raise the modes below this value (cm-1); 100 is the usual choice |
| `--energy ENERGY` | E_DFT in eV, to give G(T) |
| `-o, --outdir` | output folder |

**Physics:** [`olla-dft teoria thermochem`](THEORY.md)

### `md`

analyze a molecular dynamics trajectory: g(r), diffusion and vibrational spectrum

**Usage:** `olla-dft md [-h] [-o OUTDIR] [--skip SKIP] [--rmax RMAX] [--bins BINS] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] path`

**Arguments:**

- `path` — pw.x output with calculation='md', or its folder

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--skip SKIP` | initial steps to discard (equilibration) |
| `--rmax RMAX` | g(r) cutoff in Å; by default, half a cell edge, which is as far as the normalization is valid |
| `--bins BINS` | number of bins in the g(r) histogram (default: 200) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria md`](THEORY.md)

### `derived`

Debye, sound velocities and Slack from the Cij

**Usage:** `olla-dft derived [-h] [--cij CIJ] [--temp TEMP] [-o OUTDIR] file`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `--cij` | file with the elastic matrix (default: `ELASTIC_C.dat`) |
| `--temp TEMP` | working temperature in K (default: `300.0`) |
| `-o, --outdir` | folder where DERIVED.dat is written (default: `.`) |

**Physics:** [`olla-dft teoria derived`](THEORY.md)

## Mechanics and stability

### `converge`

convergence tests for cutoffs and k mesh

**Usage:** `olla-dft converge [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {ecutwfc,ecutrho,kmesh}] [--values VALUES] [--threshold THRESHOLD] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `convergencia`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `-k, --kind {ecutwfc,ecutrho,kmesh}` | which parameter is swept (default: ecutwfc) |
| `--values` | comma-separated values; for kmesh it accepts 8x8x8 |
| `--threshold THRESHOLD` | convergence threshold in meV/atom (default: 1) |

**Physics:** [`olla-dft teoria converge`](THEORY.md)

### `eos`

E–V equation of state and bulk modulus

**Usage:** `olla-dft eos [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--npoints NPOINTS] [--scale SCALE] [--span SPAN] [--equation {birch-murnaghan,murnaghan,vinet}] [--relax-ions] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `eos`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--npoints NPOINTS` | number of volumes (default: 9) |
| `--scale SCALE` | linear factor to center the sweep on (returned by 'olla-dft mlip scan') (default: `1.0`) |
| `--span SPAN` | relative volume variation on each side (default: 0.10) |
| `--equation {birch-murnaghan,murnaghan,vinet}` | equation to plot (default: `birch-murnaghan`) |
| `--relax-ions` | relax internal positions at each volume |

**Physics:** [`olla-dft teoria eos`](THEORY.md)

### `elastic`

elastic constants and mechanical properties

**Usage:** `olla-dft elastic [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--delta DELTA] [--npoints NPOINTS] [--2d] [--thickness A] [--ion-mode {auto,relax,fixed}] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `elastic`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--delta DELTA` | maximum applied strain (default: 0.010 = 1 %%) |
| `--npoints NPOINTS` | nonzero strains per component, even number (default: 4) |
| `--2d` | sheet: constants in N/m (not GPa), only ε1, ε2 and ε6, and 2D Born criteria |
| `--thickness A` | assumed thickness in Å to also give the GPa equivalent (a convention, not a measurement) |
| `--ion-mode {auto,relax,fixed}` | internal positions: auto = fixed for normal strains and relaxed for shears (recommended); relax = relax all; fixed = clamped-ion (default: `auto`) |

**Physics:** [`olla-dft teoria elastic`](THEORY.md)

### `strain`

strain sweep: gap, energy and moment as a function of the applied strain

**Usage:** `olla-dft strain [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}] [-r MIN:MAX:N] [--fixed-ions] [--relax-perp] [--nspin {1,2}] [--mag MAG] [--hubbard EL=U] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `strain`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `-m, --mode {biaxial,cizalla,hidrostatica,uniaxial-a,uniaxial-b,uniaxial-c}` | what gets strained (default: biaxial) |
| `-r, --range MIN:MAX:N` | range in PERCENT, e.g. -5:5:11 (from -5 %% to +5 %% in 11 points) (default: `-5:5:11`) |
| `--fixed-ions` | do not relax the internal positions at each strain (faster and less realistic) |
| `--relax-perp` | leave the axis perpendicular to the strained plane free (Poisson relaxation); essential for sheets |
| `--nspin {1,2}` | 2 turns on spin polarization (default: `1`) |
| `--mag` | starting magnetization (implies --nspin 2) |
| `--hubbard EL=U` | Hubbard U in eV per element |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion correction |

**Physics:** [`olla-dft teoria strain`](THEORY.md)

### `layers`

detect layers, basal spacing and interlayer gap

**Usage:** `olla-dft layers [-h] [--tol TOL] [--wavelength WAVELENGTH] [--slab ARCHIVO] [--vacuum VACUUM] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `--tol TOL` | bond tolerance on top of covalent radii (Å) (default: `0.45`) |
| `--wavelength` | radiation for the basal reflections (default CuKa) |
| `--slab ARCHIVO` | also write the monolayer with vacuum to this file |
| `--vacuum VACUUM` | monolayer vacuum in Å (default 20) |

**Physics:** [`olla-dft teoria layers`](THEORY.md)

### `xrd`

simulated powder diffractogram

**Usage:** `olla-dft xrd [-h] [-o OUTDIR] [--suite] [--basis {conventional,input}] [--wavelength WAVELENGTH] [--tt-min TT_MIN] [--tt-max TT_MAX] [--fwhm FWHM] [--size SIZE] [--biso BISO] [--exp EXP] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size-preset {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--suite` | also export an interchange JSON for the other apps in the suite |
| `--basis {conventional,input}` | cell in which the hkl are indexed: 'conventional' (default, the indices of the PDF cards) or 'input' (the cell of the file as is) |
| `--wavelength` | radiation: AgKa, CoKa, CrKa, CuKa, CuKa1, FeKa, MoKa or λ in Å (default: `CuKa`) |
| `--tt-min TT_MIN` | minimum 2θ (°) (default: `5.0`) |
| `--tt-max TT_MAX` | maximum 2θ (°) (default: `70.0`) |
| `--fwhm FWHM` | instrumental width (° 2θ, default 0.15) |
| `--size SIZE` | crystallite size in nm (Scherrer broadening) |
| `--biso BISO` | overall temperature factor B (Å²) |
| `--exp` | experimental diffractogram (2θ, I) to overlay |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--size-preset {paper,poster,presentation}` | type scale of the figure |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria xrd`](THEORY.md)

### `exfoliate`

exfoliation energy (bulk vs monolayer)

**Usage:** `olla-dft exfoliate [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--vacuum VACUUM] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--tol TOL] [--relax-slab] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `exfoliacion`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--vacuum VACUUM` | monolayer vacuum in Å (default 20) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion correction for both calculations |
| `--tol TOL` | bond tolerance for detecting the layers (Å) (default: `0.45`) |
| `--relax-slab` | relax the monolayer positions |

**Physics:** [`olla-dft teoria exfoliate`](THEORY.md)

### `gamma`

surface and cleavage energy by the Fiorentini–Methfessel linear fit, with convergence against slab thickness

**Usage:** `olla-dft gamma [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix N] [--relax] [--no-bulk] [--no-reduce] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `gamma`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `-m, --miller` | Miller indices of the facet, e.g. '1 1 1' (default: `1 0 0`) |
| `-l, --layers` | thicknesses to compute, comma-separated (default: 3,4,5,6). At least two are needed |
| `--vacuum VACUUM` | vacuum in Å (default: 20) |
| `--fix N` | freeze N bottom layers when relaxing |
| `--relax` | relax the positions (γ drops by 5 to 20 %%) |
| `--no-bulk` | do not compute the bulk separately; only the linear fit E_slab(N) = 2γA + N·E_bulk |
| `--no-reduce` | do not reduce the surface cell to the minimal one (it is reduced by default: same γ, much cheaper) |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion (van der Waals) correction: grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm or mbd |
| `--dipole` | dipole correction, for polar slabs |
| `--nspin {1,2}` | 2 turns on spin polarization (default: 1) |
| `--mag` | starting magnetization (implies --nspin 2) |

**Physics:** [`olla-dft teoria gamma`](THEORY.md)

## Surfaces, defects and chemistry

### `surface`

cut an (hkl) surface with vacuum

**Usage:** `olla-dft surface [-h] [-m MILLER] [-l LAYERS] [--vacuum VACUUM] [--fix FIX] [-o OUTPUT] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-m, --miller` | Miller indices, e.g. '1 1 1' or 1,1,1 (default: `1 0 0`) |
| `-l, --layers LAYERS` | number of atomic layers in the slab (default: 6) |
| `--vacuum VACUUM` | total vacuum in Å (default: `15.0`) |
| `--fix FIX` | bottom atomic planes to freeze |
| `-o, --output` | output file (CIF/POSCAR) |

**Physics:** [`olla-dft teoria surface`](THEORY.md)

### `defect`

create a point defect

**Usage:** `olla-dft defect [-h] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--supercell SUPERCELL] [--position POSITION] [-o OUTDIR] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-k, --kind {vacancy,substitution,interstitial}` | defect type: vacancy, substitution or interstitial (default: vacancy) |
| `--site SITE` | index of the affected atom (0-based) |
| `--new-element` | incoming species |
| `--supercell` | e.g. 3x3x3 |
| `--position` | fractional x,y,z (interstitial) |
| `-o, --outdir` | output folder (default: `defecto`) |

**Physics:** [`olla-dft teoria defect`](THEORY.md)

### `interface`

heterostructure: stack two materials with the smallest possible lattice strain

**Usage:** `olla-dft interface [-h] [-o OUTDIR] [--name NAME] [--max-index MAX_INDEX] [--tol TOL] [--max-atoms MAX_ATOMS] [--index INDEX] [--top TOP] [--list] [--separation SEPARATION] [--vacuum VACUUM] [--strain {first,second,both}] [--shift SHIFT] file1 file2`

**Arguments:**

- `file1` — bottom material (the substrate)
- `file2` — top material

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--name` | base name of the heterostructure files (default: heteroestructura) |
| `--max-index MAX_INDEX` | largest integer coefficient of the supercell; raising it finds more rotated cells but takes much longer (default: `4`) |
| `--tol TOL` | maximum accepted strain (0.05 = 5 %%) (default: `0.05`) |
| `--max-atoms MAX_ATOMS` | maximum number of atoms allowed in the interface supercell (default: 200) |
| `--index INDEX` | which of the candidates to build |
| `--top TOP` | how many candidates to list, from lowest to highest strain (default: 10) |
| `--list` | only list the candidates, without building anything |
| `--separation SEPARATION` | initial distance between layers in Å; by default, from the van der Waals radii |
| `--vacuum VACUUM` | vacuum above the heterostructure in Å (default: 20) |
| `--strain {first,second,both}` | who gets strained: the bottom one, the top one, or both halfway (default: `second`) |
| `--shift` | lateral shift of the top material, in fractions of the common cell |

**Physics:** [`olla-dft teoria interface`](THEORY.md)

### `adsorb`

adsorption sites on a slab and their energy

**Usage:** `olla-dft adsorb [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] --mol MOLECULA [--sites SITES] [--height HEIGHT] [--face {top,bottom}] [--rotations ROTATIONS] [--anchor ANCHOR] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--dipole] [--nspin {1,2}] [--mag MAG] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `adsorb`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--mol MOLECULA` | adsorbate: a name from the ASE database (CO2, H2O, CO, NH3, O2...) or a file with the molecule |
| `--sites` | site types to try (default: all three) |
| `--height HEIGHT` | initial height of the adsorbate above the site, in Å (default: 2.0) |
| `--face {top,bottom}` | side of the slab to adsorb on (default: `top`) |
| `--rotations ROTATIONS` | orientations to try by rotating about the normal (default: 1) |
| `--anchor ANCHOR` | atom of the molecule that sits on the site (0-based index; default: 0) |
| `--fixed-ions` | do not relax: scf only at the initial geometry |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion correction (almost mandatory for physisorption) |
| `--dipole` | dipole correction: a slab with an adsorbate on one side only is polar |
| `--nspin {1,2}` | 2 turns on spin polarization (default: 1) |
| `--mag` | starting magnetization (implies --nspin 2) |

**Physics:** [`olla-dft teoria adsorb`](THEORY.md)

### `eform`

formation energy of charged defects, transition levels and E_f vs ε_F diagram

**Usage:** `olla-dft eform [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [-k {vacancy,substitution,interstitial}] [--site SITE] [--new-element NEW_ELEMENT] [--position POSITION] [--supercell SUPERCELL] [-q CHARGES] [--epsilon EPSILON] [--correction {ninguna,makov-payne,lany-zunger}] [--mu EL=eV] [--align POT_DEF POT_PERF] [--dv DV] [--fixed-ions] [--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}] [--nspin {1,2}] [--mag MAG] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `formacion`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `-k, --kind {vacancy,substitution,interstitial}` | defect type (default: `vacancy`) |
| `--site SITE` | index of the affected atom in the supercell (0-based) |
| `--new-element` | incoming species |
| `--position` | fractional x,y,z (interstitial) |
| `--supercell` | supercell size (default: 2x2x2) |
| `-q, --charges` | comma-separated charge states, e.g. -2,-1,0,1,2 (default: `0`) |
| `--epsilon EPSILON` | dielectric constant of the material, to screen the image-charge correction |
| `--correction {ninguna,makov-payne,lany-zunger}` | finite-size correction scheme (default: `lany-zunger`) |
| `--mu EL=eV` | chemical potential per element, in eV per atom. Can be repeated |
| `--align ('POT_DEF', 'POT_PERF')` | two electrostatic-potential cube files (defective and pristine) for the ΔV term |
| `--dv DV` | ΔV alignment in eV, if you already have it computed |
| `--fixed-ions` | do not relax the defect in each charge state |
| `--vdw {grimme-d2,grimme-d3,DFT-D,ts-vdw,xdm,mbd}` | dispersion (van der Waals) correction: grimme-d2, grimme-d3, DFT-D, ts-vdw, xdm or mbd |
| `--nspin {1,2}` | 2 turns on spin polarization (default: 1) |
| `--mag` | starting magnetization (implies --nspin 2) |

**Physics:** [`olla-dft teoria eform`](THEORY.md)

### `align`

band alignment between two materials: ΔE_v, ΔE_c offsets and type I/II/III

**Usage:** `olla-dft align [-h] [--interface CARPETA] [--names NAMES] [--axis AXIS] [--window A] [--rerun] [-o OUTDIR] [--pw-cmd PW_CMD] [--nproc NPROC] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] a b`

**Arguments:**

- `a` — first material's calculation folder
- `b` — second material's calculation folder

**Options:**

| Option | Description |
|---|---|
| `--interface CARPETA` | interface folder; enables the rigorous Van de Walle-Martin method |
| `--names` | names for the report, comma-separated (by default, the folder names) |
| `--axis` | axis of the planar profile (default: `c`) |
| `--window A` | macroscopic-average window in Å (by default, one eighth of the cell) |
| `--rerun` | rerun pp.x even if the cube already exists |
| `-o, --outdir` | output folder (default: `alineamiento`) |
| `--pw-cmd` | pw.x executable for --run; the other QE binaries are located from its path |
| `--nproc NPROC` | number of MPI processes for the calculations launched with --run |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria align`](THEORY.md)

### `esm`

charged surfaces with effective screening medium: work function, capacitance and potential of zero charge

**Usage:** `olla-dft esm [-h] [-o OUTDIR] [--run] [--collect] [--pw-cmd PW_CMD] [--nproc NPROC] [-j N] [--redo] [--max-time T] [--estimate] [--timeout TIMEOUT] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--insulator] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--size {paper,poster,presentation}] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--aspect ASPECT] [--mono] [--bc {bc1,bc2,bc3}] [--charge CHARGE] [--field FIELD] [--esm-w WIDTH_ESM] [--nfit NFIT] file`

**Arguments:**

- `file` — structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | sweep folder (default: `esm`) |
| `--run` | run the calculations now, one after another |
| `--collect` | only analyze calculations that have already run |
| `--pw-cmd` | pw.x executable (overrides the configuration) |
| `--nproc NPROC` | MPI processes per calculation |
| `-j, --jobs N` | simultaneous calculations (default: 1). Without --nproc, the machine's threads are shared among them |
| `--redo` | also redo the calculations that were already finished |
| `--max-time T` | TOTAL time budget: 90m, 2h, 3600. Once it runs out no more are launched and the sweep can be resumed |
| `--estimate` | estimate how long the sweep will take and exit, using the history from 'olla-dft db' |
| `--timeout TIMEOUT` | limit in seconds per calculation |
| `--pseudo-dir` | pseudopotential folder |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff (Ry) |
| `--ecutrho ECUTRHO` | density cutoff (Ry) |
| `--kspacing KSPACING` | k spacing in Å^-1 |
| `--insulator` | occupations='fixed' |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure |
| `--size {paper,poster,presentation}` | figure size: paper, presentation or poster |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--aspect ASPECT` | height/width ratio of the figure |
| `--mono` | grayscale version: black ink and line patterns |
| `--bc {bc1,bc2,bc3}` | bc1 vacuum/vacuum (neutral slabs), bc2 metal/metal (capacitor), bc3 vacuum/metal (electrode, the only one besides bc2 that allows a net charge) (default: `bc1`) |
| `--charge` | net charges in e, comma-separated: -0.2,0,0.2 (default: `0`) |
| `--field FIELD` | applied field in Ry/a.u. (bc2 only) |
| `--esm-w WIDTH_ESM` | shift of the ESM boundary in a.u. |
| `--nfit NFIT` | fitting points for the potential at the boundary (default: `4`) |

**Physics:** [`olla-dft teoria esm`](THEORY.md)

### `echem`

computational hydrogen electrode: HER, OER, limiting potential and overpotential

**Usage:** `olla-dft echem [-h] [--her E_ads] [--oer OH=..,O=..,OOH=..] [--corrections X=eV] [-U POTENTIAL] [--ph PH] [-T TEMPERATURE] [-o OUTDIR] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono]`

**Options:**

| Option | Description |
|---|---|
| `--her E_ads` | H adsorption energy in eV (HER reaction) |
| `--oer OH=..,O=..,OOH=..` | adsorption energies of the three OER intermediates, in eV and referenced to water |
| `--corrections X=eV` | ZPE−TΔS thermal corrections per intermediate; without them the standard literature values are used |
| `-U, --potential POTENTIAL` | applied potential in V vs SHE (at pH 0 it is the same as vs RHE; the pH converts it) |
| `--ph PH` | pH |
| `-T, --temperature TEMPERATURE` | temperature in K (default: 298.15) |
| `-o, --outdir` | output folder (default: `echem`) |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria echem`](THEORY.md)

### `neb`

reaction path and activation barrier (neb.x)

**Usage:** `olla-dft neb [-h] [-o OUTDIR] [--images IMAGES] [--no-ci] [--path-thr PATH_THR] [--nstep NSTEP] [--fix FIX] [--prefix PREFIX] [--collect] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [--ecutwfc ECUTWFC] [--ecutrho ECUTRHO] [--kspacing KSPACING] [--metal] [--nspin {1,2}] [--mag MAG] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] file [final]`

**Arguments:**

- `file` — initial structure (reactant)
- `final` — final structure (product)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `neb`) |
| `--images IMAGES` | number of images in the chain (default: `7`) |
| `--no-ci` | no climbing image; the barrier will come out UNDERESTIMATED |
| `--path-thr PATH_THR` | path force threshold in eV/Å (default: `0.05`) |
| `--nstep NSTEP` | maximum path optimization steps in neb.x (default: 50) |
| `--fix` | indices of atoms to freeze (0-based) |
| `--prefix` | calculation prefix (detected automatically) |
| `--collect` | read the results of a calculation that has already run instead of preparing the inputs |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |
| `--ecutwfc ECUTWFC` | wavefunction cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--ecutrho ECUTRHO` | density cutoff in Ry (if not given, the one recommended by the UPF files) |
| `--kspacing KSPACING` | k-mesh spacing in Å^-1 |
| `--metal` | metallic system: smeared occupations instead of fixed ones |
| `--nspin {1,2}` | 2 turns on spin polarization (default: 1) |
| `--mag` | starting magnetization: a number (0.5) or per element (Fe=0.7,O=0). Implies --nspin 2 |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria neb`](THEORY.md)

### `amorphous`

amorphous solid by melt and quench with a machine-learned potential

**Usage:** `olla-dft amorphous [-h] [-n UNITS] -d G_CM3 [--melt K] [--final K] [--melt-steps MELT_STEPS] [--quench-steps QUENCH_STEPS] [--anneal-steps ANNEAL_STEPS] [--dt FS] [--model MODEL] [--min-dist F] [--seed SEED] [--pack-only] [-o OUTDIR] formula`

**Arguments:**

- `formula` — formula unit, e.g. SiO2

**Options:**

| Option | Description |
|---|---|
| `-n, --units UNITS` | formula units in the cell (default: 8) |
| `-d, --density G_CM3` | target density in g/cm³ |
| `--melt K` | melt temperature (default: 3000 K) |
| `--final K` | final temperature (default: 300 K) |
| `--melt-steps MELT_STEPS` | dynamics steps in the melt phase (default: 500) |
| `--quench-steps QUENCH_STEPS` | quench steps: these set the rate. The default (1000) is an exploratory quench at ~3e15 K/s and the report warns about it; 27000 brings it down to 1e14 K/s |
| `--anneal-steps ANNEAL_STEPS` | dynamics steps of the anneal at the final temperature (default: 200) |
| `--dt FS` | step: integration step in fs (amorphous) or temperature step in K (qha) (default: `1.0`) |
| `--model` | interatomic potential (default: `mace`) |
| `--min-dist F` | factor on the sum of covalent radii when packing (default: 0.75) |
| `--seed SEED` | seed; change it to generate another realization |
| `--pack-only` | pack only, no dynamics |
| `-o, --outdir` | output folder (default: `amorfo`) |

**Physics:** [`olla-dft teoria amorphous`](THEORY.md)

## Automation and quality

### `doctor`

diagnose a calculation: convergence, forces and why it does not converge

**Usage:** `olla-dft doctor [-h] [--system] [--project PROJECT] [--json] [--prefix PREFIX] [-o OUTDIR] [--no-plot] [--dpi DPI] [--format FORMAT] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] [path]`

**Arguments:**

- `path` — calculation folder or output file

**Options:**

| Option | Description |
|---|---|
| `--system` | check installation, resources, QE and pseudopotentials |
| `--project` | also check this project's quality gate |
| `--json` | print the diagnosis as JSON |
| `--prefix` | calculation prefix (detected automatically) |
| `-o, --outdir` | output folder (default: `.`) |
| `--no-plot` | only export the data, without generating the figure |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria doctor`](THEORY.md)

### `audit`

verify that a set of calculations is comparable before subtracting energies

**Usage:** `olla-dft audit [-h] [--index] [--db DB] paths [paths ...]`

**Arguments:**

- `paths` — folders or XML files of the calculations

**Options:**

| Option | Description |
|---|---|
| `--index` | also register them in the database |
| `--db` | SQLite file of the calculation index (default: olla-dft.db) |

**Physics:** [`olla-dft teoria audit`](THEORY.md)

### `crosscheck`

cross-check the same quantity through independent routes

**Usage:** `olla-dft crosscheck [-h] [-f FILE] [--gap-bandas GAP_BANDAS] [--gap-tauc GAP_TAUC] [project]`

**Arguments:**

- `project` — project folder

**Options:**

| Option | Description |
|---|---|
| `-f, --file` | structure (for masses and volume) |
| `--gap-bandas GAP_BANDAS` | band-structure gap in eV, to cross-check against the Tauc one |
| `--gap-tauc GAP_TAUC` | Tauc-extrapolation gap in eV, to cross-check against the band-structure one |

**Physics:** [`olla-dft teoria crosscheck`](THEORY.md)

### `cost`

what Olla-DFT knows about your machine's speed

**Usage:** `olla-dft cost [-h] [--db DB]`

**Options:**

| Option | Description |
|---|---|
| `--db` | calculation database (default: `olla-dft.db`) |

**Physics:** [`olla-dft teoria cost`](THEORY.md)

### `db`

local index of calculations

**Usage:** `olla-dft db [-h] [--db DB] [-q QUERY] [--export EXPORT] [--formula FORMULA] [--calculation CALCULATION] [--gap-min GAP_MIN] [--gap-max GAP_MAX] [--limit LIMIT] [paths ...]`

**Arguments:**

- `paths` — folders to register

**Options:**

| Option | Description |
|---|---|
| `--db` | SQLite file of the calculation index (default: olla-dft.db) |
| `-q, --query` | SQL query (SELECT only) |
| `--export` | export everything to a JSON file |
| `--formula` | filter by formula, e.g. Si |
| `--calculation` | filter by type: scf, relax, nscf... |
| `--gap-min GAP_MIN` | minimum gap in eV |
| `--gap-max GAP_MAX` | maximum gap in eV |
| `--limit LIMIT` | maximum number of filtered rows (default: `100`) |

**Physics:** [`olla-dft teoria db`](THEORY.md)

### `hull`

formation energies and convex hull

**Usage:** `olla-dft hull [-h] [-o OUTDIR] [--elements ELEMENTS] [--threshold THRESHOLD] [--force] [--dpi DPI] [--format FORMAT] [--no-plot] [-t TEMPLATE] [--font {sans,serif,latex}] [--usetex] [--palette PALETTE] [--background BACKGROUND] [--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}] [--width WIDTH] [--mono] paths [paths ...]`

**Arguments:**

- `paths` — folders or XML files of the calculations

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--elements` | element order, e.g. Zn,Al |
| `--threshold THRESHOLD` | metastability threshold in eV/atom (default: `0.025`) |
| `--force` | build the hull even if the audit fails |
| `--dpi DPI` | resolution of bitmap figures in dots per inch (default: 600) |
| `--format` | comma-separated figure formats: pdf,png,svg,eps,tif (default: pdf,png) |
| `--no-plot` | only export the data, without generating the figure |
| `-t, --template` | visual template of the figure (see the list with 'olla-dft templates list') |
| `--font {sans,serif,latex}` | font family: sans, serif or latex (Computer Modern) |
| `--usetex` | typeset the figure text with real LaTeX (needs a LaTeX installation) |
| `--palette` | color palette: grayscale, okabe-ito, okabe-ito-dark or comma-separated hex colors |
| `--background` | figure background color, e.g. '#FFFFFF' or 'none' |
| `--journal {acs,aps,elsevier,generic,iop,nature,rsc,wiley}` | figure width according to the journal (default: generic) |
| `--width` | figure width: single, onehalf or double, or a number in millimeters |
| `--mono` | grayscale version: black ink and line patterns |

**Physics:** [`olla-dft teoria hull`](THEORY.md)

### `mlip`

machine-learned potential: pre-relax and screen before spending DFT

**Usage:** `olla-dft mlip [-h] [-o OUTPUT] [--model {mace,chgnet,m3gnet}] [--size SIZE] [--device DEVICE] [--fmax FMAX] [--steps STEPS] [--fixed-cell] [--span SPAN] [--npoints NPOINTS] [--supercell SUPERCELL] {relax,scan,phonons} file`

**Arguments:**

- `action` {relax,scan,phonons} — action to perform (see the list above)
- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output structure (relax) |
| `--model {mace,chgnet,m3gnet}` | machine-learned potential: mace, chgnet or m3gnet (default: mace) |
| `--size` | MACE model size (small/medium/large) (default: `small`) |
| `--device` | device the potential runs on: cpu or cuda (default: cpu) |
| `--fmax FMAX` | target force in eV/Å (default: `0.01`) |
| `--steps STEPS` | maximum relaxation steps (default: 300) |
| `--fixed-cell` | do not relax the cell, only the positions |
| `--span SPAN` | range of the volume sweep (scan) (default: `0.1`) |
| `--npoints NPOINTS` | points in the volume sweep (default: 15) |
| `--supercell` | supercell for the screening, e.g. 2x2x2 |

**Physics:** [`olla-dft teoria mlip`](THEORY.md)

### `suggest`

suggest parameters from your previous calculations

**Usage:** `olla-dft suggest [-h] [--db DB] file`

**Arguments:**

- `file` — input structure (CIF, POSCAR, pw.x input...)

**Options:**

| Option | Description |
|---|---|
| `--db` | SQLite file of the calculation index (default: olla-dft.db) |

**Physics:** [`olla-dft teoria suggest`](THEORY.md)

### `datasheet`

material datasheet and methods paragraph

**Usage:** `olla-dft datasheet [-h] [-o OUTDIR] [--name NAME] [--methods] [project]`

**Arguments:**

- `project` — project folder (default: .)

**Options:**

| Option | Description |
|---|---|
| `-o, --outdir` | output folder (default: `.`) |
| `--name` | base name of the files |
| `--methods` | only the methods paragraph and the citations |

### `report`

local log of failures and confusions

**Usage:** `olla-dft report [-h] [--show SHOW] [--close CLOSE] [--note NOTE] [--stats] [--export EXPORT] [--only-open] [--attach ATTACH] [description ...]`

**Arguments:**

- `description` — what happened (if omitted, lists the incidents)

**Options:**

| Option | Description |
|---|---|
| `--show` | view an incident by its id |
| `--close` | mark an incident as resolved |
| `--note` | note when closing |
| `--stats` | which subcommands fail most |
| `--export` | package everything into a JSON file |
| `--only-open` | list only the open incidents |
| `--attach` | attach a file (it is copied to the local log) |

### `compare`

compare runs without subtracting incompatible energies

**Usage:** `olla-dft compare [-h] [--reference REFERENCE] [-o OUTPUT] paths [paths ...]`

**Arguments:**

- `paths` — folders or XML files of the runs

**Options:**

| Option | Description |
|---|---|
| `--reference REFERENCE` | index of the reference run (default: 0) |
| `-o, --output` | save the comparison as JSON |

### `tune`

recommend the next point of a convergence test

**Usage:** `olla-dft tune [-h] [--threshold THRESHOLD] [-o OUTPUT] file`

**Arguments:**

- `file` — CONVERGENCIA.dat

**Options:**

| Option | Description |
|---|---|
| `--threshold THRESHOLD` | threshold in meV/atom (default: 1) |
| `-o, --output` | save the recommendation as JSON |

**Physics:** [`olla-dft teoria tune`](THEORY.md)

### `results`

ingest, query and export normalized project results

**Usage:** `olla-dft results [-h] [--project PROJECT] [--db DB] [--tag TAG] [--formula FORMULA] [--calculation CALCULATION] [--status {invalid,not_converged,parsed_no_energy,parsed,converged}] [--review-status {unreviewed,accepted,rejected}] [--note NOTE] [--limit LIMIT] [--json] [-o OUTPUT] {ingest,list,show,review,export,explore} [target] [extra_paths ...]`

**Arguments:**

- `action` {ingest,list,show,review,export,explore} — action to perform (see the list above)
- `target` — input path for ingest, or id for show
- `extra_paths` — more folders/XML files for ingest

**Options:**

| Option | Description |
|---|---|
| `--project` | project folder (default: .) |
| `--db` | alternative SQLite; by default .qekit/results.sqlite3 |
| `--tag` | provenance tag for ingest |
| `--formula` | filter by formula |
| `--calculation` | filter by calculation type |
| `--status {invalid,not_converged,parsed_no_energy,parsed,converged}` | filter by status: invalid, not_converged, parsed_no_energy, parsed or converged |
| `--review-status {unreviewed,accepted,rejected}` | in review, status of the human review |
| `--note` | in review, note that goes with the decision |
| `--limit LIMIT` | maximum records: list=100, explore=10000 |
| `--json` | in list, print JSON |
| `-o, --output` | output file: export=JSON, explore=interactive HTML |

### `campaign`

create reproducible matrices of parameterized tasks

**Usage:** `olla-dft campaign [-h] [--project PROJECT] [--command CAMPAIGN_COMMAND] [--axis AXIS] [--goal GOAL] [--convergence-file CONVERGENCE_FILE] [--adaptive] [--threshold THRESHOLD] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [-o OUTPUT] {create,list,status,export,run,extend} [target]`

**Arguments:**

- `action` {create,list,status,export,run,extend} — action to perform (see the list above)
- `target` — campaign name or id

**Options:**

| Option | Description |
|---|---|
| `--project` | project folder (default: .) |
| `--command` | olla-dft command template; fields: {eje}, {index}, {id}, {structure} |
| `--axis` | axis name=v1,v2; can be repeated (default: `[]`) |
| `--goal` | scientific goal of the campaign |
| `--convergence-file` | CONVERGENCIA.dat to take a recommendation from |
| `--adaptive` | add the next recommended value to the convergence axis |
| `--threshold THRESHOLD` | convergence threshold when extending (meV/atom) |
| `--execute` | in run, execute the selected points |
| `--force` | in run, ignore the task cache |
| `--parallel PARALLEL` | in run, simultaneous independent points (default: 1) |
| `--retries RETRIES` | retries per failed point (default: 0) |
| `--timeout TIMEOUT` | maximum time per attempt, in seconds |
| `--cancel-file` | custom cooperative-cancellation marker |
| `-o, --output` | JSON file for export |

### `pseudos`

compare the available pseudopotentials and choose on merit, not alphabetically

**Usage:** `olla-dft pseudos [-h] [--element ELEMENT] [--task TASK] [--functional FUNCTIONAL] [--cheap] [--pseudo-dir PSEUDO_DIR] [--pseudo EL=UPF] [file]`

**Arguments:**

- `file` — structure

**Options:**

| Option | Description |
|---|---|
| `--element` | comma-separated elements |
| `--task` | what it is for: general, optics, soc, xanes, hubbard, fonones. Each task discards the ones that are not suitable (default: `general`) |
| `--functional` | require a specific functional (PBE, PZ, PBEsol...) |
| `--cheap` | prefer ultrasoft/PAW, which need fewer plane waves |
| `--pseudo-dir` | folder with the UPF pseudopotentials (if not given, the one from 'olla-dft config') |
| `--pseudo EL=UPF` | force a specific pseudopotential, e.g. Fe=Fe.rel-pbe.UPF. Can be repeated. Without it, Olla-DFT chooses with 'olla-dft pseudos' |

**Physics:** [`olla-dft teoria pseudos`](THEORY.md)

## Project

### `project`

manage a reproducible project: sources, workflow, quality and dashboard

**Usage:** `olla-dft project [-h] [--project PROJECT] [--name NAME] [--command TASK_COMMANDS] [--execute] [--force] [--parallel PARALLEL] [--retries RETRIES] [--timeout TIMEOUT] [--cancel-file CANCEL_FILE] [--reason REASON] [--selftest] [--advanced] [-o OUTPUT] [--pdf] [--theme {auto,light,dark}] [--language {es,en}] [--both] [--verify-environment] [--other OTHER] [--json] {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} [target]`

**Arguments:**

- `action` {init,add,plan,show,status,validate,run,dashboard,report,export,ingest,environment,diff,cancel,resume} — action on the project
- `target` — directory, file, goal, profile or task depending on the action

**Options:**

| Option | Description |
|---|---|
| `--project` | project to work from (default: .) |
| `--name` | name when initializing |
| `--command` | custom olla-dft task; can be repeated with plan |
| `--execute` | execute run/submit; by default it only simulates or writes |
| `--force` | in run, ignore the cache and prepare all tasks again |
| `--parallel PARALLEL` | in run, simultaneous independent tasks (default: 1) |
| `--retries RETRIES` | retries per failed task (default: 0) |
| `--timeout TIMEOUT` | maximum time per attempt, in seconds |
| `--cancel-file` | custom cooperative-cancellation marker |
| `--reason` | in cancel, optional reason |
| `--selftest` | in validate, run the quick validation against physical references |
| `--advanced` | in validate, check structure, commands, units and collisions |
| `-o, --output` | output for dashboard, report or export |
| `--pdf` | in report, generate a self-contained PDF report |
| `--theme {auto,light,dark}` | dashboard theme (default: auto) |
| `--language {es,en}` | dashboard language (default: es) |
| `--both` | generate Spanish and English dashboards in separate files |
| `--verify-environment` | in environment, check the saved lock |
| `--other` | in diff, snapshot or project to compare against |
| `--json` | in diff, print JSON |

### `resilient`

recoverable QE jobs for interrupted servers

**Usage:** `olla-dft resilient [-h] [--state STATE] [--pw-cmd PW_CMD] [--runtime-id RUNTIME_ID] [--checkpoint-seconds CHECKPOINT_SECONDS] [--grace-seconds GRACE_SECONDS] [--max-failures MAX_FAILURES] [--threads THREADS] [--keep KEEP] [--max-segments MAX_SEGMENTS] [--resume] [--user USER] [-o OUTPUT] {init,run,status,pause,service} target`

**Arguments:**

- `action` {init,run,status,pause,service} — action to perform (see the list above)
- `target` — input file for init; durable state directory otherwise

**Options:**

| Option | Description |
|---|---|
| `--state` | new job directory on a retained persistent disk |
| `--pw-cmd` | QE/launcher command, with fixed MPI flags (default: `pw.x`) |
| `--runtime-id` | immutable runtime image identifier |
| `--checkpoint-seconds CHECKPOINT_SECONDS` |  (default: `900`) |
| `--grace-seconds GRACE_SECONDS` |  (default: `300`) |
| `--max-failures MAX_FAILURES` |  (default: `3`) |
| `--threads THREADS` |  (default: `1`) |
| `--keep KEEP` | intact checkpoint generations to retain (minimum 2) (default: `2`) |
| `--max-segments MAX_SEGMENTS` | stop after this many saved segments; 0 means unlimited |
| `--resume` | clear an explicit pause before running |
| `--user` | unprivileged user for the generated systemd service |
| `-o, --output` | generated service file; installation is separate |

## Appearance and configuration

### `templates`

list, view or export templates

**Usage:** `olla-dft templates [-h] [-o OUTPUT] [{list,show,export}] [name]`

**Arguments:**

- `action` {list,show,export} — list (default), show or export
- `name` — template name

**Options:**

| Option | Description |
|---|---|
| `-o, --output` | output JSON file (export) |

### `config`

view or change the configuration

**Usage:** `olla-dft config [-h] [{show,set}] [key] [value]`

**Arguments:**

- `action` {show,set} — action to perform (see the list above)
- `key` — configuration key, e.g. pseudo_dir, nproc or language
- `value` — value to assign to the key

---

*Olla-DFT 1.4.0*
