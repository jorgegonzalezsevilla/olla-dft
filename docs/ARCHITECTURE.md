# Architecture

How the Olla-DFT code base is organised, the rules every module follows, and
what it takes to add a command.

The installed command is `olla-dft`; the Python package keeps its original
name, `qekit`, so that existing scripts, projects and configuration files keep
working.

## Module tree

```
qekit/
├── __init__.py            # version, product name, command name, author
├── __main__.py            # python -m qekit
├── cli.py                 # interactive menu + flat subcommands (argparse), dispatch, exit codes
├── config.py              # persistent configuration (config.ini in the per-OS folder), legacy migration
├── core/                  # infrastructure shared by every module
│   ├── atomconf.py        # atomic electronic configurations for generating pseudopotentials
│   ├── compat.py          # compatibility across versions of the dependencies (numpy, ASE)
│   ├── consola.py         # console output that never crashes: UTF-8 first, ASCII transliteration otherwise
│   ├── errors.py          # ErrorDeUso (usage error, exit 2) vs. program failure (exit 1)
│   ├── i18n.py            # interface language (es/en): --language, OLLA_DFT_LANG, config; JSON tables
│   ├── kpoints.py         # Γ-centred k-meshes from a spacing and high-symmetry k-paths (seekpath)
│   ├── layers.py          # layer detection in layered materials by periodic connectivity
│   ├── plataforma.py      # per-OS folders, binary names (.x/.exe), MPI launcher, portable run scripts
│   ├── provenance.py      # version, command line and parameters written into every output
│   ├── pseudo.py          # UPF handling: search, cutoffs and valence from the header
│   ├── qeout.py           # reading Quantum ESPRESSO output (pw.x XML, high-symmetry labels)
│   ├── runner.py          # batch execution of pw.x (--run): resumable, with translated QE failures
│   ├── structure.py       # structure reading and symmetry (ASE + spglib), primitive/conventional cells
│   ├── style.py           # figure style for publication: physical sizes, typography, palette
│   ├── themes.py          # visual templates and verified palettes
│   └── wfc.py             # reader of the binary wavefunction files of QE
├── data/
│   ├── atomic_scattering_params.json   # X-ray form factors (from pymatgen, MIT) with attribution notice
│   ├── i18n/              # cli_es/en, menu_es/en, onboarding_es/en, dashboard_es/en, recipes_en, wizard_en, docs_en
│   └── theory/            # electronica/mecanica/espectros .es.md/.en.md: the physics behind every command
└── modules/               # one file per task; each exposes the physics through its docstring
    ├── adsorb.py          # adsorption sites on a slab: enumerate, build and compare
    ├── align.py           # band alignment: where the VBM of one material sits relative to the other
    ├── amorphous.py       # amorphous solids by melt and quench (MACE)
    ├── audit.py           # consistency audit across calculations, and local database
    ├── ballistic.py       # ballistic transport: Landauer conductance with pwcond.x
    ├── bands.py           # band-structure post-processing, gap analysis, export and figure
    ├── berry.py           # electric polarisation by the Berry phase, and what follows from it (Z*)
    ├── builder.py         # structure builders: surfaces by Miller indices and point defects
    ├── campaign.py        # reproducible campaigns of parametrised calculations
    ├── charges.py         # atomic charges (Bader, Löwdin) and density difference
    ├── combined.py        # combined bands + DOS figure
    ├── compare.py         # safe comparison of Quantum ESPRESSO runs
    ├── converge.py        # convergence tests: ecutwfc, ecutrho and k mesh
    ├── corehole.py        # core-hole pseudopotentials with ld1.x, for XPS and XANES
    ├── cost.py            # how long this is going to take, before launching it
    ├── crosscheck.py      # cross-validation: the same quantity by independent routes
    ├── dashboard.py       # self-contained HTML dashboard of a project
    ├── datasheet.py       # material datasheet and methods paragraph
    ├── defects.py         # formation energy of charged defects (Madelung correction)
    ├── derived.py         # quantities derived from elastic constants and phonons: Debye, sound velocities, Slack
    ├── diagnose.py        # diagnosis of a pw.x calculation: is it usable, and if not, why
    ├── docs.py            # browsable HTML reference generated from the code itself
    ├── dos.py             # DOS and projected DOS post-processing
    ├── dynamics.py        # molecular-dynamics trajectories: g(r), MSD/diffusion, VDOS
    ├── echem.py           # computational hydrogen electrode: HER, OER and the Pourbaix diagram
    ├── effmass.py         # effective mass by parabolic fit of the bands
    ├── elastic.py         # elastic constants by the stress-strain method
    ├── elph.py            # electron-phonon coupling: lambda, alpha²F, Tc and a real tau
    ├── environment.py     # lightweight environment lock for local reproducibility
    ├── eos.py             # E-V equation of state: equilibrium volume and bulk modulus
    ├── esm.py             # charged surfaces: the effective screening medium (ESM)
    ├── exfoliate.py       # exfoliation energy of a layered material
    ├── feedback.py        # local incident log: failures, confusions and slipped errors; no telemetry
    ├── fields.py          # charge density, electrostatic potential and work function (pp.x)
    ├── health.py          # diagnosis of the installation and the execution environment
    ├── hubbard.py         # Hubbard U by linear response, with hp.x
    ├── inputgen.py        # input generator for pw.x and post-processing codes; run.sh and run.py
    ├── interface.py       # heterostructures: stacking two materials with the least strain
    ├── interop.py         # interchange format with the rest of the suite
    ├── kappa.py           # lattice thermal conductivity: phonon-phonon scattering (phono3py)
    ├── mlip.py            # machine-learned interatomic potentials: pre-relaxation and screening
    ├── neb.py             # reaction paths and activation barriers with neb.x
    ├── onboarding.py      # guided start for someone who has never used a scientific CLI
    ├── optics.py          # optical properties with epsilon.x: ε(ω), n, k, absorption, Tauc
    ├── phonons.py         # phonons by DFPT: dispersion, DOS, thermodynamics and IR
    ├── project.py         # Project Hub: projects, resumable workflows and provenance
    ├── pseudos.py         # choosing a pseudopotential with criteria, not alphabetically
    ├── qha.py             # quasi-harmonic approximation: thermal expansion and a(T)
    ├── quality.py         # scientific quality gate for a project
    ├── recipes.py         # recipes: complete sessions from structure to result
    ├── recommend.py       # recommendations from YOUR own calculation history
    ├── report.py          # compact PDF report of a project
    ├── results.py         # local engine of normalised, traceable results
    ├── selftest.py        # checks against known physics, not against itself
    ├── strain.py          # strain sweep: properties as a function of applied deformation
    ├── surfen.py          # surface energy: cut a crystal and see what it costs
    ├── sweep.py           # shared infrastructure of every sweep (prepare / run / collect)
    ├── tddft.py           # optical absorption with TDDFPT: the part epsilon.x does not see
    ├── theory.py          # the physics behind each command, readable from the terminal (olla-dft teoria)
    ├── thermo.py          # formation energies, convex hull and phase stability
    ├── thermochem.py      # thermochemistry: from a DFT energy to a comparable free energy
    ├── topology.py        # topological invariants of a Wannier Hamiltonian: Chern, Wilson loops
    ├── tphonons.py        # phonons at electronic temperature: does the imaginary mode stabilise?
    ├── transport.py       # electronic transport in the constant relaxation-time approximation
    ├── tuning.py          # adaptive recommendation from a convergence series
    ├── uncertainty.py     # small utilities to propagate experimental or numerical uncertainties
    ├── unfold.py          # band unfolding: recovering the dispersion of a supercell
    ├── validation.py      # structural and integrity validations before spending CPU
    ├── wannier.py         # Wannier functions: bringing the band structure down to a small model
    ├── wizard.py          # guided wizard: from what you want to KNOW to the files to run
    ├── xanes.py           # XANES / NEXAFS with xspectra.x
    ├── xps.py             # core-level shifts (XPS) in the initial-state approximation
    └── xrd.py             # simulated powder diffraction from the crystal structure
```

`tools/build_docs.py` generates `docs/COMMANDS.md` and `docs/THEORY.md`
from the argparse tree and from `qekit/data/theory/`; those two files are
never edited by hand (`--all` also writes the Spanish pair).

## Design rules

**Flat subcommands, grouped help.** Every command is a flat `olla-dft
<command>` (easy to script), but `olla-dft --help` and the interactive menu
present them by task: the `COMMAND_GROUPS` table in `cli.py` lists each
command exactly once and a test protects it against omissions. Three Spanish
command names have English aliases (`sistema`/`system`, `recetas`/`recipes`,
`teoria`/`theory`).

**prepare / --run / --collect.** Every sweep (convergence, EOS, elastic,
strain, surface energy, adsorption, defects, phonons, ESM, Berry, kappa…)
follows the same cycle, shared through `modules/sweep.py`: `prepare` writes
one folder per calculation plus `run.sh` and `run.py`; `--run` executes them
here through `core/runner.py` (resumable, `-j N` parallel points, `--max-time`
budget, `--timeout`, `--redo`); without `--run` the command explains how to
launch them; `--collect` reads finished calculations and writes the report
without rewriting the inputs, so a run the user edited by hand is described
as it actually ran. `--estimate` predicts the cost from the local history and
exits. Every such command shares the same option groups (`ejecución`,
`parámetros DFT`, `figura`), which argparse lists after the command's own
options.

**Exit codes.** `core/errors.py` separates errors by who has to fix them.
`ErrorDeUso` (a subclass of `ValueError`) means the command or its data do not
fit and the message already says what to do: the program did the right thing,
exits with **2** (like argparse) and is logged as type "usage" without a
trace. Any other exception is a program failure: it exits with **1** and is
archived in full — command, trace, versions — by `modules/feedback.py`.
Commands that ran correctly and *found* a problem (`doctor`, `crosscheck`,
`audit`, `selftest`, `project quality`, a `campaign` with a failed task) also
return **1**, so that a script can stop
on them; `tests/barrido_cli.sh` declares the expected code on every line.
`Ctrl-C` returns 130; a broken pipe (`| head`) closes silently with 0.

**Provenance everywhere.** `core/provenance.py` writes the Olla-DFT version,
the UTC timestamp, the exact command line and the parameters as `#` comments
at the top of every `.dat`, and as metadata in every PDF/PNG (readable with
`pdfinfo`, `exiftool` or `olla-dft info --figura`). Months later a number in a
thesis can be traced to the calculation that produced it.

**Interface language layer.** `core/i18n.py` decides the language in this
order: the global `--language en` flag (accepted anywhere on the command
line), the `OLLA_DFT_LANG` variable, the `language` key of the configuration,
Spanish. What is translated is the *interface*: the help of every command and
flag (`data/i18n/cli_en.json`), the interactive menu (`menu_*.json`), the
guided start (`onboarding_*.json`), the recipes and the wizard (`recipes_en`,
`wizard_en`, applied to the Spanish catalogue with `translate_data`), the
dashboard (`dashboard_*.json`), the HTML reference (`docs_en.json`) and the
theory (`data/theory/*.en.md`). What stays in Spanish is the scientific
report each command prints and writes: the module docstrings, the `.dat`
headers and the analysis text are not translated at run time.

**Console output that never dies.** `core/consola.py` switches the console to
UTF-8 when possible and otherwise transliterates every non-ASCII symbol to an
ASCII equivalent that reads the same (`Å → A`, `α → alpha`, `→ → ->`); `--ascii`
forces it. A test walks every real report and fails if a symbol is missing
from the table.

**No telemetry.** `modules/feedback.py` keeps the incident log in the local
configuration folder; nothing leaves the machine unless the user exports it
with `olla-dft report --export`.

## Adding a command

1. Write the physics in a new module `qekit/modules/<name>.py`. Its docstring
   is not a courtesy: `olla-dft docs` and the HTML reference quote it as "the
   physics behind" the command, so it must explain what the command answers
   and when the result is not valid. A sweep exposes `prepare`, `collect`,
   `report`, `export` and `plot` and builds its jobs with
   `sweep.prepare_common` / `sweep.write_scf_job`.
2. In `qekit/cli.py`: add the parser (inline in `build_parser`, using
   `_calc_opts` for a sweep so it inherits `--run/--collect/-j/--pseudo-dir…`),
   write `_cmd_<name>(args) -> int` (raise `ErrorDeUso` for bad input; use
   `_run_or_explain` for the run step), and register it in `_DISPATCH` and in
   `COMMAND_GROUPS`. Every option needs a Spanish `help=` string.
3. In `qekit/modules/docs.py`: add the command to a group in `GRUPOS` and, if
   its physics lives in a module with a different name, to `MODULO_DE`
   (`tests/test_docs.py` fails on orphans).
4. In `qekit/data/i18n/docs_en.json`: add the one-line summary under
   `command_summaries`. In `qekit/data/i18n/cli_en.json`: add the English
   translation of every new help string under `help` (keyed by the Spanish
   text). The Spanish catalogues (`recipes`, `wizard`) get their `{es: en}`
   entries in the matching `_en.json` file.
5. In `qekit/data/theory/<area>.es.md` and `<area>.en.md`: add a section
   headed ``### `olla-dft <name>` — title`` with the mandatory parts (*Qué
   responde / What it answers*, *Fundamento para no expertos*, *Cómo lo
   calcula Olla-DFT*, *Límites y trampas*, *Referencias*); `tests/test_teoria.py`
   requires one for every scientific command and checks es/en parity.
6. Add a test in `tests/` (with a real QE output in `tests/datos/` if the
   command reads one, and a frozen reference in `tests/referencias.py` if it
   produces a number validated against experiment).
7. Run `python tools/build_docs.py` to regenerate `docs/COMMANDS.md`
   and `docs/THEORY.md`, then
   `python -m pytest -q` and `python -m pyflakes qekit tests`.

## Tests

- `tests/` — 977 pytest tests that run without Quantum ESPRESSO in under a
  minute (`python -m pytest -q`). They cover the functions, the argparse tree
  (every command in a group, every command with a theory section, every
  README and recipe command valid), the three platforms (simulated
  `sys.platform`, forced cp1252 output), the i18n tables (the English table
  covers the whole catalogue) and the licence.
- `tests/datos/` — real Quantum ESPRESSO outputs the tests read: Si bands,
  DOS and effective mass, Si phonons, `epsilon.x` spectra, XANES, TDDFPT of
  ethylene, `pwcond.x` of an Al wire, electron-phonon of Al, an MD
  trajectory, an unfolding file and a core-hole UPF.
- `tests/referencias.py` — frozen reference values, each validated once
  against experiment, a PDF card or another implementation; from then on they
  are regression detectors, not documentation.
- `tests/barrido_cli.sh` — command-level regression sweep over already
  computed QE outputs (`OLLA_DFT_REG=/path bash tests/barrido_cli.sh`), with
  numpy `RuntimeWarning`s turned into errors so that a silent NaN does not
  pass, and the expected exit code (0, 1 or 2) declared on every line.
- Tests that need QE binaries carry the `qe` marker.
