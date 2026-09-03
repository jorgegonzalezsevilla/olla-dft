<h1 align="center">Olla-DFT</h1>

<p align="center">
<b>A command-line toolkit for Quantum ESPRESSO: from a CIF to publication-ready
band structures, phonons, elastic constants, optics and more, with the physics
behind every number written down.</b>
</p>

<p align="center">
<a href="https://github.com/jorgegonzalezsevilla/olla-dft/actions/workflows/ci.yml"><img src="https://github.com/jorgegonzalezsevilla/olla-dft/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
<a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue.svg" alt="License: GPL v3"></a>
<a href="https://doi.org/10.5281/zenodo.22263121"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.22263121.svg" alt="DOI"></a>
<img src="https://img.shields.io/badge/python-3.9%20–%203.13-blue" alt="Python 3.9 to 3.13">
<img src="https://img.shields.io/badge/tests-977-brightgreen" alt="977 tests">
<img src="https://img.shields.io/badge/commands-78-orange" alt="78 commands">
</p>

<p align="center">
<a href="https://github.com/jorgegonzalezsevilla/olla-dft-esp">Versión en español</a> ·
<a href="#install">Install</a> ·
<a href="#five-minute-tour">Five-minute tour</a> ·
<a href="docs/THEORY.md">Theory</a> ·
<a href="docs/COMMANDS.md">Commands</a> ·
<a href="examples/">Examples</a>
</p>

<p align="center">
<img src="examples/demo_Si/Si_bandas_dos.png" alt="Silicon band structure and DOS produced by olla-dft plot" width="720">
<br>
<sub>Silicon bands and DOS from a real Quantum ESPRESSO run, one command: <code>olla-dft plot .</code></sub>
</p>

Olla-DFT plays the role that VASPKIT plays for VASP, but for
[Quantum ESPRESSO](https://www.quantum-espresso.org): it reads your structure
(CIF, POSCAR or a `pw.x` input), detects the symmetry, finds your
pseudopotentials, proposes cutoffs and k-meshes, builds the high-symmetry path
and writes every input file ready to run. After the calculation it closes the
loop: band gap, bands, DOS/PDOS, equation of state, elastic constants,
phonons, optics, work function, charges, defects, surfaces, transport,
Wannier functions and a long list of derived properties, each one exported
as a table with its provenance and drawn as a journal-quality figure.

## About this project

Olla-DFT is a **one-person hobby project**. It is written, tested and
maintained by Jorge Enrique González Sevilla in his free time, out of a
personal need: a single tool that takes a Quantum ESPRESSO workflow from the
structure file to the figure without a pile of ad-hoc scripts, and that
explains the physics it applies instead of hiding it. It is not affiliated
with, funded by, or endorsed by any university, research institute or
company, and it is not an official part of Quantum ESPRESSO.

Because it is one person's spare-time work, a few things follow:

- **Releases come when they come.** There is no roadmap with dates.
- **The code is not open to outside commits.** Pull requests are not
  accepted; the author is and will remain the only one touching the code.
  What is very welcome is your feedback: bugs, wrong numbers and feature
  ideas as [issues](https://github.com/jorgegonzalezsevilla/olla-dft/issues).
- **Everything is free software.** GPL-3.0, no telemetry, nothing is ever
  sent anywhere. Clone it, read it, run it, fork it under the licence.

The name keeps the coffee wink of Quantum ESPRESSO and gives it a local
identity: *olla* is the clay pot in which *café de olla* is brewed. It is made
with love in Guadalajara, Jalisco, Mexico.

## Contents

- [About this project](#about-this-project)
- [Install](#install)
- [Five-minute tour](#five-minute-tour)
- [Two languages](#two-languages)
- [What it does](#what-it-does)
- [The physics, explained](#the-physics-explained)
- [Reproducibility and quality control](#reproducibility-and-quality-control)
- [Figures for publication](#figures-for-publication)
- [Validation](#validation)
- [Documentation](#documentation)
- [Requirements and platforms](#requirements-and-platforms)
- [Tests and contributing](#tests-and-contributing)
- [Citing](#citing)
- [License](#license)

## Install

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft.git
cd olla-dft
python3 -m venv .venv && source .venv/bin/activate   # optional but recommended
pip install .
```

That installs the `olla-dft` command and its dependencies (numpy, scipy,
matplotlib, ASE, spglib, seekpath). Quantum ESPRESSO itself is installed
separately (`apt install quantum-espresso`, `brew install quantum-espresso`,
conda, or from source); Olla-DFT only needs it to *run* calculations — preparing
inputs and post-processing outputs brought from another machine works without
a single QE binary on your laptop.

Then tell it where your pseudopotentials are (the
[SSSP](https://www.materialscloud.org/discover/sssp) library is recommended):

```bash
olla-dft config set pseudo_dir ~/pseudos/SSSP_efficiency
olla-dft sistema        # what Olla-DFT sees on this machine: QE binaries, MPI, encoding
```

Optional extras: `pip install "olla-dft[mlip]"` for machine-learned
potentials (MACE, ~1.2 GB) and `pip install "olla-dft[kappa]"` for lattice
thermal conductivity (phono3py).

## Five-minute tour

```bash
olla-dft                          # interactive menu, no flags to remember
olla-dft start --structure Si.cif # guided project for people new to the CLI
olla-dft recetas primero          # "I just installed it, where do I start?"
```

The direct, scriptable way:

```bash
olla-dft info Si.cif                        # symmetry, space group, sites
olla-dft gen Si.cif -p all -o si --insulator  # scf + nscf + bands + DOS inputs, run.sh and run.py
cd si && ./run.sh                           # run Quantum ESPRESSO
olla-dft bands . --journal aps              # gap report + bands.pdf/.png + BANDAS.dat
olla-dft dos . --mode element               # DOS/PDOS by element
olla-dft plot .                             # bands + DOS in one figure
```

Every sweep works the same way — *prepare* by default, `--run` to execute
now, `--collect` to analyse a folder that already ran:

```bash
olla-dft converge Si.cif --run              # cutoff and k-mesh convergence
olla-dft eos Si.cif --run                   # Birch–Murnaghan: V0, B0, B0', a0
olla-dft elastic Si.cif --run               # Cij, bulk/shear/Young moduli, Born stability
olla-dft phonons Si.cif --qgrid 2x2x2 --run # DFPT dispersion, DOS, F/S/Cv
olla-dft derived elastic/ELASTIC_C.dat      # Debye temperature, sound velocities, Slack κ
```

Not sure which command answers your question? Ask in plain words:

```bash
olla-dft wizard --ask "does it absorb visible light"
olla-dft teoria eos                         # the physics behind a command
```

## Two languages

This repository is the **English version** of Olla-DFT: README, documentation,
examples and the default interface (help of every command, menu, guided start,
recipes, wizard, dashboard, HTML reference and theory). The Spanish version,
with the same code and the same tests, lives at
[https://github.com/jorgegonzalezsevilla/olla-dft-esp](https://github.com/jorgegonzalezsevilla/olla-dft-esp).

The interface can also be switched to Spanish without changing repository:

```bash
olla-dft --language es bands --help
olla-dft config set language es             # make it permanent
export OLLA_DFT_LANG=es                     # or per shell
```

The Spanish-named commands have English aliases: `recipes` (`recetas`),
`theory` (`teoria`), `system` (`sistema`). The scientific reports printed by
the analysis commands are written in Spanish (the code base's source
language); every quantity in them is documented in English in
[docs/THEORY.md](docs/THEORY.md).

## What it does

78 subcommands, grouped by task. `olla-dft --help` shows the catalogue and
[docs/COMMANDS.md](docs/COMMANDS.md) lists every option.

| Area | Commands |
|---|---|
| Getting started | `start`, `wizard`, `recetas`, `teoria`, `docs`, `sistema`, `selftest` |
| Structures and inputs | `gen`, `info`, `kpath`, `prim`, `conv`, `supercell`, `convert` |
| Electronic structure | `bands`, `dos`, `plot`, `gap`, `fermi`, `effmass`, `wannier`, `unfold`, `topology`, `hubbard` |
| Spectra and response | `optics`, `tddft`, `xanes`, `xps`, `corehole`, `charge`, `charges`, `wf`, `berry` |
| Phonons, transport and temperature | `phonons`, `elph`, `transport`, `ballistic`, `kappa`, `qha`, `thermochem`, `md`, `derived` |
| Mechanics and stability | `converge`, `eos`, `elastic`, `strain`, `layers`, `xrd`, `exfoliate`, `gamma` |
| Surfaces, defects and chemistry | `surface`, `defect`, `interface`, `adsorb`, `eform`, `align`, `esm`, `echem`, `neb`, `amorphous` |
| Automation and quality | `doctor`, `audit`, `crosscheck`, `cost`, `db`, `hull`, `mlip`, `suggest`, `datasheet`, `report`, `compare`, `tune`, `results`, `campaign`, `pseudos` |
| Project | `project` |
| Appearance and configuration | `templates`, `config` |

Some highlights:

- **Inputs that run unedited.** Presets for scf, relax, vc-relax, nscf, bands,
  DOS, MD; cutoffs read from the UPF headers; spin polarisation, DFT+U (both
  QE syntaxes), spin–orbit with pseudopotential checks, hybrid functionals,
  dipole corrections.
- **Post-processing of what QE already left behind.** Gap (direct/indirect,
  VBM/CBM, per spin channel), fat bands, d-band centre, effective masses,
  Fermi surface (BXSF), band unfolding from the wavefunction files.
- **Layered materials.** Layer detection by connectivity, simulated powder
  XRD compared with an experimental pattern, exfoliation energy.
- **Wannier functions without wannier90.** Projection, Marzari–Vanderbilt
  localisation and Souza–Marzari–Vanderbilt disentanglement done in Python
  from the `pw2wannier90.x` overlaps; interpolated bands and DOS; Chern
  numbers and Wilson loops; Berry-phase polarisation and Born charges.
- **Spectroscopy.** ε(ω), n, k, α, R, Tauc gap and scissor with
  Kramers–Kronig; TDDFPT excitons; core-hole pseudopotentials, XANES and
  XPS core-level shifts; Raman activities.
- **Thermal and transport.** Harmonic and quasi-harmonic thermodynamics,
  Debye temperature and Slack κ from Cij, phonon BTE thermal conductivity,
  electron–phonon λ and Allen–Dynes Tc, Boltzmann transport (Seebeck, σ/τ,
  Lorenz number, spin-resolved), Landauer conductance.
- **Surfaces and chemistry.** Slabs, surface energy, work function, ESM
  charged surfaces, adsorption sites, charged defects with Madelung
  correction and transition levels, band alignment, NEB barriers,
  computational hydrogen electrode (HER/OER), gas-phase thermochemistry.
- **Machine-learned potentials (optional).** Pre-relax, bracket the EOS and
  screen dynamical stability with MACE before spending DFT time; melt–quench
  amorphous solids. An MLIP energy is never mixed with a DFT one without the
  provenance saying so.

## The physics, explained

Every scientific command has a section written for non-experts that says
what it answers, the formulas the code actually implements (with every
variable defined), the step-by-step procedure with the Python function and
QE binary responsible for each step, a table of *where each number comes
from*, the limits and pitfalls, and the references.

```bash
olla-dft teoria                 # index
olla-dft teoria elastic         # one command
olla-dft theory --all -o theory.md --language en
```

The same text is published as [docs/THEORY.md](docs/THEORY.md).

## Reproducibility and quality control

- **Provenance.** Every `.dat` table and every figure carries the Olla-DFT
  version, date, exact command line and parameters (also in the PDF/PNG
  metadata).
- **`doctor`** tells charge sloshing apart from slow convergence — they need
  opposite remedies — and refuses to guess with too few iterations.
- **`audit`** detects when two runs are not comparable (functional,
  pseudopotentials, cutoffs, occupations, MLIP vs DFT), the most expensive
  silent error in DFT.
- **`crosscheck`** computes the same quantity by two independent routes
  (e.g. B0 from the EOS against the trace of Cij, the Berry phase against the
  Wannier centres) and reports the discrepancy.
- **`selftest`** checks the code against published values, not against
  itself (Madelung constant, Sackur–Tetrode entropy, Allen–Dynes Tc of Al,
  Chern number of the QWZ model, …).
- **Recipes** (`olla-dft recetas`) are complete sessions that show which file
  each step leaves and which later step reads it; a test validates every
  command in them against the real argument parser, so an example cannot
  go stale.
- **Errors are honest.** A usage error (exit code 2) prints what to fix and no
  traceback; a program failure (exit code 1) is recorded locally with the
  command, trace and versions (`olla-dft report`), and QE failures are
  translated into a probable cause plus the tail of the log. Nothing is ever
  sent anywhere.

## Figures for publication

Vector output at the exact column width of the journal (`--journal aps`,
`acs`, `nature`, `elsevier`, …), interchangeable visual templates
(`journal`, `latex`, `latex-true`, `minimal`, `dark`, `slides`, `poster`,
`mono`), LaTeX typography with or without a TeX installation, colour-blind
safe palettes validated in OKLab, and a monochrome mode for journals that
charge for colour. See `olla-dft templates list` and the gallery in
[examples/plantillas](examples/plantillas).

<p align="center">
<img src="examples/plantillas/galeria_plantillas.png" alt="The same figure rendered in every visual template" width="760">
<br>
<sub>The same silicon figure in the <code>journal</code>, <code>latex</code>, <code>minimal</code>, <code>dark</code>, <code>slides</code>, <code>poster</code> and <code>mono</code> templates.</sub>
</p>

## Validation

The whole cycle (generate → run QE → post-process) was validated end to end
with Quantum ESPRESSO 6.6 on silicon, aluminium and bcc iron, and each module
against experiment or literature: Si phonons within 1–6 % of neutron data,
graphene work function 4.54 eV (exp. 4.6), Al(111) work function 4.24 eV
with ESM (exp. 4.24–4.26), Born charge of cubic BN 1.94 e (lit. 1.92),
Si lattice thermal conductivity 101 W/m·K in RTA, f-sum rule of ε₂ fulfilled
to 0.1 %, XRD peak positions within 0.05° of the PDF cards. The full list,
with references, is in [docs/VALIDATION.md](docs/VALIDATION.md). The folder
[examples/](examples/) contains real outputs and figures, not mock-ups.

<p align="center">
<img src="examples/demo_propiedades/fonones_Si.png" alt="Silicon phonon dispersion and DOS" width="360">
<img src="examples/demo_calculo/elastic.png" alt="Elastic constants of silicon" width="360">
<br>
<sub>Left: DFPT phonons of silicon (<code>olla-dft phonons</code>). Right: elastic constants (<code>olla-dft elastic</code>).</sub>
</p>

## Documentation

| Document | Content |
|---|---|
| [docs/THEORY.md](docs/THEORY.md) | The physics behind each command, for non-experts |
| [docs/COMMANDS.md](docs/COMMANDS.md) | Every command and option |
| [docs/VALIDATION.md](docs/VALIDATION.md) | Results against experiment and literature |
| [docs/PLATFORMS.md](docs/PLATFORMS.md) | Linux, macOS, Windows; requirements; configuration folders |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | How the code is organised and how to add a command |
| [examples/README.md](examples/README.md) | Worked examples with real data |
| `olla-dft docs -o reference.html` | Browsable HTML reference generated from the code |

## Requirements and platforms

Python ≥ 3.9 on Linux, macOS or Windows (native or WSL2). Quantum ESPRESSO
(`pw.x`, `ph.x`, `dos.x`, `projwfc.x`, `bands.x`, `pp.x`, `epsilon.x`,
`q2r.x`, `matdyn.x`, `dynmat.x`) to run calculations; some modules need
binaries that QE does not build by default (`make ld1 xspectra hp neb tddfpt
pwcond`) and Olla-DFT tells you the exact `make` target when one is missing.
Details, including the Windows console encoding story and `--ascii`, in
[docs/PLATFORMS.md](docs/PLATFORMS.md).

## Tests and contributing

```bash
pip install -e ".[test]"
python -m pytest -q          # 977 tests, ~30 s, no QE needed (real QE outputs live in tests/datos/)
python -m pyflakes qekit tests
```

Olla-DFT is written and maintained by a single author and does not accept
pull requests. Bug reports, questions and feature requests are very welcome
as [issues](https://github.com/jorgegonzalezsevilla/olla-dft/issues) — see
[CONTRIBUTING.md](CONTRIBUTING.md). A bug report is most useful with the
output of `olla-dft report --export incidencias.json`.

## Citing

If Olla-DFT helps your work, please cite it (see [CITATION.cff](CITATION.cff)):

> J. E. González Sevilla, *Olla-DFT: a command-line toolkit for Quantum
> ESPRESSO*, version 1.0.0 (2026). Zenodo. https://doi.org/10.5281/zenodo.22263122

and cite Quantum ESPRESSO, and the pseudopotential library you used, as their
authors request.

## License

Olla-DFT is free software under the
[GNU General Public License v3.0](LICENSE).
Copyright © 2026 Jorge Enrique González Sevilla.

It depends on numpy, scipy, matplotlib, ASE, spglib and seekpath, bundles the
atomic scattering-factor table from pymatgen (MIT), and drives Quantum
ESPRESSO as a separate process; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

---

<p align="center">
Made with love in Guadalajara, Jalisco, Mexico.<br>
<sub>A one-person project, on free time, for the pleasure of getting the physics right.</sub>
</p>
