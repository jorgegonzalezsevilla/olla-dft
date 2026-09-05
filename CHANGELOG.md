# Changelog

All notable changes to Olla-DFT. Dates are ISO 8601.

## 1.4.0 — 2026-09-05

**Interactive behavior change:** no-command invocations with redirected input or output now print help and exit (status 0). Use direct commands instead of driving the menu with a heredoc. Interactive startup offers a language choice unless `--language` is supplied.

- Unify English and Spanish in one package and canonical repository, with English as the fallback language.
- Add a bilingual selector on interactive menu startup, remember the choice, and switch languages from the menu with `l`.
- Pass the selected language through nested menu commands and catalog help, including when an environment override is set.
- Keep direct commands, help and non-terminal invocations free of language prompts; `--language` bypasses the startup selector.
- Include the Spanish README, reference, theory, platform guide and example guides alongside English documentation.
- Use one update source and Zenodo software identity; provide migration instructions for the former Spanish edition.
- Read explorer translations and assets explicitly as UTF-8 on all platforms.
- Scientific algorithms, calculation input formats and checkpoint formats are unchanged. Scientific report text is not fully translated.

## 1.3.1 — 2026-09-04

- Load Matplotlib and font discovery only when drawing figures; structure commands and help retain the same options without initializing plotting.
- Keep EOS fitting independent of structure preparation and plotting imports.
- Add fresh-process checks for English/Spanish structure commands, all three EOS fits, and figure export after deferred imports.
- Scientific formulas, tolerances, generated inputs and benchmark competitors are unchanged.

## 1.3.0 — 2026-09-04

- License project-owned software under AGPL-3.0-or-later; previous GPL releases remain available under their original terms.
- Align source notices, package metadata, citation, GitHub and Zenodo metadata; preserve third-party licenses and existing scientific example licenses.
- Include the full software license and a versioned source link in exported interactive HTML.
- Clarify the project-generated core-hole pseudopotential test fixture in the third-party inventory.
- No changes to scientific algorithms, calculation parameters or benchmark results.

## 1.2.0 — 2026-09-04

- Add an offline result explorer with numeric axes, filtering, record selection and customizable figures.
- Export scoped SVG, PNG, CSV, JSON and self-contained HTML; preserve original units, precision and uncertainties in data exports.
- Replace the dashboard's misleading connected energy series with the interactive explorer.
- Record k-grid and calculation parameters, warn about mixed methods, and declare snapshot/export limits.

- Honor stop requests during restore, reap QE after PID-record failures, and clear stale terminal PIDs.

- Add `resilient init/run/status/pause/service` for recoverable pw.x jobs.
- Preserve two complete verified checkpoints and restore private workspaces after interruption.
- Freeze input, UPFs, MPI command, threads, libraries and architecture; bound consecutive failures.
- Generate a Linux service for automatic worker restart on a retained persistent disk.
- Record restart/copy/compute overhead for Olla-Lungo cost comparisons without changing physical inputs.
- Validate local QE 7.4 SCF, relax and vc-relax recovery after SIGTERM and SIGKILL; physical power loss and disk loss recovery remain unmeasured.

## 1.1.1 — 2026-09-04

- Fix gap analysis: occupied-only calculations report insufficient bands, never a false metal.
- Report unconverged calculations and return nonzero status for an unavailable or unvalidated gap.
- Reject nonpositive k-grids/cutoffs and nonfinite sampling parameters before writing inputs.
- Add regression tests for the September benchmark audit.

## 1.1.0 — 2026-09-03

Changes driven by the first benchmark runs (olla-dft-bench):

- Start-up is about 8× cheaper: `import qekit.cli` went from ~0.6 s to ~0.07 s.
  seekpath, ase.io, matplotlib, strain and defects are now imported on first
  use instead of on every invocation.
- `gen --kgrid N N N`: explicit scf/relax k-grid (overrides `--kspacing` and
  `--klevel`). Until now only a spacing could be given.
- `mixing_beta` is 0.7 (QE's default) with fixed occupations and stays 0.4
  with smearing. On the benchmark's Si cell this takes the scf from 14 to
  about 7 iterations with the same energy.

## 1.0.1 — 2026-09-03

- New `olla-dft update` (alias `actualizar`): checks the latest published
  release on GitHub, shows what is new and the exact commands it would run,
  and installs it only after confirmation (`--check` to only look, `--yes` to
  skip the question, `--version TAG` for a specific release). Works for
  installs made with pip from GitHub and for local clones. Olla-DFT never
  checks for updates on its own.
- Links to the reproducible benchmark against ASE, pymatgen and seekpath
  (olla-dft-bench).

## 1.0.0 — 2026-09-02

First public release (GPL-3.0).

- Single command `olla-dft`; the Python package keeps the name `qekit` so
  existing scripts and project folders (`.qekit/`) keep working. Configuration
  and data now live in an `olla-dft` folder (`~/.config/olla-dft` on Linux),
  migrated automatically from the previous `qekit`/`QEkit` folders.
- Bilingual interface (Spanish default, English with `--language en`,
  `OLLA_DFT_LANG` or `config set language en`): help of all 78 commands and
  1 300 options, interactive menu, guided start, recipes, wizard, dashboard,
  HTML reference. English aliases `recipes`, `theory`, `system`.
- New `olla-dft teoria` / `theory`: the physics behind every scientific
  command — what it answers, formulas actually implemented, procedure with
  the responsible function and QE binary, where each number comes from,
  limits and references — also published as `docs/THEORY.md`.
- Help screens grouped into *options* / *execution* / *DFT parameters* /
  *figure*; every option now has a help string.
- Removed the experimental platform layer that added no physics: local LLM
  assistant, web server, plugins, HPC monitor/submission, release/SBOM
  preflight and external database connectors.
- Audit of the formulas against the code fixed 30+ defects, among them:
  Bader charges integrated with the wrong volume unit (≈6.75× too large);
  defect potential alignment added in Ry instead of eV; `adsorb --dipole`
  not writing the dipole correction; `surface --fix` not reaching
  `ATOMIC_POSITIONS`; `tddft --compare` reading reflectivity instead of
  absorption; `--raman` without `--gamma` crashing after the whole DFPT
  chain; `kappa` forcing fixed occupations on metals; QHA lattice parameter
  wrong for 1-atom primitive cells; band unfolding mixing spin channels in
  lsda runs; Allen–Dynes f₂ shape factor never applied; band-alignment
  figure drawing the wrong CBM; folded Berry-phase markers using the wrong
  modulus; `gen --soc` not verifying relativistic pseudopotentials;
  `doctor` mixing SCF cycles of a relaxation; work-function flatness
  measured outside the vacuum; `transport --spin-resolved` and `--kspacing`
  not wired; Tc column always empty in `elph`; XANES accepting edges
  xspectra.x cannot compute; and several misleading messages.
- New flags: `hubbard --hubbard-style`, `unfold --spin`, `kappa --metal`,
  `qha --structure`, `transport --nspin/--mag`, `tddft --scissor`,
  `charges --pseudo-dir`; `ballistic --ikind 2` removed (never implemented).
- Examples renamed by topic with bilingual READMEs whose commands are
  validated against the parser by a test; tests renamed by topic; 977 tests.
- Repository files for GitHub: GPL-3.0 licence, third-party notices,
  CONTRIBUTING, CITATION.cff, CI on Python 3.9–3.13, `.gitignore`.

## Before 1.0

The project grew privately as *QEkit* (0.1–0.34) and *Olla-DFT* (0.35). Milestones, kept here for reference:

| Version | Added |
|---|---|
| 0.1–0.4 | pw.x input generation, structure tools, bands/DOS/gap post-processing, publication styles and templates |
| 0.5–0.7 | convergence, EOS, elastic constants; layered materials (layers, XRD, exfoliation); optics, charge density, work function, DFPT phonons |
| 0.8–0.12 | effective mass, provenance, Raman, XPS, transport, Bader/Löwdin, surfaces and defects, SOC and DFT+U, doctor/audit/db/hull, MLIP, incident log, crosscheck, derived, QHA, datasheet |
| 0.13–0.14 | core-hole pseudopotentials, XANES, self-consistent Hubbard U, electron–phonon, unfolding, NEB, thermochemistry, MD analysis, interfaces, wizard, pseudopotential selection, TDDFPT, ballistic transport |
| 0.15–0.20 | strain, adsorption, 2D elastics, dipole correction, charged defects, parallel sweeps with time budget, cost estimator, surface energy, band alignment, fat bands, Hubbard V, phonons at electronic temperature, hybrids, Lorenz number, spin transport, selftest, CHE, HTML reference |
| 0.21–0.28 | amorphous solids, Wannier functions and disentanglement, Berry phase, lattice thermal conductivity, ESM, recipes, portability (Linux/macOS/Windows, ASCII output) |
| 0.29–0.35 | topology, project workflow, results database, campaigns, dashboard, bilingual guided start, rename to Olla-DFT |
