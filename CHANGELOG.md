# Changelog

All notable changes to Olla-DFT. Dates are ISO 8601.

## Unreleased

**Scientific values change.** Three correctness fixes alter numbers reported by
earlier versions. Results produced before this release should be recomputed for
the affected commands.

- `transport`: fix the transformation of dE/dk from fractional to Cartesian
  coordinates (the transpose was inverted). Band velocities, and therefore σ/τ,
  κ_e/τ and the Seebeck tensor, were wrong for every cell whose reciprocal
  lattice is not symmetric — hexagonal, trigonal, monoclinic and triclinic. The
  Seebeck coefficient of cubic, tetragonal and orthorhombic cells is unaffected.
- `transport`: carry the spin degeneracy in the k-point weights, as Quantum
  ESPRESSO does. Without spin polarization σ/τ, κ_e/τ and the power factor were
  reported at half their value; the Seebeck coefficient and the Lorenz number
  are ratios and were already correct.
- `echem --her`: apply the potential and pH terms with the sign of a reduction.
  The HER consumes H⁺+e⁻, so ΔG = ΔG₀ + eU: its limiting potential is negative
  versus RHE, and the ΔG(U,pH) column and the energy ladder now move in the
  right direction. The ΔG_H* descriptor and the overpotential magnitude are
  unchanged.
- `kappa`: weight the cumulative curve with phono3py's own per-mode `mode_kappa`
  instead of reconstructing C·v²·τ/3 by hand, and build the mean free path from
  the effective linewidth Γ_ph-ph + Γ_iso + Γ_boundary, the same sum phono3py
  uses for κ. The curve now adds up to the reported κ by construction and
  matches phono3py's own `kaccum`. Before, `--isotopes` and `--grain` changed κ
  but left the mean free paths at their pure-infinite-crystal values: a silicon
  run with a 1 µm grain reported Λ₉₀ = 14 µm, longer than the grain itself.
- `kappa`: warn that Λ₉₀ is not converged below a 25³ q-grid. Measured in
  silicon with a Stillinger-Weber potential on a 3×3×3 supercell, Λ₉₀ grows from
  4.4 to 21.9 µm between 11³ and 31³ while Λ₅₀ moves only 21 % (0.68 → 0.82 µm).
- `crosscheck`: stop applying a kbar→GPa factor to a pressure column that is
  already in GPa. The third route to B₀ was reported ten times too small and
  disagreed with the equation of state even when both agreed.

- `kappa --grain`: pass the grain size to phono3py in micrometres, the unit it
  documents, instead of converting it to Angstrom. The conversion made boundary
  scattering 10,000 times weaker than requested, so `--grain` had no effect on κ
  while the report stated it had been applied.
- `kappa`: use phono3py's own lifetime convention, τ = 1/(2·2π·Γ). The missing
  2π — which converts cyclic to angular frequency — made every reported mean
  free path 2π times too long. κ, its temperature dependence and the shape of
  the cumulative curve are unaffected (the factor cancels on normalisation);
  only the Λ axis moves. **The Λ₅₀ = 1.0 µm figure recorded in the validation
  documents was obtained with the old convention and is pending recomputation.**

Other fixes, with no effect on scientific values:

- Do not split an executable path on spaces when building the pw.x command.
  `C:\Program Files\QE\bin\pw.exe` — the default location on Windows —
  failed with «not found: C:\Program».
- Read `KPATH.txt` as UTF-8: the Γ, Δ, Σ and Λ tick labels of band figures came
  out as mojibake on Windows, or raised `UnicodeDecodeError`.
- Validate numeric configuration values in `config set`, instead of accepting
  them and failing later in another command.
- Report bad `--position`, `--miller` and `--fix` values as usage errors rather
  than as program incidents with a traceback; `--miller` now checks it got
  three indices.
- Write the `kappa` run script with POSIX line endings and the execute bit, as
  the other generated scripts already did.
- Sort the volumes in `qha` before fitting, instead of assuming the input file
  is ordered.
- Return no accumulation in `kappa` when no mode contributes, instead of
  raising `IndexError`.
- Open documentation with a valid `file://` URI on Windows.
- Write every generated file as UTF-8 explicitly, and read back as UTF-8 the
  files the program itself wrote. 51 writes across 33 modules used the machine's
  locale encoding: on Windows, exporting any report containing Δ, Å, κ or ⁻¹
  raised `UnicodeEncodeError`, and JSON state written as UTF-8 came back as
  mojibake. A test now enforces the convention for the whole package.
- `exfoliate`: warn when the layers in the cell are not equivalent. E_exf
  divides E(bulk) by the number of layers, which is meaningless for a
  heterostructure; the number came out silently.
- `converge`: reject malformed `--values` k-meshes (`4x4` raised `IndexError`),
  and require two finished calculations before drawing the curve (one raised a
  matplotlib error). The "not converged" branch of the report was unreachable —
  the densest point is its own reference and always qualifies — so its advice
  never reached anyone; the reachable branch now carries it.
- `doctor --system`: measure available memory on macOS and Windows instead of
  only Linux, look for the platform's binary names (`pw.exe` on Windows), and
  stop reporting a Quantum ESPRESSO installation as fine when pw.x is missing.
- `kappa`: build the report warnings without mutating the run, so they are no
  longer duplicated in `KAPPA.txt`; and say so when phono3py does not expose the
  per-mode data instead of dropping a section, a file and a figure in silence.
- Run the lattice-thermal-conductivity tests in CI. Their `phono3py` skip was at
  module level, so all 30 were skipped and three of them had been failing since
  the 1.5.0 English translation without anyone seeing it.
- Add the resilience scope caveat to `README.de.md`, which did not mention
  `resilient` at all, so German readers never saw that recovery after a physical
  power outage or disk loss has not been demonstrated.
- Write Hubbard cards, thermochemistry reports and exported themes as UTF-8.

## 1.5.0 — 2026-09-05

- Write scientific reports, diagnostics and new plot labels in English in every interface locale. Text output changes; structured identifiers, scientific values and calculation logic remain compatible.
- Add Deutsch to the startup selector, saved preferences and `--language de`.
- Correct the energy-axis labels of band and combined band/DOS figures, including unshifted energies.
- Include German command help, menus, guided setup, recipes, wizard and explorer catalogs, a German README and generated command reference. Detailed theory for the German interface is shared in English.
- Generate the initial explorer labels in the selected language before JavaScript starts.
- Add `--all-languages` for reference/dashboard generation; retain the existing English/Spanish `--both` option. Standalone dashboards no longer link to language files that were not generated.
- Use English as the default for standalone documentation and visualization APIs.

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
