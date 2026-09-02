# Platforms and requirements

## Where it runs

Olla-DFT is pure Python: it runs the same on Linux, macOS and Windows. What is
not portable is Quantum ESPRESSO, and that changes the advice per system.
Before anything else:

```
olla-dft sistema
```

(`olla-dft system` also works) prints what Olla-DFT sees of your machine — the
console encoding, where it stores its configuration, which QE binaries it
finds and where, whether there is an MPI launcher, how many cores are visible
— and what to do about whatever is missing. It is the first command to run on
a new system.

Olla-DFT is useful without a single QE binary on the machine: it generates the
inputs and post-processes outputs brought from elsewhere. Only `--run` needs
the binaries locally.

## Installing Quantum ESPRESSO per system

**Linux.** Quantum ESPRESSO is packaged in all major distributions:
`sudo apt install quantum-espresso` on Debian and Ubuntu, `sudo dnf install
quantum-espresso` on Fedora and RHEL, `sudo zypper install quantum-espresso`
on openSUSE, `sudo pacman -S quantum-espresso` on Arch, or `conda install -c
conda-forge qe`. Compiling from source at https://www.quantum-espresso.org
gives better performance. If you compiled it by hand and it is not in the
PATH: `olla-dft config set pw_cmd /path/to/bin/pw.x`. Olla-DFT also looks in
`/usr/bin`, `/usr/local/bin`, `/opt/qe/bin`, `~/q-e/bin` and
`/usr/lib64/openmpi/bin` without any configuration.

**macOS.** `brew install quantum-espresso` or `sudo port install
quantum-espresso`. On Apple Silicon Homebrew installs into
`/opt/homebrew/bin`, which is not always in the PATH of a non-interactive
shell; Olla-DFT looks there anyway (and in `/usr/local/bin`, `/opt/local/bin`
for MacPorts, and `~/q-e/bin`). The generated scripts count cores with
`sysctl -n hw.ncpu` as well as with `nproc`, which does not exist on macOS.

**Windows.** There are three routes, from most to least recommended:

1. **WSL2.** `wsl --install`, and inside Ubuntu everything behaves as on
   Linux (`sudo apt install quantum-espresso python3-pip`, then install
   Olla-DFT with pip). It is the best-tested route and the one with the
   fewest surprises.
2. **Native QE binaries.** They are called `pw.exe`, not `pw.x`; Olla-DFT
   tries both endings and also looks in `C:\Program Files\QE\bin`,
   `C:\Program Files\quantum-espresso\bin` and `C:\qe\bin`. If they are
   elsewhere: `olla-dft config set pw_cmd C:\path\to\pw.exe`. Sweeps are
   launched with `python run.py` instead of `run.sh`: **every working folder
   carries both scripts**, the `.sh` for POSIX and the `.py`, which does the
   same without needing bash or xargs and detects on its own whether
   `mpiexec` is available.
3. **Olla-DFT here, QE on a cluster.** Generate the inputs on the laptop, run
   them wherever, bring the outputs back: all the post-processing works
   without a single QE binary on your machine.

Generated `.sh` scripts are always written with POSIX line endings, even from
Windows, so they do not fail with `bad interpreter: /bin/bash^M` when run in
WSL or on a cluster.

## MPI

The launcher is detected in this order: `mpirun -np N` (OpenMPI), `mpiexec -n
N` (MPICH, MS-MPI), `srun -n N` (Slurm). If none exists the calculations run
serially, which is the right thing on a laptop; there, running several points
at once with `-j N` uses 2 or 4 cores better than a single `pw.x` under MPI.

## The console encoding and `--ascii`

The reports carry Å, α, ε, →, ①, ✓, and the legacy Windows code page (cp1252)
cannot write them: `print` raises `UnicodeEncodeError` and the command **dies
halfway through its output** with exit code 1, so it looks as if the
calculation failed. It was verified: `info`, `selftest` and `recetas` all
three crashed. Olla-DFT solves it in two steps — first it tries to switch the
console to UTF-8 (Windows 10 and later support it), and if it cannot, it
transliterates to ASCII preserving the meaning:

```
  │ a0 = 5.402 Å   κ_L → 100.7 W/m·K   ✓          (UTF-8)
  | a0 = 5.402 A   kappa_L -> 100.7 W/m.K   ok    (--ascii)
```

It never substitutes `?`: a report that says "the lattice parameter is
5.43 ?" is worse than printing nothing, because it looks like corrupt data.
`--ascii` forces the ASCII output on any system, and it is what to use when
redirecting to a file or pasting the output into an e-mail; like
`--language`, it is accepted anywhere on the command line. On Windows you can
also enable UTF-8 with `chcp 65001` or by setting `PYTHONUTF8=1`.

There are tests for the three platforms: they simulate `sys.platform`, launch
the whole program with the output forced to cp1252, and one walks through
**every report Olla-DFT really generates** checking that each non-ASCII
character can be transliterated. Adding a new symbol to a report without
teaching it to the table breaks pytest before anyone suffers it — that is how
four were found (`Λ`, `⇌`, `∝` and the typographic minus `−`, which looks
like a hyphen and is not) and a hard-wired `mpirun` in the XPS module.

## Configuration and data folders

Configuration (`config.ini`, templates) follows the convention of each
system; large, replaceable data (models, history) go to a separate folder so
that backups and desktop conventions are respected.

| | Configuration | Data |
|---|---|---|
| Linux | `~/.config/olla-dft` (respects `XDG_CONFIG_HOME`) | `~/.local/share/olla-dft` (respects `XDG_DATA_HOME`) |
| macOS | `~/Library/Application Support/olla-dft` | `~/Library/Application Support/olla-dft` |
| Windows | `%APPDATA%\olla-dft` | `%LOCALAPPDATA%\olla-dft` |

Environment variables override the defaults:

| Variable | Effect |
|---|---|
| `OLLA_DFT_CONFIG_DIR` | configuration folder (run from a USB stick, or on a cluster with a full HOME; the tests use it so they never touch the real configuration) |
| `OLLA_DFT_DATA_DIR` | data folder (put models and history on a disk with more space) |
| `OLLA_DFT_LANG` | interface language, `es` or `en` (below `--language`, above `olla-dft config set language`) |

The older names `QEKIT_CONFIG_DIR` and `QEKIT_DATA_DIR` are still honoured.
A configuration left by a previous version — in `~/.config/qekit` or in the
`QEkit` folder of each system — is copied automatically the first time
Olla-DFT runs, and `olla-dft sistema` says so when it finds one.

The local incident log (`olla-dft report`) lives inside the configuration
folder. Nothing is ever sent anywhere: there is no telemetry.

## Requirements

- Python ≥ 3.9 with numpy ≥ 1.20, scipy ≥ 1.8, matplotlib ≥ 3.5, ASE ≥ 3.22,
  spglib ≥ 2.0 and seekpath ≥ 2.0 (installed automatically by `pip install .`
  or `pip install -e .`).
- Quantum ESPRESSO to run the calculations: `pw.x`, `dos.x`, `projwfc.x`,
  `bands.x`, `pp.x`, `epsilon.x`, `ph.x`, `q2r.x`, `matdyn.x`, `dynmat.x`, and
  `pw2wannier90.x` for `wannier`. All of these are part of a standard QE
  build.
- Some modules need binaries that QE **does not build by default**. From the
  QE source folder:

  ```bash
  make ld1        # olla-dft corehole  (core-hole pseudopotentials)
  make xspectra   # olla-dft xanes
  make hp         # olla-dft hubbard
  make neb        # olla-dft neb
  make tddfpt     # olla-dft tddft     (turbo_lanczos.x, turbo_spectrum.x)
  make pwcond     # olla-dft ballistic
  ```

  Olla-DFT warns with the exact command when one is missing, instead of
  failing with a "command not found".
- Optional extras:

  | Install | Adds | Used by |
  |---|---|---|
  | `pip install "olla-dft[mlip]"` | torch, mace-torch (about 1.2 GB) | `olla-dft mlip`, `olla-dft amorphous`, MACE forces in `kappa`, `selftest --mlip` |
  | `pip install "olla-dft[kappa]"` | phono3py ≥ 3 | `olla-dft kappa` (phonon Boltzmann equation) |
  | `pip install "olla-dft[test]"` | pytest ≥ 7, pyflakes ≥ 3 | the test suite |
