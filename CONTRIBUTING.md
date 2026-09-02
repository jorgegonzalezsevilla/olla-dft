# Contributing to Olla-DFT

## How this project is developed

Olla-DFT is written and maintained by a single author, Jorge Enrique González
Sevilla, and it will stay that way. **Pull requests are not accepted**: any
PR opened on this repository will be closed without review. This is not a
comment on your work; it is simply how the project is run.

What is very welcome is your feedback. If you find a bug, get a wrong number,
or would like a feature, **open an issue** at
https://github.com/jorgegonzalezsevilla/olla-dft/issues and describe it. The
author reads every issue and decides what goes into the next version.

## Reporting bugs

When a command fails unexpectedly, Olla-DFT records the incident on your
machine — exact command, traceback, versions of Python and the dependencies,
whether QE was available — and prints its id. Package everything with

```bash
olla-dft report --export incidencias.json
```

and attach that file to the issue, together with the structure file if it is
not confidential (`--attach archivo.cif` copies it into the log).
`olla-dft report "what happened"` records something that did not crash but
was confusing; `olla-dft report --stats` shows which commands fail most.
Nothing is ever sent automatically: there is no telemetry, the log lives in
your configuration folder and you decide whether to share it.

For a wrong number rather than a crash, say which reference you compared
against and where it comes from; `olla-dft crosscheck` and `olla-dft selftest`
output helps.

## Requesting features

Open an issue with:

- what you are trying to compute or automate, and why;
- what Olla-DFT does today and where it falls short;
- if possible, a small example (structure file, command line, expected
  result).

Clear use cases are what most often make it into a release.

## Running the test suite yourself

You are free to clone, read, and run the code under the GPL. If you want to
check that everything works on your machine:

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft
cd olla-dft
python -m venv .venv && source .venv/bin/activate   # optional
pip install -e ".[test]"
python -m pytest -q                 # no QE needed, under a minute
```

Quantum ESPRESSO is not needed to run the test suite; `olla-dft sistema`
tells you what your machine has. Requirements and per-OS notes are in
[docs/PLATFORMS.md](docs/PLATFORMS.md).

## Licence

Olla-DFT is free software under the GNU General Public License, version 3
(see [LICENSE](LICENSE)). Third-party components and data are listed in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Copyright © 2026 Jorge Enrique González Sevilla.
