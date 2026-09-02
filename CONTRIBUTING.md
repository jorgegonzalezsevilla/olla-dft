# Contributing to Olla-DFT

## Set up

```bash
git clone https://github.com/jorgegonzalezsevilla/olla-dft-en
cd olla-dft-en
python -m venv .venv && source .venv/bin/activate   # optional
pip install -e ".[test]"
```

That installs the `olla-dft` command in editable mode plus pytest and
pyflakes. Quantum ESPRESSO is not needed to develop or to run the test suite;
`olla-dft sistema` tells you what your machine has. Requirements and per-OS
notes are in [docs/PLATFORMS.md](docs/PLATFORMS.md).

## Before opening a pull request

```bash
python -m pytest -q                 # the whole suite, no QE needed, under a minute
python -m pyflakes qekit tests      # must print nothing
python tools/build_docs.py          # regenerates docs/COMMANDS.md and docs/THEORY.md
```

Run `build_docs.py` whenever you add, rename or change the options of a
command, or touch `qekit/data/theory/`: `tests/test_docs.py` and
`tests/test_teoria.py` fail if the generated files are stale. Do not edit those
four files by hand. Optional but appreciated: `olla-dft selftest` (seconds,
no QE) and, if you have `pw.x`, `olla-dft selftest --full --pseudo-dir
/path/to/upf` (about ten minutes).

## Style conventions

- **Spanish identifiers and comments are the house style.** Functions,
  variables, docstrings, help strings and the scientific reports are written
  in Spanish, as the existing code is (`ErrorDeUso`, `preparar`, `informe`).
  Do not translate existing names.
- **English goes in the i18n tables.** Every user-facing interface string has
  its English counterpart in `qekit/data/i18n/`: help strings in
  `cli_en.json` (keyed by the Spanish text), one-line command summaries in
  `docs_en.json`, menu/onboarding/dashboard labels in their `_en.json`, and
  the theory in `qekit/data/theory/*.en.md`. The scientific reports stay in
  Spanish. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the layer.
- **Every new option needs a help string** (`help=` in Spanish, translation in
  `cli_en.json`). Options that repeat across commands take theirs from the
  `defaults` table of `cli_es.json`.
- **Every physics change needs a test and a note in the theory docs.** A new
  formula, constant, default or correction gets a test in `tests/` (with a
  real QE output in `tests/datos/` if it reads one, and a frozen value in
  `tests/referencias.py` if it was validated against experiment) and an update
  of the matching section in `qekit/data/theory/<area>.es.md` **and**
  `<area>.en.md` (the test suite checks es/en parity and the mandatory
  subsections). Never change a value in `tests/referencias.py` just to make a
  test pass: it is a regression detector, and it is updated only when the new
  value was validated again against the external source.
- Every `.dat` and figure goes through `core/provenance.py`; every sweep goes
  through `modules/sweep.py` and follows prepare / `--run` / `--collect`.
- Usage errors raise `ErrorDeUso` (exit 2, clean message, no trace). Anything
  else is a program failure (exit 1) and is logged locally.
- Output must survive a cp1252 console: if you print a new non-ASCII symbol,
  add it to `TRANSLITERACION` in `core/consola.py` or `test_portabilidad.py`
  will tell you.
- Keep `COMMAND_GROUPS` (cli.py), `docs.GRUPOS`/`MODULO_DE` (docs.py) and
  the theory sections in sync when adding a command; the step-by-step list is
  in [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md#adding-a-command).
- No emojis in code, output or documentation. Line width about 79 columns.

## Reporting bugs

When a command fails unexpectedly, Olla-DFT records the incident on your
machine — exact command, traceback, versions of Python and the dependencies,
whether QE was available — and prints its id. Package everything with

```bash
olla-dft report --export incidencias.json
```

and attach that file to the issue at
https://github.com/jorgegonzalezsevilla/olla-dft-en/issues, together with the
structure file if it is not confidential (`--attach archivo.cif` copies it
into the log). `olla-dft report "what happened"` records something that did
not crash but was confusing; `olla-dft report --stats` shows which commands
fail most. Nothing is ever sent automatically: there is no telemetry, the
log lives in your configuration folder and you decide whether to share it.

For a wrong number rather than a crash, say which reference you compared
against and where it comes from; `olla-dft crosscheck` and `olla-dft selftest`
output helps.

## Licence

Olla-DFT is free software under the GNU General Public License, version 3
(see [LICENSE](LICENSE)). By contributing you agree that your contribution is
distributed under the same licence. Third-party components and data are
listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Copyright © 2026 Jorge Enrique González Sevilla.
