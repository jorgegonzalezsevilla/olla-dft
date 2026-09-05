# One application, two interface languages

[English README](../README.md) · [README en español](../README.es.md)

From **1.4.0**, `olla-dft` is the only maintained repository and package.
English and Spanish share the same commands, scientific implementation and releases.
The former Spanish repository is retained as a read-only historical archive.

## Choose a language

Run `olla-dft` in a terminal. Choose **1 — English** or **2 — Español**.
Press Enter to keep the suggested language. The menu remembers the choice;
press **l** to change it without restarting. Choosing a language changes only
interface text, not units, calculation settings or result data. If preferences
cannot be saved, the selected language still works for the current session.

| Use | Command |
|---|---|
| Open directly in Spanish | `olla-dft --language es` |
| Open directly in English | `olla-dft --language en` |
| Save Spanish as the default | `olla-dft config set language es` |
| Help in Spanish | `olla-dft eos --help --language es` |
| Guided project in English | `olla-dft start --language en` |

For direct commands, precedence is `--language` → `OLLA_DFT_LANG` → saved
preference → English. The menu selection applies explicitly to its child
commands, even when the environment has a different language. An explicit
`--language` bypasses the startup selector and does not save a preference.
With no command and no interactive terminal, Olla-DFT prints help and exits.
Direct commands, `--help` and `--version` never show the selector.
This changes older heredoc-driven menu sessions: replace them with direct
commands in batch scripts. Language codes are `en` and `es`; the environment
also accepts uppercase codes. Other environment values fall back to the saved
preference (or English); unsupported `--language` values are usage errors.
Saving preferences is atomic, but simultaneous edits by multiple sessions
are not merged. The configuration file is written with owner-only permissions.

Translated surfaces include command help, menus, guided setup, recipes,
theory, the dashboard, reference HTML and the result explorer. **Some scientific
reports and error messages remain in Spanish.** Raw Quantum ESPRESSO output,
command names, file formats and scientific identifiers are preserved.
The explorer also has its own language selector; it does not change the
terminal preference.

## Migrar desde la edición española

El paquete ya se llamaba `olla-dft` en ambas ediciones: no instales dos copias
en el mismo entorno. Activa el entorno de Python que quieres actualizar e
instala la versión unificada desde el repositorio principal:

```bash
python -m pip install --upgrade "olla-dft @ git+https://github.com/jorgegonzalezsevilla/olla-dft.git@v1.4.0"
olla-dft config set language es
olla-dft --version
olla-dft
```

Esto sustituye la instalación de Python, incluidas las instalaciones editables;
no borra el clon anterior ni tus proyectos. Si prefieres seguir desarrollando
con una instalación editable, clona el repositorio principal y ejecuta
`python -m pip install -e .` dentro de ese nuevo clon. Conserva aparte cualquier
cambio local que tengas en el repositorio anterior.

La carpeta de configuración sigue siendo la misma, incluidas las rutas de
pseudopotenciales y parámetros guardados. Si antes usabas español sin una
preferencia guardada, elige Español una vez. Los comandos, alias en español,
`python -m qekit`, proyectos y formatos de resultados se conservan. Cambiar
la instalación usada por un trabajo de QE en curso requiere una validación
separada; esta unificación no demuestra la migración de un trabajo activo.

El actualizador de las versiones españolas antiguas consulta su repositorio
anterior: **haz esta primera migración con pip**. A partir de 1.4.0 ambos idiomas
consultan las versiones de `olla-dft`.

La referencia y la teoría se mantienen juntas en [español](COMANDOS.md) y
[inglés](COMMANDS.md); los ejemplos incluyen `README.es.md` y `README.md`.
Se cita una sola aplicación mediante [CITATION.cff](../CITATION.cff). Los DOI
históricos permanecen disponibles y no se crean nuevas ediciones españolas
en Zenodo.
