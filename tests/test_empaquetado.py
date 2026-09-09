"""Que el paquete se pueda publicar y que su ficha de PyPI se vea bien.

La ficha de PyPI se genera a partir de README.md, pero PyPI NO resuelve los
enlaces relativos contra el repositorio: un `](docs/COMMANDS.md)` queda como
un 404 en la primera página que ve alguien que acaba de descubrir el
proyecto. Estas pruebas fijan las dos condiciones que hacen falta: que los
enlaces sean absolutos y que apunten a algo que existe.
"""

import re
from pathlib import Path

import pytest

from tests.marcadores import solo_repositorio

# tomllib entró en Python 3.11 y el paquete soporta desde 3.9, así que las dos
# pruebas que leen pyproject.toml se saltan en las versiones viejas: el archivo
# es el mismo, y en CI lo cubren las matrices 3.11+.
try:
    import tomllib
except ModuleNotFoundError:                     # pragma: no cover
    tomllib = None

necesita_tomllib = pytest.mark.skipif(
    tomllib is None, reason="tomllib necesita Python 3.11 o superior")

RAIZ = Path(__file__).resolve().parent.parent
READMES = ("README.md", "README.es.md", "README.de.md")
BASE = "https://github.com/jorgegonzalezsevilla/olla-dft/blob/main/"
ARBOL = "https://github.com/jorgegonzalezsevilla/olla-dft/tree/main/"


def _enlaces(nombre):
    texto = (RAIZ / nombre).read_text(encoding="utf-8")
    return re.findall(r"\]\(([^)]+)\)", texto)


@pytest.mark.parametrize("nombre", READMES)
def test_los_readme_no_llevan_enlaces_relativos(nombre):
    relativos = [d for d in _enlaces(nombre)
                 if not d.startswith(("http://", "https://", "#", "mailto:"))]
    assert not relativos, (
        f"{nombre} tiene enlaces relativos; en PyPI serían 404:\n  "
        + "\n  ".join(relativos))


@solo_repositorio
@pytest.mark.parametrize("nombre", READMES)
def test_los_enlaces_al_repositorio_apuntan_a_algo_que_existe(nombre):
    rotos = []
    for destino in _enlaces(nombre):
        for prefijo in (BASE, ARBOL):
            if destino.startswith(prefijo):
                rel = destino[len(prefijo):].split("#")[0].rstrip("/")
                if rel and not (RAIZ / rel).exists():
                    rotos.append(rel)
    assert not rotos, f"{nombre} apunta a archivos que no existen: {rotos}"


@necesita_tomllib
def test_los_metadatos_de_pypi_estan_completos():
    cfg = tomllib.loads((RAIZ / "pyproject.toml").read_text(encoding="utf-8"))
    p = cfg["project"]
    assert p["name"] == "olla-dft"
    assert p["readme"] == "README.md"
    assert p["license"] == "AGPL-3.0-or-later", "la licencia debe ir en SPDX (PEP 639)"
    assert p["license-files"], "sin license-files la rueda no lleva la licencia"
    # El clasificador de licencia está obsoleto desde PEP 639 y duplica el campo
    assert not [c for c in p["classifiers"] if c.startswith("License ::")], \
        "quita el clasificador de licencia: lo sustituye la expresión SPDX"
    for clave in ("Homepage", "Documentation", "Source", "Changelog", "Issues"):
        assert clave in p["urls"], f"falta el enlace '{clave}' en la ficha"
    assert p["description"] and len(p["description"]) < 200
    assert cfg["build-system"]["requires"] == ["setuptools>=77"], \
        "PEP 639 (licencia SPDX) necesita setuptools >= 77 al construir"


@necesita_tomllib
def test_la_version_es_la_misma_en_todas_partes():
    """Cinco archivos declaran la versión; en PyPI solo se sube una vez."""
    import qekit
    cfg = tomllib.loads((RAIZ / "pyproject.toml").read_text(encoding="utf-8"))
    v = cfg["project"]["version"]
    assert qekit.__version__ == v
    assert f"version: {v}" in (RAIZ / "CITATION.cff").read_text(encoding="utf-8")
    assert f'"version": "{v}"' in (RAIZ / ".zenodo.json").read_text(encoding="utf-8")
    assert f"## {v} " in (RAIZ / "CHANGELOG.md").read_text(encoding="utf-8")
