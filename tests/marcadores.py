"""Marcadores compartidos por la suite.

Un sdist NO lleva la galería ni las figuras de los ejemplos: son varios MB que
solo sirven en GitHub. Las pruebas que comprueban que la documentación del
REPOSITORIO cita archivos existentes necesitan un clon completo, y con este
marcador se saltan solas al correr la suite desde el paquete distribuido —
que es justamente la comprobación de que el sdist es autoconsistente.
"""

from pathlib import Path

import pytest

RAIZ = Path(__file__).resolve().parent.parent
ARBOL_COMPLETO = (RAIZ / "docs" / "gallery" / "manifest.json").exists()

solo_repositorio = pytest.mark.skipif(
    not ARBOL_COMPLETO,
    reason="necesita un clon completo del repositorio (el sdist no lleva la "
           "galería ni las figuras de los ejemplos)")
