# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Elegir pseudopotencial con criterio, no por orden alfabético.

EL PROBLEMA
-----------
En una carpeta de pseudopotenciales suele haber varios para el mismo
elemento: `Fe.pbe-nd-rrkjus.UPF`, `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`,
`Fe.rel-pbe-spn-kjpaw_psl.1.0.0.UPF`... Hasta ahora Olla-DFT tomaba el
primero por orden alfabético. Eso funciona por casualidad, y cuando falla
lo hace en silencio: un pseudo escalar-relativista con `lspinorb` da un
desdoblamiento espín-órbita de CERO sin ningún error, y un ultrasuave con
`epsilon.x` da un espectro entero que está mal.

La elección depende de lo que vayas a hacer, y esa información Olla-DFT ya la
tiene. Este módulo la usa.

LO QUE MIRA
-----------
- **El funcional.** Mezclar funcionales entre elementos de la misma
  estructura invalida la energía total. Es lo primero que se comprueba.
- **El tipo.** Norma conservada (NC), ultrasuave (US) o PAW. Los NC son
  caros pero los únicos que sirven para `epsilon.x`; los US son baratos y
  necesitan un `ecutrho` mucho mayor; PAW da mejores densidades.
- **Lo relativista.** Solo un pseudo 'full' sirve para espín-órbita.
- **Los electrones de valencia.** Más valencia (semicore) es más caro y
  más transferible; para metales de transición y DFT+U suele hacer falta.
- **El cutoff sugerido.** El que declara el propio archivo. Un pseudo que
  pide 90 Ry cuesta el doble que uno que pide 45.

CÓMO PUNTÚA
-----------
No hay un "mejor pseudopotencial": hay uno adecuado para lo que vas a
hacer. El módulo aplica REQUISITOS DUROS (que descartan) y PREFERENCIAS
(que ordenan), las dos declaradas en una tabla que se lee, y explica por
qué quedó fuera cada uno. La última palabra es siempre del usuario.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path

from qekit.core import pseudo as ps
from qekit.core.errors import ErrorDeUso

#: Requisitos duros por tarea: descartan un pseudo, no lo penalizan.
#: Cada entrada es (nombre legible, comprobación, por qué).
TAREAS = {
    "optics": {
        "nombre": "optics with epsilon.x",
        "tipo": ("NC",),
        "razon_tipo": "epsilon.x ONLY works with norm-conserving pseudopotentials. With "
                      "ultrasoft or PAW it returns a spectrum without complaint, "
                      "and it is wrong.",
    },
    "soc": {
        "nombre": "spin-orbit",
        "relativista": ("full",),
        "razon_rel": "spin-orbit coupling needs a FULLY relativistic "
                     "pseudopotential. A scalar one has already averaged the "
                     "SOC and would give a zero splitting disguised as a "
                     "result.",
    },
    "xanes": {
        "nombre": "XANES with xspectra.x",
        "gipaw": True,
        "razon_gipaw": "xspectra.x reconstructs the all-electron wave "
                       "function with GIPAW: the pseudopotential must carry "
                       "that information.",
    },
    "hubbard": {
        "nombre": "DFT+U",
        "prefiere_semicore": True,
    },
    "fonones": {
        "nombre": "phonons (DFPT)",
        "prefiere_tipo": ("NC", "US"),
        "razon_pref": "DFPT with PAW is more fragile and more expensive; with "
                      "norm-conserving or ultrasoft it works better.",
    },
    "general": {"nombre": "general calculation"},
}



#: Por debajo de este Z el espin-orbita es despreciable en la practica.
Z_SOC = 19


def _es_ligero(elemento: str) -> bool:
    from qekit.core.atomconf import Z_DE
    return Z_DE.get((elemento or "").capitalize(), 999) < Z_SOC

@dataclass
class Candidato:
    ruta: str = ""
    nombre: str = ""
    elemento: str = ""
    tipo: str = None
    funcional: str = None
    relativista: str = None
    z_valence: float = None
    ecutwfc: float = None
    ecutrho: float = None
    gipaw: bool = False
    tamano_kb: int = 0
    descartado: str = ""          # motivo, si lo hay
    puntos: float = 0.0
    notas: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.descartado


_RE_FUNC = re.compile(r'functional\s*=\s*"([^"]*)"', re.I)
_RE_FUNC_V1 = re.compile(r"^\s*(.+?)\s+Exchange-Correlation functional",
                         re.M)


#: Nombre corto <- las cuatro piezas con que QE escribe el mismo funcional.
NOMBRE_CORTO = {
    "SLA PW PBX PBC": "PBE",
    "SLA PZ NOGX NOGC": "PZ",
    "SLA PW PSX PSC": "PBESOL",
    "SLA B88 LYP BLYP": "BLYP",
    "SLA PW RPB PBC": "REVPBE",
    "SLA PW PBE": "PBE",
}


def _funcional(head: str):
    """El funcional, normalizado al nombre corto cuando se reconoce.

    QE escribe unas veces 'PBE' y otras 'SLA PW PBX PBC': son lo mismo.
    Mostrar las dos formas en la misma tabla hace pensar que dos pseudos
    son incompatibles cuando no lo son.
    """
    m = _RE_FUNC.search(head) or _RE_FUNC_V1.search(head)
    if not m:
        return None
    piezas = m.group(1).upper().split()
    if not piezas:
        return None
    crudo = " ".join(piezas)
    if crudo in NOMBRE_CORTO:
        return NOMBRE_CORTO[crudo]
    # UPF v1 escribe "SLA PW PBX PBC   PBE  Exchange-Correlation functional":
    # las cuatro piezas Y el nombre corto al final. Si la ultima pieza ya es
    # un nombre reconocible, esa es la respuesta.
    if piezas[-1] in NOMBRE_CORTO.values():
        return piezas[-1]
    if len(piezas) > 4 and " ".join(piezas[:4]) in NOMBRE_CORTO:
        return NOMBRE_CORTO[" ".join(piezas[:4])]
    return crudo


def leer(ruta) -> Candidato:
    """Lo que se puede saber de un UPF sin cargarlo entero."""
    p = Path(ruta)
    head = p.read_text(errors="ignore")[:30000]
    c = Candidato(ruta=str(p), nombre=p.name,
                  tipo=ps.pseudo_type(p), funcional=_funcional(head),
                  relativista=ps.relativistic(p), z_valence=ps.z_valence(p),
                  gipaw="PP_GIPAW" in p.read_text(errors="ignore"),
                  tamano_kb=int(p.stat().st_size / 1024))
    c.ecutwfc, c.ecutrho = ps.suggested_cutoffs(p)
    m = re.search(r'element\s*=\s*"?\s*(\w+)', head, re.I) or \
        re.search(r"^\s*(\w+)\s+Element", head, re.M)
    c.elemento = (m.group(1).strip() if m else p.name.split(".")[0])
    return c


def candidatos(elemento: str, pseudo_dir: str) -> list:
    """Todos los UPF de un elemento en la carpeta, leídos."""
    pdir = Path(pseudo_dir).expanduser()
    if not pdir.is_dir():
        raise ErrorDeUso(
            f"the pseudopotential folder '{pseudo_dir}' does not exist.\n"
            "It is set with:  olla-dft config set pseudo_dir /path/to/your/pseudos")
    return [leer(f) for f in ps.find_for_element(elemento, pdir)]


# ----------------------------------------------------------------------
# Selección
# ----------------------------------------------------------------------
def evaluar(cands: list, tarea: str = "general",
            funcional: str = None, prefiere_ligero: bool = False) -> list:
    """Descarta los que no sirven y ordena el resto. No elige por ti."""
    if tarea not in TAREAS:
        raise ErrorDeUso(
            f"unknown task '{tarea}'. Options: "
            + ", ".join(sorted(TAREAS)))
    reglas = TAREAS[tarea]

    for c in cands:
        c.descartado, c.puntos, c.notas = "", 0.0, []

        tipos = reglas.get("tipo")
        if tipos and c.tipo and c.tipo not in tipos:
            c.descartado = (f"it is {c.tipo} and "
                            f"{'/'.join(tipos)} is required: " + reglas["razon_tipo"])
            continue
        rel = reglas.get("relativista")
        if rel and c.relativista and c.relativista not in rel:
            if _es_ligero(c.elemento):
                # En un elemento ligero el espin-orbita es despreciable y
                # no existen pseudos relativistas para casi ninguno.
                # Descartarlo dejaria sin opciones a un calculo que es
                # perfectamente valido.
                c.notas.append(
                    f"it is {c.relativista}, but in {c.elemento} (small Z) "
                    "spin-orbit is\n    negligible: it can be used "
                    "together with relativistic pseudopotentials of the\n    heavy "
                    "elements.")
                c.puntos -= 0.5
            else:
                c.descartado = (f"it is {c.relativista}: " + reglas["razon_rel"])
                continue
        if reglas.get("gipaw") and not c.gipaw:
            c.descartado = "no GIPAW data: " + reglas["razon_gipaw"]
            continue
        if funcional and c.funcional and \
                not _mismo_funcional(c.funcional, funcional):
            c.descartado = (f"its functional is {c.funcional} and the requested one is "
                            f"{funcional}: mixing functionals invalidates the "
                            "total energy.")
            continue

        # preferencias
        pref = reglas.get("prefiere_tipo")
        if pref and c.tipo in pref:
            c.puntos += 2.0
            c.notas.append(f"type {c.tipo}: " + reglas.get("razon_pref", ""))
        if reglas.get("prefiere_semicore") and c.z_valence:
            c.puntos += 0.15 * c.z_valence
            c.notas.append(f"{c.z_valence:g} valence electrons: more "
                           "semicore is more transferable for DFT+U")
        if c.ecutwfc:
            # un cutoff bajo es dinero; se premia, pero poco
            c.puntos += max(0.0, (90.0 - c.ecutwfc) / 30.0)
        else:
            c.notas.append("no suggested cutoff declared: it will have to be "
                           "converged blindly")
            c.puntos -= 0.5
        if prefiere_ligero and c.tipo in ("US", "PAW"):
            c.puntos += 1.0
            c.notas.append("ultrasoft/PAW: fewer plane waves, cheaper")
        if c.gipaw:
            c.puntos += 0.3
            c.notas.append("has GIPAW data: also usable for XANES and NMR")
        if c.relativista == "full":
            c.puntos += 0.2

    return sorted(cands, key=lambda c: (not c.ok, -c.puntos, c.nombre))


def elegir(elemento: str, pseudo_dir: str, tarea: str = "general",
           funcional: str = None, prefiere_ligero: bool = False) -> tuple:
    """(el mejor candidato, la lista entera evaluada)."""
    cands = candidatos(elemento, pseudo_dir)
    if not cands:
        raise ErrorDeUso(
            f"there is no pseudopotential for {elemento} in "
            f"'{pseudo_dir}'.\nThey can be downloaded from pseudo-dojo.org or from "
            "quantum-espresso.org/pseudopotentials.")
    ev = evaluar(cands, tarea, funcional, prefiere_ligero)
    buenos = [c for c in ev if c.ok]
    if not buenos:
        raise ErrorDeUso(
            f"there are {len(ev)} pseudopotential(s) for {elemento} but none "
            f"is suitable for {TAREAS[tarea]['nombre']}:\n" +
            "\n".join(f"  {c.nombre}: {c.descartado}" for c in ev))
    return buenos[0], ev


def _mismo_funcional(a: str, b: str) -> bool:
    """Compara funcionales tolerando las distintas formas de escribirlos.

    'PBE' y 'SLA PW PBX PBC' son el mismo funcional escrito de dos formas:
    la corta y la lista de los cuatro trozos. Compararlos como cadenas
    descartaría pseudos perfectamente válidos.
    """
    ALIAS = {
        "PBE": {"PBE", "SLA PW PBX PBC"},
        "PZ": {"PZ", "SLA PZ NOGX NOGC", "LDA"},
        "PBESOL": {"PBESOL", "SLA PW PSX PSC"},
        "BLYP": {"BLYP", "SLA B88 LYP BLYP"},
        "REVPBE": {"REVPBE", "SLA PW RPB PBC"},
    }
    na = " ".join(a.upper().split())
    nb = " ".join(b.upper().split())
    if na == nb:
        return True
    for grupo in ALIAS.values():
        if na in grupo and nb in grupo:
            return True
    return False


def coherencia(elegidos: dict) -> list:
    """Avisos si los pseudos elegidos para distintos elementos no casan."""
    avisos = []
    funcs = {c.funcional for c in elegidos.values() if c.funcional}
    if len(funcs) > 1:
        # puede que sean el mismo escrito distinto
        base = list(funcs)[0]
        if not all(_mismo_funcional(base, f) for f in funcs):
            avisos.append(
                "DIFFERENT FUNCTIONALS between elements: "
                + "; ".join(f"{k}={v.funcional}" for k, v in elegidos.items()
                            if v.funcional)
                + ".\nThe total energy of a structure with pseudopotentials of "
                  "different functionals means\nnothing. This gives no "
                  "error in QE: a perfectly plausible and\nwrong number "
                  "comes out.")
    tipos = {c.tipo for c in elegidos.values() if c.tipo}
    if "NC" in tipos and tipos & {"US", "PAW"}:
        avisos.append(
            "Norm-conserving and ultrasoft/PAW are mixed. QE allows it, "
            "but the\necutrho is dictated by the most demanding one: use the dual of the "
            "ultrasoft (8-12x) for\nall of them, not that of the NC (4x).")
    cut = [c.ecutwfc for c in elegidos.values() if c.ecutwfc]
    if cut and max(cut) / min(cut) > 2.5:
        peor = max(elegidos.values(), key=lambda c: c.ecutwfc or 0)
        avisos.append(
            f"The suggested cutoffs range from {min(cut):.0f} to {max(cut):.0f} "
            f"Ry. The largest one rules ({peor.nombre}),\nso that element "
            "decides the cost of the whole calculation. If there is another softer\n"
            "pseudopotential for that element, it can save a lot of time.")
    return avisos


# ----------------------------------------------------------------------
# Reporte
# ----------------------------------------------------------------------
def report(elemento: str, evaluados: list, tarea: str = "general") -> str:
    lines = [f"--- Pseudopotentials for {elemento} ---",
             f"For: {TAREAS[tarea]['nombre']}", ""]
    buenos = [c for c in evaluados if c.ok]
    malos = [c for c in evaluados if not c.ok]

    if buenos:
        lines.append(f"{'':2s} {'file':<42s} {'type':>5s} {'func':>12s} "
                     f"{'rel':>7s} {'zval':>5s} {'ecut':>6s} {'rho':>6s}")
        for i, c in enumerate(buenos):
            marca = "->" if i == 0 else "  "
            lines.append(
                f"{marca} {c.nombre[:42]:<42s} {(c.tipo or '?'):>5s} "
                f"{(c.funcional or '?')[:12]:>12s} "
                f"{(c.relativista or '?')[:7]:>7s} "
                f"{(f'{c.z_valence:g}' if c.z_valence else '?'):>5s} "
                f"{(f'{c.ecutwfc:.0f}' if c.ecutwfc else '-'):>6s} "
                f"{(f'{c.ecutrho:.0f}' if c.ecutrho else '-'):>6s}")
        lines += ["", f"Recommended: {buenos[0].nombre}"]
        for n in buenos[0].notas:
            if n:
                lines.append(f"  - {n}")

    if malos:
        lines += ["", "Discarded:"]
        for c in malos:
            lines.append(f"  {c.nombre}")
            lines.append(f"      {c.descartado}")

    lines += ["",
              "This is a recommendation, not a truth: the appropriate "
              "pseudopotential depends\non the system. A specific one is forced with "
              "--pseudo " + elemento + "=file.UPF, and whichever\none you "
              "choose, the cutoff must be converged with 'olla-dft converge'."]
    return "\n".join(lines)


def report_coherencia(elegidos: dict) -> str:
    avisos = coherencia(elegidos)
    if not avisos:
        return ("The chosen pseudopotentials are mutually consistent "
                "(same functional,\ncompatible types, cutoffs of the same "
                "order).")
    return "\n\n".join(avisos)
