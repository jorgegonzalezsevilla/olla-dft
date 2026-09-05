# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Registro local de incidencias: fallas, confusiones y errores colados.

QUÉ ES Y QUÉ NO ES
------------------
Esto NO manda nada a ningún lado. Olla-DFT no lleva telemetría: todo se
guarda en tu máquina, en la carpeta de configuración de tu sistema
(`olla-dft sistema` dice cuál es), y tú decides si alguna vez lo compartes.

Tampoco "aprende solo". Un programa no se arregla a sí mismo. Lo que hace
este módulo es cerrar el ciclo de la única forma que funciona de verdad:

1. cuando algo falla, captura AUTOMÁTICAMENTE todo lo necesario para
   reproducirlo — versión, comando exacto, traza del error, versiones de
   las dependencias, si Quantum ESPRESSO estaba disponible;
2. te deja registrar a mano lo que no revienta pero estorba: una salida
   confusa, una bandera que no hace lo que parece, un número sospechoso;
3. resume cuáles comandos fallan más, que es la señal de dónde está mal
   la interfaz;
4. lo empaqueta en UN archivo para entregárselo a quien vaya a arreglarlo
   (quien mantenga el programa, o tú mismo dentro de seis meses).

El paso 1 es el que más vale. Un reporte que dice "falló optics" no sirve
para nada; uno que trae el comando, la traza y las versiones se arregla en
minutos.

SOBRE LO QUE SE GUARDA
----------------------
Por defecto se registran rutas de archivos, no su contenido. Los archivos
solo se copian si los adjuntas explícitamente con --attach, porque pueden
ser estructuras que no quieras mover de sitio ni siquiera dentro de tu
propia máquina.
"""

import json
import platform
import shutil
import subprocess
import sys
import traceback
import uuid
from dataclasses import dataclass, field
from pathlib import Path

from qekit import __version__
from qekit import config as qcfg
from qekit.core import provenance

DIR = qcfg.CONFIG_DIR / "incidencias"
DEPENDENCIAS = ("numpy", "scipy", "ase", "spglib", "seekpath", "matplotlib",
                "torch", "mace-torch")


def _versiones() -> dict:
    import importlib.metadata as md
    d = {"python": sys.version.split()[0],
         "sistema": f"{platform.system()} {platform.release()}",
         "qekit": __version__}
    for p in DEPENDENCIAS:
        try:
            d[p] = md.version(p)
        except Exception:                              # noqa: BLE001
            d[p] = None
    return d


def _quantum_espresso() -> dict:
    """¿Está pw.x, y de qué versión? Muchos fallos son de versión de QE."""
    cfg = qcfg.load()
    exe = shutil.which(cfg.get("pw_cmd", "pw.x")) or shutil.which("pw.x")
    if not exe:
        return {"disponible": False}
    ver = None
    try:
        out = subprocess.run([exe], input="", capture_output=True,
                             text=True, timeout=15).stdout
        for linea in out.splitlines():
            if "Program PWSCF" in linea:
                ver = linea.split("starts")[0].replace("Program PWSCF", "")
                ver = ver.strip()
                break
    except Exception:                                  # noqa: BLE001
        pass
    return {"disponible": True, "ruta": exe, "version": ver}


@dataclass
class Incidencia:
    id: str = ""
    fecha: str = ""
    tipo: str = "manual"          # "manual" | "error" | "uso"
    comando: str = ""
    descripcion: str = ""
    traceback: str = ""
    excepcion: str = ""
    cwd: str = ""
    versiones: dict = field(default_factory=dict)
    qe: dict = field(default_factory=dict)
    adjuntos: list = field(default_factory=list)
    estado: str = "abierta"       # "abierta" | "cerrada"
    nota: str = ""


def _nueva(tipo: str) -> Incidencia:
    return Incidencia(
        id=uuid.uuid4().hex[:8],
        fecha=provenance.fields()["generado"],
        tipo=tipo,
        comando=provenance.command_line(),
        cwd=str(Path.cwd()),
        versiones=_versiones(),
        qe=_quantum_espresso(),
    )


def registrar(descripcion: str = "", exc: BaseException = None,
              adjuntos=None, dir_=None, tipo: str = None) -> Incidencia:
    """Crea y guarda una incidencia. Con `exc`, captura la traza.

    `tipo="uso"` es para los errores de uso: se guarda el mensaje y el
    comando, pero NO la traza —no aporta nada cuando el programa hizo lo
    correcto al rechazar la entrada— y quedan aparte en las estadísticas.
    """
    if tipo is None:
        tipo = "error" if exc is not None else "manual"
    inc = _nueva(tipo)
    inc.descripcion = descripcion or ""
    if exc is not None:
        inc.excepcion = f"{type(exc).__name__}: {exc}"
        if tipo != "uso":
            inc.traceback = "".join(traceback.format_exception(
                type(exc), exc, exc.__traceback__))

    carpeta = Path(dir_ or DIR) / inc.id
    carpeta.mkdir(parents=True, exist_ok=True)
    for a in (adjuntos or []):
        origen = Path(a)
        if not origen.is_file():
            continue
        try:
            shutil.copy2(origen, carpeta / origen.name)
            inc.adjuntos.append(origen.name)
        except OSError:
            continue

    (carpeta / "incidencia.json").write_text(
        json.dumps(inc.__dict__, ensure_ascii=False, indent=2) + "\n")
    return inc


def listar(dir_=None) -> list:
    base = Path(dir_ or DIR)
    if not base.is_dir():
        return []
    out = []
    for f in sorted(base.glob("*/incidencia.json")):
        try:
            out.append(Incidencia(**json.loads(f.read_text())))
        except Exception:                              # noqa: BLE001
            continue
    return sorted(out, key=lambda i: i.fecha, reverse=True)


def cerrar(ident: str, nota: str = "", dir_=None) -> bool:
    base = Path(dir_ or DIR)
    f = base / ident / "incidencia.json"
    if not f.exists():
        return False
    d = json.loads(f.read_text())
    d["estado"] = "cerrada"
    if nota:
        d["nota"] = nota
    f.write_text(json.dumps(d, ensure_ascii=False, indent=2) + "\n")
    return True


def estadisticas(dir_=None) -> dict:
    """Qué comandos fallan más: la señal de dónde está mal la interfaz."""
    incs = listar(dir_)
    por_comando, por_excepcion, uso_por_comando = {}, {}, {}
    for i in incs:
        partes = i.comando.split()
        sub = partes[1] if len(partes) > 1 else "(no subcommand)"
        if i.tipo == "uso":
            uso_por_comando[sub] = uso_por_comando.get(sub, 0) + 1
            continue
        por_comando[sub] = por_comando.get(sub, 0) + 1
        if i.excepcion:
            clave = i.excepcion.split(":")[0]
            por_excepcion[clave] = por_excepcion.get(clave, 0) + 1
    return {
        "total": len(incs),
        "abiertas": sum(1 for i in incs if i.estado == "abierta"),
        "errores": sum(1 for i in incs if i.tipo == "error"),
        "uso": sum(1 for i in incs if i.tipo == "uso"),
        "por_comando": dict(sorted(por_comando.items(),
                                   key=lambda t: -t[1])),
        "por_excepcion": dict(sorted(por_excepcion.items(),
                                     key=lambda t: -t[1])),
        "uso_por_comando": dict(sorted(uso_por_comando.items(),
                                       key=lambda t: -t[1])),
    }


def exportar(destino="incidencias_qekit.json", dir_=None,
             solo_abiertas: bool = False) -> str:
    """Empaqueta todo en UN archivo para entregárselo a quien lo arregle."""
    incs = listar(dir_)
    if solo_abiertas:
        incs = [i for i in incs if i.estado == "abierta"]
    doc = {
        "que_es": (
            "Olla-DFT incident log exported for review. "
            "Each entry carries the exact command, the error traceback if there "
            "was one, and the versions of Olla-DFT, Python, the dependencies and "
            "Quantum ESPRESSO. With that the failure can be reproduced without "
            "asking anything further."),
        "qekit_version": __version__,
        "generado": provenance.fields()["generado"],
        "estadisticas": estadisticas(dir_),
        "incidencias": [i.__dict__ for i in incs],
    }
    Path(destino).write_text(
        json.dumps(doc, ensure_ascii=False, indent=2) + "\n")
    return str(destino)


def report_lista(incs: list) -> str:
    if not incs:
        return ("No incidents registered.\n\n"
                "They are registered automatically when a command fails. To note "
                "something that does not crash\nbut gets in the way —a confusing output, "
                "a flag that does not do what it seems—:\n"
                "    olla-dft report \"what happened\"")
    lines = [f"--- Incidents ({len(incs)}) ---",
             f"{'id':>9s} {'date':>17s} {'type':>7s} {'state':>8s}  "
             "command"]
    for i in incs:
        cmd = i.comando or "(no command)"
        lines.append(f"{i.id:>9s} {i.fecha[:16]:>17s} {i.tipo:>7s} "
                     f"{i.estado:>8s}  {cmd[:60]}")
    lines += ["", "Details:  olla-dft report --show <id>",
              "Close:    olla-dft report --close <id>",
              "Export:   olla-dft report --export incidencias.json"]
    return "\n".join(lines)


def report_detalle(inc: Incidencia) -> str:
    lines = [f"--- Incident {inc.id} ({inc.estado}) ---",
             f"Date: {inc.fecha}   |   type: {inc.tipo}",
             f"Command: {inc.comando or '(no command)'}",
             f"Directory: {inc.cwd}"]
    if inc.descripcion:
        lines += ["", "Description:", f"  {inc.descripcion}"]
    if inc.excepcion:
        lines += ["", f"Exception: {inc.excepcion}"]
    v = inc.versiones
    lines += ["", "Environment:",
              f"  Olla-DFT {v.get('qekit')}  |  Python {v.get('python')}  |  "
              f"{v.get('sistema')}"]
    deps = "  ".join(f"{k}={v[k]}" for k in DEPENDENCIAS
                     if v.get(k) is not None)
    if deps:
        lines.append(f"  {deps}")
    if inc.qe.get("disponible"):
        lines.append(f"  Quantum ESPRESSO: {inc.qe.get('version') or '?'}")
    else:
        lines.append("  Quantum ESPRESSO: not found in PATH")
    if inc.adjuntos:
        lines.append(f"  Attachments: {', '.join(inc.adjuntos)}")
    if inc.traceback:
        lines += ["", "Traceback:"]
        lines += [f"  {l}" for l in inc.traceback.rstrip().splitlines()]
    if inc.nota:
        lines += ["", f"Closing note: {inc.nota}"]
    return "\n".join(lines)


def report_estadisticas(st: dict) -> str:
    if not st["total"]:
        return "No incidents registered."
    lines = ["--- Incident summary ---",
             f"Total: {st['total']}  |  open: {st['abiertas']}  |  "
             f"errors: {st['errores']}  |  usage: {st.get('uso', 0)}", ""]
    if st["por_comando"]:
        lines.append("By subcommand (where the interface fails most):")
        for k, n in st["por_comando"].items():
            lines.append(f"  {n:4d}  {k}")
    if st["por_excepcion"]:
        lines += ["", "By exception type:"]
        for k, n in st["por_excepcion"].items():
            lines.append(f"  {n:4d}  {k}")
    if st.get("uso_por_comando"):
        lines += ["", "USAGE errors by subcommand (the program warned "
                  "correctly; the interface is confusing):"]
        for k, n in st["uso_por_comando"].items():
            lines.append(f"  {n:4d}  {k}")
    lines += ["",
              "A subcommand that accumulates failures is not bad luck: it is a "
              "confusing interface\nor an unhandled case. That is what "
              "must be fixed first.",
              "USAGE ones are not program failures, but if one repeats "
              "often the flag\nis badly named or badly documented."]
    return "\n".join(lines)
