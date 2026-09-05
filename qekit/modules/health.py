# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Diagnóstico de la instalación y del entorno de ejecución de Olla-DFT.

Este diagnóstico es deliberadamente independiente de ``doctor`` para una
salida de pw.x: una instalación puede estar sana aunque todavía no exista un
cálculo, y un cálculo puede fallar aunque Python esté correctamente instalado.
"""

from __future__ import annotations

import importlib.metadata
import json
import platform
import shutil
import sys
from pathlib import Path

from qekit import __version__
from qekit import config


DEPENDENCIES = ("numpy", "ase", "spglib", "seekpath", "matplotlib", "scipy")
QE_BINARIES = ("pw.x", "ph.x", "q2r.x", "matdyn.x", "bands.x", "dos.x",
               "projwfc.x", "epsilon.x", "pp.x")


def _item(code, title, level, detail, evidence=""):
    return {"code": code, "title": title, "level": level,
            "detail": detail, "evidence": evidence}


def _memory_available_gb():
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return float(line.split()[1]) / 1024 / 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def check(path=".", project_path=None) -> dict:
    """Comprueba instalación, recursos, QE y un proyecto opcional."""
    checks = [_item("qekit.version", "Olla-DFT", "ok", __version__,
                    str(Path(__file__).resolve()))]
    version_ok = sys.version_info >= (3, 9)
    checks.append(_item(
        "python.version", "Python", "ok" if version_ok else "fail",
        platform.python_version(), "Python >= 3.9 is required"))

    missing = []
    installed = {}
    for name in DEPENDENCIES:
        try:
            installed[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            missing.append(name)
    checks.append(_item(
        "python.dependencies", "Python dependencies",
        "fail" if missing else "ok",
        "missing: " + ", ".join(missing) if missing else
        "all declared dependencies are installed",
        json.dumps(installed, ensure_ascii=False, sort_keys=True)))

    found_qe = {name: shutil.which(name) for name in QE_BINARIES}
    found_qe = {name: value for name, value in found_qe.items() if value}
    if found_qe:
        detail = f"{len(found_qe)}/{len(QE_BINARIES)} binaries found"
        level = "ok"
    else:
        detail = "pw.x not found; Olla-DFT can still prepare and analyze files"
        level = "warn"
    checks.append(_item("qe.binaries", "Quantum ESPRESSO", level, detail,
                        ", ".join(f"{key}={value}" for key, value in found_qe.items())))

    mpi = shutil.which("mpirun") or shutil.which("mpiexec")
    checks.append(_item(
        "mpi.launcher", "MPI", "ok" if mpi else "warn",
        f"launcher found: {mpi}" if mpi else
        "no mpirun/mpiexec; a single process can be used",
        "optional for preparation and post-processing"))

    try:
        values = config.load()
        pseudo = Path(values["pseudo_dir"]).expanduser()
        checks.append(_item(
            "pseudos.directory", "Pseudopotentials", "ok" if pseudo.is_dir() else "warn",
            str(pseudo) if pseudo.is_dir() else
            f"the configured folder does not exist: {pseudo}",
            "configure with olla-dft config set pseudo_dir PATH"))
    except Exception as exc:  # noqa: BLE001
        checks.append(_item("pseudos.config", "Pseudopotentials", "warn",
                            "could not read the configuration", str(exc)))

    target = Path(path or ".").expanduser()
    if not target.exists():
        target = Path.cwd()
    try:
        free = shutil.disk_usage(target).free / 1024 ** 3
        level = "fail" if free < 0.5 else "warn" if free < 2.0 else "ok"
        checks.append(_item("resources.disk", "Available disk space", level,
                            f"{free:.2f} GiB free", str(target.resolve())))
    except OSError as exc:
        checks.append(_item("resources.disk", "Available disk space", "warn",
                            "could not query the disk", str(exc)))

    memory = _memory_available_gb()
    if memory is not None:
        level = "fail" if memory < 0.5 else "warn" if memory < 2.0 else "ok"
        checks.append(_item("resources.memory", "Available memory", level,
                            f"{memory:.2f} GiB available", "/proc/meminfo"))
    else:
        checks.append(_item("resources.memory", "Available memory", "warn",
                            "could not be measured on this platform"))

    if project_path:
        try:
            from qekit.modules import project, quality
            root, data = project.load(project_path)
            gate = quality.evaluate(root, data)
            level = "fail" if gate["fails"] else "warn" if gate["warnings"] else "ok"
            checks.append(_item(
                "project.quality", "Project", level,
                f"{data['name']}: {gate['verdict']} ({gate['score']}/100)",
                str(root)))
        except Exception as exc:  # noqa: BLE001
            checks.append(_item("project.load", "Project", "fail",
                                "could not open the project", str(exc)))

    fails = sum(item["level"] == "fail" for item in checks)
    warnings = sum(item["level"] == "warn" for item in checks)
    return {"qekit_version": __version__, "python": platform.python_version(),
            "checks": checks, "fails": fails, "warnings": warnings,
            "ok": not fails}


def report(result: dict) -> str:
    lines = ["--- Olla-DFT installation diagnostics ---",
             f"Status: {'READY' if result['ok'] else 'BLOCKED'} · "
             f"warnings={result['warnings']} failures={result['fails']}"]
    marks = {"ok": "OK", "warn": "WARN", "fail": "FAIL"}
    for item in result["checks"]:
        lines.append(f"  [{marks[item['level']]:5s}] {item['title']}: {item['detail']}")
        if item.get("evidence"):
            lines.append(f"         {item['evidence']}")
    lines.append("\nTo repair: olla-dft doctor --help")
    return "\n".join(lines)
