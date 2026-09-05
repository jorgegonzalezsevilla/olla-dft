# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Puerta de calidad científica para un Project Hub de Olla-DFT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from qekit.modules import project


@dataclass(frozen=True)
class Check:
    code: str
    title: str
    level: str  # ok, warn, fail
    detail: str
    evidence: str = ""


#: English display names of the verdict identifiers (the identifiers
#: themselves are kept as stored in the result dict).
VERDICT_EN = {"bloqueado": "BLOCKED", "revisar": "REVIEW", "listo": "READY"}


def evaluate(root: Path, data: dict) -> dict:
    checks = []
    state = project.status(root, data)
    if data.get("sources"):
        if state["changed_sources"]:
            checks.append(Check(
                "sources.changed", "Reproducible sources", "fail",
                "there are inputs that no longer match the recorded SHA-256",
                ", ".join(state["changed_sources"])))
        else:
            checks.append(Check(
                "sources.locked", "Reproducible sources", "ok",
                "the recorded sources keep their size and SHA-256"))
    else:
        checks.append(Check("sources.missing", "Reproducible sources", "fail",
                            "the project does not record any source"))

    # El bloqueo no sustituye a un contenedor ni a una receta de instalación,
    # pero sí hace visible si el equipo que reabre el proyecto cambió de
    # dependencias, Python o binarios de Quantum ESPRESSO.
    lock_path = root / project.PROJECT_DIR / "environment.lock.json"
    if lock_path.is_file():
        try:
            from qekit.modules import environment
            locked = environment.verify(root)
            checks.append(Check(
                "environment.locked", "Reproducible environment",
                "ok" if locked["ok"] else "warn",
                "Python, dependencies and binaries match" if locked["ok"] else
                "the current environment differs from the saved lock",
                ", ".join(locked.get("changed", []))))
        except Exception as exc:  # noqa: BLE001
            checks.append(Check("environment.unreadable", "Reproducible environment",
                                "warn", "the lock could not be verified",
                                str(exc)))
    elif data.get("sources"):
        checks.append(Check(
            "environment.missing", "Reproducible environment", "warn",
            "environment.lock.json is missing; create it before sharing or publishing",
            "olla-dft project environment"))

    advanced = data.get("metadata", {}).get("advanced_validation")
    if advanced:
        advanced_level = ("fail" if not advanced.get("passed") else
                          "warn" if advanced.get("warnings") else "ok")
        checks.append(Check(
            "validation.advanced", "Advanced validation",
            advanced_level,
            "structure, commands, units and outputs reviewed" if advanced_level == "ok"
            else "the advanced validation found failures" if advanced_level == "fail"
            else "advanced validation completed with warnings",
            str(advanced.get("at", ""))))

    # La presencia del índice no es una prueba de validez física, pero sí
    # evita que un proyecto con tareas terminadas pierda silenciosamente sus
    # salidas normalizadas antes de publicarse.
    try:
        from qekit.modules import results
        indexed = results.summary(results.project_db(root))
    except Exception as exc:  # noqa: BLE001
        checks.append(Check(
            "results.unreadable", "Traceable results", "fail",
            "the results index exists but cannot be read",
            str(exc).splitlines()[0]))
        indexed = {"count": 0, "by_status": {}}
    if indexed.get("count"):
        invalid = indexed.get("by_status", {}).get("invalid", 0)
        level = "warn" if invalid else "ok"
        detail = (f"{indexed['count']} result(s) keep metrics and hashes"
                  + (f"; {invalid} could not be interpreted" if invalid else ""))
        checks.append(Check("results.indexed", "Traceable results", level, detail,
                            str(results.project_db(root))))
    elif any(x.get("status") == "succeeded" for x in data.get("tasks", [])):
        checks.append(Check(
            "results.missing", "Traceable results", "warn",
            "there are finished tasks but no XML has been ingested yet"))

    counts = state["counts"]
    if counts.get("failed"):
        checks.append(Check(
            "tasks.failed", "Workflow executed", "fail",
            f"{counts['failed']} task(s) finished with an error"))
    elif counts.get("cancelled"):
        checks.append(Check(
            "tasks.cancelled", "Workflow executed", "warn",
            f"{counts['cancelled']} task(s) were left cancelled; resume and review before publishing"))
    elif data.get("tasks"):
        checks.append(Check(
            "tasks.state", "Workflow executed",
            "ok" if not counts.get("pending") and not counts.get("blocked")
            else "warn",
            "all tasks finished" if not counts.get("pending")
            and not counts.get("blocked") else
            "there are planned tasks that have not been executed yet"))
    else:
        checks.append(Check("tasks.empty", "Workflow executed", "warn",
                            "there are no tasks in the project yet"))

    if any(x.get("status") == "succeeded" for x in data.get("tasks", [])):
        checks.append(Check(
            "provenance.logs", "Logs and provenance", "ok",
            "the executed tasks have persistent state; export the "
            "snapshot to archive it"))
    else:
        checks.append(Check(
            "provenance.pending", "Logs and provenance", "warn",
            "there are no executed tasks to audit yet"))

    # No se presenta como aprobación de publicación: solo constata si el
    # proyecto ha pasado explícitamente la suite independiente.
    selftest = data.get("metadata", {}).get("selftest", {})
    if selftest.get("passed"):
        checks.append(Check("selftest.passed", "Independent formulas", "ok",
                            "the project records a passed selftest",
                            str(selftest.get("at"))))
    else:
        checks.append(Check(
            "selftest.missing", "Independent formulas", "warn",
            "no passed run of olla-dft selftest is recorded; this does not "
            "invalidate an exploration, but independent evidence is missing"))

    fails = sum(c.level == "fail" for c in checks)
    warns = sum(c.level == "warn" for c in checks)
    score = max(0, 100 - 40 * fails - 10 * warns)
    verdict = "bloqueado" if fails else "revisar" if warns else "listo"
    return {"checks": checks, "fails": fails, "warnings": warns,
            "score": score, "verdict": verdict}


def report(result: dict) -> str:
    lines = ["--- Scientific quality gate ---",
             f"Verdict: {VERDICT_EN.get(result['verdict'], result['verdict'].upper())}  |  indicative score: "
             f"{result['score']}/100"]
    for check in result["checks"]:
        mark = {"ok": "OK", "warn": "WARN", "fail": "FAIL"}[check.level]
        lines.append(f"  [{mark:5s}] {check.title}: {check.detail}")
        if check.evidence:
            lines.append(f"         evidence: {check.evidence}")
    lines.append("\nThis gate organizes evidence; it does not replace scientific review nor "
                 "authorize automatic publication.")
    return "\n".join(lines)
