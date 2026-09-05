# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Recomendación adaptativa a partir de una serie de convergencia."""

from __future__ import annotations

import json
import math
from pathlib import Path

from qekit.core.errors import ErrorDeUso


def read(path) -> list:
    """Lee ``CONVERGENCIA.dat`` sin asumir que todos los puntos terminaron."""
    rows = []
    for number, line in enumerate(Path(path).read_text(encoding="utf-8",
                                                         errors="replace").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) < 3:
            continue
        try:
            value, energy_ry, delta = map(float, tokens[:3])
        except ValueError:
            continue
        if not all(math.isfinite(v) for v in (value, energy_ry, delta)):
            continue
        rows.append({"line": number, "value": value, "energy_ry": energy_ry,
                     "delta_mev_atom": abs(delta)})
    if not rows:
        raise ErrorDeUso(f"'{path}' contains no numeric convergence rows.")
    return rows


def analyze(path, threshold=None) -> dict:
    rows = read(path)
    if threshold is None:
        threshold = 1.0
    if threshold <= 0:
        raise ErrorDeUso("the threshold must be positive.")
    index = None
    for i in range(len(rows)):
        if all(row["delta_mev_atom"] <= threshold for row in rows[i:]):
            index = i
            break
    values = [row["value"] for row in rows]
    if index is None:
        status = "extend"
        recommendation = _next_value(values)
        reason = "no point keeps the whole tail within the threshold"
    elif index == len(rows) - 1:
        status = "confirm"
        recommendation = _next_value(values)
        reason = "only the last point complies; one more point is needed to confirm"
    else:
        status = "ready"
        recommendation = rows[index]["value"]
        reason = "from this point on the whole tail stays within the threshold"
    return {"file": str(Path(path).resolve()), "threshold": float(threshold),
            "rows": rows, "converged_index": index, "status": status,
            "recommended_value": recommendation, "reason": reason}


def _next_value(values):
    if len(values) < 2:
        return values[-1] * 1.25 if values[-1] > 0 else values[-1] + 1.0
    diffs = [b - a for a, b in zip(values, values[1:]) if b > a]
    if diffs:
        step = sorted(diffs)[len(diffs) // 2]
        return values[-1] + max(step, abs(values[-1]) * 0.10)
    return values[-1] * 1.25 if values[-1] > 0 else values[-1] + 1.0


def report(result: dict) -> str:
    lines = ["--- Adaptive convergence recommendation ---",
             f"File: {result['file']}",
             f"Threshold: {result['threshold']:g} meV/atom"]
    index = result["converged_index"]
    if index is None:
        lines.append("Status: EXTEND — the series does not converge yet.")
    elif result["status"] == "confirm":
        lines.append("Status: CONFIRM — the last point is not enough as evidence.")
    else:
        lines.append(f"Status: READY — use from point {index + 1} of the series.")
    lines.append(f"Recommendation: try value {result['recommended_value']:g}.")
    lines.append(f"Reason: {result['reason']}.")
    lines.append("The energy may converge before forces, phonons or tensors.")
    return "\n".join(lines)


def export(result: dict, destination="CONVERGENCIA_RECOMENDACION.json") -> Path:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target
