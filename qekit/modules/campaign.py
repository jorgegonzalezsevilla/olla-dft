# Olla-DFT — command-line toolkit for Quantum ESPRESSO
# Copyright (C) 2026 Jorge Enrique González Sevilla
# SPDX-License-Identifier: AGPL-3.0-or-later
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. See the LICENSE file for details.

"""Campañas reproducibles de cálculos parametrizados.

Una campaña es una matriz explícita de parámetros que se convierte en tareas
independientes dentro del Project Hub. El usuario revisa los comandos y luego
decide si ejecuta el workflow; crear una campaña nunca lanza Quantum
ESPRESSO. Puede tomar el siguiente valor recomendado por ``olla-dft tune`` para
extender una serie de convergencia de forma controlada.
"""

from __future__ import annotations

import itertools
import json
import re
import shlex
from pathlib import Path

from qekit.core.errors import ErrorDeUso
from qekit.modules import project


_AXIS_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_-]{0,31}$")


def _value(text: str):
    text = str(text).strip()
    if not text:
        raise ErrorDeUso("a campaign axis cannot have empty values.")
    try:
        if re.fullmatch(r"[+-]?\d+", text):
            return int(text)
        if re.fullmatch(r"[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][+-]?\d+)?", text):
            return float(text)
    except ValueError:
        pass
    return text


def parse_axes(specs) -> dict:
    axes = {}
    for spec in specs or []:
        if "=" not in spec:
            raise ErrorDeUso(f"invalid axis '{spec}'; use name=v1,v2,v3.")
        name, raw_values = spec.split("=", 1)
        name = name.strip()
        if not _AXIS_NAME.fullmatch(name):
            raise ErrorDeUso(f"invalid axis name '{name}'.")
        values = [_value(item) for item in raw_values.split(",")]
        if not values or len(values) > 256:
            raise ErrorDeUso(f"axis '{name}' must have between 1 and 256 values.")
        if name in axes:
            raise ErrorDeUso(f"axis '{name}' appears more than once.")
        axes[name] = values
    if not axes:
        raise ErrorDeUso("campaign create needs at least one --axis name=v1,v2.")
    return axes


def _format(template: str, parameters: dict, index: int, total: int,
            structure: str) -> str:
    values = dict(parameters)
    values.update(index=index, id=f"point-{index:03d}", total=total,
                 structure=shlex.quote(structure))
    try:
        command = template.format(**values)
    except KeyError as exc:
        raise ErrorDeUso(f"the template uses the unknown field '{exc.args[0]}'.") from None
    except ValueError as exc:
        raise ErrorDeUso(f"invalid campaign template: {exc}") from None
    # La validación final la hace el mismo contrato que el Project Hub.
    project._valid_command(command)
    return command


def create(root, data, name, command, axis_specs, goal=None,
           convergence_file=None, adaptive=False) -> dict:
    if not name or not str(name).strip():
        raise ErrorDeUso("campaign create needs a name.")
    if not command or not str(command).strip():
        raise ErrorDeUso("campaign create needs --command.")
    axes = parse_axes(axis_specs)
    campaign_id = project._slug(name)
    if any(item.get("id") == campaign_id for item in data.get("campaigns", [])):
        raise ErrorDeUso(f"campaign '{campaign_id}' already exists.")
    recommendation = None
    if convergence_file and adaptive:
        from qekit.modules import tuning
        recommendation = tuning.analyze(convergence_file)["recommended_value"]
        candidate = next((key for key in axes
                          if key.lower() in ("ecutwfc", "ecutrho", "kmesh")), None)
        if candidate is None:
            raise ErrorDeUso("--adaptive needs an ecutwfc, ecutrho or kmesh axis.")
        if recommendation not in axes[candidate]:
            axes[candidate].append(recommendation)

    source = project._relative(root, project.source_path(root, data))
    keys = list(axes)
    combinations = [dict(zip(keys, values))
                    for values in itertools.product(*(axes[key] for key in keys))]
    if len(combinations) > 1024:
        raise ErrorDeUso("the campaign would produce more than 1024 points; reduce the axes.")
    commands = [_format(command, params, index, len(combinations), source)
                for index, params in enumerate(combinations, 1)]
    tasks = project.plan(root, data, goal or f"campaign:{campaign_id}", commands,
                         task_prefix=campaign_id)
    ids = [task["id"] for task in tasks]
    # Los puntos son independientes: una campaña no debe bloquearse porque
    # falló otro punto de la matriz.
    by_id = {task["id"]: task for task in data["tasks"]}
    for index, (task_id, params) in enumerate(zip(ids, combinations), 1):
        task = by_id[task_id]
        task["depends_on"] = []
        task["campaign_id"] = campaign_id
        task["parameters"] = params
        task["label"] = f"{name} · point {index}/{len(combinations)}"
        # Si la plantilla declara un destino estándar, registrarlo para que
        # Project Hub pueda calcular hashes e ingerir XML al terminar.
        tokens = project._valid_command(task["command"])
        for flag in ("--outdir", "-o"):
            if flag in tokens:
                position = tokens.index(flag)
                if position + 1 < len(tokens):
                    task["outputs"] = [tokens[position + 1]]
                break
    record = {
        "id": campaign_id, "name": str(name), "goal": goal or "",
        "command_template": str(command), "axes": axes,
        "tasks": ids, "points": len(ids), "created": project._now(),
        "adaptive": bool(adaptive),
        "output_isolation": ("explicit" if len(ids) == 1 or
                              "{id}" in str(command) or "{index}" in str(command)
                              else "review_required"),
        "adaptive_recommendation": recommendation,
        "convergence_file": (str(Path(convergence_file).resolve())
                              if convergence_file else None),
    }
    data.setdefault("campaigns", []).append(record)
    data.setdefault("metadata", {})["last_campaign"] = campaign_id
    return record


def get(data: dict, campaign_id: str) -> dict:
    wanted = str(campaign_id).strip()
    for item in data.get("campaigns", []):
        if item.get("id") == wanted:
            return item
    raise ErrorDeUso(f"campaign '{campaign_id}' not found.")


def run(root, data, campaign_id, execute=False, force=False, parallel=1,
        retries=0, timeout=None, cancel_file=None) -> list:
    """Corre solo los puntos de una campaña, nunca el resto del proyecto."""
    item = get(data, campaign_id)
    return project.run(root, data, execute=execute,
                       task_ids=item.get("tasks", []), force=force,
                       parallel=parallel, retries=retries, timeout=timeout,
                       cancel_file=cancel_file)


def extend(root, data, campaign_id, convergence_file, threshold=None) -> dict:
    """Añade la recomendación de ``tune`` como nuevos puntos de la campaña."""
    from qekit.modules import tuning

    item = get(data, campaign_id)
    analysis = tuning.analyze(convergence_file, threshold=threshold)
    recommendation = analysis["recommended_value"]
    candidate = next((key for key in item.get("axes", {})
                      if key.lower() in ("ecutwfc", "ecutrho", "kmesh")), None)
    if candidate is None:
        raise ErrorDeUso("the campaign has no ecutwfc, ecutrho or kmesh axis.")
    values = item["axes"][candidate]
    if any(abs(float(recommendation) - float(value)) < 1e-12 for value in values):
        return {"extended": False, "reason": "the recommended value is already in the campaign",
                "recommended_value": recommendation, "campaign": item}
    values.append(recommendation)
    keys = list(item["axes"])
    combinations = [dict(zip(keys, values_))
                    for values_ in itertools.product(*(item["axes"][key] for key in keys))]
    existing = {tuple(sorted((task.get("parameters") or {}).items()))
                for task in data.get("tasks", [])
                if task.get("campaign_id") == campaign_id}
    new_combinations = [params for params in combinations
                        if tuple(sorted(params.items())) not in existing]
    offset = len(item.get("tasks", []))
    total = offset + len(new_combinations)
    structure = project._relative(root, project.source_path(root, data))
    prefix = f"{campaign_id}-ext{offset}"
    commands = [_format(item["command_template"], params, offset + index,
                        total, structure)
                for index, params in enumerate(new_combinations, 1)]
    tasks = project.plan(root, data, item.get("goal") or f"campaign:{campaign_id}",
                         commands, task_prefix=prefix) if commands else []
    by_id = {task["id"]: task for task in data["tasks"]}
    for index, (task, params) in enumerate(zip(tasks, new_combinations), offset + 1):
        record = by_id[task["id"]]
        record["depends_on"] = []
        record["campaign_id"] = campaign_id
        record["parameters"] = params
        record["label"] = f"{item['name']} · point {index}/{total}"
        tokens = project._valid_command(record["command"])
        for flag in ("--outdir", "-o"):
            if flag in tokens and tokens.index(flag) + 1 < len(tokens):
                record["outputs"] = [tokens[tokens.index(flag) + 1]]
                break
    item["tasks"].extend(task["id"] for task in tasks)
    item["points"] = len(item["tasks"])
    item["adaptive_recommendation"] = recommendation
    item.setdefault("extensions", []).append({
        "at": project._now(), "file": str(Path(convergence_file).resolve()),
        "value": recommendation, "points_added": len(tasks),
    })
    return {"extended": bool(tasks), "recommended_value": recommendation,
            "points_added": len(tasks), "campaign": item}


def status(data: dict, campaign_id: str) -> dict:
    item = get(data, campaign_id)
    tasks = {task.get("id"): task for task in data.get("tasks", [])}
    counts = {key: 0 for key in ("pending", "running", "succeeded", "failed", "blocked")}
    for task_id in item.get("tasks", []):
        state = tasks.get(task_id, {}).get("status", "missing")
        counts[state] = counts.get(state, 0) + 1
    return {"campaign": item, "counts": counts,
            "missing_tasks": [x for x in item.get("tasks", []) if x not in tasks]}


def report(data: dict, campaign_id=None) -> str:
    campaigns = data.get("campaigns", [])
    if campaign_id:
        items = [get(data, campaign_id)]
    else:
        items = campaigns
    if not items:
        return "--- Campaigns ---\nNo campaigns registered."
    lines = ["--- Reproducible campaigns ---"]
    tasks = {task.get("id"): task for task in data.get("tasks", [])}
    for item in items:
        counts = {}
        for task_id in item.get("tasks", []):
            state = tasks.get(task_id, {}).get("status", "missing")
            counts[state] = counts.get(state, 0) + 1
        state = ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))
        lines.append(f"{item['id']}: {item['name']} · points={item['points']} · {state}")
        lines.append(f"  template: {item['command_template']}")
        if item.get("output_isolation") == "review_required":
            lines.append("  WARNING: the template does not include {id} or {index}; check for possible output collisions.")
        if item.get("adaptive_recommendation") is not None:
            lines.append(f"  adaptive: last suggested value {item['adaptive_recommendation']:g}")
        lines.append("  axes: " + ", ".join(
            f"{key}={','.join(map(str, values))}" for key, values in item["axes"].items()))
    return "\n".join(lines)


def export(data: dict, destination, campaign_id=None) -> Path:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {"campaigns": ([get(data, campaign_id)] if campaign_id
                              else data.get("campaigns", [])),
               "tasks": data.get("tasks", [])}
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return target
