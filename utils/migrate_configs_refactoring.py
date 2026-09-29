#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: migrate_configs_refactoring.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  migrate_configs_refactoring.py -- Migrate JSON configs to EfficiencyMC
#  naming conventions (Refactoring_Changelog.md)
# ============================================================================
"""Migrate JSON configs to match Refactoring_Changelog.md conventions."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CONFIG_DIRS = [REPO / "configs", REPO / "outputs"]


# ----------------------------------------------------------------------------
# flatten_cluster_to_efficiency
#   Copy of a legacy cluster_mc block with its nested diffusion and binning sub-blocks flattened into the top level, as efficiency_mc expects.
# ----------------------------------------------------------------------------
def flatten_cluster_to_efficiency(cluster: dict) -> dict:
    out = deepcopy(cluster)
    if "diffusion" in out:
        diff = out.pop("diffusion")
        if isinstance(diff, dict):
            out.update(diff)
    if "binning" in out:
        binning = out.pop("binning")
        if isinstance(binning, dict):
            out.update(binning)
    return out


# ----------------------------------------------------------------------------
# migrate_response
#   Rename legacy response keys in place (mode cluster_mc -> pattern, pattern_mc -> efficiency_mc, cluster_mc block -> efficiency_mc) and return the list of changes made.
# ----------------------------------------------------------------------------
def migrate_response(response: dict) -> list[str]:
    changes: list[str] = []

    if response.get("mode") == "cluster_mc":
        response["mode"] = "pattern"
        changes.append("mode: cluster_mc -> pattern")

    if "pattern_mc" in response and "efficiency_mc" not in response:
        response["efficiency_mc"] = response.pop("pattern_mc")
        changes.append("pattern_mc -> efficiency_mc")

    if "cluster_mc" in response:
        cmc = response["cluster_mc"]
        if "efficiency_mc" not in response:
            response["efficiency_mc"] = flatten_cluster_to_efficiency(cmc)
            changes.append("cluster_mc -> efficiency_mc (flattened)")
        del response["cluster_mc"]
        changes.append("removed cluster_mc block")

    if response.get("mode") == "pcd" and "pcd" in response:
        pcd = response["pcd"]
        for key, default in (("sigma_res_e", 0.21), ("Dqmin", 0.5), ("Dqmax", 0.5)):
            if key not in pcd:
                pcd[key] = default
                changes.append(f"pcd: added default {key}")

    return changes


# ----------------------------------------------------------------------------
# migrate_file
#   Migrate the response block of one config file and rewrite it if anything changed; returns the changes (empty for non-config files).
# ----------------------------------------------------------------------------
def migrate_file(path: Path) -> list[str]:
    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return []

    if not isinstance(data, dict) or "response" not in data:
        return []

    changes = migrate_response(data["response"])
    if not changes:
        return []

    path.write_text(json.dumps(data, indent=2) + "\n")
    return changes


# ----------------------------------------------------------------------------
# main
#   Migrate every JSON config in the configured directories and print which files were updated.
# ----------------------------------------------------------------------------
def main() -> None:
    updated = []
    for cfg_dir in CONFIG_DIRS:
        if not cfg_dir.is_dir():
            continue
        for path in sorted(cfg_dir.rglob("*.json")):
            changes = migrate_file(path)
            if changes:
                rel = path.relative_to(REPO)
                updated.append((str(rel), changes))

    print(f"Updated {len(updated)} JSON file(s):\n")
    for rel, changes in updated:
        print(f"  {rel}")
        for c in changes:
            print(f"    - {c}")
        print()


if __name__ == "__main__":
    main()
