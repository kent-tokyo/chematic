"""The comparator platform policy (#754) and the workflow that applies it."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
POLICY = json.loads((ROOT / "validation/rdkit-comparator-platform-policy-v1.json").read_text())
WORKFLOW = yaml.safe_load((ROOT / ".github/workflows/published-wheel-chemistry-gates.yml").read_text())


def runners(job: str) -> set[str]:
    return {entry["os"] for entry in WORKFLOW["jobs"][job]["strategy"]["matrix"]["include"]}


def test_same_process_lane_runs_exactly_where_a_pinned_build_exists():
    pinned = {p["runner"] for p in POLICY["platforms"].values() if p["same_process_lanes"] and p["runner"]}
    assert runners("autoconf") == pinned
    absent = {p["runner"] for p in POLICY["platforms"].values() if not p["same_process_lanes"]}
    assert not runners("autoconf") & absent


def test_platforms_without_a_build_still_run_chematics_side():
    for key, entry in POLICY["platforms"].items():
        if not entry["wheel"]:
            assert entry["runner"] in runners("chematic"), key


def test_every_pinned_wheel_names_its_platform_and_a_digest():
    tags = {"linux-x86_64": "manylinux_2_28_x86_64", "linux-aarch64": "manylinux_2_28_aarch64",
            "macos-arm64": "macosx_11_0_arm64", "windows-amd64": "win_amd64"}
    for key, tag in tags.items():
        entry = POLICY["platforms"][key]
        assert entry["wheel"].endswith(f"-{POLICY['comparator']['python']}-{POLICY['comparator']['python']}-{tag}.whl")
        assert len(entry["sha256"]) == 64


def test_installer_refuses_a_platform_without_a_build():
    script = ROOT / "scripts/install_pinned_rdkit.py"
    result = subprocess.run([sys.executable, str(script), "--platform", "macos-x86_64"],
                            capture_output=True, text=True)
    assert result.returncode == 3
    result = subprocess.run([sys.executable, str(script), "--platform", "macos-arm64", "--check-only"],
                            capture_output=True, text=True)
    assert result.returncode == 0
    assert result.stdout.strip() == POLICY["platforms"]["macos-arm64"]["wheel"]
