#!/usr/bin/env python3
"""Validate the checked-in bounded Indigo workflow evidence packet."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "benchmarks/2026-10-12-indigo-workflows-source.json"
EXPECTED_VERSIONS = {
    "chematic": "1.1.0",
    "indigo": "1.46.0",
    "rdkit": "2026.03.6",
}
EXPECTED_COUNTS = {
    "reaction": {"chematic": (15, 15), "indigo": (13, 15), "timed": 13},
    "standardization": {"chematic": (39, 39), "indigo": (38, 39), "timed": 38},
    "rgroup": {"chematic": (16, 16), "indigo": (16, 16), "timed": 16},
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"Indigo workflow evidence invalid: {message}")


def check_corpus(entry: dict[str, object]) -> None:
    path = ROOT / str(entry["path"])
    if not path.is_file():
        fail(f"missing corpus {path}")
    if entry.get("sha256") != sha256(path):
        fail(f"corpus hash drift for {path}")


def main() -> None:
    data = json.loads(REPORT.read_text())
    if data.get("schema_version") != 1:
        fail("schema_version must be 1")
    if data.get("engines") != EXPECTED_VERSIONS:
        fail(f"engine versions differ: {data.get('engines')!r}")
    if data.get("blocks") != 21:
        fail("expected 21 alternating blocks")
    if not re.fullmatch(r"[0-9a-f]{40}", str(data.get("source_commit", ""))):
        fail("source_commit is not a full Git SHA")

    corpora = data.get("corpora", {})
    check_corpus(corpora["reaction"])
    check_corpus(corpora["rgroup"])
    for entry in corpora["standardization"]:
        check_corpus(entry)

    workflows = data.get("workflows", {})
    if set(workflows) != set(EXPECTED_COUNTS):
        fail(f"unexpected workflow keys: {sorted(workflows)}")
    for name, expected in EXPECTED_COUNTS.items():
        workflow = workflows[name]
        for engine in ("chematic", "indigo"):
            actual = workflow["accuracy"][engine]
            if (actual["match"], actual["total"]) != expected[engine]:
                fail(f"{name}/{engine} counts changed: {actual}")
        timing = workflow["timing"]
        if timing.get("rows") != expected["timed"]:
            fail(f"{name} timed-row count changed")
        comparison = timing["comparison"]
        if comparison.get("blocks") != 21 or comparison.get("winning_blocks") != 21:
            fail(f"{name} no longer wins every recorded block")
        if comparison.get("ci95", [0])[0] <= 1.0:
            fail(f"{name} paired interval no longer favors CheMatic")
        for engine in ("chematic", "indigo"):
            samples = timing["samples_ns"][engine]
            if len(samples) != 21 or any(not isinstance(value, int) or value <= 0 for value in samples):
                fail(f"{name}/{engine} timing samples are invalid")

    print("Indigo workflow evidence: OK")


if __name__ == "__main__":
    main()
