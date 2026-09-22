#!/usr/bin/env python3
"""Execute one argv-safe command from the RDKit rebaseline catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "rdkit_rebaseline_execution.json"
PLACEHOLDER = re.compile(r"\$\{([A-Z][A-Z0-9_]*)\}")


def parse_values(items: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in items:
        key, separator, value = item.partition("=")
        if not separator or not key or not value or not re.fullmatch(r"[A-Z][A-Z0-9_]*", key):
            raise ValueError(f"invalid --set value {item!r}; expected NAME=value")
        result[key] = value
    return result


def expand(argv: list[str], values: dict[str, str]) -> list[str]:
    expanded: list[str] = []
    missing: set[str] = set()
    for token in argv:
        def replace(match: re.Match[str]) -> str:
            name = match.group(1)
            if name not in values:
                missing.add(name)
                return match.group(0)
            return values[name]

        expanded.append(PLACEHOLDER.sub(replace, token))
    if missing:
        raise ValueError(f"missing placeholders: {', '.join(sorted(missing))}")
    return expanded


def find_command(lane_id: str, command_id: str) -> tuple[dict, dict]:
    document = json.loads(MANIFEST.read_text(encoding="utf-8"))
    lane = next((item for item in document["lanes"] if item["id"] == lane_id), None)
    if lane is None:
        raise ValueError(f"unknown lane {lane_id!r}")
    command = next((item for item in lane["commands"] if item["id"] == command_id), None)
    if command is None:
        raise ValueError(f"unknown command {command_id!r} for lane {lane_id!r}")
    return lane, command


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lane", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--set", action="append", default=[])
    parser.add_argument("--execution-record", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        lane, command = find_command(args.lane, args.command)
        argv = expand(command["argv"], parse_values(args.set))
    except (OSError, json.JSONDecodeError, KeyError, ValueError) as exc:
        print(f"RDKit rebaseline command error: {exc}", file=sys.stderr)
        return 2
    if args.dry_run:
        print(json.dumps({"lane": lane["id"], "command": command["id"], "argv": argv}, indent=2))
        return 0
    started_at = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter_ns()
    process = subprocess.run(argv, cwd=ROOT, text=True, capture_output=True)
    elapsed_ns = time.perf_counter_ns() - started
    sys.stdout.write(process.stdout)
    sys.stderr.write(process.stderr)
    if args.execution_record is not None:
        record = {
            "schema_version": 1,
            "lane": lane["id"],
            "command": command["id"],
            "covers": command["covers"],
            "argv": argv,
            "started_at_utc": started_at,
            "elapsed_ns": elapsed_ns,
            "returncode": process.returncode,
            "stdout": {
                "bytes": len(process.stdout.encode()),
                "sha256": hashlib.sha256(process.stdout.encode()).hexdigest(),
                "tail": process.stdout[-4096:],
            },
            "stderr": {
                "bytes": len(process.stderr.encode()),
                "sha256": hashlib.sha256(process.stderr.encode()).hexdigest(),
                "tail": process.stderr[-4096:],
            },
        }
        args.execution_record.parent.mkdir(parents=True, exist_ok=True)
        args.execution_record.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return process.returncode


if __name__ == "__main__":
    raise SystemExit(main())
