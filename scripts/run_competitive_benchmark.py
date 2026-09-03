#!/usr/bin/env python3
"""Run the competitive benchmark protocol with resumable per-operation state."""

from __future__ import annotations

import argparse
import json
import platform
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = ROOT / "validation" / "competitive_benchmark_manifest.json"


def now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def atomic_write(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def load_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read {path}: {exc}") from exc


def installed_chematic_version() -> str:
    probe = subprocess.run(
        [sys.executable, "-c", "import chematic; print(getattr(chematic, '__version__', 'unknown'))"],
        cwd=ROOT, text=True, capture_output=True, check=False,
    )
    if probe.returncode != 0:
        raise SystemExit("chematic is not importable in the benchmark Python environment")
    return probe.stdout.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--result-dir", type=Path,
                    default=ROOT / "validation" / "results" / "competitive-benchmark")
    ap.add_argument("--state", type=Path, default=None,
                    help="State JSON; defaults to <result-dir>/state.json")
    ap.add_argument("--resume", action="store_true",
                    help="Resume from an existing state and skip measured operations")
    ap.add_argument("--operation", action="append", dest="operations",
                    help="Run only this operation id; may be repeated")
    ap.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    args = ap.parse_args()

    manifest = load_json(args.manifest)
    expected_version = manifest["target_version"]
    actual_version = installed_chematic_version()
    if actual_version != expected_version:
        raise SystemExit(
            f"chematic version mismatch: expected {expected_version}, found {actual_version}; "
            "build/install the current workspace package before measuring"
        )
    result_dir = args.result_dir if args.result_dir.is_absolute() else ROOT / args.result_dir
    result_dir.mkdir(parents=True, exist_ok=True)
    state_path = args.state or result_dir / "state.json"
    if not state_path.is_absolute():
        state_path = ROOT / state_path

    state = load_json(state_path) if args.resume and state_path.exists() else {
        "schema_version": 1,
        "protocol": "competitive-benchmark",
        "target_version": manifest["target_version"],
        "environment": {"chematic_version": actual_version},
        "status": "not_started",
        "started_at_utc": now(),
        "host": {"os": platform.platform(), "machine": platform.machine(),
                 "python": platform.python_version()},
        "operations": {}
    }
    state["status"] = "dry_run" if args.dry_run else "in_progress"
    atomic_write(state_path, state)

    selected = set(args.operations or [item["id"] for item in manifest["operations"]])
    unknown = selected - {item["id"] for item in manifest["operations"]}
    if unknown:
        raise SystemExit(f"unknown operation(s): {', '.join(sorted(unknown))}")

    try:
        for item in manifest["operations"]:
            operation_id = item["id"]
            if operation_id not in selected:
                continue
            previous = state["operations"].get(operation_id, {})
            if args.resume and previous.get("status") == "measured":
                print(f"SKIP {operation_id}: already measured")
                continue

            command = [part.format(result_dir=str(result_dir), operation_id=operation_id)
                       for part in item["command"]]
            log_path = result_dir / f"{operation_id}.log"
            record = {"status": "running", "command": command,
                      "command_text": shlex.join(command), "log": str(log_path),
                      "started_at_utc": now()}
            state["operations"][operation_id] = record
            atomic_write(state_path, state)
            print(f"RUN {operation_id}: {record['command_text']}")
            if args.dry_run:
                record["status"] = "not_run"
                record["finished_at_utc"] = now()
                atomic_write(state_path, state)
                continue

            started = time.monotonic()
            try:
                completed = subprocess.run(command, cwd=ROOT, text=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                           timeout=3600, check=False)
                log_path.write_text(completed.stdout)
                record["exit_code"] = completed.returncode
                record["duration_seconds"] = round(time.monotonic() - started, 3)
                record["finished_at_utc"] = now()
                record["status"] = "measured" if completed.returncode == 0 else "failed"
            except (OSError, subprocess.TimeoutExpired, KeyboardInterrupt) as exc:
                record["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
                record["error"] = str(exc)
                record["finished_at_utc"] = now()
                atomic_write(state_path, state)
                raise
            atomic_write(state_path, state)
            if record["status"] != "measured":
                state["status"] = "failed"
                atomic_write(state_path, state)
                return 1
    except KeyboardInterrupt:
        state["status"] = "interrupted"
        atomic_write(state_path, state)
        print(f"Interrupted; resume with --resume --state {state_path}", file=sys.stderr)
        return 130

    state["status"] = "not_run" if args.dry_run else "complete"
    state["finished_at_utc"] = now()
    atomic_write(state_path, state)
    print(f"State written to {state_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
