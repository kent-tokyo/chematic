#!/usr/bin/env python3
"""Execute sourced parser-security cases in isolated, bounded child processes.

The corpus is intentionally external to this runner: each case must retain its
source URL, license, affected version, cause and byte digest.  This prevents a
large but unattributed generated corpus from being presented as competitive
security evidence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORMATS = ("smiles", "smarts", "mol_v2000", "mol_v3000", "sdf")
REQUIRED = ("id", "format", "payload", "source_url", "license", "affected_version", "cause", "sha256", "expected_status")
DEFAULT_RUNNER = [str(ROOT / "target/release/examples/parser_security_case")]


def digest(payload: str) -> str:
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def preexec(memory_bytes: int) -> None:
    # macOS does not permit reducing RLIMIT_AS from a child preexec hook.
    # Do not pretend that an unenforced local setting is a memory gate: CI
    # uses Linux/RLIMIT_AS, while macOS requires --allow-unenforced-memory for
    # parser-functional diagnostics only.
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    target_soft = memory_bytes if hard in (-1, resource.RLIM_INFINITY) else min(memory_bytes, hard)
    resource.setrlimit(resource.RLIMIT_AS, (target_soft, hard))
    os.environ.pop("http_proxy", None)
    os.environ.pop("https_proxy", None)
    os.environ.pop("HTTP_PROXY", None)
    os.environ.pop("HTTPS_PROXY", None)


def isolated_command(runner: list[str], network_enforced: bool, fmt: str, input_path: Path) -> list[str]:
    """Return the actual child command without treating proxy cleanup as isolation.

    Linux CI creates a fresh network namespace under the GitHub runner's
    passwordless `sudo`. This avoids relying on unprivileged user namespaces,
    which GitHub-hosted kernels can disable. If that boundary cannot be made,
    the case remains a runner error instead of an untruthful pass. macOS
    diagnostics have no equivalent standard primitive here and remain
    explicitly unmeasured.
    """
    prefix = ["sudo", "--non-interactive", "unshare", "--net", "--"] if network_enforced else []
    return [*prefix, *runner, "--format", fmt, "--input", str(input_path)]


def load_cases(path: Path) -> list[dict[str, str]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    cases = document.get("cases")
    if not isinstance(cases, list):
        raise ValueError("corpus must contain a cases array")
    counts = {fmt: 0 for fmt in FORMATS}
    seen: set[str] = set()
    for case in cases:
        if not isinstance(case, dict) or any(not isinstance(case.get(key), str) or not case[key] for key in REQUIRED):
            raise ValueError(f"case missing required string field: {case!r}")
        if case["format"] not in counts:
            raise ValueError(f"unsupported format {case['format']}")
        if case["id"] in seen:
            raise ValueError(f"duplicate case id {case['id']}")
        if case["expected_status"] not in {"accepted", "rejected"}:
            raise ValueError(f"invalid expected_status for {case['id']}")
        if case["sha256"] != digest(case["payload"]):
            raise ValueError(f"payload digest mismatch for {case['id']}")
        seen.add(case["id"])
        counts[case["format"]] += 1
    if any(counts[fmt] < 20 for fmt in FORMATS):
        raise ValueError(f"need at least 20 sourced cases per format, got {counts}")
    if not all(any(case["format"] == fmt and case["expected_status"] == "accepted" for case in cases) for fmt in FORMATS):
        raise ValueError("each format requires an accepted valid-input control")
    return cases


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", type=Path)
    parser.add_argument("--runner", nargs="+", default=DEFAULT_RUNNER)
    parser.add_argument("--timeout-seconds", type=float, default=2.0)
    parser.add_argument("--memory-mib", type=int, default=256)
    parser.add_argument("--allow-unenforced-memory", action="store_true", help="run macOS parser diagnostics without claiming a memory gate")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.memory_mib <= 0:
        parser.error("timeout and memory budget must be positive")
    try:
        cases = load_cases(args.corpus.resolve())
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    memory_bytes = args.memory_mib * 1024 * 1024
    if not Path(args.runner[0]).is_file():
        parser.error(
            f"runner not found: {args.runner[0]}; build it before timing cases "
            "(cargo build --release --offline -p chematic-cli --example parser_security_case)"
        )
    memory_enforced = sys.platform == "linux"
    network_enforced = sys.platform == "linux"
    if not memory_enforced and not args.allow_unenforced_memory:
        parser.error("macOS cannot enforce RLIMIT_AS here; rerun on Linux or pass --allow-unenforced-memory for non-gate diagnostics")
    results: list[dict[str, object]] = []
    for case in cases:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", suffix=".input", delete=False) as handle:
            handle.write(case["payload"])
            input_path = Path(handle.name)
        started = time.monotonic()
        try:
            completed = subprocess.run(
                isolated_command(args.runner, network_enforced, case["format"], input_path),
                cwd=ROOT,
                capture_output=True,
                text=True,
                timeout=args.timeout_seconds,
                preexec_fn=(lambda: preexec(memory_bytes)) if memory_enforced else None,
            )
            elapsed_ms = round((time.monotonic() - started) * 1000, 3)
            try:
                output = json.loads(completed.stdout)
            except json.JSONDecodeError:
                output = None
            status = output.get("status") if isinstance(output, dict) else None
            passed = completed.returncode == 0 and status == case["expected_status"]
            peak_rss_kib = output.get("peak_rss_kib") if isinstance(output, dict) else None
            results.append({"id": case["id"], "format": case["format"], "status": status or "runner_error", "exit_code": completed.returncode, "signal": -completed.returncode if completed.returncode < 0 else None, "wall_ms": elapsed_ms, "peak_rss_kib": peak_rss_kib if isinstance(peak_rss_kib, int) else None, "passed": passed})
        except subprocess.TimeoutExpired:
            results.append({"id": case["id"], "format": case["format"], "status": "timeout", "exit_code": None, "signal": None, "wall_ms": round((time.monotonic() - started) * 1000, 3), "peak_rss_kib": None, "passed": False})
        finally:
            input_path.unlink(missing_ok=True)
    failures = [row for row in results if not row["passed"]]
    report = {"schema": "chematic.isolated-parser-security.v1", "execution": {"network": "Linux user+network namespace via unshare" if network_enforced else "not enforced locally; macOS diagnostic subprocess inherits host network policy", "network_enforced": network_enforced, "timeout_seconds": args.timeout_seconds, "memory_mib": args.memory_mib, "memory_enforced": memory_enforced, "peak_rss_kib": "Linux runner self-reports VmHWM; unavailable on macOS"}, "cases": results, "status": "pass" if not failures and memory_enforced and network_enforced else "not_measured" if not failures else "fail", "failure_count": len(failures)}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        output = args.json_out if args.json_out.is_absolute() else ROOT / args.json_out
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
