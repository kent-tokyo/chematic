#!/usr/bin/env python3
"""Enforce the checked-in refactoring debt ratchet.

The repository already contains large, chemistry-sensitive modules.  A hard
limit would either break CI immediately or encourage unsafe mechanical splits.
Instead this checker records the current over-budget items and enforces three
rules:

1. new files and functions must stay below the soft limits;
2. existing over-budget items may shrink, but may not grow past the baseline;
3. resolved baseline entries must be removed, so the allowance only decreases.

This is a lexical Rust scanner, not a parser.  It deliberately uses a small,
documented subset that is stable for this workspace: ordinary free functions
and impl methods with brace bodies.  Comments and string/character literals are
masked before matching, and inline ``#[cfg(test)] mod tests`` bodies are ignored
for function budgets.  File budgets still count inline tests because moving a
large test module out of production code is itself useful refactoring work.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from lib.refactor_budget_scanner import (
    collect,
    function_metrics,
    mask_rust,
    rust_sources,
    validate_excluded_paths,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / "validation" / "refactor_budget_v1.json"

def generated_baseline(config: dict[str, object]) -> dict[str, object]:
    files, functions = collect(ROOT, config)
    limits = config["limits"]
    max_file = int(limits["max_file_lines"])
    max_function = int(limits["max_function_lines"])
    max_decisions = int(limits["max_function_decisions"])
    return {
        "schema_version": 1,
        "limits": limits,
        "excluded_paths": config.get("excluded_paths", []),
        "file_allowances": {
            path: lines for path, lines in files.items() if lines > max_file
        },
        "function_allowances": {
            key: {
                "lines": metric["lines"],
                "decisions": metric["decisions"],
            }
            for key, metric in functions.items()
            if metric["lines"] > max_function or metric["decisions"] > max_decisions
        },
    }


def check(config: dict[str, object]) -> list[str]:
    files, functions = collect(ROOT, config)
    limits = config["limits"]
    max_file = int(limits["max_file_lines"])
    max_function = int(limits["max_function_lines"])
    max_decisions = int(limits["max_function_decisions"])
    file_allowances = config.get("file_allowances", {})
    function_allowances = config.get("function_allowances", {})
    errors = validate_excluded_paths(ROOT, config.get("excluded_paths", []))

    for path, lines in files.items():
        allowed = file_allowances.get(path)
        if lines > max_file and allowed is None:
            errors.append(f"new oversized file: {path} has {lines} lines (limit {max_file})")
        elif allowed is not None and lines > int(allowed):
            errors.append(f"file grew: {path} has {lines} lines (baseline {allowed})")
    for path in sorted(set(file_allowances) - set(files)):
        errors.append(f"stale file allowance: {path}")
    for path, allowed in file_allowances.items():
        if path in files and files[path] <= max_file:
            errors.append(
                f"resolved file allowance must be removed: {path} now has {files[path]} lines"
            )

    for key, metric in functions.items():
        over = metric["lines"] > max_function or metric["decisions"] > max_decisions
        allowed = function_allowances.get(key)
        if over and allowed is None:
            errors.append(
                f"new oversized function: {key} at line {metric['start_line']} has "
                f"{metric['lines']} lines/{metric['decisions']} decisions "
                f"(limits {max_function}/{max_decisions})"
            )
        elif allowed is not None:
            if metric["lines"] > int(allowed["lines"]):
                errors.append(
                    f"function grew: {key} has {metric['lines']} lines "
                    f"(baseline {allowed['lines']})"
                )
            if metric["decisions"] > int(allowed["decisions"]):
                errors.append(
                    f"function complexity grew: {key} has {metric['decisions']} decisions "
                    f"(baseline {allowed['decisions']})"
                )
    for key in sorted(set(function_allowances) - set(functions)):
        errors.append(f"stale function allowance: {key}")
    for key, allowed in function_allowances.items():
        if key not in functions:
            continue
        metric = functions[key]
        if metric["lines"] <= max_function and metric["decisions"] <= max_decisions:
            errors.append(
                f"resolved function allowance must be removed: {key} now has "
                f"{metric['lines']} lines/{metric['decisions']} decisions"
            )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument(
        "--print-baseline",
        action="store_true",
        help="print a refreshed baseline to stdout without writing files",
    )
    args = parser.parse_args()

    try:
        config = json.loads(args.baseline.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"refactor budget check failed: {exc}", file=sys.stderr)
        return 2

    if config.get("schema_version") != 1:
        print("refactor budget check failed: unsupported schema_version", file=sys.stderr)
        return 2

    if args.print_baseline:
        errors = validate_excluded_paths(ROOT, config.get("excluded_paths", []))
        if errors:
            print("Refactor budget failures:", file=sys.stderr)
            print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
            return 1
        print(json.dumps(generated_baseline(config), indent=2, sort_keys=True))
        return 0

    errors = check(config)
    if errors:
        print("Refactor budget failures:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1

    files, functions = collect(ROOT, config)
    print(
        "Refactor budget OK: "
        f"{len(files)} production Rust files, {len(functions)} functions, "
        f"{len(config.get('file_allowances', {}))} file allowances, "
        f"{len(config.get('function_allowances', {}))} function allowances"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
