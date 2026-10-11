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
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASELINE = ROOT / "validation" / "refactor_budget_v1.json"

FUNCTION_RE = re.compile(
    r"(?m)^[ \t]*"
    r"(?:pub(?:\([^\n)]*\))?\s+)?"
    r"(?:(?:const|async|unsafe)\s+)*"
    r"(?:extern\s+\"[^\"]+\"\s+)?"
    r"fn\s+([A-Za-z_][A-Za-z0-9_]*)\b"
)
TEST_MODULE_RE = re.compile(
    r"(?ms)#\s*\[\s*cfg\s*\(\s*test\s*\)\s*\]\s*"
    r"(?:#\s*\[[^\]]*\]\s*)*mod\s+tests\s*\{"
)
DECISION_RE = re.compile(r"\b(?:if|match|for|while|loop)\b|&&|\|\||\?")


@dataclass(frozen=True)
class FunctionMetric:
    key: str
    path: str
    name: str
    occurrence: int
    start_line: int
    lines: int
    decisions: int


def mask_rust(text: str) -> str:
    """Mask comments and literals while preserving bytes and newlines."""

    chars = list(text)
    out = list(text)
    i = 0
    block_depth = 0
    state = "code"
    raw_hashes = 0

    def blank(index: int) -> None:
        if out[index] != "\n":
            out[index] = " "

    while i < len(chars):
        if state == "line_comment":
            blank(i)
            if chars[i] == "\n":
                state = "code"
            i += 1
            continue

        if state == "block_comment":
            blank(i)
            if i + 1 < len(chars) and chars[i] == "/" and chars[i + 1] == "*":
                blank(i + 1)
                block_depth += 1
                i += 2
            elif i + 1 < len(chars) and chars[i] == "*" and chars[i + 1] == "/":
                blank(i + 1)
                block_depth -= 1
                i += 2
                if block_depth == 0:
                    state = "code"
            else:
                i += 1
            continue

        if state in {"string", "char"}:
            blank(i)
            delimiter = '"' if state == "string" else "'"
            if chars[i] == "\\" and i + 1 < len(chars):
                blank(i + 1)
                i += 2
            elif chars[i] == delimiter:
                state = "code"
                i += 1
            else:
                i += 1
            continue

        if state == "raw":
            blank(i)
            if chars[i] == '"':
                suffix = "#" * raw_hashes
                if text.startswith(suffix, i + 1):
                    for j in range(i + 1, i + 1 + raw_hashes):
                        blank(j)
                    i += 1 + raw_hashes
                    state = "code"
                    continue
            i += 1
            continue

        if i + 1 < len(chars) and chars[i] == "/" and chars[i + 1] == "/":
            blank(i)
            blank(i + 1)
            state = "line_comment"
            i += 2
            continue
        if i + 1 < len(chars) and chars[i] == "/" and chars[i + 1] == "*":
            blank(i)
            blank(i + 1)
            state = "block_comment"
            block_depth = 1
            i += 2
            continue
        if chars[i] in {"b", "r"}:
            match = re.match(r"(?:br|r)(#*)\"", text[i:])
            if match:
                raw_hashes = len(match.group(1))
                for j in range(i, i + len(match.group(0))):
                    blank(j)
                i += len(match.group(0))
                state = "raw"
                continue
        if chars[i] == "b" and i + 1 < len(chars) and chars[i + 1] == "'":
            blank(i)
            blank(i + 1)
            state = "char"
            i += 2
            continue
        if chars[i] == '"':
            blank(i)
            state = "string"
            i += 1
            continue
        if chars[i] == "'":
            # Lifetimes (for example, 'a) are not character literals.
            if (
                i + 1 < len(chars)
                and re.match(r"[A-Za-z_]", chars[i + 1])
                and (i + 2 >= len(chars) or chars[i + 2] != "'")
            ):
                i += 1
                continue
            blank(i)
            state = "char"
            i += 1
            continue
        i += 1

    return "".join(out)


def matching_brace(masked: str, opening: int) -> int | None:
    depth = 0
    for index in range(opening, len(masked)):
        char = masked[index]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index
    return None


def test_module_ranges(masked: str) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    for match in TEST_MODULE_RE.finditer(masked):
        opening = masked.find("{", match.start(), match.end())
        closing = matching_brace(masked, opening)
        if closing is not None:
            ranges.append((match.start(), closing + 1))
    return ranges


def in_ranges(position: int, ranges: Iterable[tuple[int, int]]) -> bool:
    return any(start <= position < end for start, end in ranges)


def find_body_start(masked: str, signature_end: int) -> int | None:
    paren = bracket = angle = 0
    index = signature_end
    while index < len(masked):
        char = masked[index]
        if char == "(":
            paren += 1
        elif char == ")":
            paren = max(0, paren - 1)
        elif char == "[":
            bracket += 1
        elif char == "]":
            bracket = max(0, bracket - 1)
        elif char == "<":
            angle += 1
        elif char == ">":
            angle = max(0, angle - 1)
        elif char == ";" and paren == bracket == angle == 0:
            return None
        elif char == "{" and paren == bracket == 0:
            return index
        index += 1
    return None


def function_metrics(path: Path, relative: str) -> list[FunctionMetric]:
    text = path.read_text(encoding="utf-8")
    masked = mask_rust(text)
    ignored = test_module_ranges(masked)
    occurrences: Counter[str] = Counter()
    metrics: list[FunctionMetric] = []

    for match in FUNCTION_RE.finditer(masked):
        if in_ranges(match.start(), ignored):
            continue
        name = match.group(1)
        body_start = find_body_start(masked, match.end())
        if body_start is None:
            continue
        body_end = matching_brace(masked, body_start)
        if body_end is None:
            continue
        occurrences[name] += 1
        occurrence = occurrences[name]
        start_line = text.count("\n", 0, match.start()) + 1
        end_line = text.count("\n", 0, body_end) + 1
        body = masked[body_start : body_end + 1]
        metrics.append(
            FunctionMetric(
                key=f"{relative}::{name}#{occurrence}",
                path=relative,
                name=name,
                occurrence=occurrence,
                start_line=start_line,
                lines=end_line - start_line + 1,
                decisions=len(DECISION_RE.findall(body)),
            )
        )
    return metrics


def rust_sources(root: Path, excluded: set[str]) -> list[Path]:
    sources: list[Path] = []
    for path in sorted((root / "crates").glob("*/src/**/*.rs")):
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        if "/tests/" in relative or relative.endswith("/tests.rs"):
            continue
        sources.append(path)
    return sources


def collect(root: Path, config: dict[str, object]) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
    excluded = set(config.get("excluded_paths", []))
    file_lines: dict[str, int] = {}
    functions: dict[str, dict[str, int]] = {}
    for path in rust_sources(root, excluded):
        relative = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8")
        file_lines[relative] = len(text.splitlines())
        for metric in function_metrics(path, relative):
            functions[metric.key] = {
                "lines": metric.lines,
                "decisions": metric.decisions,
                "start_line": metric.start_line,
            }
    return file_lines, functions


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
    errors: list[str] = []

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
