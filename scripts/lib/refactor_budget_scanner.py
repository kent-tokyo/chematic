"""Lexical Rust scanning primitives for the refactoring debt ratchet."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


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
GENERATED_MARKER_RE = re.compile(r"\b(?:auto[- ]?)?generated\b", re.IGNORECASE)


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
    index = 0
    block_depth = 0
    state = "code"
    raw_hashes = 0

    def blank(position: int) -> None:
        if out[position] != "\n":
            out[position] = " "

    while index < len(chars):
        if state == "line_comment":
            blank(index)
            if chars[index] == "\n":
                state = "code"
            index += 1
            continue

        if state == "block_comment":
            blank(index)
            if index + 1 < len(chars) and chars[index : index + 2] == ["/", "*"]:
                blank(index + 1)
                block_depth += 1
                index += 2
            elif index + 1 < len(chars) and chars[index : index + 2] == ["*", "/"]:
                blank(index + 1)
                block_depth -= 1
                index += 2
                if block_depth == 0:
                    state = "code"
            else:
                index += 1
            continue

        if state in {"string", "char"}:
            blank(index)
            delimiter = '"' if state == "string" else "'"
            if chars[index] == "\\" and index + 1 < len(chars):
                blank(index + 1)
                index += 2
            elif chars[index] == delimiter:
                state = "code"
                index += 1
            else:
                index += 1
            continue

        if state == "raw":
            blank(index)
            if chars[index] == '"':
                suffix = "#" * raw_hashes
                if text.startswith(suffix, index + 1):
                    for position in range(index + 1, index + 1 + raw_hashes):
                        blank(position)
                    index += 1 + raw_hashes
                    state = "code"
                    continue
            index += 1
            continue

        if index + 1 < len(chars) and chars[index : index + 2] == ["/", "/"]:
            blank(index)
            blank(index + 1)
            state = "line_comment"
            index += 2
            continue
        if index + 1 < len(chars) and chars[index : index + 2] == ["/", "*"]:
            blank(index)
            blank(index + 1)
            state = "block_comment"
            block_depth = 1
            index += 2
            continue
        if chars[index] in {"b", "r"}:
            match = re.match(r'(?:br|r)(#*)"', text[index:])
            if match:
                raw_hashes = len(match.group(1))
                for position in range(index, index + len(match.group(0))):
                    blank(position)
                index += len(match.group(0))
                state = "raw"
                continue
        if chars[index] == "b" and index + 1 < len(chars) and chars[index + 1] == "'":
            blank(index)
            blank(index + 1)
            state = "char"
            index += 2
            continue
        if chars[index] == '"':
            blank(index)
            state = "string"
            index += 1
            continue
        if chars[index] == "'":
            if (
                index + 1 < len(chars)
                and re.match(r"[A-Za-z_]", chars[index + 1])
                and (index + 2 >= len(chars) or chars[index + 2] != "'")
            ):
                index += 1
                continue
            blank(index)
            state = "char"
            index += 1
            continue
        index += 1

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
        if (
            "/tests/" in relative
            or relative.endswith("/tests.rs")
            or relative.endswith("_tests.rs")
        ):
            continue
        sources.append(path)
    return sources


def validate_excluded_paths(root: Path, configured: object) -> list[str]:
    """Require every budget exclusion to remain an explicit generated source."""

    if not isinstance(configured, list):
        return ["excluded_paths must be a list"]

    errors: list[str] = []
    paths = [path for path in configured if isinstance(path, str)]
    if len(paths) != len(configured):
        errors.append("excluded_paths entries must be strings")
    duplicates = sorted(path for path, count in Counter(paths).items() if count > 1)
    errors.extend(f"duplicate excluded path: {path}" for path in duplicates)

    for relative in paths:
        candidate = Path(relative)
        if (
            candidate.is_absolute()
            or ".." in candidate.parts
            or not relative.startswith("crates/")
            or "/src/" not in relative
            or candidate.suffix != ".rs"
        ):
            errors.append(f"excluded path is outside crate production Rust: {relative}")
            continue
        path = root / candidate
        if not path.is_file():
            errors.append(f"excluded generated file does not exist: {relative}")
            continue
        header = "\n".join(path.read_text(encoding="utf-8").splitlines()[:20])
        if not GENERATED_MARKER_RE.search(header):
            errors.append(f"excluded file lacks a generated header: {relative}")
    return errors


def collect(
    root: Path, config: dict[str, object]
) -> tuple[dict[str, int], dict[str, dict[str, int]]]:
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
