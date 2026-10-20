"""Shared, side-effect-free I/O primitives for committed evidence gates."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


def read_json(path: Path) -> dict[str, Any]:
    """Read one UTF-8 JSON object."""
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_bytes(data: bytes) -> str:
    """Return the lowercase SHA-256 digest of ``data``."""
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    """Raise the stable gate error used by fail-fast evidence checkers."""
    if not condition:
        raise ValueError(message)


def collect_error(condition: bool, message: str, errors: list[str]) -> None:
    """Append an error for checkers that report every failure together."""
    if not condition:
        errors.append(message)


def read_gzip_jsonl(path: Path) -> tuple[bytes, list[dict[str, Any]]]:
    """Return decompressed bytes and decoded JSONL rows."""
    body = gzip.decompress(path.read_bytes())
    return body, [json.loads(line) for line in body.splitlines()]


def read_gzip_jsonl_with_sha(path: Path) -> tuple[list[dict[str, Any]], str]:
    """Return decoded rows and the digest of their uncompressed bytes."""
    body, rows = read_gzip_jsonl(path)
    return rows, sha256_bytes(body)
