#!/usr/bin/env python3
"""Validate a checked-in, version-specific release-channel verification record."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_CHANNELS = {"github_release", "npm", "pypi", "crates_io", "docs_rs", "site"}
STATUSES = {"verified", "mismatch", "unreachable", "not_measured"}


def main() -> int:
    version = re.search(r'^version\s*=\s*"([^"]+)"\s*$', (ROOT / "Cargo.toml").read_text(encoding="utf-8"), re.MULTILINE)
    if version is None:
        print("release channel evidence invalid: workspace version missing", file=sys.stderr)
        return 1
    path = ROOT / "validation" / "results" / f"release-channel-verification-v{version.group(1)}.json"
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"release channel evidence invalid: {exc}", file=sys.stderr)
        return 1
    errors: list[str] = []
    release = document.get("release", {})
    if document.get("schema_version") != 1 or release.get("version") != version.group(1) or release.get("tag") != f"v{version.group(1)}":
        errors.append("release identity does not match workspace version")
    channels = document.get("channels")
    if not isinstance(channels, list):
        errors.append("channels must be an array")
        channels = []
    ids = [entry.get("id") for entry in channels if isinstance(entry, dict)]
    if set(ids) != EXPECTED_CHANNELS or len(ids) != len(EXPECTED_CHANNELS):
        errors.append("channels must contain each required channel exactly once")
    verified = True
    for entry in channels:
        if not isinstance(entry, dict) or entry.get("status") not in STATUSES:
            errors.append("channel has invalid status")
            continue
        if not isinstance(entry.get("url"), str) or not entry["url"].startswith("https://"):
            errors.append(f"{entry.get('id')}: URL is required")
        if not isinstance(entry.get("detail"), str) or not entry["detail"]:
            errors.append(f"{entry.get('id')}: detail is required")
        if entry["status"] != "verified":
            verified = False
        elif entry.get("observed_version") != version.group(1):
            errors.append(f"{entry.get('id')}: verified channel has wrong observed version")
    if document.get("release_ready") is not verified:
        errors.append("release_ready must exactly reflect all channels being verified")
    if errors:
        print("release channel evidence invalid:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    status = "READY" if verified else "PARTIAL"
    print(f"release channel evidence {status}: {len(channels)} channels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
