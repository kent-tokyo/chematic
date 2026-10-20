#!/usr/bin/env python3
"""Validate that a release tag names the checked-out, versioned main commit."""

from __future__ import annotations

import argparse
import re
import subprocess
import tomllib
from pathlib import Path

TAG = re.compile(r"v(?P<version>[0-9]+\.[0-9]+\.[0-9]+)(?:-[0-9A-Za-z][0-9A-Za-z.-]*)?\Z")


def version_from_tag(tag: str, workspace_version: str) -> str:
    match = TAG.fullmatch(tag)
    if not match:
        raise ValueError(f"invalid release tag: {tag}")
    version = match.group("version")
    if version != workspace_version:
        raise ValueError(f"tag {tag} does not match workspace version {workspace_version}")
    return version


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--main-ref", default="origin/main")
    args = parser.parse_args()
    manifest = tomllib.loads(Path("Cargo.toml").read_text())
    version = version_from_tag(args.tag, manifest["workspace"]["package"]["version"])
    tag_commit = subprocess.check_output(
        ["git", "rev-parse", f"refs/tags/{args.tag}^{{commit}}"], text=True
    ).strip()
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    if head != tag_commit:
        raise SystemExit(f"checkout {head} differs from tag {args.tag}: {tag_commit}")
    if subprocess.run([
        "git", "merge-base", "--is-ancestor", tag_commit, args.main_ref
    ]).returncode != 0:
        raise SystemExit(f"tag {args.tag} is not on {args.main_ref}")
    print(version)


if __name__ == "__main__":
    main()
