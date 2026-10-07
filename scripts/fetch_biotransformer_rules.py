#!/usr/bin/env python3
"""Fetch the three pinned BioTransformer rule tables without a checkout.

A partial clone with no working tree, then ``git show <commit>:<path>`` for
each table, checked against the SHA-256 the gate driver pins. The earlier
workflow step looked the tables up as ``database_<table>.json`` while the
repository stores them as ``database/[ENVMICRO/]<table>.json``, so it found
nothing and failed (the Windows runs stopped there; macOS stopped earlier).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from published_wheel_chemistry_gates import RULE_FILES, sha256  # noqa: E402

REPO = "https://github.com/Wishartlab-openscience/Biotransformer.git"
COMMIT = "7149f7ec"


def git(*args: str, capture: bool = True) -> bytes:
    return subprocess.run(["git", *args], check=True, capture_output=capture).stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--clone-dir", type=Path, default=Path("biotransformer-src"))
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--commit", default=COMMIT)
    args = ap.parse_args()
    if not args.clone_dir.exists():
        git("clone", "--quiet", "--filter=blob:none", "--no-checkout", args.repo, str(args.clone_dir),
            capture=False)
    paths = git("-C", str(args.clone_dir), "ls-tree", "-r", "--name-only", args.commit).decode().splitlines()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, digest in RULE_FILES.items():
        # The tables are stored as database/[ENVMICRO/]<table>.json; the
        # driver names them by that path with "_" for "/".
        match = [p for p in paths if p == name.replace("_", "/") or p.endswith("/" + name)]
        if not match:
            raise SystemExit(f"{name} not found at {args.commit}")
        target = args.out_dir / name
        target.write_bytes(git("-C", str(args.clone_dir), "show", f"{args.commit}:{match[0]}"))
        if sha256(target) != digest:
            raise SystemExit(f"{match[0]} at {args.commit} does not hash to the pinned {digest}")
        print(f"{name}: {match[0]} ({digest[:12]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
