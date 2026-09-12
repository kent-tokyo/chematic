"""Dependency-free provenance records for validation reports."""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from pathlib import Path


def _git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=root, text=True, capture_output=True, check=True
    )
    return result.stdout


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _untracked_source_digest(root: Path) -> tuple[str, int]:
    paths = [
        Path(line)
        for line in _git(root, "ls-files", "--others", "--exclude-standard").splitlines()
        if line and Path(line).suffix in {".rs", ".py", ".js", ".mjs", ".ts", ".tsx"}
    ]
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(str(path).encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest(), len(paths)


def collect(root: Path) -> dict[str, object]:
    root = root.resolve()
    status = _git(root, "status", "--porcelain")
    diff = _git(root, "diff", "--binary").encode("utf-8")
    untracked_digest, untracked_count = _untracked_source_digest(root)
    lockfile = root / "Cargo.lock"
    return {
        "source_commit": _git(root, "rev-parse", "HEAD").strip(),
        "dirty_worktree": bool(status.strip()),
        "tracked_diff_sha256": _sha256(diff),
        "untracked_source_sha256": untracked_digest,
        "untracked_source_files": untracked_count,
        "lockfile_sha256": _sha256(lockfile.read_bytes()) if lockfile.is_file() else None,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
    }
