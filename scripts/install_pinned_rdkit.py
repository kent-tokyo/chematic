#!/usr/bin/env python3
"""Install the RDKit build the comparator policy pins for this platform (#754).

Reads ``validation/rdkit-comparator-platform-policy-v1.json``, downloads the
platform's RDKit wheel from PyPI as a binary, checks its SHA-256 against the
policy and installs it. A platform the policy lists without a wheel exits
with status 3, so a workflow can skip a same-process lane there instead of
installing another RDKit version.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "validation/rdkit-comparator-platform-policy-v1.json"


def platform_key() -> str:
    system = {"Linux": "linux", "Darwin": "macos", "Windows": "windows"}[platform.system()]
    machine = platform.machine().lower()
    machine = {"amd64": "amd64" if system == "windows" else "x86_64", "x86_64": "x86_64",
               "arm64": "arm64" if system == "macos" else "aarch64", "aarch64": "aarch64"}[machine]
    if system == "windows" and machine == "x86_64":
        machine = "amd64"
    return f"{system}-{machine}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--policy", type=Path, default=POLICY)
    ap.add_argument("--platform", default=None, help="policy platform key (default: this host)")
    ap.add_argument("--check-only", action="store_true", help="print the wheel and exit")
    args = ap.parse_args()
    policy = json.loads(args.policy.read_text(encoding="utf-8"))
    key = args.platform or platform_key()
    entry = policy["platforms"].get(key)
    if entry is None:
        print(f"platform {key} is not in {args.policy.name}", file=sys.stderr)
        return 2
    if not entry["wheel"]:
        print(f"{key}: no pinned RDKit build ({entry.get('note', '')})")
        return 3
    if args.check_only:
        print(entry["wheel"])
        return 0
    comparator = policy["comparator"]
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run([sys.executable, "-m", "pip", "download", "--disable-pip-version-check",
                        "--only-binary=:all:", "--no-deps", "-d", tmp,
                        f"{comparator['package']}=={comparator['version']}"], check=True)
        wheel = Path(tmp) / entry["wheel"]
        if not wheel.is_file():
            print(f"pip did not download {entry['wheel']}: {sorted(p.name for p in Path(tmp).iterdir())}",
                  file=sys.stderr)
            return 1
        digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
        if digest != entry["sha256"]:
            print(f"{wheel.name}: sha256 {digest} != pinned {entry['sha256']}", file=sys.stderr)
            return 1
        # RDKit's own dependencies (numpy, Pillow) come from PyPI as usual.
        subprocess.run([sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
                        str(wheel)], check=True)
    from importlib import metadata

    print(f"{key}: installed {entry['wheel']} ({metadata.version('rdkit')})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
