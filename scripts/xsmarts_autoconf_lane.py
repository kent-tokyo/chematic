#!/usr/bin/env python3
"""Pinned xsmarts-autoconf lane (#754): chematic's SMARTS/SMIRKS behaviour flags.

Runs ``xsmarts-autoconf update chematic`` at the pinned commit against the
installed ``chematic`` and compares every flag with the tool's RDKit
2026.03.6 config. The flags that may differ, with their values and classes,
are listed in ``validation/xsmarts_autoconf_expected.json``; any other
difference fails, and so does an expected difference that disappeared (update
the expectations when a flag is aligned).

Install the tool at the pinned commit and RDKit 2026.03.6 (the tool uses RDKit
to generate atom reorderings) first, e.g.::

    pip install 'rdkit==2026.3.6' \
        "xsmarts-autoconf @ git+https://github.com/swamidasslab/xsmarts-autoconf@<commit>"
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def installed_commit() -> str | None:
    """The VCS commit pip recorded for the installed tool, if any."""
    dist = importlib.metadata.distribution("xsmarts-autoconf")
    text = dist.read_text("direct_url.json")
    if not text:
        return None
    info = json.loads(text)
    return info.get("vcs_info", {}).get("commit_id")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--expected", type=Path, default=ROOT / "validation/xsmarts_autoconf_expected.json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--allow-unpinned", action="store_true",
                        help="accept a local (editable) install whose commit pip did not record")
    args = parser.parse_args()
    expected = json.loads(args.expected.read_text(encoding="utf-8"))

    import xsmarts_autoconf

    commit = installed_commit()
    if commit != expected["autoconf_commit"] and not (args.allow_unpinned and commit is None):
        parser.error(f"xsmarts-autoconf commit {commit} != pinned {expected['autoconf_commit']}")
    shipped = Path(xsmarts_autoconf.__file__).resolve().parent / "data" / "configs"
    ref = expected["reference"]
    reference = json.loads((shipped / ref["adapter"] / f"{ref['version']}.json").read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory() as tmp:
        config_dir = Path(tmp) / "configs"
        shutil.copytree(shipped, config_dir)
        for stale in (config_dir / "chematic").glob("*.json"):
            stale.unlink()
        env = {**os.environ, "XSMARTS_CONFIG_DIR": str(config_dir)}
        subprocess.run([sys.executable, "-m", "xsmarts_autoconf", "update", "chematic"],
                       check=True, env=env, stdout=subprocess.DEVNULL)
        written = list((config_dir / "chematic").glob("*.json"))
        if len(written) != 1:
            raise ValueError(f"expected one chematic config, got {written}")
        candidate = json.loads(written[0].read_text(encoding="utf-8"))

    if candidate["rules_digest"] != expected["rules_digest"]:
        raise ValueError(f"rules digest {candidate['rules_digest']} != {expected['rules_digest']}")
    flags = reference["flags"]
    if len(flags) != expected["flag_count"] or set(flags) != set(candidate["flags"]):
        raise ValueError("flag set changed")
    differing = {
        name: {"rdkit": flags[name]["value"], "chematic": candidate["flags"][name]["value"]}
        for name in sorted(flags)
        if candidate["flags"][name]["value"] != flags[name]["value"]
    }
    want = {name: {"rdkit": row["rdkit"], "chematic": row["chematic"]}
            for name, row in expected["differing"].items()}
    new = {k: v for k, v in differing.items() if want.get(k) != v}
    aligned = sorted(set(want) - set(differing))
    report = {"schema": "xsmarts-autoconf-lane/v1", "autoconf_commit": commit,
              "rules_digest": candidate["rules_digest"],
              "chematic_version": candidate["version"],
              "reference": f"{ref['adapter']} {ref['version']}",
              "flags": len(flags), "same_as_reference": len(flags) - len(differing),
              "differing": differing, "unexpected": new, "newly_aligned": aligned}
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(f"{report['same_as_reference']}/{len(flags)} flags give {report['reference']}'s value; "
          f"{len(new)} unexpected, {len(aligned)} newly aligned")
    if new or aligned:
        print(rendered)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
