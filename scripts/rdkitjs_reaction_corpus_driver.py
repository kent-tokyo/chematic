#!/usr/bin/env python3
"""Run ``rdkitjs_reaction_corpus.mjs`` one rule per Node process.

RDKit.js (MinimalLib) can abort its WebAssembly instance on some reaction
inputs. Each rule runs in its own process; when one aborts, the row it was
running is read from the progress file, recorded as ``"error": "abort"``
and skipped on a rerun of that rule. Output: one JSONL with a header line,
every rule's rows, and the aborted rows.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rdkit", required=True, help="RDKit.js dist directory")
    ap.add_argument("--rules", type=Path, required=True)
    ap.add_argument("--reactants", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--node", default="node")
    ap.add_argument("--resume", action="store_true",
                    help="keep the rules already in --output and run the rest")
    args = ap.parse_args()
    script = Path(__file__).with_name("rdkitjs_reaction_corpus.mjs")
    rules = json.loads(args.rules.read_text())
    started = time.time()
    header = None
    aborts = 0
    done_rules: set[str] = set()
    if args.resume and args.output.exists():
        kept = []
        current: list[str] = []
        for line in args.output.read_text().splitlines():
            row = json.loads(line)
            if header is None:
                header = row
                continue
            if "seconds" in row:
                continue
            if row.get("start"):
                kept.extend(current)
                current = [line]
                done_rules.add(row["rule"])
            elif current:
                current.append(line)
            else:
                kept.append(line)
        # The last rule may have been cut off mid-copy: run it again.
        if current:
            done_rules.discard(json.loads(current[0])["rule"])
        args.output.write_text("".join(l + "\n" for l in [json.dumps(header)] + kept)
                               if header else "")
    mode = "a" if header is not None else "w"
    with args.output.open(mode) as out, tempfile.TemporaryDirectory() as tmp:
        part, progress = Path(tmp) / "part.jsonl", Path(tmp) / "progress"
        for i in range(len(rules)):
            if rules[i]["id"] in done_rules:
                continue
            skip: list[str] = []
            while True:
                cmd = [args.node, str(script), "--rdkit", args.rdkit, "--rules", str(args.rules),
                       "--reactants", str(args.reactants), "--output", str(part),
                       "--rule-index", str(i), "--progress", str(progress)]
                if skip:
                    cmd += ["--skip", ",".join(skip)]
                done = subprocess.run(cmd, capture_output=True, text=True)
                if done.returncode == 0:
                    break
                row = progress.read_text().strip()
                if not row or row in skip:
                    raise SystemExit(f"rule {i}: abort without progress\n{done.stderr[-2000:]}")
                skip.append(row)
                aborts += 1
            lines = part.read_text().splitlines()
            if header is None:
                header = json.loads(lines[0])
                header["rules"] = len(rules)
                out.write(json.dumps(header) + "\n")
            for line in lines[1:]:
                if '"seconds"' not in line:
                    out.write(line + "\n")
        out.write(json.dumps({"seconds": round(time.time() - started, 1), "aborted_rows": aborts}) + "\n")
    print(args.output, "aborted rows:", aborts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
