#!/usr/bin/env python3
"""Dump one engine's answers for every row of a SMILES corpus.

Writes JSONL: a header line, then one line per input row with each
operation's status (``ok``, ``refused``, ``unsupported``, ``error``) and value.
Rows are processed in shards, each in a fresh interpreter, so a native crash
costs one shard (recorded as ``crash`` rows) instead of the whole run.

    python run_corpus.py --engine cosmolkit --corpus CORPUS.smi --output out.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def run_shard(engine_name: str, rows: list[str], start: int) -> list[dict]:
    from corpus_engines import ENGINES, RefusedError, UnsupportedError

    engine = ENGINES[engine_name]()
    out = []
    for offset, smiles in enumerate(rows):
        rec = {"index": start + offset, "smiles": smiles, "ops": {}}
        try:
            mol = engine["parse"](smiles)
        except Exception as exc:  # noqa: BLE001 - every parse failure is recorded
            rec["parse"] = {"status": "parse_error", "error": f"{type(exc).__name__}: {exc}"[:300]}
            out.append(rec)
            continue
        rec["parse"] = {"status": "ok"}
        for name, fn in engine["ops"].items():
            try:
                rec["ops"][name] = {"status": "ok", "value": fn(mol)}
            except RefusedError as exc:
                rec["ops"][name] = {"status": "refused", "error": str(exc)[:200]}
            except UnsupportedError as exc:
                rec["ops"][name] = {"status": "unsupported", "error": str(exc)[:200]}
            except Exception as exc:  # noqa: BLE001
                rec["ops"][name] = {"status": "error", "error": f"{type(exc).__name__}: {exc}"[:300]}
        out.append(rec)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", choices=["rdkit", "chematic", "cosmolkit"], required=True)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--shard-size", type=int, default=500)
    ap.add_argument("--limit", type=int, default=0, help="first N rows only (0 = all)")
    ap.add_argument("--_shard", nargs=2, type=int, help=argparse.SUPPRESS)
    args = ap.parse_args()

    data = args.corpus.read_bytes()
    rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
    if args.limit:
        rows = rows[: args.limit]

    if args._shard:  # child process
        start, stop = args._shard
        for rec in run_shard(args.engine, rows[start:stop], start):
            sys.stdout.write(json.dumps(rec, sort_keys=True) + "\n")
        return 0

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from corpus_engines import API_NOTES, ENGINES  # noqa: E402

    version = ENGINES[args.engine]()["version"]
    notes_key = args.engine
    if args.engine == "cosmolkit" and not str(version).startswith("0.3"):
        notes_key = "cosmolkit-0.5"
    header = {
        "schema": "cosmolkit-corpus-comparison/v1",
        "engine": args.engine,
        "engine_version": version,
        "api_notes": API_NOTES.get(notes_key, {}),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "corpus": args.corpus.name,
        "corpus_sha256": hashlib.sha256(data).hexdigest(),
        "rows": len(rows),
    }
    tmp = args.output.with_name(args.output.name + ".partial")
    t0 = time.time()
    with tmp.open("w") as out:
        out.write(json.dumps(header, sort_keys=True) + "\n")
        for start in range(0, len(rows), args.shard_size):
            stop = min(start + args.shard_size, len(rows))
            cmd = [sys.executable, __file__, "--engine", args.engine, "--corpus", str(args.corpus),
                   "--output", str(args.output), "--_shard", str(start), str(stop)]
            if args.limit:
                cmd += ["--limit", str(args.limit)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            got = [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]
            done = {r["index"] for r in got}
            for rec in got:
                out.write(json.dumps(rec, sort_keys=True) + "\n")
            if proc.returncode != 0:
                # The shard died (native crash or uncaught error); keep what it
                # wrote, retry the rest one row at a time to isolate the culprit.
                for i in range(start, stop):
                    if i in done:
                        continue
                    one = subprocess.run(cmd[:-2] + [str(i), str(i + 1)], capture_output=True, text=True)
                    lines = [l for l in one.stdout.splitlines() if l.strip()]
                    if one.returncode == 0 and lines:
                        out.write(lines[0] + "\n")
                    else:
                        out.write(json.dumps({"index": i, "smiles": rows[i], "ops": {},
                                              "parse": {"status": "crash",
                                                        "error": (one.stderr or "")[-300:]}},
                                             sort_keys=True) + "\n")
            print(f"{args.engine}: {stop}/{len(rows)} rows ({time.time() - t0:.0f}s)", file=sys.stderr)
    tmp.replace(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
