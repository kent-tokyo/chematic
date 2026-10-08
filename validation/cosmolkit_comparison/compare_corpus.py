#!/usr/bin/env python3
"""Score engine dumps from ``run_corpus.py`` against a reference dump.

For each operation and engine: ``match`` (identical value), ``mismatch``,
``refused`` (typed refusal), ``unsupported``, ``error`` and ``uncomparable``
(the reference or the engine did not produce a value). Float values are
compared exactly; ``match_1e-9`` additionally counts relative differences
below 1e-9. The 31 SMARTS queries are reported per query and as one
``smarts_31`` group. Unsupported, refused and failed cells are never counted
as matches.

    python compare_corpus.py --reference rdkit=ref.jsonl \
        --engine chematic=chem.jsonl --engine cosmolkit=ck.jsonl --output summary.json
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def load(path: Path):
    lines = path.read_text().splitlines()
    header = json.loads(lines[0])
    rows = {r["index"]: r for r in map(json.loads, lines[1:])}
    return header, rows


def close(a, b) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            return False
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
    return a == b


def score(ref_rows, rows, max_examples=8):
    ops = sorted({op for r in ref_rows.values() for op in r["ops"]} - {"cip_abstain", "input_canonical"})
    counts = defaultdict(Counter)
    examples = defaultdict(list)
    parse = Counter()
    for idx, ref in ref_rows.items():
        row = rows.get(idx)
        ref_ok = ref["parse"]["status"] == "ok"
        eng_ok = row is not None and row["parse"]["status"] == "ok"
        parse[f"ref_{'ok' if ref_ok else 'fail'}__engine_{'ok' if eng_ok else (row['parse']['status'] if row else 'missing')}"] += 1
        for op in ops:
            # round trips are judged against the input molecule, not the
            # reference engine's own (possibly lossy) writer
            r = ref["ops"].get("input_canonical" if op.endswith("_rt") else op) if ref_ok else None
            e = row["ops"].get(op) if eng_ok else None
            if r is None or r["status"] != "ok" or r.get("value") is None:
                counts[op]["uncomparable"] += 1
                continue
            if e is None:
                counts[op]["engine_parse_failed" if not eng_ok else "unsupported"] += 1
                continue
            if e["status"] != "ok":
                counts[op][e["status"]] += 1
                continue
            if e["value"] == r["value"] and type(e["value"]) is type(r["value"]):
                counts[op]["match"] += 1
                counts[op]["match_1e-9"] += 1
            elif close(e["value"], r["value"]):
                counts[op]["match_1e-9"] += 1
                counts[op]["float_last_bits"] += 1
            else:
                counts[op]["mismatch"] += 1
                if len(examples[op]) < max_examples:
                    examples[op].append({"index": idx, "smiles": ref["smiles"],
                                         "reference": r["value"], "engine": e["value"]})
    group = Counter()
    for op, c in counts.items():
        if op.startswith("smarts:"):
            group.update(c)
    out = {op: dict(sorted(c.items())) for op, c in sorted(counts.items())}
    out["smarts_31"] = dict(sorted(group.items()))
    return {"parse": dict(sorted(parse.items())), "operations": out,
            "examples": {k: v for k, v in sorted(examples.items())}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, metavar="NAME=PATH")
    ap.add_argument("--engine", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    ref_name, ref_path = args.reference.split("=", 1)
    ref_header, ref_rows = load(Path(ref_path))
    report = {"schema": "cosmolkit-corpus-comparison-summary/v1",
              "corpus": ref_header["corpus"], "corpus_sha256": ref_header["corpus_sha256"],
              "rows": ref_header["rows"],
              "reference": {"engine": ref_name, "version": ref_header["engine_version"]},
              "engines": {}}
    for spec in args.engine:
        name, path = spec.split("=", 1)
        header, rows = load(Path(path))
        if header["corpus_sha256"] != ref_header["corpus_sha256"]:
            raise SystemExit(f"{name}: corpus hash differs from reference")
        result = score(ref_rows, rows)
        result["version"] = header["engine_version"]
        result["api_notes"] = header.get("api_notes", {})
        abst = sum(len(r["ops"].get("cip_abstain", {}).get("value") or []) for r in rows.values())
        if abst:
            result["cip_abstained_centres"] = abst
        report["engines"][name] = result
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    # compact console table
    names = list(report["engines"])
    ops = list(next(iter(report["engines"].values()))["operations"])
    print(f"{'operation':24}" + "".join(f"{n:>34}" for n in names))
    for op in ops:
        if op.startswith("smarts:"):
            continue
        cells = []
        for n in names:
            c = report["engines"][n]["operations"].get(op, {})
            tot = sum(v for k, v in c.items() if k not in ("match_1e-9", "uncomparable"))
            extra = "".join(f" {k[:5]}={v}" for k, v in c.items()
                            if k in ("refused", "unsupported", "error", "engine_parse_failed") and v)
            cells.append(f"{c.get('match', 0)}/{tot} (~{c.get('match_1e-9', 0)}){extra}")
        print(f"{op:24}" + "".join(f"{x:>34}" for x in cells))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
