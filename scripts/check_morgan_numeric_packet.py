#!/usr/bin/env python3
"""Recompute the path-free Morgan packet's timing counts, statistics and gate.

This verifies numeric consistency. Chemistry equality depends on the original
full-output local audit; hashes in this export commit to those output journals.
"""
import argparse
import itertools
import gzip
import json
from pathlib import Path
import random
import statistics

ARMS = ("baseline", "candidate", "rdkit", "cosmolkit")
PRIMARY = ("morgan2_2048", "morgan2_chiral", "morgan3_2048")
CONTROLS = ("tpsa", "labute_asa", "num_rings", "maccs")
MODES = ("first_use", "repeated_use", "parse_inclusive")


def stats(a, b, seed):
    ratios = [x / y for x, y in zip(a, b, strict=True)]
    rng = random.Random(seed)
    boot = sorted(
        statistics.median(ratios[rng.randrange(24)] for _ in range(24))
        for _ in range(10000)
    )
    return statistics.median(ratios), boot[249], boot[9749], sum(x > 1 for x in ratios)


def check(packet):
    assert packet["schema"] == "morgan-numeric-packet/v1"
    assert len(packet["profiles"]) == 2 and len(packet["summary"]) == 78
    assert packet["protocol"]["primary"] == [
        "0/morgan2_chiral/first_use",
        "5000/morgan2_chiral/first_use",
    ]
    saved = {
        (r["offset"], r["operation"], r["mode"], r["comparator"]): r
        for r in packet["summary"]
    }
    assert len(saved) == 78
    cells = batches = rows = journals = 0
    source = []
    primary = []
    for profile in packet["profiles"]:
        offset = profile["offset"]
        assert offset in (0, 5000) and profile["rows"] == 1000
        orders = list(itertools.permutations(ARMS))
        random.Random(20261017 + offset).shuffle(orders)
        assert [b["order"] for b in profile["blocks"]] == [list(o) for o in orders]
        hashes = {}
        for block in profile["blocks"]:
            assert set(block["arms"]) == set(ARMS)
            for arm, ops in block["arms"].items():
                journals += 1
                assert set(ops) == set(
                    PRIMARY + CONTROLS if arm in ARMS[:2] else PRIMARY
                )
                for op, modes in ops.items():
                    assert set(modes) == set(MODES)
                    for mode, d in modes.items():
                        budget = 50000000 if op in PRIMARY else 20000000
                        ts = d["batch_ns"]
                        assert ts and all(type(t) is int and t > 0 for t in ts)
                        assert (
                            d["mode"] == mode and d["minimum_accumulated_ns"] == budget
                        )
                        assert d["elapsed_ns_raw"] == sum(ts) >= budget
                        assert d["batches"] == len(ts) and d[
                            "processed_rows"
                        ] == 1000 * len(ts)
                        assert d["elapsed_ns_per_batch"] == sum(ts) // len(ts)
                        hs = d["output_sha256_per_batch"]
                        assert len(hs) == len(ts) and len(set(hs)) == 1
                        assert len(hs[0]) == 64 and all(
                            c in "0123456789abcdef" for c in hs[0]
                        )
                        assert hashes.setdefault((arm, op), hs[0]) == hs[0]
                        cells += 1
                        batches += len(ts)
                        rows += d["processed_rows"]
        for op in PRIMARY + CONTROLS:
            assert hashes["baseline", op] == hashes["candidate", op]
            if op in PRIMARY:
                assert (
                    hashes["rdkit", op]
                    == hashes["cosmolkit", op]
                    == hashes["candidate", op]
                )
            for mode in MODES:
                b = [
                    block["arms"]["candidate"][op][mode]["elapsed_ns_per_batch"]
                    for block in profile["blocks"]
                ]
                for seed, arm in enumerate(("baseline", "rdkit", "cosmolkit")):
                    if op not in PRIMARY and arm != "baseline":
                        continue
                    a = [
                        block["arms"][arm][op][mode]["elapsed_ns_per_batch"]
                        for block in profile["blocks"]
                    ]
                    row = saved[offset, op, mode, arm]
                    m, lo, hi, wins = stats(a, b, seed)
                    assert (
                        row["paired_median"],
                        row["ci95_low"],
                        row["ci95_high"],
                        row["winning_blocks"],
                    ) == (m, lo, hi, wins)
                    assert row["comparator_median_ns"] == statistics.median(a) and row[
                        "candidate_median_ns"
                    ] == statistics.median(b)
                    assert row["blocks"] == 24 and row["strict_external_win"] == (
                        arm != "baseline" and lo > 1 and wins == 24
                    )
                    if arm == "baseline":
                        source.append(row)
                        if op == "morgan2_chiral" and mode == "first_use":
                            primary.append(row)
    assert {p["offset"] for p in packet["profiles"]} == {0, 5000}
    assert len(source) == 42 and len(primary) == 2 and journals == 192 and cells == 2880
    broad = packet["broad"]
    assert (
        broad["source_rows_equal"]
        and broad["source_exact_rows_including_refusals"] == 14072
    )
    assert len(broad["comparisons"]["candidate"]) == 8
    assert all(
        v == {"equal": 10000 if k.startswith("rdkit-js") else 4072}
        for k, v in broad["comparisons"]["candidate"].items()
    )
    adopted = all(r["ci95_low"] > 1 for r in primary) and all(
        r["paired_median"] >= 0.99 for r in source
    )
    d = packet["decision"]
    assert d == {
        "source_adoption_passed": adopted,
        "journals": journals,
        "accounting_cells": cells,
        "total_batches": batches,
        "total_processed_rows_in_timed_batches": rows,
        "broad_source_and_rdkit_equal_cells": 56288,
    }
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    args = parser.parse_args()
    raw = args.packet.read_bytes()
    if args.packet.suffix == ".gz":
        raw = gzip.decompress(raw)
    print(json.dumps(check(json.loads(raw)), indent=2))


if __name__ == "__main__":
    main()
