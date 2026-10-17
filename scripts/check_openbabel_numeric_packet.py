#!/usr/bin/env python3
"""Audit the small, path-free Open Babel numeric packet without chemistry SDKs."""
import argparse
import collections
import itertools
import json
from pathlib import Path
import random
import statistics


def stats(a, b, seed):
    ratios = [x / y for x, y in zip(a, b, strict=True)]
    rng = random.Random(seed)
    boot = sorted(
        statistics.median(ratios[rng.randrange(24)] for _ in range(24))
        for _ in range(10000)
    )
    return statistics.median(ratios), boot[249], boot[9749], sum(x > 1 for x in ratios)


def check(packet):
    assert packet["schema"] == "openbabel-numeric-packet/v1"
    assert len(packet["profiles"]) == 10 and len(packet["summary"]) == 36
    saved = {(r["profile"], r["mode"], r["comparator"]): r for r in packet["summary"]}
    assert len(saved) == 36
    source, primary, external = [], [], []
    calls = batches = records = 0
    for profile in packet["profiles"]:
        blocks, count = profile["blocks"], profile["records"]
        assert len(blocks) == 24
        arms = (
            {"baseline", "candidate", "openbabel"}
            if profile["kind"] == "sdf"
            else {"baseline", "candidate"}
        )
        orders = collections.Counter(tuple(b["order"]) for b in blocks)
        assert set(orders) == set(itertools.permutations(sorted(arms)))
        assert all(n == 24 // len(orders) for n in orders.values())
        assert (
            profile["snapshot_sha256"]["baseline"]
            == profile["snapshot_sha256"]["candidate"]
        )
        assert (
            profile["rewrite_sha256"]["baseline"]
            == profile["rewrite_sha256"]["candidate"]
        )
        units = {}
        for block in blocks:
            assert set(block) == arms | {"order"}
            for arm in arms:
                assert set(block[arm]) == {"parse", "write", "roundtrip"}
                for mode, row in block[arm].items():
                    assert all(type(v) is int and v > 0 for v in row.values())
                    assert row["minimum_repeats"] == 3 and row["repeats"] >= 3
                    assert (
                        row["records_per_batch"] == count
                        and row["records"] == count * row["repeats"]
                    )
                    assert (
                        row["minimum_warmup_ns"] == 500000000
                        and row["warmup_elapsed_ns"] >= 500000000
                    )
                    assert (
                        row["minimum_timed_ns"] == 100000000
                        and row["elapsed_ns_raw"] >= 100000000
                    )
                    assert (
                        row["elapsed_ns"]
                        == row["elapsed_ns_per_batch"]
                        == row["elapsed_ns_raw"] // row["repeats"]
                    )
                    assert row["output_units"] % row["repeats"] == 0
                    unit = row["output_units"] // row["repeats"]
                    assert units.setdefault((arm, mode), unit) == unit
                    calls += 1
                    batches += row["repeats"]
                    records += row["records"]
        for mode in ("parse", "write", "roundtrip"):
            assert units["baseline", mode] == units["candidate", mode]
            for arm, seed in [("baseline", 0)] + (
                [("openbabel", 1)] if "openbabel" in arms else []
            ):
                row = saved[profile["id"], mode, arm]
                m, lo, hi, wins = stats(
                    [b[arm][mode]["elapsed_ns"] for b in blocks],
                    [b["candidate"][mode]["elapsed_ns"] for b in blocks],
                    seed,
                )
                assert (m, lo, hi, wins) == (
                    row["paired_median"],
                    row["ci95_low"],
                    row["ci95_high"],
                    row["winning_blocks"],
                )
                assert row["blocks"] == 24
                assert row["strict_external_win"] == (
                    arm == "openbabel" and row["ranked"] and lo > 1 and wins == 24
                )
                if arm == "baseline":
                    assert row["ranked"]
                    source.append(row)
                    if profile["id"] + "/" + mode in packet["protocol"]["primary"]:
                        primary.append(row)
                else:
                    external.append(row)
    assert (
        calls == 1584 and len(source) == 30 and len(primary) == 2 and len(external) == 6
    )
    unranked = [r for r in external if not r["ranked"]]
    assert [(r["profile"], r["mode"]) for r in unranked] == [("sdf-0", "parse")]
    passed = all(r["ci95_low"] > 1 for r in primary) and all(
        r["paired_median"] >= 0.99 for r in source
    )
    assert passed == packet["gate_passed"]
    return {
        "gate_passed": passed,
        "native_calls": calls,
        "timed_batches": batches,
        "processed_records": records,
        "strict_external_wins": sum(r["strict_external_win"] for r in external),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    args = parser.parse_args()
    print(json.dumps(check(json.loads(args.packet.read_text())), indent=2))
