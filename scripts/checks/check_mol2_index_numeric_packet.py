#!/usr/bin/env python3
"""Recompute the fixed MOL2 experiment's complete numeric adoption gate."""
from __future__ import annotations

import argparse
import collections
import gzip
import itertools
import json
import random
import statistics
from pathlib import Path

PROFILES = (
    "ethanol",
    "chain-32",
    "chain-64",
    "chain-256",
    "chain-1024",
    "chain-4096",
    "sparse-chain-1024",
    "reversed-chain-1024",
    "mol2-stereo-fixtures",
)
MODES = ("parse", "write", "roundtrip")
READ_PRIMARY = [
    "chain-1024/parse",
    "chain-4096/parse",
    "sparse-chain-1024/parse",
    "reversed-chain-1024/parse",
]
WRITE_PRIMARY = ["ethanol/write", "chain-1024/write", "chain-4096/write"]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def interval(values: list[float], seed: int) -> dict:
    randomizer = random.Random(seed)
    count = len(values)
    bootstrap = sorted(
        statistics.median(randomizer.choices(values, k=count)) for _ in range(10_000)
    )
    return {
        "median": statistics.median(values),
        "ci95": [bootstrap[249], bootstrap[9749]],
        "winning_blocks": sum(x > 1 for x in values),
        "blocks": count,
    }


def read_packet(path: Path) -> dict:
    raw = (
        gzip.decompress(path.read_bytes())
        if path.suffix == ".gz"
        else path.read_bytes()
    )
    return json.loads(raw)


def check(packet: dict) -> dict:
    require(packet["schema"] == "chematic-mol2-index-numeric/v1", "schema")
    require(packet["label"] in ("index-only", "index-and-fixed4"), "epoch")
    require(set(packet["profiles"]) == set(PROFILES), "all profiles required")
    preflight = packet["preflight"]
    require(
        preflight["records"] == 32 and preflight["formal_samples"] == 0,
        "preflight denominator",
    )
    source_gate = packet["source_gate"]
    primary = READ_PRIMARY + (
        WRITE_PRIMARY if packet["label"] == "index-and-fixed4" else []
    )
    require(
        source_gate["primary"] == primary
        and source_gate["all_lane_median_floor"] == 0.99,
        "preregistered gate",
    )
    decision = packet["decision"]
    require(
        decision["complete"] and not decision["universal_superiority"], "claim boundary"
    )
    totals = dict(blocks=0, invocations=0, batches=0, processed_records=0)
    lanes = {}
    for name in PROFILES:
        report = packet["profiles"][name]
        info = preflight["profiles"][name]
        require(report["complete"] and report["profile"] == name, "profile completion")
        require(
            info["records"] == (24 if name == "mol2-stereo-fixtures" else 1),
            "input count",
        )
        require(
            info["source_snapshots_equal"] and info["source_rewrites_equal"],
            "source outputs",
        )
        eligible = info["openbabel_read_success"] and not info["external_errors"]
        require(info["external_eligible"] == eligible, "semantic eligibility")
        blocks = report["blocks"]
        require(len(blocks) == 24, "all blocks required")
        orders = [row["order"] for row in blocks]
        require(orders == info["orders"], "fixed order")
        require(
            packet["external_timing"] == (packet["label"] == "index-only"),
            "epoch timing scope",
        )
        arms = ["baseline", "candidate"] + (
            ["openbabel"]
            if packet["external_timing"] and info["openbabel_read_success"]
            else []
        )
        expected = collections.Counter(
            {
                order: 24 // len(list(itertools.permutations(arms)))
                for order in itertools.permutations(arms)
            }
        )
        require(collections.Counter(map(tuple, orders)) == expected, "balanced orders")
        units = collections.defaultdict(set)
        for number, block in enumerate(blocks):
            require(block["index"] == number, "block index")
            for arm in block["order"]:
                for mode in MODES:
                    row = block[arm][mode]
                    for field in (
                        "repeats",
                        "records",
                        "output_units",
                        "elapsed_ns",
                        "warmup_repeats",
                        "warmup_elapsed_ns",
                        "minimum_warmup_ns",
                        "minimum_timed_ns",
                    ):
                        require(
                            type(row[field]) is int and row[field] > 0,
                            "positive integer accounting",
                        )
                    repetitions = row["repeats"]
                    require(row["mode"] == mode and repetitions >= 3, "operation")
                    require(
                        row["engine"]
                        == ("openbabel" if arm == "openbabel" else "chematic"),
                        "engine",
                    )
                    require(
                        row["records"] == info["records"] * repetitions,
                        "record accounting",
                    )
                    require(
                        row["minimum_warmup_ns"] == 250_000_000
                        and row["warmup_elapsed_ns"] >= 250_000_000,
                        "warmup",
                    )
                    require(
                        row["minimum_timed_ns"] == 100_000_000
                        and row["elapsed_ns"] >= 100_000_000,
                        "timed budget",
                    )
                    require(row["output_units"] % repetitions == 0, "unit accounting")
                    units[arm, mode].add(row["output_units"] // repetitions)
                    totals["invocations"] += 1
                    totals["batches"] += repetitions
                    totals["processed_records"] += row["records"]
        require(all(len(v) == 1 for v in units.values()), "stable work units")
        for mode in MODES:
            require(
                units["baseline", mode] == units["candidate", mode], "equal source work"
            )

            def speed(competitor: str, seed: int) -> dict:
                values = []
                for block in blocks:
                    a, b = block[competitor][mode], block["candidate"][mode]
                    values.append(
                        (a["elapsed_ns"] / a["repeats"])
                        / (b["elapsed_ns"] / b["repeats"])
                    )
                return interval(values, seed)

            source = speed("baseline", 0)
            external = speed("openbabel", 1) if "openbabel" in arms else None
            summary = {
                "source": source,
                "openbabel": external,
                "external_eligible": eligible,
                "strict_external_win": bool(
                    eligible
                    and external
                    and external["ci95"][0] > 1
                    and external["winning_blocks"] == 24
                ),
            }
            require(summary == report["summary"][mode], "recomputed profile interval")
            require(
                summary == decision["lanes"][name + "/" + mode], "decision interval"
            )
            lanes[name + "/" + mode] = summary
        totals["blocks"] += len(blocks)
    floor = [name for name, row in lanes.items() if row["source"]["median"] < 0.99]
    negative = [name for name, row in lanes.items() if row["source"]["ci95"][1] < 1]
    passed = not floor and all(
        lanes[name]["source"]["ci95"][0] > 1
        and lanes[name]["source"]["winning_blocks"] == 24
        for name in primary
    )
    require(
        decision["primary"] == primary and decision["source_floor_failures"] == floor,
        "floor decision",
    )
    require(decision["source_negative_intervals"] == negative, "negative intervals")
    require(decision["source_gate_passed"] == passed, "adoption decision")
    for name, value in totals.items():
        require(packet["audit"][name] == value, "independent totals")
    require(packet["audit"]["source_gate_passed"] == passed, "audit gate")
    return dict(passed=True, source_gate_passed=passed, **totals)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("packet", type=Path)
    args = parser.parse_args()
    print(json.dumps(check(read_packet(args.packet)), sort_keys=True))


if __name__ == "__main__":
    main()
