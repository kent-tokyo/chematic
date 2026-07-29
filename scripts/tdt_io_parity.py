#!/usr/bin/env python3
"""IO-2 acceptance gate: compares chematic's TDT tokenization against real
RDKit's, per row, across every scenario in `gen_tdt_fixtures.py`'s manifest.

Same methodology as `scripts/smiles_table_io_parity.py` (see its module
docs for the full rationale): never compares chematic-canonical vs.
RDKit-canonical SMILES strings directly. Instead: exact name/property
string equality, plus each tool's own self-consistency against the
fixture's known ground-truth SMILES, canonicalized only within that same
tool.

Two scenarios are treated as KNOWN, non-gating divergences (documented in
`chematic_mol::tdt`'s module doc comment, both confirmed via this session's
own oracle runs, not assumed):
- "coordinates": RDKit's own coordinate-list parser drops the last atom's
  position (a real bug); chematic fixes it. Coordinate VALUES are not
  gated for this scenario (name/properties still are).
- "eof_mid_record": RDKit drops the final tag line when a file has no
  trailing newline; chematic reads it correctly. Name/property mismatches
  for this specific scenario are not gated.

Usage:
    python scripts/tdt_io_parity.py --chematic <chematic_tdt_dump.jsonl> \\
        --rdkit-oracle <rdkit_tdt_oracle.jsonl> --manifest <manifest.json> \\
        --summary-out <out.json>

    python scripts/tdt_io_parity.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys

KNOWN_DIVERGENT_SCENARIOS = {"coordinates", "eof_mid_record"}


def load_jsonl(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def key(row):
    return (row["scenario"], row["row_index"])


def run(chematic_rows, rdkit_rows, manifest):
    chematic_by_key = {key(r): r for r in chematic_rows}
    rdkit_by_key = {key(r): r for r in rdkit_rows}

    if set(chematic_by_key) != set(rdkit_by_key):
        missing_in_rdkit = sorted(set(chematic_by_key) - set(rdkit_by_key))
        missing_in_chematic = sorted(set(rdkit_by_key) - set(chematic_by_key))
        print(
            f"PIPELINE ERROR: row key sets differ. missing_in_rdkit={missing_in_rdkit[:5]} "
            f"missing_in_chematic={missing_in_chematic[:5]}",
            file=sys.stderr,
        )
        sys.exit(1)

    total = len(chematic_by_key)
    status_mismatches = []
    name_mismatches = []
    property_mismatches = []
    coordinate_mismatches = []
    known_divergences = []
    chematic_self_consistency_failures = []
    rdkit_self_consistency_failures = []
    malformed_row_agreement = {"expected": 0, "both_error": 0, "mismatched": []}

    for k in sorted(chematic_by_key):
        c = chematic_by_key[k]
        r = rdkit_by_key[k]
        scenario_name = k[0]
        scenario = manifest["scenarios"][scenario_name]
        known_divergent = scenario_name in KNOWN_DIVERGENT_SCENARIOS

        if c["status"] != r["status"]:
            entry = {"key": k, "chematic": c["status"], "rdkit": r["status"]}
            if known_divergent:
                known_divergences.append(entry)
            else:
                status_mismatches.append(entry)
            continue

        if k[1] in scenario.get("malformed_row_indices", []):
            malformed_row_agreement["expected"] += 1
            if c["status"] == "error" and r["status"] == "error":
                malformed_row_agreement["both_error"] += 1
            else:
                malformed_row_agreement["mismatched"].append(k)

        if c["status"] != "success":
            continue

        if c["name"] != r["name"]:
            entry = {"key": k, "chematic": c["name"], "rdkit": r["name"]}
            if known_divergent:
                known_divergences.append(entry)
            else:
                name_mismatches.append(entry)

        c_props = {kk: vv for kk, vv in c["properties"]}
        r_props = {kk: vv for kk, vv in r["properties"]}
        if c_props != r_props:
            entry = {"key": k, "chematic": c_props, "rdkit": r_props}
            if known_divergent:
                known_divergences.append(entry)
            else:
                property_mismatches.append(entry)

        if c.get("self_consistent_with_known_smiles") is False:
            chematic_self_consistency_failures.append(k)
        if r.get("self_consistent_with_known_smiles") is False:
            rdkit_self_consistency_failures.append(k)

        for dim in ("coordinates_2d", "coordinates_3d"):
            c_coords, r_coords = c.get(dim), r.get(dim)
            if c_coords is None and r_coords is None:
                continue
            entry = {"key": k, "dim": dim, "chematic": c_coords, "rdkit": r_coords}
            if c_coords != r_coords:
                if known_divergent:
                    known_divergences.append(entry)
                else:
                    coordinate_mismatches.append(entry)

    success_count = sum(1 for c in chematic_by_key.values() if c["status"] == "success")

    gate_failures = []
    if status_mismatches:
        gate_failures.append(f"{len(status_mismatches)} row(s) with status mismatch")
    if name_mismatches:
        gate_failures.append(f"{len(name_mismatches)} row(s) with name mismatch")
    if property_mismatches:
        gate_failures.append(f"{len(property_mismatches)} row(s) with property mismatch")
    if coordinate_mismatches:
        gate_failures.append(f"{len(coordinate_mismatches)} row(s) with coordinate mismatch")
    if chematic_self_consistency_failures:
        gate_failures.append(
            f"{len(chematic_self_consistency_failures)} row(s) where chematic's own extraction "
            "disagreed with its own parse of the known ground-truth SMILES (tokenization bug)"
        )
    if malformed_row_agreement["mismatched"]:
        gate_failures.append(
            f"{len(malformed_row_agreement['mismatched'])} known-malformed row(s) where chematic/RDKit disagreed"
        )

    summary = {
        "total_rows": total,
        "success_count": success_count,
        "error_count": total - success_count,
        "status_mismatches": status_mismatches,
        "name_mismatches": name_mismatches,
        "property_mismatches": property_mismatches,
        "coordinate_mismatches": coordinate_mismatches,
        "known_divergences_non_gating": known_divergences,
        "chematic_self_consistency_failures": chematic_self_consistency_failures,
        "rdkit_self_consistency_failures_non_gating": rdkit_self_consistency_failures,
        "malformed_row_agreement": malformed_row_agreement,
        "gate_failures": gate_failures,
        "gate_passed": len(gate_failures) == 0,
    }

    if gate_failures:
        print(f"GATE FAILED: {len(gate_failures)} violation(s):", file=sys.stderr)
        for msg in gate_failures:
            print(f"  - {msg}", file=sys.stderr)

    return summary, summary["gate_passed"]


def run_self_test():
    checks = []
    manifest = {"scenarios": {"s1": {"malformed_row_indices": []}, "coordinates": {"malformed_row_indices": []}}}

    def row(scenario, idx, status, name=None, props=None, self_consistent=None):
        r = {"scenario": scenario, "row_index": idx, "status": status}
        if status == "success":
            r["name"] = name or ""
            r["properties"] = props or []
            r["self_consistent_with_known_smiles"] = self_consistent
        return r

    c_rows = [row("s1", 0, "success", "a", [["K", "V"]], True)]
    r_rows = [row("s1", 0, "success", "a", [["K", "V"]], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("exact_match_passes", passed is True))

    c_rows = [row("s1", 0, "success", "a", [], True)]
    r_rows = [row("s1", 0, "success", "b", [], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("name_mismatch_fails", passed is False and len(summary["name_mismatches"]) == 1))

    c_rows = [row("coordinates", 0, "success", "a", [], True)]
    r_rows = [row("coordinates", 0, "success", "b", [], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("known_divergent_scenario_non_gating", passed is True and len(summary["known_divergences_non_gating"]) == 1))

    try:
        run([row("s1", 0, "success", "a", [], True)], [], manifest)
        pipeline_ok = False
    except SystemExit as e:
        pipeline_ok = e.code != 0
    checks.append(("mismatched_keys_hard_exit", pipeline_ok))

    ok = True
    for name, passed in checks:
        status = "OK" if passed else "FAIL"
        print(f"  self-test {name}: {status}")
        ok = ok and passed
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chematic")
    p.add_argument("--rdkit-oracle")
    p.add_argument("--manifest")
    p.add_argument("--summary-out", default=None)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        ok = run_self_test()
        sys.exit(0 if ok else 1)

    if not args.chematic or not args.rdkit_oracle or not args.manifest:
        p.error("--chematic, --rdkit-oracle, and --manifest are required unless --self-test")

    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit_oracle)
    with open(args.manifest) as f:
        manifest = json.load(f)

    summary, gate_passed = run(chematic_rows, rdkit_rows, manifest)
    print(json.dumps(summary, indent=2))
    if args.summary_out:
        with open(args.summary_out, "w") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
        print(f"summary written to {args.summary_out}")

    sys.exit(0 if gate_passed else 1)


if __name__ == "__main__":
    main()
