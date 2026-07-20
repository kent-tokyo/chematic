#!/usr/bin/env python3
"""IO-1 acceptance gate: compares chematic's SMILES-table tokenization
against real RDKit's, per row, across every scenario in
`gen_smiles_table_fixtures.py`'s manifest.

**Why this never compares chematic-canonical vs. RDKit-canonical SMILES
strings directly:** chematic/RDKit canonical-SMILES divergence is a known,
already-tracked, separate issue in this project (see the project's own
Morgan/ECFP-parity notes -- canonicalization is an open roadmap item, not
something this I/O milestone is scoped to fix or measure). A naive
string-equality gate here would report a large "parity gap" that has
nothing to do with whether the SMILES-table *tokenizer* correctly split
each line into (smiles, name, properties) columns -- the actual thing this
gate needs to prove.

Instead, three independent, tool-neutral checks are combined:

1. **Status parity** -- did both tools agree a row is `success` vs.
   unparseable, per row.
2. **Name/property parity** -- exact string equality of the extracted
   `name` and (key, value) property pairs. Pure tokenization, no chemistry,
   no canonicalizer involved on either side.
3. **Self-consistency, per tool, against the manifest's own known ground
   truth** -- did each tool's *own* canonicalizer agree that what it
   extracted from the SMILES column represents the same molecule as the
   manifest's authored SMILES string for that row, when parsed and
   re-canonicalized entirely within that same tool? This validates each
   tool's own tokenizer independently (proves the right substring was
   extracted) without ever comparing chematic's canonical form against
   RDKit's.

Usage:
    python scripts/smiles_table_io_parity.py \\
        --chematic <chematic_smiles_table_dump.jsonl> \\
        --rdkit-oracle <rdkit_smiles_table_oracle.jsonl> \\
        --manifest <manifest.json> \\
        --summary-out <out.json>

    python scripts/smiles_table_io_parity.py --self-test
"""

from __future__ import annotations

import argparse
import json
import sys


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


# Scenarios where chematic's own, documented behavior deliberately diverges
# from RDKit's -- reported separately below, never counted as a gate
# failure. See `smiles_table.rs`'s module doc comment for the source-cited
# justification of each.
KNOWN_NAME_DIVERGENCE_SCENARIOS = {
    # name_column=None: RDKit falls back to the physical line number as
    # `_Name`; chematic's MoleculeRecord::name is simply empty.
    "no_name_column",
}

KNOWN_PROPERTY_DIVERGENCE_SCENARIOS = {
    # RDKit's SmilesMolSupplier has NO CSV-quote-awareness for its comma
    # delimiter mode at all (confirmed via this oracle, not just inferred):
    # a quoted field with an embedded comma/quote is split into extra raw
    # columns by RDKit's literal splitting. Chematic's RFC 4180-subset
    # quoting support is a deliberate, beneficial divergence.
    "csv_quoted",
}


def normalize_no_title_line_properties(props_list, extra_keys, name_column, prefix):
    """When a scenario has no title line, each tool names extra columns
    with its own fallback convention ("column_N" for chematic, "Column_N"
    for RDKit, both 0-indexed over the full row) rather than the manifest's
    abstract key names. Remap `{prefix}_N` -> the manifest's own extra_keys[i]
    so both sides compare on the same semantic key set -- the naming
    convention difference is cosmetic, not a real divergence to gate on."""
    base = 1 if name_column is None else 2
    index_to_abstract = {base + i: k for i, k in enumerate(extra_keys)}
    out = {}
    for k, v in props_list:
        if k.startswith(prefix):
            try:
                idx = int(k[len(prefix) :])
            except ValueError:
                out[k] = v
                continue
            out[index_to_abstract.get(idx, k)] = v
        else:
            out[k] = v
    return out


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
    known_name_divergences = []
    known_property_divergences = []
    chematic_self_consistency_failures = []
    rdkit_self_consistency_failures = []
    malformed_row_agreement = {"expected": 0, "both_error": 0, "mismatched": []}

    for k in sorted(chematic_by_key):
        c = chematic_by_key[k]
        r = rdkit_by_key[k]
        scenario_name = k[0]
        scenario = manifest["scenarios"][scenario_name]

        if c["status"] != r["status"]:
            status_mismatches.append({"key": k, "chematic": c["status"], "rdkit": r["status"]})
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
            if scenario_name in KNOWN_NAME_DIVERGENCE_SCENARIOS:
                known_name_divergences.append({"key": k, "chematic": c["name"], "rdkit": r["name"]})
            else:
                name_mismatches.append({"key": k, "chematic": c["name"], "rdkit": r["name"]})

        opts = scenario["options"]
        if opts["title_line"]:
            c_props = {kk: vv for kk, vv in c["properties"]}
            r_props = {kk: vv for kk, vv in r["properties"]}
        else:
            extra_keys = scenario["extra_property_keys"]
            c_props = normalize_no_title_line_properties(c["properties"], extra_keys, opts["name_column"], "column_")
            r_props = normalize_no_title_line_properties(r["properties"], extra_keys, opts["name_column"], "Column_")
        if c_props != r_props:
            entry = {"key": k, "chematic": c_props, "rdkit": r_props}
            if scenario_name in KNOWN_PROPERTY_DIVERGENCE_SCENARIOS:
                known_property_divergences.append(entry)
            else:
                property_mismatches.append(entry)

        if c.get("self_consistent_with_known_smiles") is False:
            chematic_self_consistency_failures.append(k)
        if r.get("self_consistent_with_known_smiles") is False:
            rdkit_self_consistency_failures.append(k)

    success_count = sum(1 for c in chematic_by_key.values() if c["status"] == "success")

    gate_failures = []
    if status_mismatches:
        gate_failures.append(f"{len(status_mismatches)} row(s) with status mismatch (asymmetric parse failure)")
    if name_mismatches:
        gate_failures.append(f"{len(name_mismatches)} row(s) with name mismatch")
    if property_mismatches:
        gate_failures.append(f"{len(property_mismatches)} row(s) with property mismatch")
    if chematic_self_consistency_failures:
        gate_failures.append(
            f"{len(chematic_self_consistency_failures)} row(s) where chematic's own extraction+canonicalization "
            "disagreed with its own parse of the known ground-truth SMILES (tokenization bug)"
        )
    if malformed_row_agreement["mismatched"]:
        gate_failures.append(
            f"{len(malformed_row_agreement['mismatched'])} known-malformed row(s) where chematic/RDKit disagreed "
            "on whether the row is parseable"
        )

    summary = {
        "total_rows": total,
        "success_count": success_count,
        "error_count": total - success_count,
        "status_mismatches": status_mismatches,
        "name_mismatches": name_mismatches,
        "known_name_divergences_non_gating": known_name_divergences,
        "property_mismatches": property_mismatches,
        "known_property_divergences_non_gating": known_property_divergences,
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

    manifest = {
        "scenarios": {
            "s1": {"malformed_row_indices": []},
        }
    }

    def row(scenario, idx, status, name=None, props=None, self_consistent=None):
        r = {"scenario": scenario, "row_index": idx, "status": status}
        if status == "success":
            r["name"] = name or ""
            r["properties"] = props or []
            r["self_consistent_with_known_smiles"] = self_consistent
        return r

    # exact match
    c_rows = [row("s1", 0, "success", "a", [["K", "V"]], True)]
    r_rows = [row("s1", 0, "success", "a", [["K", "V"]], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("exact_match_passes", passed is True and summary["name_mismatches"] == []))

    # name mismatch
    c_rows = [row("s1", 0, "success", "a", [], True)]
    r_rows = [row("s1", 0, "success", "b", [], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("name_mismatch_fails", passed is False and len(summary["name_mismatches"]) == 1))

    # property mismatch
    c_rows = [row("s1", 0, "success", "a", [["K", "V1"]], True)]
    r_rows = [row("s1", 0, "success", "a", [["K", "V2"]], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("property_mismatch_fails", passed is False and len(summary["property_mismatches"]) == 1))

    # status mismatch
    c_rows = [row("s1", 0, "success", "a", [], True)]
    r_rows = [row("s1", 0, "error")]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("status_mismatch_fails", passed is False and len(summary["status_mismatches"]) == 1))

    # chematic self-consistency failure gates
    c_rows = [row("s1", 0, "success", "a", [], False)]
    r_rows = [row("s1", 0, "success", "a", [], True)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("chematic_self_consistency_gates", passed is False and len(summary["chematic_self_consistency_failures"]) == 1))

    # rdkit self-consistency failure is non-gating
    c_rows = [row("s1", 0, "success", "a", [], True)]
    r_rows = [row("s1", 0, "success", "a", [], False)]
    summary, passed = run(c_rows, r_rows, manifest)
    checks.append(("rdkit_self_consistency_non_gating", passed is True and len(summary["rdkit_self_consistency_failures_non_gating"]) == 1))

    # both-error on a known-malformed row is fine
    manifest2 = {"scenarios": {"s1": {"malformed_row_indices": [0]}}}
    c_rows = [row("s1", 0, "error")]
    r_rows = [row("s1", 0, "error")]
    summary, passed = run(c_rows, r_rows, manifest2)
    checks.append(("known_malformed_both_error_ok", passed is True and summary["malformed_row_agreement"]["both_error"] == 1))

    # pipeline error: mismatched key sets -> hard exit
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
