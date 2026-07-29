#!/usr/bin/env python3
"""Generate SMILES-table I/O fixtures for the IO-1 oracle gate.

Produces a set of scenario files (each with a small manifest of KNOWN ground
truth this script itself authors -- not derived from either chematic or
RDKit) covering the categories required by the IO-1 acceptance gate:
space/tab/comma delimiters, header/no-header, name/no-name column, extra
properties, quoted CSV, blank lines, comments, malformed SMILES, isotopes,
charges, disconnected structures, stereochemistry.

The "general" scenario draws real, chemically diverse SMILES from the
Morgan M4-A0 corpus (`/tmp/m4a0/combined_input.csv`, already validated
elsewhere in this project) rather than hand-authoring hundreds of molecules.

Usage:
    python scripts/gen_smiles_table_fixtures.py --out-dir <dir> \\
        --corpus /tmp/m4a0/combined_input.csv --manifest-out <manifest.json>
"""

from __future__ import annotations

import argparse
import json
import random
import sys

try:
    from rdkit import Chem
except ImportError:
    print("rdkit is required to author fixtures (validates parseability while generating)", file=sys.stderr)
    raise


def load_corpus(path, limit=None):
    smis = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                smis.append(s)
    if limit:
        smis = smis[:limit]
    return smis


def rdkit_parses(smi):
    try:
        return Chem.MolFromSmiles(smi) is not None
    except Exception:
        return False


def build_scenarios(corpus):
    random.seed(42)
    valid_corpus = [s for s in corpus if rdkit_parses(s)]
    scenarios = {}

    # --- general: space-delimited, title line, 2 extra property columns ---
    general_smis = valid_corpus[:150]
    rows = []
    for i, smi in enumerate(general_smis):
        rows.append(
            {
                "smiles": smi,
                "name": f"mol_{i}",
                "properties": {"Activity": f"{(i % 10) + 0.1:.2f}", "MW": f"{100 + i}.0"},
            }
        )
    scenarios["general_space_titleline"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": rows,
        "malformed_row_indices": [],
    }

    # --- tab-delimited, no collapsing of empty fields ---
    tab_rows = []
    for i, smi in enumerate(valid_corpus[150:165]):
        name = f"tab_{i}" if i % 3 != 0 else ""  # some rows have an empty name field
        tab_rows.append({"smiles": smi, "name": name, "properties": {"Note": f"n{i}"}})
    scenarios["tab_delimited"] = {
        "options": {"delimiter": "\t", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": tab_rows,
        "malformed_row_indices": [],
    }

    # --- comma CSV with quoted embedded commas/quotes ---
    csv_rows = []
    notes = [
        "plain note",
        "has, a comma",
        'has "quotes" inside',
        "has, both \"kinds\" of, trouble",
        "trailing comma,",
    ]
    for i, smi in enumerate(valid_corpus[165:180]):
        csv_rows.append(
            {"smiles": smi, "name": f"csv_{i}", "properties": {"Note": notes[i % len(notes)]}}
        )
    scenarios["csv_quoted"] = {
        "options": {"delimiter": ",", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": csv_rows,
        "malformed_row_indices": [],
    }

    # --- no title line: stable fallback column naming ---
    no_title_rows = []
    for i, smi in enumerate(valid_corpus[180:195]):
        no_title_rows.append(
            {"smiles": smi, "name": f"nt_{i}", "properties": {"__pos2__": f"v{i}"}}
        )
    scenarios["no_title_line"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": 1, "title_line": False},
        "rows": no_title_rows,
        "malformed_row_indices": [],
    }

    # --- no name column at all ---
    no_name_rows = []
    for i, smi in enumerate(valid_corpus[195:205]):
        no_name_rows.append({"smiles": smi, "name": "", "properties": {"__pos1__": f"e{i}"}})
    scenarios["no_name_column"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": None, "title_line": False},
        "rows": no_name_rows,
        "malformed_row_indices": [],
    }

    # --- malformed SMILES interspersed with valid ones ---
    malformed_smis = ["C1CC", "not(a(smiles", "[Xx]", "C(", ")C"]
    mal_rows = []
    mal_indices = []
    good_pool = iter(valid_corpus[205:215])
    for i in range(10):
        if i % 2 == 0:
            mal_rows.append({"smiles": malformed_smis[(i // 2) % len(malformed_smis)], "name": f"bad_{i}", "properties": {}})
            mal_indices.append(i)
        else:
            mal_rows.append({"smiles": next(good_pool), "name": f"good_{i}", "properties": {}})
    scenarios["malformed_smiles"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": mal_rows,
        "malformed_row_indices": mal_indices,
    }

    # --- blank lines and comment lines interspersed ---
    blank_comment_smis = valid_corpus[215:225]
    scenarios["blank_and_comment_lines"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": [{"smiles": s, "name": f"bc_{i}", "properties": {}} for i, s in enumerate(blank_comment_smis)],
        "malformed_row_indices": [],
        "inject_blanks_and_comments": True,
    }

    # --- isotopes, charges, stereochemistry, disconnected fragments ---
    special = [
        ("[13CH4]", "carbon_13"),
        ("CC(=O)[O-]", "acetate_anion"),
        ("[NH4+]", "ammonium"),
        ("C/C=C/C", "trans_2_butene"),
        ("C/C=C\\C", "cis_2_butene"),
        ("N[C@@H](C)C(=O)O", "l_alanine"),
        ("N[C@H](C)C(=O)O", "d_alanine"),
        ("[Na+].[Cl-]", "sodium_chloride"),
        ("CCO.CC(=O)O", "disconnected_ethanol_acetic_acid"),
        ("[2H]C([2H])([2H])O", "deuterated_methanol"),
    ]
    scenarios["special_chemistry"] = {
        "options": {"delimiter": " ", "smiles_column": 0, "name_column": 1, "title_line": True},
        "rows": [{"smiles": s, "name": n, "properties": {"Category": "special"}} for s, n in special],
        "malformed_row_indices": [],
    }

    return scenarios


def write_scenario_file(path, scenario, delimiter):
    opts = scenario["options"]
    inject = scenario.get("inject_blanks_and_comments", False)
    lines = []
    header_cols = ["SMILES"]
    if opts["name_column"] is not None:
        header_cols.append("Name")
    extra_keys = sorted({k for r in scenario["rows"] for k in r["properties"]})
    header_cols += extra_keys
    if opts["title_line"]:
        lines.append(delimiter.join(header_cols))

    for i, row in enumerate(scenario["rows"]):
        if inject and i == 2:
            lines.append("# a comment line")
        if inject and i == 5:
            lines.append("")
        fields = [row["smiles"]]
        if opts["name_column"] is not None:
            fields.append(row["name"])
        for k in extra_keys:
            v = row["properties"].get(k, "")
            if delimiter == "," and ("," in v or '"' in v):
                v = '"' + v.replace('"', '""') + '"'
            fields.append(v)
        lines.append(delimiter.join(fields))

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return header_cols, extra_keys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--corpus", required=True)
    p.add_argument("--manifest-out", required=True)
    p.add_argument("--corpus-limit", type=int, default=1000)
    args = p.parse_args()

    import os

    os.makedirs(args.out_dir, exist_ok=True)
    corpus = load_corpus(args.corpus, args.corpus_limit)
    scenarios = build_scenarios(corpus)

    manifest = {"scenarios": {}}
    total_rows = 0
    for name, scenario in scenarios.items():
        delim = scenario["options"]["delimiter"]
        path = os.path.join(args.out_dir, f"{name}.txt")
        header_cols, extra_keys = write_scenario_file(path, scenario, delim)
        manifest["scenarios"][name] = {
            "file": os.path.basename(path),
            "options": scenario["options"],
            "extra_property_keys": extra_keys,
            "rows": scenario["rows"],
            "malformed_row_indices": scenario["malformed_row_indices"],
        }
        total_rows += len(scenario["rows"])

    manifest["total_rows"] = total_rows
    with open(args.manifest_out, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"wrote {len(scenarios)} scenario files, {total_rows} total data rows, manifest -> {args.manifest_out}")


if __name__ == "__main__":
    main()
