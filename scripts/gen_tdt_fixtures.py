#!/usr/bin/env python3
"""Generate TDT (Daylight Tagged Data) I/O fixtures for the IO-2 oracle gate.

Produces scenario `.tdt` files + a manifest of KNOWN ground truth this
script itself authors, covering: SMI record, NAME record, arbitrary
properties, `|` terminator, empty property, repeated property, unknown tag,
malformed tag (recovery), EOF mid-record, 2D/3D coordinate tags, isotopes,
charges, stereochemistry, disconnected fragments.

The "general" scenario draws real, chemically diverse SMILES from the
Morgan M4-A0 corpus rather than hand-authoring hundreds of molecules.

Usage:
    python scripts/gen_tdt_fixtures.py --out-dir <dir> --corpus <SMILES.csv> \\
        --manifest-out <manifest.json>
"""

from __future__ import annotations

import argparse
import json
import sys

try:
    from rdkit import Chem
except ImportError:
    print("rdkit is required to author fixtures", file=sys.stderr)
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
    valid_corpus = [s for s in corpus if rdkit_parses(s)]
    scenarios = {}

    general = valid_corpus[:170]
    scenarios["general"] = {
        "rows": [
            {"smiles": s, "name": f"mol_{i}", "properties": {"ACTIVITY": f"{(i % 10) + 0.1:.2f}"}}
            for i, s in enumerate(general)
        ],
        "malformed_row_indices": [],
    }

    empty_prop = valid_corpus[170:180]
    scenarios["empty_property"] = {
        "rows": [{"smiles": s, "name": f"ep_{i}", "properties": {"NOTE": ""}} for i, s in enumerate(empty_prop)],
        "malformed_row_indices": [],
    }

    unknown_tags = valid_corpus[180:190]
    scenarios["unknown_tags"] = {
        "rows": [
            {"smiles": s, "name": f"ut_{i}", "properties": {"MFCD": f"{1000+i}", "CAS": f"{64+i}-17-5"}}
            for i, s in enumerate(unknown_tags)
        ],
        "malformed_row_indices": [],
    }

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
        "rows": [{"smiles": s, "name": n, "properties": {"CATEGORY": "special"}} for s, n in special],
        "malformed_row_indices": [],
    }

    return scenarios


def write_tdt_scenario(path, scenario):
    lines = []
    for row in scenario["rows"]:
        lines.append(f"$SMI<{row['smiles']}>")
        lines.append(f"NAME<{row['name']}>")
        for k, v in row["properties"].items():
            lines.append(f"{k}<{v}>")
        lines.append("|")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_repeated_tag_fixture(path):
    # FOO appears twice -- last value wins, same position.
    content = "$SMI<CC>\nNAME<ethane>\nFOO<first>\nBAR<x>\nFOO<second>\n|\n"
    with open(path, "w") as f:
        f.write(content)


def write_malformed_tag_fixture(path):
    # Row 0: malformed (missing closing '>'). Row 1: valid. Tests recovery.
    content = "$SMI<CC>\nBROKEN<no_close\n|\n$SMI<CCO>\nNAME<ethanol>\n|\n"
    with open(path, "w") as f:
        f.write(content)


def write_eof_mid_record_fixture(path):
    # No trailing '|' at all -- EOF mid-record.
    content = "$SMI<CC>\nNAME<ethane>"
    with open(path, "w") as f:
        f.write(content)


def write_coordinate_fixture(path):
    # 3-atom molecule (CCO), 2D and 3D coordinate tags -- exercises the
    # last-atom-drop bug real RDKit has and chematic fixes.
    content = "$SMI<CCO>\nNAME<ethanol_with_coords>\n2D<0.0,0.0,1.0,0.0,2.0,1.0;>\n3D<0.0,0.0,0.0,1.0,0.0,0.0,2.0,1.0,0.5;>\n|\n"
    with open(path, "w") as f:
        f.write(content)


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
        path = os.path.join(args.out_dir, f"{name}.tdt")
        write_tdt_scenario(path, scenario)
        manifest["scenarios"][name] = {
            "file": os.path.basename(path),
            "rows": scenario["rows"],
            "malformed_row_indices": scenario["malformed_row_indices"],
        }
        total_rows += len(scenario["rows"])

    # Special single-file scenarios, not corpus-derived, each with its own
    # exact ground truth described inline (not enumerable as "rows").
    repeated_path = os.path.join(args.out_dir, "repeated_tag.tdt")
    write_repeated_tag_fixture(repeated_path)
    manifest["scenarios"]["repeated_tag"] = {
        "file": "repeated_tag.tdt",
        "rows": [{"smiles": "CC", "name": "ethane", "properties": {"FOO": "second", "BAR": "x"}}],
        "malformed_row_indices": [],
    }
    total_rows += 1

    malformed_path = os.path.join(args.out_dir, "malformed_tag.tdt")
    write_malformed_tag_fixture(malformed_path)
    manifest["scenarios"]["malformed_tag"] = {
        "file": "malformed_tag.tdt",
        "rows": [
            {"smiles": None, "name": None, "properties": {}},  # row 0: malformed, no ground truth
            {"smiles": "CCO", "name": "ethanol", "properties": {}},
        ],
        "malformed_row_indices": [0],
    }
    total_rows += 2

    eof_path = os.path.join(args.out_dir, "eof_mid_record.tdt")
    write_eof_mid_record_fixture(eof_path)
    manifest["scenarios"]["eof_mid_record"] = {
        "file": "eof_mid_record.tdt",
        "rows": [{"smiles": "CC", "name": "ethane", "properties": {}}],
        "malformed_row_indices": [],
    }
    total_rows += 1

    coord_path = os.path.join(args.out_dir, "coordinates.tdt")
    write_coordinate_fixture(coord_path)
    manifest["scenarios"]["coordinates"] = {
        "file": "coordinates.tdt",
        "rows": [{"smiles": "CCO", "name": "ethanol_with_coords", "properties": {}}],
        "malformed_row_indices": [],
        "known_divergent": "RDKit's own coordinate-list parser drops the last atom's position (confirmed bug); chematic fixes it. Coordinate values are NOT gated in the comparator for this scenario.",
    }
    total_rows += 1

    manifest["total_rows"] = total_rows
    with open(args.manifest_out, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"wrote {len(manifest['scenarios'])} scenario files, {total_rows} total rows, manifest -> {args.manifest_out}")


if __name__ == "__main__":
    main()
