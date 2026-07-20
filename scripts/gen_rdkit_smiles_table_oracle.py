#!/usr/bin/env python3
"""Real RDKit oracle for the IO-1 (SMILES table) acceptance gate.

Runs RDKit's actual `Chem.SmilesMolSupplier` over every scenario file named
in `gen_smiles_table_fixtures.py`'s manifest, using that scenario's own
configured options, and dumps each row's extracted
`(status, name, properties, rdkit_canonical_smiles)` as JSONL, in the same
row order as `smiles_table_dump.rs`'s chematic-side dump.

Also computes, for every successfully-parsed row, a *self-consistency*
check against the manifest's own known ground-truth SMILES for that row:
does `Chem.MolToSmiles(Chem.MolFromSmiles(known_smiles))` equal
`Chem.MolToSmiles(mol_from_supplier)`? This validates RDKit's own
tokenization extracted the right substring, using ONLY RDKit's own
canonicalizer -- never compared against chematic's canonical output (see
`scripts/smiles_table_io_parity.py`'s module docs for why: chematic/RDKit
canonical-SMILES divergence is a known, separately-tracked, unrelated
issue -- this oracle is not evidence about it, and must not be misread as
such).

Usage:
    python scripts/gen_rdkit_smiles_table_oracle.py --manifest <manifest.json> \\
        --fixtures-dir <dir> --out <out.jsonl>
"""

from __future__ import annotations

import argparse
import json

from rdkit import Chem


def rdkit_style_split(line, delimiter):
    """RDKit's own tokenization rule (see the IO-1 source audit): a
    space/tab delimiter string is a CHARACTER CLASS (any char in it is an
    independent separator, `boost::keep_empty_tokens` -- consecutive
    delimiters yield empty tokens), a comma or other single character is
    a literal split with no quoting logic in RDKit's own tokenizer (RDKit's
    SmilesMolSupplier has no CSV-quote awareness at all)."""
    if all(c in " \t" for c in delimiter):
        out, cur = [], ""
        for ch in line:
            if ch in " \t":
                out.append(cur)
                cur = ""
            else:
                cur += ch
        out.append(cur)
        return out
    return line.split(delimiter)


def run_scenario(name, scenario, fixtures_dir):
    opts = scenario["options"]
    delimiter = opts["delimiter"]
    smiles_column = opts["smiles_column"]
    name_column = opts["name_column"]
    title_line = opts["title_line"]

    path = f"{fixtures_dir}/{scenario['file']}"

    # RDKit's SmilesMolSupplier delimiter kwarg: pass through verbatim; a
    # multi-char whitespace-only class (unused here, kept single-char to
    # avoid RDKit's own cross-entry-point delimiter-default inconsistency
    # noted in the IO-1 source audit) works the same as a plain " ".
    sup = Chem.SmilesMolSupplier(
        path,
        delimiter=delimiter,
        smilesColumn=smiles_column,
        nameColumn=name_column if name_column is not None else -1,
        titleLine=title_line,
        sanitize=True,
    )

    known_rows = scenario["rows"]
    rows_out = []
    for row_index, mol in enumerate(sup):
        if mol is None:
            rows_out.append({"scenario": name, "row_index": row_index, "status": "error", "error": "rdkit_returned_none"})
            continue

        # Dump RDKit's raw property dict as-is, under whatever key RDKit
        # itself used (the manifest's own authored key when a title line
        # named the column, "Column_N" 0-indexed-over-the-full-row when not
        # -- reconciling the two tools' differing fallback-name CONVENTIONS
        # is the comparator's job, done once in one place rather than
        # duplicated across both oracle-generating scripts). `GetProp` (not
        # `GetPropsAsDict`, which opportunistically coerces numeric-looking
        # strings to float/int) is used deliberately -- chematic's own
        # property store never does that type inference (source-confirmed:
        # SmilesMolSupplier stores every extra column as a plain string),
        # so comparing against RDKit's *raw* string avoids a false mismatch
        # that's really just a Python-API-layer type-coercion difference.
        props = [
            [k, mol.GetProp(k)]
            for k in mol.GetPropNames(includePrivate=False, includeComputed=False)
            if k != "_Name" and not k.startswith("__")
        ]
        props.sort()

        rdkit_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        rdkit_canonical = Chem.MolToSmiles(mol)

        self_consistent = None
        if row_index < len(known_rows):
            known_smi = known_rows[row_index]["smiles"]
            known_mol = Chem.MolFromSmiles(known_smi)
            if known_mol is not None:
                self_consistent = Chem.MolToSmiles(known_mol) == rdkit_canonical

        rows_out.append(
            {
                "scenario": name,
                "row_index": row_index,
                "status": "success",
                "name": rdkit_name,
                "properties": props,
                "rdkit_canonical_smiles": rdkit_canonical,
                "self_consistent_with_known_smiles": self_consistent,
            }
        )
    return rows_out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--fixtures-dir", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with open(args.manifest) as f:
        manifest = json.load(f)

    total = 0
    with open(args.out, "w") as out:
        for name in sorted(manifest["scenarios"].keys()):
            scenario = manifest["scenarios"][name]
            for row in run_scenario(name, scenario, args.fixtures_dir):
                out.write(json.dumps(row) + "\n")
                total += 1

    print(f"total_rows={total} rdkit_version={Chem.rdBase.rdkitVersion} out={args.out}")


if __name__ == "__main__":
    main()
