#!/usr/bin/env python3
"""Real RDKit oracle for the IO-2 (TDT) acceptance gate.

Runs RDKit's actual `Chem.TDTMolSupplier` over every scenario file named in
`gen_tdt_fixtures.py`'s manifest, using `nameRecord="NAME"` explicitly
(matching what the fixtures actually write, and what chematic's own
TdtReaderOptions default now recognizes -- see `chematic_mol::tdt`'s module
doc comment for why RDKit's OWN default of `nameRecord=""` would silently
fail to recognize it). Dumps each row's extracted
`(status, name, properties, coordinates)` as JSONL, same row order as
`tdt_dump.rs`'s chematic-side dump.

Also computes each successfully-parsed row's self-consistency against the
manifest's own known ground-truth SMILES -- see
`scripts/smiles_table_io_parity.py`'s module docs for why this project never
compares chematic-canonical vs. RDKit-canonical SMILES directly.

Usage:
    python scripts/gen_rdkit_tdt_oracle.py --manifest <manifest.json> \\
        --fixtures-dir <dir> --out <out.jsonl>
"""

from __future__ import annotations

import argparse
import json

from rdkit import Chem


def run_scenario(name, scenario, fixtures_dir):
    path = f"{fixtures_dir}/{scenario['file']}"
    read_coords = name == "coordinates"

    # Distinct conformer IDs for 2D vs 3D -- using the SAME id for both would
    # make RDKit's own addConformer(id=0, assignId=False) calls collide.
    sup = Chem.TDTMolSupplier(
        path,
        nameRecord="NAME",
        confId2D=0 if read_coords else -1,
        confId3D=1 if read_coords else -1,
        sanitize=True,
    )

    known_rows = scenario["rows"]
    rows_out = []
    # Index-based access, not a plain `for mol in sup` loop: RDKit's own
    # TDTMolSupplier has a confirmed recovery-hazard bug (see
    # chematic_mol::tdt's module doc comment and this session's source
    # audit) -- a malformed generic tag throws an exception that is NOT
    # caught inside next()'s own position-advance bookkeeping, so naively
    # retrying .next()/iterating re-throws on the SAME record indefinitely.
    # Explicit index access (sup[idx]) sidesteps this: a failure at index i
    # never blocks reaching index i+1.
    for row_index in range(len(known_rows)):
        try:
            mol = sup[row_index]
        except Exception as e:
            rows_out.append({"scenario": name, "row_index": row_index, "status": "error", "error": str(e)})
            continue
        if mol is None:
            rows_out.append({"scenario": name, "row_index": row_index, "status": "error", "error": "rdkit_returned_none"})
            continue

        props = [
            [k, mol.GetProp(k)]
            for k in mol.GetPropNames(includePrivate=False, includeComputed=False)
            if k != "_Name"
        ]
        props.sort()

        rdkit_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""

        coords_2d = None
        coords_3d = None
        if read_coords:
            try:
                conf2d = mol.GetConformer(0)
                coords_2d = [[conf2d.GetAtomPosition(i).x, conf2d.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())]
            except ValueError:
                pass
            try:
                conf3d = mol.GetConformer(1)
                coords_3d = [
                    [conf3d.GetAtomPosition(i).x, conf3d.GetAtomPosition(i).y, conf3d.GetAtomPosition(i).z]
                    for i in range(mol.GetNumAtoms())
                ]
            except ValueError:
                pass

        self_consistent = None
        if row_index < len(known_rows) and known_rows[row_index].get("smiles"):
            known_mol = Chem.MolFromSmiles(known_rows[row_index]["smiles"])
            if known_mol is not None:
                self_consistent = Chem.MolToSmiles(known_mol) == Chem.MolToSmiles(mol)

        rows_out.append(
            {
                "scenario": name,
                "row_index": row_index,
                "status": "success",
                "name": rdkit_name,
                "properties": props,
                "coordinates_2d": coords_2d,
                "coordinates_3d": coords_3d,
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
