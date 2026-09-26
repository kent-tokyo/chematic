#!/usr/bin/env python3
"""Classify every CIP and SMARTS residual of an RDKit rebaseline rows file.

Input is the JSONL written by ``run_rdkit_python_chemistry_lane.py``. Every
differing CIP row and every differing SMARTS cell gets exactly one
root-cause family (issues #634 and #635); nothing is dropped. Families are
decided from atom identities, bond endpoints and ring facts — never from
engine-local bond indices — and a family never adopts the RDKit answer by
itself: it only records why the two engines differ.

CIP families (per differing label):

* ``phosphorus_oracle_unstable`` — CheMatic abstains with the typed reason
  ``oracle_unstable`` on a phosphorus centre (RDKit's labels for these flip
  under neutral Kekulé respellings; see the P boundary in docs/validation.md).
* ``trivalent_nitrogen_unsupported`` — RDKit labels a three-coordinate
  nitrogen (a bridgehead amine whose inversion is locked); CheMatic does not
  model a lone pair as a fourth ligand.
* ``adjudication_required`` — both engines label the centre and disagree;
  needs independent chemical adjudication before either answer is adopted.
* ``other`` — anything else (must be zero for a closed classification).

SMARTS families (per differing cell), for queries whose only difference is
ring semantics:

* ``ring_semantics:sssr_choice`` — the molecule's SSSR is not unique and the
  two engines select different rings of equal total size;
* ``ring_semantics:symmetrized_rings`` — RDKit's ring information holds more
  rings than the cycle rank (symmetrized SSSR), which ``[Rn]`` counts;
* ``ring_semantics:organometallic`` — rings through a metal atom;
* ``other`` — anything else (must be zero for a closed classification).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

RING_QUERY = re.compile(r"^(\[(R\d?|R0|k\d|x\d|r\d)\]|\*@\*|\*!@\*)$")
METALS = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 72, 73,
    74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
}


def open_rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ring_facts(smiles: str) -> dict[str, object]:
    import chematic
    from rdkit import Chem

    rd = Chem.MolFromSmiles(smiles)
    ch = chematic.from_smiles(smiles)
    rd_rings = sorted(len(r) for r in rd.GetRingInfo().AtomRings())
    cycle_rank = (
        rd.GetNumBonds() - rd.GetNumAtoms() + len(Chem.GetMolFrags(rd))
    )
    ch_rings = sorted(len(r) for r in ch.sssr_atom_rings)
    has_metal = any(a.GetAtomicNum() in METALS for a in rd.GetAtoms())
    return {
        "rdkit_ring_sizes": rd_rings,
        "chematic_sssr_sizes": ch_rings,
        "cycle_rank": cycle_rank,
        "has_metal": has_metal,
    }


def smarts_family(query: str, facts: dict[str, object]) -> str:
    if not RING_QUERY.match(query):
        return "other"
    if facts["has_metal"]:
        return "ring_semantics:organometallic"
    if len(facts["rdkit_ring_sizes"]) > facts["cycle_rank"]:
        return "ring_semantics:symmetrized_rings"
    if facts["rdkit_ring_sizes"] != facts["chematic_sssr_sizes"] or sum(
        facts["rdkit_ring_sizes"]
    ) == sum(facts["chematic_sssr_sizes"]):
        return "ring_semantics:sssr_choice"
    return "other"


def cip_families(row: dict) -> list[dict[str, object]]:
    from rdkit import Chem

    cip = row["cip"]
    rd_mol = Chem.MolFromSmiles(row["smiles"])
    out = []
    for key in sorted(set(cip["rdkit_atoms"]) | set(cip["chematic_atoms"]), key=int):
        rd_label = cip["rdkit_atoms"].get(key)
        ch_label = cip["chematic_atoms"].get(key)
        if rd_label == ch_label:
            continue
        atom = rd_mol.GetAtomWithIdx(int(key))
        reason = cip["chematic_unresolved"].get(key)
        if reason == "oracle_unstable" and atom.GetAtomicNum() == 15:
            family = "phosphorus_oracle_unstable"
        elif ch_label is None and atom.GetAtomicNum() == 7 and atom.GetDegree() == 3:
            family = "trivalent_nitrogen_unsupported"
        elif rd_label is not None and ch_label is not None:
            family = "adjudication_required"
        else:
            family = "other"
        out.append(
            {
                "kind": "atom",
                "atom": int(key),
                "element": atom.GetSymbol(),
                "rdkit": rd_label,
                "chematic": ch_label,
                "chematic_unresolved": reason,
                "family": family,
            }
        )
    for key in sorted(set(cip["rdkit_bonds"]) | set(cip["chematic_bonds"])):
        rd_label = cip["rdkit_bonds"].get(key)
        ch_label = cip["chematic_bonds"].get(key)
        if rd_label != ch_label:
            out.append(
                {
                    "kind": "bond",
                    "bond_endpoints": key,
                    "rdkit": rd_label,
                    "chematic": ch_label,
                    "family": "other",
                }
            )
    if not out and not all(cip["index_correspondence"].values()):
        out.append({"kind": "correspondence", "family": "other"})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import chematic
    from rdkit import rdBase

    cip_counts: Counter[str] = Counter()
    smarts_counts: Counter[str] = Counter()
    smarts_by_query: dict[str, Counter[str]] = defaultdict(Counter)
    cip_rows: list[dict[str, object]] = []
    smarts_rows: list[dict[str, object]] = []
    total_rows = 0
    total_cells = 0
    for row in open_rows(args.rows):
        total_rows += 1
        if row.get("status") != "completed":
            continue
        total_cells += row["smarts"]["query_count"]
        if not row["cip"]["exact"]:
            labels = cip_families(row)
            for item in labels:
                cip_counts[item["family"]] += 1
            cip_rows.append(
                {"input_index": row["input_index"], "smiles": row["smiles"], "labels": labels}
            )
        diffs = row["smarts"]["differences"]
        if diffs:
            facts = ring_facts(row["smiles"])
            cells = []
            for cell in diffs:
                family = smarts_family(cell["query"], facts)
                smarts_counts[family] += 1
                smarts_by_query[cell["query"]][family] += 1
                cells.append({"query": cell["query"], "family": family})
            smarts_rows.append(
                {
                    "input_index": row["input_index"],
                    "smiles": row["smiles"],
                    "ring_facts": facts,
                    "cells": cells,
                }
            )

    summary = {
        "schema_version": 1,
        "profile": "rdkit_rebaseline_residual_classification_v1",
        "rows_file": {"path": str(args.rows), "sha256": sha256(args.rows)},
        "chematic_version": getattr(chematic, "__version__", None),
        "rdkit_version": rdBase.rdkitVersion,
        "rows": total_rows,
        "smarts_cells": total_cells,
        "cip": {
            "differing_rows": len(cip_rows),
            "differing_labels_by_family": dict(sorted(cip_counts.items())),
            "unclassified_labels": cip_counts.get("other", 0),
            "rows": cip_rows,
        },
        "smarts": {
            "differing_cells": sum(smarts_counts.values()),
            "differing_rows": len(smarts_rows),
            "cells_by_family": dict(sorted(smarts_counts.items())),
            "cells_by_query": {
                query: dict(sorted(counts.items()))
                for query, counts in sorted(smarts_by_query.items())
            },
            "unclassified_cells": smarts_counts.get("other", 0),
            "rows": smarts_rows,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k]["differing_labels_by_family"] if k == "cip" else summary[k]["cells_by_family"] for k in ("cip", "smarts")}, indent=1))
    return 0 if summary["cip"]["unclassified_labels"] == 0 and summary["smarts"]["unclassified_cells"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
