#!/usr/bin/env python3
"""Dump chematic's answers for the #734/#754 chemistry checks, one JSON line
per corpus row, without importing RDKit.

Published chematic Linux wheels are CPython 3.9 builds and RDKit 2026.03.6
has no CPython 3.9 wheel, so the published artifact is run here and
``compare_chemistry_dump_rdkit.py`` compares the dump with RDKit in another
interpreter. Per row: accurate CIP labels and abstentions, RDKit-model
hybridization, MMFF94 numeric types (after ``add_hydrogens``), the V2000 MOL
block and whether ``to_mol_block(strict=True)`` refuses it.
Python 3.9 syntax, standard library and chematic only.
"""

import argparse
import hashlib
import json
import sys

import chematic


def row(smiles):
    out = {"smiles": smiles}
    try:
        mol = chematic.from_smiles(smiles)
    except ValueError as exc:
        out["error"] = str(exc)
        return out
    try:
        out["cip"] = [
            {k: d[k] for k in ("atom_idx", "descriptor", "kekule_dependent") if k in d}
            for d in mol.cip_stereo(mode="accurate")
            if "bond_idx" not in d
        ]
        out["cip_unresolved"] = mol.cip_stereo_unresolved()
    except ValueError as exc:
        out["cip_error"] = str(exc)
    out["hybridization"] = list(mol.hybridization_per_atom())
    try:
        out["mmff"] = list(mol.add_hydrogens().mmff94_numeric_atom_types())
    except Exception as exc:  # noqa: BLE001 - a typed refusal is recorded
        out["mmff_error"] = str(exc)
    if any(c in smiles for c in "@/\\"):
        out["molblock"] = mol.to_mol_block()
        try:
            mol.to_mol_block(strict=True)
            out["mol_strict_refused"] = False
        except ValueError:
            out["mol_strict_refused"] = True
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    data = open(args.corpus, "rb").read()
    lines = [l.split()[0] for l in data.decode().splitlines() if l.strip()]
    with open(args.output, "w") as out:
        out.write(json.dumps({
            "chematic": chematic.__version__,
            "module": chematic.__file__,
            "python": sys.version.split()[0],
            "corpus": args.corpus,
            "corpus_sha256": hashlib.sha256(data).hexdigest(),
        }) + "\n")
        for smiles in lines:
            out.write(json.dumps(row(smiles)) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
