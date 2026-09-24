#!/usr/bin/env python3
"""Row-level output agreement between chematic and RDKit, per corpus.

Counts, for each corpus row both libraries parse, whether chematic's output
equals RDKit's for a fixed set of RDKit-defined operations (bit-identical
fingerprints, descriptor values within a stated tolerance, identical
canonical scaffolds). No timing. Operations the installed chematic build does
not expose are recorded as unavailable rather than failed, so older wheels can
be measured with the same script.

Writes one JSON record with the chematic revision/wheel hash, the RDKit
version and each corpus's SHA-256.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone

import chematic
from rdkit import Chem, RDLogger
from rdkit.Chem import QED, Crippen, Descriptors, MACCSkeys, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator as rfg
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def on_bits(data: bytes, n: int) -> set[int]:
    return {i for i in range(n) if data[i // 8] >> (i % 8) & 1}


def canon(mol) -> str | None:
    if isinstance(mol, chematic.Mol):
        mol = Chem.MolFromSmiles(mol.smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None


def catalog(which) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(which)
    return FilterCatalog(params)


MORGAN = rfg.GetMorganGenerator(radius=2, fpSize=2048)
ATOM_PAIR = rfg.GetAtomPairGenerator(fpSize=2048)
TORSION = rfg.GetTopologicalTorsionGenerator(fpSize=2048)
PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
BRENK = catalog(FilterCatalogParams.FilterCatalogs.BRENK)


def maccs_equal(m, r) -> bool:
    data = m.maccs()
    # chematic stores MACCS key k at bit k-1; RDKit at bit k (bit 0 unused).
    return {i + 1 for i in on_bits(data, 166)} == set(MACCSkeys.GenMACCSKeys(r).GetOnBits()) - {0}


def close(tol):
    return lambda a, b: abs(float(a) - float(b)) <= tol


# name -> (chematic attribute or callable, RDKit callable, comparison, note)
OPERATIONS = {
    "morgan_r2_2048(rdkit-compatible)": (
        lambda m: on_bits(m.rdkit_ecfp4(), 2048), lambda r: set(MORGAN.GetFingerprint(r).GetOnBits()),
        lambda a, b: a == b, "bit-identical"),
    "atom_pair(rdkit-compatible)": (
        lambda m: on_bits(m.rdkit_atom_pair_fp(), 2048), lambda r: set(ATOM_PAIR.GetFingerprint(r).GetOnBits()),
        lambda a, b: a == b, "bit-identical"),
    "torsion(rdkit-compatible)": (
        lambda m: on_bits(m.rdkit_torsion_fp(), 2048), lambda r: set(TORSION.GetFingerprint(r).GetOnBits()),
        lambda a, b: a == b, "bit-identical"),
    "maccs": (None, None, None, "all 166 keys identical"),
    "logp": (lambda m: m.logp, Crippen.MolLogP, close(1e-6), "|diff| <= 1e-6"),
    "mr": (lambda m: m.molar_refractivity, Crippen.MolMR, close(1e-2), "|diff| <= 0.01"),
    "rdkit_tpsa": (lambda m: m.rdkit_tpsa, rdMolDescriptors.CalcTPSA, close(0.1), "|diff| <= 0.1; N/O-only TPSA"),
    "tpsa(S/P included)": (lambda m: m.tpsa, lambda r: rdMolDescriptors.CalcTPSA(r, includeSandP=True),
                           close(0.1), "|diff| <= 0.1 vs includeSandP=True"),
    "qed": (lambda m: m.qed, QED.qed, close(1e-3), "|diff| <= 1e-3"),
    "rdkit_mw": (lambda m: m.rdkit_mw, Descriptors.MolWt, close(1e-3), "|diff| <= 1e-3"),
    "hbd": (lambda m: m.hbd, rdMolDescriptors.CalcNumHBD, close(0), "exact"),
    "rotatable_bonds": (lambda m: m.rotatable_bonds, rdMolDescriptors.CalcNumRotatableBonds, close(0), "exact"),
    "murcko_scaffold": (lambda m: canon(m.scaffold()), lambda r: canon(MurckoScaffold.GetScaffoldForMol(r)),
                        lambda a, b: a is not None and a == b, "same RDKit canonical SMILES"),
    "pains_passes": (lambda m: m.pains_passes, lambda r: not PAINS.HasMatch(r), lambda a, b: a == b,
                     "vs RDKit PAINS FilterCatalog"),
    "brenk_passes": (lambda m: m.brenk_passes, lambda r: not BRENK.HasMatch(r), lambda a, b: a == b,
                     "vs RDKit BRENK FilterCatalog (chematic uses its own alert list)"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", action="append", required=True)
    ap.add_argument("--limit", type=int, default=5000)
    ap.add_argument("--output-json", required=True)
    ap.add_argument("--chematic-revision", required=True)
    ap.add_argument("--chematic-artifact", default=None)
    args = ap.parse_args()

    probe = chematic.from_smiles("CCO")
    available = {}
    for name, (c_fn, *_rest) in OPERATIONS.items():
        if name == "maccs":
            available[name] = hasattr(probe, "maccs")
            continue
        try:
            c_fn(probe)
            available[name] = True
        except AttributeError:
            available[name] = False

    record = {
        "schema": "chematic-rdkit-output-agreement/v1",
        "created_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "chematic": {
            "version": chematic.__version__,
            "revision": args.chematic_revision,
            "artifact": os.path.basename(args.chematic_artifact) if args.chematic_artifact else None,
            "artifact_sha256": sha256_file(args.chematic_artifact) if args.chematic_artifact else None,
        },
        "rdkit": {"version": Chem.rdBase.rdkitVersion},
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "operations": {name: {"note": spec[3], "available": available[name]} for name, spec in OPERATIONS.items()},
        "corpora": [],
    }
    for path in args.corpus:
        counts = {name: 0 for name in OPERATIONS if available[name]}
        rows = 0
        for line in open(path, encoding="utf-8"):
            if rows >= args.limit:
                break
            if not line.strip():
                continue
            smiles = line.split()[0]
            rdmol = Chem.MolFromSmiles(smiles)
            if rdmol is None:
                continue
            try:
                mol = chematic.from_smiles(smiles)
            except Exception:
                continue
            rows += 1
            for name in counts:
                c_fn, r_fn, cmp, _ = OPERATIONS[name]
                try:
                    ok = maccs_equal(mol, rdmol) if name == "maccs" else cmp(c_fn(mol), r_fn(rdmol))
                except Exception:
                    ok = False
                counts[name] += bool(ok)
        entry = {
            "path": os.path.relpath(os.path.abspath(path), ROOT),
            "sha256": sha256_file(path),
            "rows_compared": rows,
            "agree": counts,
        }
        record["corpora"].append(entry)
        print(entry["path"], rows, json.dumps(counts), flush=True)
    with open(args.output_json, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1)
        fh.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
