#!/usr/bin/env python3
"""Depiction speed of one engine on a corpus: SVG, 2D coordinates, MOL block.

Run once per engine (``--engine chematic`` under the interpreter with the
chematic build to measure, ``--engine rdkit`` under one with RDKit); the
molecules are parsed first and only the depiction calls are timed, each the
minimum of ``--repeats`` passes over every row. The rows each pass writes
are counted and their output bytes summed, so two runs do the same work.

    chematic: Mol.svg(), Mol.depict_data(), Mol.to_mol_block()
    rdkit:    MolDraw2DSVG(300, 300) (draw + finish), rdDepictor.Compute2DCoords,
              Chem.MolToMolBlock on a copy without conformers (so it lays out)
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path


def rows(corpus: Path) -> list[str]:
    return [l.split()[0] for l in corpus.read_text(encoding="utf-8").splitlines() if l.strip()]


def chematic_ops(smiles):
    import chematic

    mols = []
    for s in smiles:
        try:
            mols.append(chematic.from_smiles(s))
        except ValueError:
            pass
    ops = {
        "svg": lambda m: len(m.svg()),
        "depict_data": lambda m: len(m.depict_data()["atoms"]),
        "mol_block": lambda m: len(m.to_mol_block()),
    }
    return mols, ops, {"chematic": chematic.__version__}


def rdkit_ops(smiles):
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    RDLogger.DisableLog("rdApp.*")
    mols = [m for m in (Chem.MolFromSmiles(s) for s in smiles) if m is not None]

    def svg(m):
        d = rdMolDraw2D.MolDraw2DSVG(300, 300)
        d.DrawMolecule(m)
        d.FinishDrawing()
        return len(d.GetDrawingText())

    def coords(m):
        rdDepictor.Compute2DCoords(m)
        return m.GetNumAtoms()

    def mol_block(m):
        # Without the coordinates the depict_data pass left, as chematic's
        # to_mol_block computes its layout.
        copy = Chem.Mol(m)
        copy.RemoveAllConformers()
        return len(Chem.MolToMolBlock(copy))

    ops = {"svg": svg, "depict_data": coords, "mol_block": mol_block}
    return mols, ops, {"rdkit": rdBase.rdkitVersion}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--engine", choices=("chematic", "rdkit"), required=True)
    ap.add_argument("--corpus", type=Path, required=True)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    smiles = rows(args.corpus)
    mols, ops, versions = (chematic_ops if args.engine == "chematic" else rdkit_ops)(smiles)
    result = {}
    for name, op in ops.items():
        best, size = float("inf"), 0
        for _ in range(args.repeats):
            t0 = time.perf_counter()
            size = sum(op(m) for m in mols)
            best = min(best, time.perf_counter() - t0)
        result[name] = {"ms": round(best * 1000, 1), "rows": len(mols), "output_size": size}
    report = {"schema": "depiction-speed/v1", "engine": args.engine, **versions,
              "corpus": args.corpus.name, "input_rows": len(smiles), "repeats": args.repeats,
              "python": platform.python_version(), "machine": platform.machine(), "ops": result}
    args.output.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(report["ops"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
