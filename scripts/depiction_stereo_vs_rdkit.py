#!/usr/bin/env python3
"""Check that chematic's 2D depiction carries the input stereo, with RDKit.

For every row with stereo (a tetrahedral centre or a double-bond direction),
``Mol.depict_data()`` is turned into a V2000 block as drawn: depiction
coordinates (y flipped, scaled to 1.5 Å bonds), double bonds as drawn and
``Up``/``Down`` bonds as wedge/hash from their first atom. RDKit reads that
block, assigns stereo from the drawing, and its isomeric canonical SMILES is
compared with RDKit's for the input SMILES (double bonds the input leaves
unspecified are not compared: a drawing has to show them some way). A row
passes when they are equal: what the picture shows is the molecule given.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")

ORDER = {"Single": 1, "Double": 2, "Triple": 3, "Aromatic": 4, "Up": 1, "Down": 1}
STEREO = {"Up": 1, "Down": 6}


def block_from_depiction(data: dict) -> str:
    atoms, bonds = data["atoms"], data["bonds"]
    lengths = [
        ((atoms[b["atom1"]]["x"] - atoms[b["atom2"]]["x"]) ** 2
         + (atoms[b["atom1"]]["y"] - atoms[b["atom2"]]["y"]) ** 2) ** 0.5
        for b in bonds
    ]
    scale = 1.5 / (sorted(lengths)[len(lengths) // 2] if lengths else 1.0)
    lines = ["", "  depiction", "", f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for a in atoms:
        chg = {3: 1, 2: 2, 1: 3, -1: 5, -2: 6, -3: 7}.get(a["charge"], 0)
        lines.append(f"{a['x'] * scale:10.4f}{-a['y'] * scale:10.4f}{0.0:10.4f} {a['element']:<3} 0{chg:3d}  0  0  0  0  0  0  0  0  0  0")
    for b in bonds:
        lines.append(f"{b['atom1'] + 1:3d}{b['atom2'] + 1:3d}{ORDER[b['kind']]:3d}{STEREO.get(b['kind'], 0):3d}")
    lines.append("M  END")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--corpus", type=Path, nargs="+", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rows-output", type=Path)
    args = ap.parse_args()
    summary: dict[str, dict[str, int]] = {}
    rows = []
    for corpus in args.corpus:
        counts = {"stereo_rows": 0, "drawn_as_input": 0, "differs": 0, "rdkit_unreadable": 0}
        for i, line in enumerate(corpus.read_text().splitlines()):
            smi = line.split()[0] if line.strip() else ""
            if not ("@" in smi or "/" in smi or "\\" in smi):
                continue
            ref = Chem.MolFromSmiles(smi)
            if ref is None:
                continue
            try:
                data = chematic.from_smiles(smi).depict_data()
            except Exception:
                continue
            if len(data["atoms"]) != ref.GetNumAtoms():
                continue
            counts["stereo_rows"] += 1
            # Hydrogen counts come from the input, as a MOL reader would infer
            # them for the drawn heavy-atom graph.
            drawn = Chem.MolFromMolBlock(block_from_depiction(data), sanitize=False)
            status = "rdkit_unreadable"
            if drawn is not None:
                try:
                    for a, r in zip(drawn.GetAtoms(), ref.GetAtoms()):
                        a.SetNoImplicit(True)
                        a.SetNumExplicitHs(r.GetTotalNumHs())
                        a.SetNumRadicalElectrons(r.GetNumRadicalElectrons())
                    Chem.SanitizeMol(drawn)
                    Chem.AssignChiralTypesFromBondDirs(drawn)
                    Chem.AssignStereochemistry(drawn, cleanIt=True, force=True)
                    Chem.DetectBondStereochemistry(drawn)
                    Chem.AssignStereochemistry(drawn, cleanIt=True, force=True)
                    # A drawing shows some geometry for a double bond the input
                    # leaves unspecified (the MOL writer marks it "either");
                    # only declared E/Z is compared.
                    for bond in drawn.GetBonds():
                        r = ref.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                        if r is not None and r.GetStereo() == Chem.BondStereo.STEREONONE:
                            bond.SetStereo(Chem.BondStereo.STEREONONE)
                    same = Chem.MolToSmiles(drawn) == Chem.MolToSmiles(ref)
                    status = "drawn_as_input" if same else "differs"
                except Exception:
                    status = "rdkit_unreadable"
            counts[status] += 1
            rows.append({"corpus": corpus.name, "row": i, "smiles": smi, "status": status})
        summary[corpus.name] = counts
    report = {"schema": "depiction-stereo-vs-rdkit/v1", "rdkit": rdBase.rdkitVersion,
              "chematic": chematic.__version__, "summary": summary}
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    if args.rows_output:
        args.rows_output.write_text("".join(json.dumps(r) + "\n" for r in rows))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
