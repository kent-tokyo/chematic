#!/usr/bin/env python3
"""2D layout quality of chematic's depiction coordinates against RDKit's.

For each corpus row, both chematic (``Mol.depict_data()``, the layout the
MOL writer and SVG use) and RDKit (``rdDepictor.Compute2DCoords``) lay the
molecule out; each layout is scaled to a median bond length of 1 and scored:

* ``clashes``: non-bonded atom pairs closer than 0.4 bond lengths;
* ``crossings``: pairs of bonds without a common atom that intersect;
* ``bond_dev``: the largest |bond length - 1|.

A row is *clean* with no clash and no crossing. Rows are split by whether
the molecule has bridgehead atoms (bridged or cage ring systems).
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdDepictor, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def score(points: list[tuple[float, float]], bonds: list[tuple[int, int]]) -> dict:
    lengths = [((points[a][0] - points[b][0]) ** 2 + (points[a][1] - points[b][1]) ** 2) ** 0.5 for a, b in bonds]
    if not lengths:
        return {"clashes": 0, "crossings": 0, "bond_dev": 0.0}
    unit = statistics.median(lengths) or 1.0
    pts = [(x / unit, y / unit) for x, y in points]
    bonded = {frozenset(b) for b in bonds}
    clashes = 0
    n = len(pts)
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset((i, j)) in bonded:
                continue
            if ((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2) ** 0.5 < 0.4:
                clashes += 1

    def cross(p, q, r, s):
        def orient(a, b, c):
            v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            return 0 if abs(v) < 1e-9 else (1 if v > 0 else -1)
        o1, o2, o3, o4 = orient(p, q, r), orient(p, q, s), orient(r, s, p), orient(r, s, q)
        return o1 * o2 < 0 and o3 * o4 < 0

    crossings = 0
    for k in range(len(bonds)):
        a, b = bonds[k]
        for m in range(k + 1, len(bonds)):
            c, d = bonds[m]
            if len({a, b, c, d}) < 4:
                continue
            if cross(pts[a], pts[b], pts[c], pts[d]):
                crossings += 1
    dev = max(abs(l / unit - 1.0) for l in lengths)
    return {"clashes": clashes, "crossings": crossings, "bond_dev": dev}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--corpus", type=Path, nargs="+", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--rows-output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for corpus in args.corpus:
        for i, line in enumerate(corpus.read_text().splitlines()):
            smi = line.split()[0] if line.strip() else ""
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            try:
                cm = chematic.from_smiles(smi)
            except Exception:
                continue
            d = cm.depict_data()
            if len(d["atoms"]) != m.GetNumAtoms():
                continue
            cpts = [(a["x"], a["y"]) for a in d["atoms"]]
            cbonds = [(b["atom1"], b["atom2"]) for b in d["bonds"]]
            rdDepictor.Compute2DCoords(m)
            conf = m.GetConformer()
            rpts = [(conf.GetAtomPosition(k).x, conf.GetAtomPosition(k).y) for k in range(m.GetNumAtoms())]
            rbonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
            rows.append({
                "corpus": corpus.name, "row": i, "smiles": smi,
                "bridged": rdMolDescriptors.CalcNumBridgeheadAtoms(m) > 0,
                "chematic": score(cpts, cbonds), "rdkit": score(rpts, rbonds),
            })
    args.rows_output.write_text("".join(json.dumps(r) + "\n" for r in rows))

    def tally(subset):
        out = {"rows": len(subset)}
        for engine in ("chematic", "rdkit"):
            s = [r[engine] for r in subset]
            out[engine] = {
                "clean": sum(x["clashes"] == 0 and x["crossings"] == 0 for x in s),
                "with_clash": sum(x["clashes"] > 0 for x in s),
                "with_crossing": sum(x["crossings"] > 0 for x in s),
                "bond_dev_over_0.2": sum(x["bond_dev"] > 0.2 for x in s),
            }
        return out

    report = {
        "schema": "layout-2d-quality-vs-rdkit/v1",
        "rdkit": rdBase.rdkitVersion,
        "chematic": chematic.__version__,
        "all": tally(rows),
        "bridged": tally([r for r in rows if r["bridged"]]),
        "not_bridged": tally([r for r in rows if not r["bridged"]]),
    }
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
