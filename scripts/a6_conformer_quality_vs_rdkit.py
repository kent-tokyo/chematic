#!/usr/bin/env python3
"""Score chematic's A6 conformers against RDKit's own conformers, with RDKit.

For every successful row of a chematic A6 run (heavy-atom coordinates, in
the SMILES atom order), RDKit 2026.03.6:

* adds hydrogens to the chematic geometry (``AddHs(addCoords=True)``),
  relaxes only the hydrogens under MMFF94 (heavy atoms fixed), and records
  that energy (``e_fixed``);
* then relaxes everything to MMFF94 convergence (``e_relaxed``) and records
  the heavy-atom RMSD from the chematic geometry to that minimum
  (``rmsd_to_own_minimum``: how far the returned geometry is from the
  nearest MMFF94 minimum);
* embeds ``--n-conformers`` ETKDGv3 conformers (``randomSeed``), relaxes
  each to convergence and records the lowest energy (``e_rdkit_best``), the
  first conformer's energy (``e_rdkit_first``, RDKit's one-conformer
  default) and the symmetry-aware heavy-atom RMSD from the chematic
  geometry to the closest RDKit conformer.

Energies are RDKit's MMFF94 with its defaults (``ignoreInterfragInteractions``).
This scores conformer quality independently of chematic's own force field;
it is not a speed comparison.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import statistics
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")
MAX_ITERS = 20_000


def minimize(mol: Chem.Mol, conf_id: int = -1, fixed: list[int] | None = None) -> tuple[float, int] | None:
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    if ff is None:
        return None
    for idx in fixed or []:
        ff.AddFixedPoint(idx)
    ff.Initialize()
    try:
        rc = ff.Minimize(maxIts=MAX_ITERS)
    except RuntimeError:
        # RDKit's BFGS asserts on a non-descent direction (e.g. nothing left
        # to move); the geometry it stopped at is kept.
        rc = -1
    return ff.CalcEnergy(), rc


def score(task: tuple[dict, int, int]) -> dict:
    row, n_conf, seed = task
    out = {"name": row["name"], "row_index": row.get("row_index"), "tier": row.get("tier")}
    mol = Chem.MolFromSmiles(row["smiles"])
    coords = row.get("coords")
    if mol is None or not coords or len(coords) != mol.GetNumAtoms():
        out["status"] = "not_comparable"
        return out
    heavy = mol.GetNumAtoms()
    conf = Chem.Conformer(heavy)
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf, assignId=True)
    chem = Chem.AddHs(mol, addCoords=True)
    if AllChem.MMFFGetMoleculeProperties(chem) is None:
        out["status"] = "rdkit_mmff_unsupported"
        return out
    fixed = minimize(chem, fixed=list(range(heavy)))
    if fixed is None:
        out["status"] = "rdkit_mmff_unsupported"
        return out
    heavy_chem = Chem.Mol(chem)
    relaxed = minimize(chem)
    out["e_fixed"], _ = fixed
    out["e_relaxed"], out["relaxed_rc"] = relaxed
    rmsd_own = rdMolAlign.AlignMol(Chem.RemoveHs(heavy_chem), Chem.RemoveHs(chem))
    out["rmsd_to_own_minimum"] = rmsd_own

    ref = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    ids = list(AllChem.EmbedMultipleConfs(ref, numConfs=n_conf, params=params))
    if not ids:
        out["status"] = "rdkit_embed_failed"
        return out
    energies = []
    for cid in ids:
        res = minimize(ref, conf_id=cid)
        energies.append((res[0] if res else float("inf"), cid))
    out["e_rdkit_first"] = energies[0][0]
    out["e_rdkit_best"] = min(e for e, _ in energies)
    ref_heavy = Chem.RemoveHs(ref)
    probe = Chem.RemoveHs(heavy_chem)
    out["rmsd_to_closest_rdkit"] = min(
        rdMolAlign.GetBestRMS(probe, ref_heavy, prbId=0, refId=cid, maxMatches=2000) for cid in ids
    )
    out["status"] = "scored"
    return out


def summary(rows: list[dict]) -> dict:
    scored = [r for r in rows if r.get("status") == "scored"]
    def q(values: list[float]) -> dict:
        values = sorted(values)
        if not values:
            return {}
        pick = lambda f: values[min(len(values) - 1, int(f * (len(values) - 1)))]
        return {"median": pick(0.5), "p90": pick(0.9), "max": values[-1]}
    d_relaxed = [r["e_relaxed"] - r["e_rdkit_best"] for r in scored]
    d_first = [r["e_rdkit_first"] - r["e_rdkit_best"] for r in scored]
    return {
        "rows": len(rows),
        "scored": len(scored),
        "statuses": {s: sum(r.get("status") == s for r in rows) for s in sorted({r.get("status") for r in rows})},
        "chematic_relaxed_minus_rdkit_best_kcal_mol": q(d_relaxed),
        "chematic_relaxed_within_1_kcal_mol_of_rdkit_best": sum(d <= 1.0 for d in d_relaxed),
        "chematic_relaxed_at_or_below_rdkit_best": sum(d <= 1e-3 for d in d_relaxed),
        "rdkit_first_minus_rdkit_best_kcal_mol": q(d_first),
        "rdkit_first_within_1_kcal_mol_of_rdkit_best": sum(d <= 1.0 for d in d_first),
        "chematic_fixed_minus_relaxed_kcal_mol": q([r["e_fixed"] - r["e_relaxed"] for r in scored]),
        "rmsd_to_own_minimum_angstrom": q([r["rmsd_to_own_minimum"] for r in scored]),
        "rmsd_to_closest_rdkit_conformer_angstrom": q([r["rmsd_to_closest_rdkit"] for r in scored]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--chematic-rows", type=Path, required=True)
    ap.add_argument("--arm", default="chematic_pipeline_v2_mmff94_strict_stereo_safe")
    ap.add_argument("--n-conformers", type=int, default=10)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True, help="summary JSON")
    ap.add_argument("--rows-output", type=Path, required=True, help="per-molecule JSONL")
    args = ap.parse_args()
    rows = [json.loads(line) for line in args.chematic_rows.read_text().splitlines() if line.strip()]
    rows = [r for r in rows if r.get("arm") == args.arm and r.get("status") == "success"]
    tasks = [(r, args.n_conformers, args.seed) for r in rows]
    with mp.Pool(args.jobs) as pool:
        scored = pool.map(score, tasks, chunksize=4)
    args.rows_output.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in scored))
    report = {
        "schema": "a6-conformer-quality-vs-rdkit/v1",
        "rdkit": rdBase.rdkitVersion,
        "arm": args.arm,
        "n_conformers": args.n_conformers,
        "seed": args.seed,
        "max_iters": MAX_ITERS,
        "summary": summary(scored),
    }
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report["summary"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
