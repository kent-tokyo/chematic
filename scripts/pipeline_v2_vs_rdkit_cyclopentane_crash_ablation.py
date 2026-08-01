#!/usr/bin/env python3
"""Wave 1 follow-up: scope the RDKit `RuntimeError: Invariant Violation --
bad direction in linearSearch` crash on cyclopentane found during the Wave 1
benchmark, rather than reporting it as a general "RDKit crashes on
cyclopentane" claim.

Ablates:
  - useSmallRingTorsions: True / False
  - enforceChirality: True / False
  - 5 seeds (the benchmark's own seed 20260801, plus 4 more)
  - stage: raw embed / +UFF / +MMFF94

Classifies the result as one of:
  - "rdkit_default_reproducible": crashes even with RDKit's own ETKDGv3 defaults
  - "nondefault_small_ring_torsion_only": only crashes with useSmallRingTorsions=True
  - "seed_specific": crash depends on the specific seed, not the config
  - "broadly_reproducible": crashes across the config/seed matrix tested

This is the minimal-repro script referenced by the benchmark report -- run
directly to reproduce: `.venv/bin/python
scripts/pipeline_v2_vs_rdkit_cyclopentane_crash_ablation.py`
"""

import json
import sys

CYCLOPENTANE_SMILES = "C1CCCC1"
BENCHMARK_SEED = 20260801
SEEDS = [BENCHMARK_SEED, 1, 2, 3, 4]


def try_stage(mol_template, use_small_ring_torsions, enforce_chirality, seed, stage):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.Mol(mol_template))
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = use_small_ring_torsions
    params.enforceChirality = enforce_chirality
    params.randomSeed = seed

    try:
        cid = AllChem.EmbedMolecule(mol, params)
    except Exception as e:
        return {"ok": False, "stage": "embed", "error": f"{type(e).__name__}: {e}"}
    if cid < 0:
        return {"ok": False, "stage": "embed", "error": "EmbedMolecule returned -1 (no exception)"}

    if stage == "raw":
        return {"ok": True, "stage": "raw"}

    try:
        if stage == "uff":
            AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=200)
        elif stage == "mmff94":
            mp = AllChem.MMFFGetMoleculeProperties(mol)
            if mp is None:
                return {"ok": False, "stage": "mmff94_setup", "error": "MMFF_parameters_unavailable"}
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
    except Exception as e:
        return {"ok": False, "stage": stage, "error": f"{type(e).__name__}: {e}"}

    return {"ok": True, "stage": stage}


def main():
    try:
        from rdkit import Chem
        import rdkit
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")

    mol_template = Chem.MolFromSmiles(CYCLOPENTANE_SMILES)
    if mol_template is None:
        sys.exit("cyclopentane SMILES failed to parse -- cannot ablate")

    print(f"rdkit_version={rdkit.__version__}", file=sys.stderr)

    results = []
    for use_small_ring_torsions in (True, False):
        for enforce_chirality in (True, False):
            for seed in SEEDS:
                for stage in ("raw", "uff", "mmff94"):
                    r = try_stage(mol_template, use_small_ring_torsions, enforce_chirality, seed, stage)
                    row = {
                        "use_small_ring_torsions": use_small_ring_torsions,
                        "enforce_chirality": enforce_chirality,
                        "seed": seed,
                        "stage": stage,
                        **r,
                    }
                    results.append(row)
                    print(json.dumps(row))

    # --- Classification ---
    crashes = [r for r in results if not r["ok"]]
    crash_configs = {(r["use_small_ring_torsions"], r["enforce_chirality"]) for r in crashes}
    crash_seeds = {r["seed"] for r in crashes}

    default_crashes = [
        r for r in crashes if r["use_small_ring_torsions"] is False and r["enforce_chirality"] is True
    ]
    small_ring_only = all(r["use_small_ring_torsions"] is True for r in crashes) if crashes else False

    if not crashes:
        classification = "not_reproducible_in_this_ablation"
    elif default_crashes:
        classification = "rdkit_default_reproducible"
    elif small_ring_only and len(crash_configs) >= 1:
        classification = "nondefault_small_ring_torsion_only"
    elif len(crash_seeds) < len(SEEDS):
        classification = "seed_specific"
    else:
        classification = "broadly_reproducible"

    summary = {
        "classification": classification,
        "n_total_trials": len(results),
        "n_crashes": len(crashes),
        "crash_configs": sorted(str(c) for c in crash_configs),
        "crash_seeds": sorted(crash_seeds),
        "default_config_crashes": len(default_crashes),
    }
    print("SUMMARY:", json.dumps(summary), file=sys.stderr)
    print(json.dumps({"_summary": summary}))


if __name__ == "__main__":
    main()
