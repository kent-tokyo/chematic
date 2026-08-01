#!/usr/bin/env python3
"""Wave 1 ("RDKit alternative" program): independent RDKit ETKDGv3 oracle for
the pipeline v2 vs RDKit ETKDGv3 benchmark.

Reads the SAME Tier A/B corpus manifests the Rust chematic dump
(`crates/chematic-3d/examples/pipeline_v2_vs_rdkit_dump.rs`) reads --
independently, never consuming chematic's own output. Neither side feeds the
other: both read the same fixed manifest files, producing two independently
generated JSONL files that the aggregation script (
`scripts/gen_pipeline_v2_vs_rdkit_report.py`) cross-references afterward.

Four arms, matching the benchmark spec:
  - etkdgv3_raw: ETKDGv3 embedding only, no force-field optimization.
  - etkdgv3_uff: ETKDGv3 + UFFOptimizeMolecule.
  - etkdgv3_mmff94: ETKDGv3 + MMFFOptimizeMolecule (RDKit's own MMFF94).
  - etkdgv3_best_of_n: N=10 ETKDGv3 conformers, each UFF-optimized, lowest
    final UFF energy kept -- a genuinely distinct arm (multi-restart), not
    folded into etkdgv3_uff's single-attempt result.

AddHs: explicit hydrogens are added before embedding (RDKit's own documented
best practice for ETKDG; chematic-3d's pipeline_v2 does not add explicit
hydrogens itself, so this is disclosed as a REAL, not incidental, difference
in operating conditions between the two engines -- not something this script
tries to normalize away). Heavy-atom mapping: `Chem.AddHs` appends new H
atoms at the end and never reorders existing atoms (documented RDKit
behavior), so heavy-atom index i in the post-AddHs molecule equals atom index
i in the pre-AddHs parse, which is the same index chematic assigns when
parsing the identical SMILES string. `heavy_atom_elements` is emitted here so
the aggregation script can independently verify this correspondence
per-molecule (element-symbol sequence match) rather than assume it.

Run: `.venv/bin/python scripts/pipeline_v2_vs_rdkit_oracle.py
  > validation/results/pipeline_v2_vs_rdkit_rdkit_rows.jsonl`
"""

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).parent.parent

RANDOM_SEED = 20260801  # same numeric seed as the Rust dump's EMBED_SEED; not a cross-engine determinism claim
MAX_ATTEMPTS = 8  # matches chematic's EmbedParameters.max_attempts for etkdgv3_raw/uff/mmff94
BEST_OF_N = 10
FORCE_FIELD_MAX_ITERATIONS = 200


def load_manifest(path):
    with open(ROOT / path) as f:
        return json.load(f)


def heavy_atom_elements(mol):
    return [atom.GetSymbol() for atom in mol.GetAtoms()]


def coords_to_json(conf, n_heavy):
    return [list(conf.GetAtomPosition(i)) for i in range(n_heavy)]


def embed_with_retry(mol, params, max_attempts):
    """Mirrors chematic's own max_attempts retry convention: each attempt
    gets its own derived seed (base_seed + attempt_index), not a single fixed
    seed retried identically. Returns the conformer id, or -1 on total failure."""
    from rdkit.Chem import AllChem

    for attempt in range(max_attempts):
        params.randomSeed = RANDOM_SEED + attempt
        cid = AllChem.EmbedMolecule(mol, params)
        if cid >= 0:
            return cid, attempt + 1
    return -1, max_attempts


def run_etkdgv3_raw(mol_noHs, n_heavy):
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(mol_noHs)
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.enforceChirality = True

    t0 = time.monotonic()
    cid, attempts_used = embed_with_retry(mol, params, MAX_ATTEMPTS)
    elapsed_ms = int((time.monotonic() - t0) * 1000)

    if cid < 0:
        return {
            "status": "typed_failure",
            "failure_cause": "EmbedMolecule_failed",
            "elapsed_ms": elapsed_ms,
        }

    conf = mol.GetConformer(cid)
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "embed_attempts_used": attempts_used,
        "coords": coords_to_json(conf, n_heavy),
        "force_field": "none",
    }


def run_etkdgv3_plus_ff(mol_noHs, n_heavy, ff_name):
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(mol_noHs)
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.enforceChirality = True

    t0 = time.monotonic()
    cid, attempts_used = embed_with_retry(mol, params, MAX_ATTEMPTS)
    if cid < 0:
        elapsed_ms = int((time.monotonic() - t0) * 1000)
        return {
            "status": "typed_failure",
            "failure_cause": "EmbedMolecule_failed",
            "elapsed_ms": elapsed_ms,
        }

    ff_ret = None
    ff_error = None
    try:
        if ff_name == "uff":
            ff_ret = AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=FORCE_FIELD_MAX_ITERATIONS)
        elif ff_name == "mmff94":
            mp = AllChem.MMFFGetMoleculeProperties(mol)
            if mp is None:
                ff_error = "MMFF_parameters_unavailable"
            else:
                ff_ret = AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=FORCE_FIELD_MAX_ITERATIONS)
    except Exception as e:  # RDKit force-field calls can raise on pathological inputs
        ff_error = f"{type(e).__name__}: {e}"

    elapsed_ms = int((time.monotonic() - t0) * 1000)

    if ff_error is not None:
        return {
            "status": "typed_failure",
            "failure_cause": ff_error,
            "elapsed_ms": elapsed_ms,
            "embed_attempts_used": attempts_used,
        }

    conf = mol.GetConformer(cid)
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "embed_attempts_used": attempts_used,
        "coords": coords_to_json(conf, n_heavy),
        "force_field": ff_name,
        # RDKit's *OptimizeMolecule returns 0 (converged), 1 (not converged),
        # or -1 (couldn't be set up, e.g. UFF-unsupported element)
        "force_field_converged": ff_ret == 0 if ff_ret is not None else None,
        "force_field_return_code": ff_ret,
    }


def run_best_of_n(mol_noHs, n_heavy):
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(mol_noHs)
    params = AllChem.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    params.useMacrocycle14config = True
    params.enforceChirality = True
    params.randomSeed = RANDOM_SEED

    t0 = time.monotonic()
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=BEST_OF_N, params=params))
    if not cids:
        elapsed_ms = int((time.monotonic() - t0) * 1000)
        return {
            "status": "typed_failure",
            "failure_cause": "EmbedMultipleConfs_failed",
            "elapsed_ms": elapsed_ms,
        }

    best_cid = None
    best_energy = None
    n_optimized = 0
    for cid in cids:
        try:
            ff_ret = AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=FORCE_FIELD_MAX_ITERATIONS)
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            if ff is None:
                continue
            energy = ff.CalcEnergy()
            n_optimized += 1
            if best_energy is None or energy < best_energy:
                best_energy = energy
                best_cid = cid
        except Exception:
            continue

    elapsed_ms = int((time.monotonic() - t0) * 1000)

    if best_cid is None:
        return {
            "status": "typed_failure",
            "failure_cause": "all_conformers_failed_uff_optimization",
            "elapsed_ms": elapsed_ms,
            "conformers_embedded": len(cids),
        }

    conf = mol.GetConformer(best_cid)
    return {
        "status": "success",
        "elapsed_ms": elapsed_ms,
        "conformers_embedded": len(cids),
        "conformers_optimized": n_optimized,
        "best_conformer_id": int(best_cid),
        "best_uff_energy": best_energy,
        "coords": coords_to_json(conf, n_heavy),
        "force_field": "uff",
    }


ARMS = {
    "rdkit_etkdgv3_raw": run_etkdgv3_raw,
    "rdkit_etkdgv3_uff": lambda m, n: run_etkdgv3_plus_ff(m, n, "uff"),
    "rdkit_etkdgv3_mmff94": lambda m, n: run_etkdgv3_plus_ff(m, n, "mmff94"),
    "rdkit_etkdgv3_best_of_n": run_best_of_n,
}


def main():
    global Chem
    try:
        from rdkit import Chem as _Chem
        from rdkit.Chem import AllChem as _AllChem
        import rdkit
    except ImportError:
        sys.exit("rdkit not installed. pip install rdkit")
    Chem = _Chem

    eprint = lambda *a: print(*a, file=sys.stderr)
    eprint(
        f"config_snapshot rdkit_version={rdkit.__version__} random_seed={RANDOM_SEED} "
        f"max_attempts={MAX_ATTEMPTS} best_of_n={BEST_OF_N} "
        f"add_hs=True enforce_chirality=True use_small_ring_torsions=True "
        f"use_macrocycle_torsions=True use_macrocycle14config=True "
        f"force_field_max_iterations={FORCE_FIELD_MAX_ITERATIONS}"
    )

    for tier, path in [
        ("A", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
        ("B", "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
    ]:
        manifest = load_manifest(path)
        for m in manifest["molecules"]:
            name = m["name"]
            smiles = m["smiles"]
            primary_category = m.get("primary_category", "unknown")

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                row = {
                    "tier": tier,
                    "name": name,
                    "smiles": smiles,
                    "primary_category": primary_category,
                    "arm": "all",
                    "status": "parse_failure",
                }
                print(json.dumps(row))
                continue

            n_heavy = mol.GetNumAtoms()
            elements = heavy_atom_elements(mol)

            for arm_name, fn in ARMS.items():
                try:
                    result = fn(mol, n_heavy)
                except Exception as e:  # never let one molecule/arm crash the whole run
                    result = {
                        "status": "internal_error",
                        "failure_cause": f"{type(e).__name__}: {e}",
                        "elapsed_ms": 0,
                    }
                row = {
                    "tier": tier,
                    "name": name,
                    "smiles": smiles,
                    "primary_category": primary_category,
                    "arm": arm_name,
                    "heavy_atom_elements": elements,
                    "heavy_atom_count": n_heavy,
                    **result,
                }
                print(json.dumps(row))


if __name__ == "__main__":
    main()
