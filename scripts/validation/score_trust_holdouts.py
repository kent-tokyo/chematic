#!/usr/bin/env python3
"""Score a frozen trust holdout once, without changing its inputs.

The output keeps success, mismatch, refusal, invalid-input and internal-error
counts separate for every observable. RDKit is a pinned compatibility oracle,
not an independent chemistry truth.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import tempfile
from collections import Counter
from pathlib import Path

import chematic
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdFMCS

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "validation/trust_holdouts/v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(name: str) -> list[dict]:
    return [json.loads(line) for line in (BASE / name).read_text().splitlines() if line]


def status_counter() -> Counter:
    return Counter({key: 0 for key in ("success", "mismatch", "refusal", "invalid_input", "internal_error")})


def score_mcs(source: list[dict]) -> tuple[dict, list[dict]]:
    columns = {key: status_counter() for key in ("parse", "completion", "atom_count", "bond_count")}
    details = []
    for row in source:
        detail = {"id": row["id"]}
        try:
            inputs = row["inputs"] if "inputs" in row else [row["left"], row["right"]]
            rd_inputs = [Chem.MolFromSmiles(item["smiles"]) for item in inputs]
            cm_inputs = [chematic.from_smiles(item["smiles"]) for item in inputs]
            if any(mol is None for mol in rd_inputs):
                raise ValueError("RDKit rejected a frozen input")
        except (ValueError, TypeError) as error:
            columns["parse"]["invalid_input"] += 1
            detail.update(status="invalid_input", error=str(error))
            details.append(detail)
            continue
        except Exception as error:  # pragma: no cover - recorded fail-closed
            columns["parse"]["internal_error"] += 1
            detail.update(status="internal_error", error=f"{type(error).__name__}: {error}")
            details.append(detail)
            continue
        columns["parse"]["success"] += 1
        atom_rd = rdFMCS.AtomCompare.CompareElements if row["atom_compare"] == "elements" else rdFMCS.AtomCompare.CompareAny
        bond_rd = rdFMCS.BondCompare.CompareOrder if row["bond_compare"] == "order" else rdFMCS.BondCompare.CompareAny
        atom_cm = "elements" if row["atom_compare"] == "elements" else "any"
        bond_cm = "order_or_aromatic" if row["bond_compare"] == "order" else "any"
        try:
            rd = rdFMCS.FindMCS(rd_inputs, atomCompare=atom_rd, bondCompare=bond_rd, timeout=10)
            cm, cm_timeout = chematic.find_mcs_checked(
                cm_inputs, atom_compare=atom_cm, bond_compare=bond_cm, timeout_ms=10_000
            )
            if rd.canceled or cm_timeout:
                columns["completion"]["refusal"] += 1
                detail.update(status="refusal", rdkit_timeout=bool(rd.canceled), chematic_timeout=cm_timeout)
            elif cm is None:
                columns["completion"]["mismatch"] += 1
                detail.update(status="mismatch", reason="chematic_returned_none")
            else:
                columns["completion"]["success"] += 1
                cm_atoms = cm.heavy_atoms
                cm_bonds = len(cm.bond_table)
                for name, actual, expected in (
                    ("atom_count", cm_atoms, rd.numAtoms),
                    ("bond_count", cm_bonds, rd.numBonds),
                ):
                    columns[name]["success" if actual == expected else "mismatch"] += 1
                detail.update(
                    status="success" if (cm_atoms, cm_bonds) == (rd.numAtoms, rd.numBonds) else "mismatch",
                    chematic_atoms=cm_atoms, rdkit_atoms=rd.numAtoms,
                    chematic_bonds=cm_bonds, rdkit_bonds=rd.numBonds,
                )
        except (ValueError, TypeError) as error:
            columns["completion"]["refusal"] += 1
            detail.update(status="refusal", error=str(error))
        except Exception as error:  # pragma: no cover - recorded fail-closed
            columns["completion"]["internal_error"] += 1
            detail.update(status="internal_error", error=f"{type(error).__name__}: {error}")
        details.append(detail)
    return {key: dict(value) for key, value in columns.items()}, details


def score_sdf(source: list[dict]) -> tuple[dict, list[dict]]:
    columns = {key: status_counter() for key in ("parse", "graph", "title", "properties", "stereo")}
    details = []
    with tempfile.TemporaryDirectory(prefix="chematic-trust-sdf-") as tmp:
        for n, row in enumerate(source):
            detail = {"id": row["id"]}
            try:
                mol = chematic.from_smiles(row["smiles"])
                expected_graph = mol.rdkit_smiles
                mol.SetName(row["title"])
                for key, value in row["fields"].items():
                    mol.SetProp(key, value)
            except (ValueError, TypeError) as error:
                columns["parse"]["invalid_input"] += 1
                detail.update(status="invalid_input", error=str(error))
                details.append(detail)
                continue
            columns["parse"]["success"] += 1
            path = Path(tmp) / f"{n}.sdf"
            try:
                with chematic.SDWriter(str(path)) as writer:
                    writer.write(mol)
                reread = next(iter(chematic.SDMolSupplier(str(path), removeHs=False)))
                if reread is None:
                    raise ValueError("reader returned None")
                observed = {
                    "graph": reread.rdkit_smiles == expected_graph,
                    "title": reread.GetName() == row["title"],
                    "properties": {
                        key: value for key, value in reread.GetPropsAsDict().items() if key != "_Name"
                    } == row["fields"],
                    "stereo": Chem.MolToSmiles(Chem.MolFromMolBlock(path.read_text().split("$$$$", 1)[0], removeHs=False), isomericSmiles=True) == expected_graph,
                }
                for key, matched in observed.items():
                    columns[key]["success" if matched else "mismatch"] += 1
                detail.update(status="success" if all(observed.values()) else "mismatch", **observed)
            except (ValueError, TypeError, StopIteration) as error:
                columns["graph"]["refusal"] += 1
                detail.update(status="refusal", error=str(error))
            except Exception as error:  # pragma: no cover - recorded fail-closed
                columns["graph"]["internal_error"] += 1
                detail.update(status="internal_error", error=f"{type(error).__name__}: {error}")
            details.append(detail)
    return {key: dict(value) for key, value in columns.items()}, details


def max_distance_delta(left, right) -> float:
    """Rigid-transform-invariant maximum pair-distance difference."""
    if len(left) != len(right):
        return math.inf
    maximum = 0.0
    for i in range(len(left)):
        for j in range(i):
            left_d = math.dist(left[i], left[j])
            right_d = math.dist(right[i], right[j])
            maximum = max(maximum, abs(left_d - right_d))
    return maximum


def score_stereo_3d(source: list[dict]) -> tuple[dict, list[dict]]:
    columns = {key: status_counter() for key in (
        "parse", "stereo_graph", "rdkit_embed", "chematic_embed", "joint_embed_status",
        "conformer_count", "coordinates",
    )}
    details = []
    for row in source:
        detail = {"id": row["id"]}
        try:
            rd_base = Chem.MolFromSmiles(row["smiles"])
            cm_base = chematic.from_smiles(row["smiles"])
            if rd_base is None:
                raise ValueError("RDKit rejected a frozen input")
            rd_smiles = Chem.MolToSmiles(rd_base, canonical=True, isomericSmiles=True)
            same_graph = cm_base.rdkit_smiles == rd_smiles
            columns["stereo_graph"]["success" if same_graph else "mismatch"] += 1
            rd = Chem.AddHs(rd_base)
            cm = cm_base.add_hydrogens()
            columns["parse"]["success"] += 1
        except (ValueError, TypeError) as error:
            columns["parse"]["invalid_input"] += 1
            detail.update(status="invalid_input", error=str(error))
            details.append(detail)
            continue
        rd_ids = None
        cm_confs = None
        rd_error = None
        cm_error = None
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = int(row["random_seed"])
            rd_ids = list(AllChem.EmbedMultipleConfs(rd, numConfs=int(row["conformers"]), params=params))
            columns["rdkit_embed"]["success"] += 1
        except (ValueError, TypeError, RuntimeError) as error:
            rd_error = f"{type(error).__name__}: {error}"
            columns["rdkit_embed"]["refusal"] += 1
        except Exception as error:  # pragma: no cover - recorded fail-closed
            rd_error = f"{type(error).__name__}: {error}"
            columns["rdkit_embed"]["internal_error"] += 1
        try:
            cm_confs = cm.rdkit_embed_multiple_confs(int(row["conformers"]), random_seed=int(row["random_seed"]))
            columns["chematic_embed"]["success"] += 1
        except (ValueError, TypeError, RuntimeError) as error:
            cm_error = f"{type(error).__name__}: {error}"
            columns["chematic_embed"]["refusal"] += 1
        except Exception as error:  # pragma: no cover - recorded fail-closed
            cm_error = f"{type(error).__name__}: {error}"
            columns["chematic_embed"]["internal_error"] += 1

        if rd_ids is None or cm_confs is None:
            if rd_ids is None and cm_confs is None:
                columns["joint_embed_status"]["success"] += 1
                detail.update(status="refusal", rdkit_error=rd_error, chematic_error=cm_error)
            else:
                columns["joint_embed_status"]["mismatch"] += 1
                detail.update(status="mismatch", rdkit_error=rd_error, chematic_error=cm_error)
            details.append(detail)
            continue

        columns["joint_embed_status"]["success"] += 1
        try:
            same_count = len(rd_ids) == len(cm_confs)
            columns["conformer_count"]["success" if same_count else "mismatch"] += 1
            deltas = []
            if same_count:
                for conf_id, cm_coords in zip(rd_ids, cm_confs, strict=True):
                    conf = rd.GetConformer(conf_id)
                    rd_coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(rd.GetNumAtoms())]
                    deltas.append(max_distance_delta(cm_coords, rd_coords))
            matched_coords = same_count and all(math.isfinite(delta) and delta <= 1e-4 for delta in deltas)
            columns["coordinates"]["success" if matched_coords else "mismatch"] += 1
            detail.update(
                status="success" if same_graph and same_count and matched_coords else "mismatch",
                stereo_graph=same_graph,
                chematic_conformers=len(cm_confs), rdkit_conformers=len(rd_ids),
                max_pair_distance_delta=max(deltas, default=0.0),
            )
        except (ValueError, TypeError, RuntimeError) as error:
            columns["coordinates"]["refusal"] += 1
            detail.update(status="refusal", error=str(error))
        except Exception as error:  # pragma: no cover - recorded fail-closed
            columns["coordinates"]["internal_error"] += 1
            detail.update(status="internal_error", error=f"{type(error).__name__}: {error}")
        details.append(detail)
    return {key: dict(value) for key, value in columns.items()}, details


def main() -> None:
    global BASE
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--holdout-version", default="v1")
    parser.add_argument("--wheel", type=Path)
    args = parser.parse_args()
    BASE = ROOT / "validation/trust_holdouts" / args.holdout_version
    manifest = json.loads((BASE / "manifest.json").read_text())
    if args.candidate_commit != manifest["candidate_frozen_commit"]:
        raise SystemExit("candidate commit does not match the frozen manifest")
    for name, contract in manifest["lanes"].items():
        if sha(BASE / name) != contract["sha256"]:
            raise SystemExit(f"frozen input hash changed: {name}")
    scorers = {
        "mcs": ("mcs-comparison.jsonl", score_mcs),
        "sdf": ("sdf-attributes.jsonl", score_sdf),
        "stereo_3d": ("stereo-3d.jsonl", score_stereo_3d),
    }
    lanes = {}
    for lane, (name, scorer) in scorers.items():
        columns, details = scorer(rows(name))
        lanes[lane] = {"input": name, "columns": columns, "rows": details}
    module_spec = importlib.util.find_spec("chematic")
    module_path = Path(module_spec.origin) if module_spec and module_spec.origin else None
    native_extensions = {
        path.name: sha(path)
        for path in sorted(module_path.parent.glob("*.so"))
    } if module_path else {}
    payload = {
        "schema": "chematic.trust-holdout-score.v2",
        "holdout_version": args.holdout_version,
        "candidate_commit": args.candidate_commit,
        "frozen_manifest_sha256": sha(BASE / "manifest.json"),
        "wheel": {
            "path": str(args.wheel) if args.wheel else None,
            "sha256": sha(args.wheel) if args.wheel else None,
        },
        "versions": {"chematic": chematic.__version__, "rdkit": rdBase.rdkitVersion, "python": platform.python_version()},
        "module": {
            "path": str(module_path),
            "sha256": sha(module_path) if module_path and module_path.is_file() else None,
            "native_extensions": native_extensions,
        },
        "policy": "One-shot scoring only; no tuning on these rows after inspection.",
        "lanes": lanes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({lane: data["columns"] for lane, data in lanes.items()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
