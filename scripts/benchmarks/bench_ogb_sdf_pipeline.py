#!/usr/bin/env python3
"""Measure SDF reading, categorical chemistry, and complete OGB arrays.

Use a source-built chematic wheel, the pinned sdfrust wheel, and RDKit in the
same Python environment. Timings are local observations, not product rankings.
The second stage includes chemistry perception and categorical graph creation;
the third also constructs the OGB NumPy arrays. All stages reopen the same SDF.

Example:
  PYTHONPATH=/path/to/chematic:/path/to/sdfrust python3 \
    scripts/benchmarks/bench_ogb_sdf_pipeline.py --output /tmp/ogb-pipeline.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import statistics
import tempfile
import time
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase

import chematic
import sdfrust


ATOM_COLUMNS = (
    list(range(1, 119)) + ["misc"],
    ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW", "CHI_OTHER", "misc"],
    list(range(11)) + ["misc"],
    list(range(-5, 6)) + ["misc"],
    list(range(9)) + ["misc"],
    list(range(5)) + ["misc"],
    ["SP", "SP2", "SP3", "SP3D", "SP3D2", "misc"],
)
BOND_TYPES = ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "misc"]
BOND_STEREO = ["STEREONONE", "STEREOZ", "STEREOE", "STEREOCIS", "STEREOTRANS", "STEREOANY"]


def index(allowed: list, value: object) -> int:
    try:
        return allowed.index(value)
    except ValueError:
        return len(allowed) - 1


def rdkit_raw_graph(mol: Chem.Mol) -> dict:
    """Apply OGB's published feature categories to an RDKit molecule.

    Category definitions and edge order: https://github.com/snap-stanford/ogb/
    blob/master/ogb/utils/features.py and ogb/utils/mol.py.
    """
    nodes = []
    for atom in mol.GetAtoms():
        nodes.append([
            index(ATOM_COLUMNS[0], atom.GetAtomicNum()),
            index(ATOM_COLUMNS[1], str(atom.GetChiralTag())),
            index(ATOM_COLUMNS[2], atom.GetTotalDegree()),
            index(ATOM_COLUMNS[3], atom.GetFormalCharge()),
            index(ATOM_COLUMNS[4], atom.GetTotalNumHs()),
            index(ATOM_COLUMNS[5], atom.GetNumRadicalElectrons()),
            index(ATOM_COLUMNS[6], str(atom.GetHybridization())),
            int(atom.GetIsAromatic()),
            int(atom.IsInRing()),
        ])
    sources, destinations, features = [], [], []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        feature = [
            index(BOND_TYPES, str(bond.GetBondType())),
            BOND_STEREO.index(str(bond.GetStereo())),
            int(bond.GetIsConjugated()),
        ]
        sources.extend((a, b))
        destinations.extend((b, a))
        features.extend((feature, feature))
    return {
        "node_feat": nodes,
        "edge_index": [sources, destinations],
        "edge_feat": features,
        "num_nodes": len(nodes),
    }


def arrays(raw: dict) -> dict:
    return {
        "node_feat": np.asarray(raw["node_feat"], dtype=np.int64).reshape(-1, 9),
        "edge_index": np.asarray(raw["edge_index"], dtype=np.int64).reshape(2, -1),
        "edge_feat": np.asarray(raw["edge_feat"], dtype=np.int64).reshape(-1, 3),
        "num_nodes": raw["num_nodes"],
    }


def sdfrust_arrays(raw: dict) -> dict:
    return arrays({
        "node_feat": raw["atom_features"],
        "edge_index": [raw["edge_src"], raw["edge_dst"]],
        "edge_feat": raw["bond_features"],
        "num_nodes": raw["num_atoms"],
    })


def read_chematic(path: Path) -> list:
    molecules = []
    for index_, record, error in chematic.iter_sdf_checked(str(path)):
        if error:
            raise ValueError(f"chematic record {index_}: {error}")
        molecules.append(record)
    return molecules


def read_sdfrust(path: Path) -> list:
    return list(sdfrust.iter_sdf_file(str(path)))


def read_sdfrust_v3000(path: Path) -> list:
    return list(sdfrust.iter_sdf_v3000_file(str(path)))


def read_rdkit(path: Path) -> list:
    molecules = list(Chem.ForwardSDMolSupplier(str(path), sanitize=False, removeHs=False))
    if any(mol is None for mol in molecules):
        raise ValueError("RDKit rejected an SDF record")
    return molecules


def read_rdkit_sanitized(path: Path) -> list:
    molecules = list(Chem.ForwardSDMolSupplier(str(path), sanitize=True, removeHs=False))
    if any(mol is None for mol in molecules):
        raise ValueError("RDKit rejected an SDF record during sanitization")
    return molecules


READERS = {"chematic": read_chematic, "sdfrust": read_sdfrust, "rdkit": read_rdkit}


def run_stage(path: Path, tool: str, stage: str) -> int:
    if tool == "chematic" and stage == "read":
        molecules = list(chematic.SDMolSupplier(str(path)))
    elif tool == "rdkit" and stage != "read":
        molecules = read_rdkit_sanitized(path)
    else:
        molecules = READERS[tool](path)
    if stage == "read":
        return len(molecules)
    count = 0
    for molecule in molecules:
        if tool == "chematic":
            graph = (molecule._ogb_graph_features_raw() if stage == "categorical"
                     else molecule.ogb_graph_features())
        elif tool == "sdfrust":
            raw = molecule.ogb_graph_features()
            graph = raw if stage == "categorical" else sdfrust_arrays(raw)
        else:
            raw = rdkit_raw_graph(molecule)
            graph = raw if stage == "categorical" else arrays(raw)
        count += graph["num_atoms"] if tool == "sdfrust" and stage == "categorical" else graph["num_nodes"]
    return count


def make_fixture(source: Path, target: Path, limit: int, retain_hs: bool,
                 v3000: bool) -> int:
    writer = Chem.SDWriter(str(target))
    writer.SetForceV3000(v3000)
    accepted = 0
    for molecule in Chem.ForwardSDMolSupplier(str(source), sanitize=True, removeHs=False):
        if molecule is None:
            raise ValueError("source SDF contains a rejected record")
        writer.write(molecule if retain_hs else Chem.RemoveHs(molecule))
        accepted += 1
        if accepted == limit:
            break
    writer.close()
    if accepted == 0:
        raise ValueError("source SDF has no molecules")
    return accepted


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def flattened_rdkit_graph(molecule: Chem.Mol) -> dict:
    """Read the same RDKit-written CTAB with 2D coordinates and wedge flags."""
    lines = Chem.MolToMolBlock(molecule).splitlines()
    if len(lines) < 4 or "V3000" in lines[3]:
        return {}
    atom_count = int(lines[3][:3])
    lines[1] = lines[1].replace("3D", "2D")
    for line in range(4, 4 + atom_count):
        lines[line] = lines[line][:20] + "    0.0000" + lines[line][30:]
    flat = Chem.MolFromMolBlock("\n".join(lines) + "\n", sanitize=True, removeHs=False)
    return arrays(rdkit_raw_graph(flat)) if flat is not None else {}


def correctness(path: Path) -> dict:
    mols = {"chematic": read_chematic(path), "sdfrust": read_sdfrust(path),
            "rdkit": read_rdkit_sanitized(path)}
    counts = {name: len(items) for name, items in mols.items()}
    if len(set(counts.values())) != 1:
        return {"record_counts": counts, "comparable": False}
    stats = {
        name: {"equal_graphs": 0, "shape_mismatch": 0, "node_cells_different": 0,
               "edge_cells_different": 0, "edge_index_cells_different": 0,
               "mismatched_records": []}
        for name in ("chematic", "sdfrust")
    }
    for index_, (cmol, smol, rmol) in enumerate(zip(mols["chematic"], mols["sdfrust"], mols["rdkit"])):
        reference = arrays(rdkit_raw_graph(rmol))
        graphs = {"chematic": cmol.ogb_graph_features(),
                  "sdfrust": sdfrust_arrays(smol.ogb_graph_features())}
        flat_reference = None
        for name, graph in graphs.items():
            if any(graph[key].shape != reference[key].shape for key in
                   ("node_feat", "edge_index", "edge_feat")):
                stats[name]["shape_mismatch"] += 1
                stats[name]["mismatched_records"].append({"index": index_,
                    "name": rmol.GetProp("_Name") if rmol.HasProp("_Name") else "",
                    "reason": "array_shape"})
                continue
            differences = {
                "node_cells_different": int(np.count_nonzero(graph["node_feat"] != reference["node_feat"])),
                "edge_cells_different": int(np.count_nonzero(graph["edge_feat"] != reference["edge_feat"])),
                "edge_index_cells_different": int(np.count_nonzero(graph["edge_index"] != reference["edge_index"])),
            }
            for key, value in differences.items():
                stats[name][key] += value
            if all(value == 0 for value in differences.values()):
                stats[name]["equal_graphs"] += 1
            else:
                if flat_reference is None:
                    flat_reference = flattened_rdkit_graph(rmol)
                changed = np.argwhere(graph["node_feat"] != reference["node_feat"])
                stats[name]["mismatched_records"].append({
                    "index": index_,
                    "name": rmol.GetProp("_Name") if rmol.HasProp("_Name") else "",
                    "node_cells": [[int(row), int(column)] for row, column in changed],
                    "matches_rdkit_2d_wedges": bool(flat_reference) and all(
                        np.array_equal(graph[key], flat_reference[key])
                        for key in ("node_feat", "edge_index", "edge_feat")),
                })
    return {"record_counts": counts, "comparable": True, "against_rdkit_ogb": stats}


def measure(path: Path, repeats: int) -> dict:
    stages = ("read", "categorical", "numpy_graph")
    tools = tuple(READERS)
    for stage in stages:
        for name in tools:
            run_stage(path, name, stage)
    samples = {stage: {name: [] for name in tools} for stage in stages}
    results = {stage: {name: None for name in tools} for stage in stages}
    for stage in stages:
        for repeat in range(repeats):
            for name in tools[repeat % len(tools):] + tools[:repeat % len(tools)]:
                started = time.perf_counter_ns()
                count = run_stage(path, name, stage)
                elapsed = time.perf_counter_ns() - started
                if results[stage][name] is None:
                    results[stage][name] = count
                elif results[stage][name] != count:
                    raise RuntimeError(f"{name} {stage} produced inconsistent counts")
                samples[stage][name].append(elapsed)
    return {
        stage: {
            name: {"count": results[stage][name],
                   "median_ms": round(statistics.median(values) / 1e6, 3),
                   "p95_ms": round(sorted(values)[math.ceil(len(values)*0.95)-1] / 1e6, 3)}
            for name, values in measurements.items()
        }
        for stage, measurements in samples.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdf", type=Path, default=None)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--retain-hs", action="store_true")
    parser.add_argument("--v3000", action="store_true")
    parser.add_argument("--candidate-wheel", type=Path)
    parser.add_argument("--sdfrust-wheel", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.limit < 1 or args.repeats < 3:
        parser.error("--limit must be positive and --repeats must be at least 3")
    source = args.sdf or Path(Chem.__file__).resolve().parents[1] / "Contrib/PBF/testData/egfr.sdf"
    source = source.resolve()
    if not source.is_file():
        parser.error(f"SDF not found: {source}")
    if args.v3000:
        READERS["sdfrust"] = read_sdfrust_v3000
    with tempfile.TemporaryDirectory(prefix="chematic-ogb-benchmark-") as directory:
        fixture = Path(directory) / "fixture.sdf"
        records = make_fixture(source, fixture, args.limit, args.retain_hs, args.v3000)
        report = {
            "schema_version": 1,
            "status": "local-source-candidate",
            "source": {"path": str(source), "sha256": sha256(source)},
            "fixture": {"sha256": sha256(fixture), "bytes": fixture.stat().st_size,
                        "records": records, "format": "V3000" if args.v3000 else "V2000",
                        "hydrogen_policy": "retained" if args.retain_hs else "RDKit RemoveHs"},
            "versions": {"chematic": chematic.__version__, "sdfrust": sdfrust.__version__,
                         "rdkit": rdBase.rdkitVersion, "python": platform.python_version()},
            "module_paths": {"chematic": chematic.__file__, "sdfrust": sdfrust.__file__},
            "wheels": {
                "chematic": ({"path": str(args.candidate_wheel.resolve()),
                              "sha256": sha256(args.candidate_wheel)} if args.candidate_wheel else None),
                "sdfrust": ({"path": str(args.sdfrust_wheel.resolve()),
                             "sha256": sha256(args.sdfrust_wheel)} if args.sdfrust_wheel else None),
            },
            "host": {"machine": platform.machine(), "platform": platform.platform()},
            "repeats": args.repeats,
            "contracts": {
                "read": "stream file into molecule objects; chematic SDMolSupplier fast path; sdfrust iter_sdf_file; RDKit sanitize=False",
                "categorical": "read plus chemistry perception and OGB category graph as Python lists; chematic uses checked SdfRecord diagnostics",
                "numpy_graph": "read plus chemistry perception and complete OGB int64 NumPy graph; chematic uses checked SdfRecord diagnostics",
                "scope": "single local fixture and Python API; no equal-chemistry speed claim unless graphs agree",
            },
            "correctness": correctness(fixture),
            "timing": measure(fixture, args.repeats),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
