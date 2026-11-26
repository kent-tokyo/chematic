#!/usr/bin/env python3
"""Compare local bond lengths and angles with development crystal structures."""

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import math
import multiprocessing
import platform
import time
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem

from platinum_development_probe import graph_matches


def _worker(pipe, engine, smiles):
    try:
        if engine == "chematic":
            molecule = chematic.from_smiles(smiles)
            coords = molecule.generate_3d()
            status = "generated"
        else:
            molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.numThreads = 1
            if AllChem.EmbedMolecule(molecule, params) != 0:
                pipe.send({"status": "embedding_failed"})
                return
            status = AllChem.MMFFOptimizeMolecule(molecule, maxIters=200)
            coords = [list(molecule.GetConformer().GetAtomPosition(i))
                      for i in range(molecule.GetNumAtoms())]
        pipe.send({"status": status, "coordinates": coords})
    except Exception as error:
        pipe.send({"status": "error", "error": repr(error)})
    finally:
        pipe.close()


def generate(engine, smiles, seconds):
    context = multiprocessing.get_context("spawn")
    parent, child = context.Pipe(False)
    process = context.Process(target=_worker, args=(child, engine, smiles))
    start = time.perf_counter()
    process.start()
    child.close()
    try:
        result = parent.recv() if parent.poll(seconds) else {"status": "timeout"}
    except EOFError:
        result = {"status": "worker_exit"}
    finally:
        if process.is_alive():
            process.terminate()
        process.join()
        parent.close()
    result["elapsed_with_startup_and_cleanup"] = time.perf_counter() - start
    return result


def local_geometry(coords, center, neighbors):
    vectors = [[coords[n][i] - coords[center][i] for i in range(3)] for n in neighbors]
    lengths = [math.sqrt(sum(value * value for value in vector)) for vector in vectors]
    if not all(math.isfinite(length) and length > 1e-8 for length in lengths):
        raise ValueError("invalid local coordinates")
    angles = []
    for i, j in itertools.combinations(range(len(neighbors)), 2):
        cosine = sum(a * b for a, b in zip(vectors[i], vectors[j])) / (lengths[i] * lengths[j])
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return sorted(lengths), sorted(angles)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdf", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--boron-limit", type=int, default=40)
    parser.add_argument("--seconds", type=float, default=5.0)
    args = parser.parse_args()
    if args.boron_limit < 1 or not math.isfinite(args.seconds) or args.seconds <= 0:
        parser.error("limits must be positive")
    split = json.loads(args.split.read_text())
    if hashlib.sha256(args.sdf.read_bytes()).hexdigest() != split["source_sha256"]:
        parser.error("source hash differs from split")
    RDLogger.DisableLog("rdApp.*")
    supplier = Chem.SDMolSupplier(str(args.sdf), removeHs=False)
    records = []
    for row in split["rows"]:
        if row.get("split") != "development":
            continue
        reference = supplier[row["index"]]
        if reference is None:
            continue
        smiles = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
        graph = Chem.MolFromSmiles(smiles)
        if graph is None:
            continue
        boron = any(atom.GetAtomicNum() == 5 and atom.GetDegree() == 3
                    and atom.GetFormalCharge() == 0 for atom in graph.GetAtoms())
        high = any(atom.GetDegree() >= 5 for atom in graph.GetAtoms())
        if boron or high:
            records.append({"index": row["index"], "name": row["name"],
                            "smiles": smiles, "boron": boron, "high": high})
    key = lambda row: hashlib.sha256(row["name"].encode()).digest()
    boron_rows = sorted((r for r in records if r["boron"]), key=key)[:args.boron_limit]
    high_rows = sorted((r for r in records if r["high"]), key=key)
    selected = {row["index"]: row for row in boron_rows + high_rows}
    output = {"source_sha256": split["source_sha256"],
              "split_sha256": hashlib.sha256(args.split.read_bytes()).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "rdkit": rdBase.rdkitVersion,
              "chematic": importlib.metadata.version("chematic"),
              "platform": platform.platform(), "seconds": args.seconds,
              "timing": "per-engine process wall time including startup and cleanup",
              "selection": "development only; first neutral degree-3 boron by SHA256(name), plus all degree>=5 centers",
              "eligible_boron": sum(r["boron"] for r in records),
              "eligible_high": sum(r["high"] for r in records), "results": []}
    for row in sorted(selected.values(), key=key):
        reference = Chem.RemoveHs(supplier[row["index"]])
        graph = Chem.MolFromSmiles(row["smiles"])
        other = chematic.from_smiles(row["smiles"])
        item = dict(row)
        try:
            if not graph_matches(other, graph):
                raise ValueError("atom graph differs")
            mapping = reference.GetSubstructMatch(graph)
            if len(mapping) != graph.GetNumAtoms():
                raise ValueError("reference graph does not map")
            ref_coords = [list(reference.GetConformer().GetAtomPosition(mapping[i]))
                          for i in range(graph.GetNumAtoms())]
            centers = [a.GetIdx() for a in graph.GetAtoms()
                       if (a.GetAtomicNum() == 5 and a.GetDegree() == 3 and a.GetFormalCharge() == 0)
                       or a.GetDegree() >= 5]
            for engine in ("chematic", "rdkit"):
                result = generate(engine, row["smiles"], args.seconds)
                coords = result.pop("coordinates", None)
                if coords is not None:
                    try:
                        result["centers"] = []
                        for center in centers:
                            neighbors = [a.GetIdx() for a in graph.GetAtomWithIdx(center).GetNeighbors()]
                            ref_lengths, ref_angles = local_geometry(ref_coords, center, neighbors)
                            lengths, angles = local_geometry(coords, center, neighbors)
                            result["centers"].append({"atom": center, "element": graph.GetAtomWithIdx(center).GetSymbol(),
                                "reference_lengths": ref_lengths, "lengths": lengths,
                                "reference_angles": ref_angles, "angles": angles})
                    except Exception as error:
                        result["scoring_error"] = repr(error)
                        result.pop("centers", None)
                item[engine] = result
        except Exception as error:
            item["error"] = repr(error)
        output["results"].append(item)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        temporary.replace(args.output)
        print(row["index"], item.get("error", "ok"), flush=True)


if __name__ == "__main__":
    main()
