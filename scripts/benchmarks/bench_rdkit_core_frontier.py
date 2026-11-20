#!/usr/bin/env python3
"""Compare core, CIP, Morgan, and search against a pinned RDKit wheel.

One run uses one named corpus. All rows keep their original line indices.
Independent gold labels and larger search corpora are separate evaluations.
"""

import argparse
from collections import Counter
import gzip
import hashlib
import heapq
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import sys
from time import perf_counter
from zipfile import ZipFile

import chematic
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import rdCIPLabeler, rdFingerprintGenerator, rdFMCS, rdSubstructLibrary
from rdkit.Chem.MolStandardize import rdMolStandardize

from bench_python_op_matrix_vs_rdkit import cpu_model, paired_median_ci

QUERIES = (
    "c1ccccc1", "[#7]", "[r5]", "[R2]", "C(=O)[OH]",
    "[NX3;H2,H1;!$(NC=O)]", "[!#6;!#1]", "[CX4]", "C[C@H](O)Cl",
)
MCS_INPUTS = (
    ("c1ccccc1", "Cc1ccccc1"),
    ("C1CC1", "CCC"),
    ("CC(O)C", "CC(=O)C", "CC(N)C"),
    ("C1CC2CCC1C2", "C1CCC2CCCCC2C1"),
    ("c1ccc2ccccc2c1", "c1ccc2c(c1)cco2"),
)


def digest(path):
    with path.open("rb") as reader:
        return hashlib.file_digest(reader, "sha256").hexdigest()


def package_identity(wheel, name, version):
    distribution = metadata.distribution(name)
    package_root = distribution.locate_file("")
    checked = 0
    with ZipFile(wheel) as archive:
        for member in archive.namelist():
            if not member.startswith(f"{name}/") or member.endswith("/"):
                continue
            installed = package_root / member
            if not installed.is_file() or installed.read_bytes() != archive.read(member):
                raise ValueError(f"installed package differs from {wheel.name}: {member}")
            checked += 1
    if not checked:
        raise ValueError(f"{wheel.name} has no {name} package files")
    return {"version": version, "wheel": wheel.name,
            "sha256": digest(wheel), "package_files_verified": checked}


def error_record(error):
    return {"status": "error", "kind": type(error).__name__,
            "code": getattr(error, "code", None), "message": str(error)}


def capture(operation, value):
    try:
        result = operation(value)
    except Exception as error:
        return error_record(error)
    if result is None:
        return {"status": "rejected"}
    return {"status": "ok", "value": result}


def clear_cip_labels(molecule):
    for item in [*molecule.GetAtoms(), *molecule.GetBonds()]:
        if item.HasProp("_CIPCode"):
            item.ClearProp("_CIPCode")



def cip_labels(molecule, *, clear=True):
    if clear:
        clear_cip_labels(molecule)
    rdCIPLabeler.AssignCIPLabels(molecule)
    return sorted(
        [["atom", a.GetIdx(), a.GetProp("_CIPCode")]
         for a in molecule.GetAtoms() if a.HasProp("_CIPCode")]
        + [["bond", *sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), b.GetProp("_CIPCode")]
           for b in molecule.GetBonds() if b.HasProp("_CIPCode")]
    )


def chematic_cip(molecule):
    values = molecule.cip_stereo(mode="accurate")
    labels = sorted(
        ["bond", *sorted(label["bond_atoms"]), label["descriptor"]]
        if "bond_idx" in label else ["atom", label["atom_idx"], label["descriptor"]]
        for label in values
    )
    return {"labels": labels, "unresolved": molecule.cip_stereo_unresolved()}


def normalized_smiles(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(molecule) if molecule is not None else None


def morgan_result(molecule, generator):
    raw = rdFingerprintGenerator.AdditionalOutput()
    folded = rdFingerprintGenerator.AdditionalOutput()
    raw.AllocateBitInfoMap()
    folded.AllocateBitInfoMap()
    counts = generator.GetSparseCountFingerprint(molecule, additionalOutput=raw)
    bits = generator.GetFingerprint(molecule, additionalOutput=folded)
    return (DataStructs.BitVectToBinaryText(bits), counts.GetNonzeroElements(),
            raw.GetBitInfoMap(), folded.GetBitInfoMap())


def normalize_morgan(values):
    bits, counts, raw, folded = values
    return [bytes(bits).hex(), counts,
            {key: sorted(tuple(env) for env in envs) for key, envs in raw.items()},
            {key: sorted(tuple(env) for env in envs) for key, envs in folded.items()}]


def compare(left, right, equal=lambda a, b: a == b):
    if left["status"] == right["status"] == "ok":
        return "match" if equal(left["value"], right["value"]) else "mismatch"
    return f'{left["status"]}/{right["status"]}'


def time_blocks(candidate, oracle, repeats):
    candidate()()
    oracle()()
    blocks = []
    for i in range(repeats):
        order = ("chematic", "rdkit") if i % 2 == 0 else ("rdkit", "chematic")
        functions = {"chematic": candidate, "rdkit": oracle}
        block = {"order": order}
        for name in order:
            run = functions[name]()
            start = perf_counter()
            result = run()
            block[name] = perf_counter() - start
            del result
        blocks.append(block)
    return {"blocks": blocks, "ratio": paired_median_ci(
        [b["chematic"] for b in blocks], [b["rdkit"] for b in blocks]
    )}


def mcs_probe(chematic_build, rdkit_build):
    def describe(smarts, timed_out, members):
        query = Chem.MolFromSmarts(smarts) if smarts else None
        return {"atoms": query.GetNumAtoms() if query else 0,
                "bonds": query.GetNumBonds() if query else 0,
                "smarts": smarts, "timed_out": timed_out,
                "matches": [m.HasSubstructMatch(query) for m in members] if query else []}

    for strings in MCS_INPUTS:
        for order in (strings, strings[::-1]):
            oracle_mols = list(map(rdkit_build, order))

            def candidate(inputs):
                mols = list(map(chematic_build, inputs))
                if hasattr(chematic, "find_mcs_smarts_checked"):
                    smarts, limited = chematic.find_mcs_smarts_checked(mols, timeout_ms=1000)
                    return describe(smarts, limited, oracle_mols)
                mol, limited = chematic.find_mcs_checked(mols, timeout_ms=1000)
                concrete = Chem.MolFromSmiles(mol.smiles) if mol is not None else None
                result = describe(Chem.MolToSmarts(concrete) if concrete else None,
                                  limited, oracle_mols)
                result["output_kind"] = "concrete_projection"
                return result

            def oracle(mols):
                result = rdFMCS.FindMCS(mols, timeout=1)
                return describe(result.smartsString, result.canceled, mols)

            yield {"input": order, "chematic": capture(candidate, order),
                   "rdkit": capture(oracle, oracle_mols)}


def search_probe(rows, joint, repeats):
    def ids(flags):
        if len(flags) != len(joint):
            raise ValueError("match result dropped a library member")
        return [joint[i][0] for i, flag in enumerate(flags) if flag]

    batch = chematic.MoleculeBatch([cm for _, cm, _ in joint])
    library = rdSubstructLibrary.SubstructLibrary(
        rdSubstructLibrary.MolHolder(), rdSubstructLibrary.PatternHolder()
    )
    for _, _, molecule in joint:
        library.AddMol(molecule)

    for smarts in QUERIES:
        use_chirality = "@" in smarts
        query = Chem.MolFromSmarts(smarts)
        compiled = capture(chematic.compile_smarts, smarts)

        def oracle(q):
            hits = library.GetMatches(q, useChirality=use_chirality, numThreads=1,
                                      maxResults=len(joint))
            return sorted(joint[i][0] for i in hits)

        def candidate(q):
            return ids(q.matches_batch(batch, use_chirality=use_chirality))

        left = capture(candidate, compiled["value"]) if compiled["status"] == "ok" else compiled
        right = capture(oracle, query)
        result = compare(left, right)
        yield {"query": smarts, "chematic": left, "rdkit": right, "comparison": result,
               "use_chirality": use_chirality,
               "timing": time_blocks(lambda: lambda: candidate(compiled["value"]),
                                     lambda: lambda: oracle(query), repeats)
               if repeats and result == "match" and len(rows) == len(joint) else None}


def similarity_probe(rows, joint, generator, repeats):
    dataset = [smiles for _, smiles in rows]
    index = chematic.PreparedFingerprintIndex.from_smiles(dataset, fp="rdkit_ecfp4")
    fingerprints = generator.GetFingerprints([rm for _, _, rm in joint], numThreads=1)
    candidate_failed = index.failed_indices()
    records = []
    for position, smiles in rows[: min(16, len(rows))]:
        def candidate(smi):
            return [(rows[i][0], value) for i, value in index.search(smi, k=10)]

        def oracle(smi):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            scores = DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(mol), fingerprints)
            return heapq.nsmallest(
                10, ((row[0], score) for row, score in zip(joint, scores, strict=True)),
                key=lambda pair: (-pair[1], pair[0])
            )

        def equal(left, right):
            return len(left) == len(right) and all(
                a == b and abs(x - y) <= 1e-12 for (a, x), (b, y) in zip(left, right)
            )

        left = capture(candidate, smiles)
        right = capture(oracle, smiles)
        result = compare(left, right, equal)
        records.append({"index": position, "chematic": left, "rdkit": right,
                        "comparison": result, "timing": time_blocks(
                            lambda: lambda: candidate(smiles), lambda: lambda: oracle(smiles), repeats)
                        if repeats and result == "match" and not candidate_failed
                        and len(rows) == len(joint) else None})
    return {"candidate_failed_rows": [rows[i][0] for i in candidate_failed],
            "oracle_library_rows": len(joint), "queries": records}


def score_rows(rows, out, generator):
    joint = []
    counts = Counter()
    for index, smiles in rows:
        cm = capture(chematic.from_smiles, smiles)
        rm = capture(Chem.MolFromSmiles, smiles)
        parsed = {"chematic": cm["status"], "rdkit": rm["status"]}
        result = {"index": index, "smiles": smiles, "parse": parsed}
        counts[f'parse:{parsed["chematic"]}/{parsed["rdkit"]}'] += 1

        lanes = {
            "cleanup": (lambda m: m.cleanup().rdkit_smiles,
                        lambda m: Chem.MolToSmiles(rdMolStandardize.Cleanup(m))),
            "cip": (chematic_cip, cip_labels),
            "morgan": (lambda m: normalize_morgan(m.rdkit_ecfp4_detail()),
                       lambda m: normalize_morgan(morgan_result(m, generator))),
        }
        for name, (candidate, oracle) in lanes.items():
            left = capture(candidate, cm["value"]) if cm["status"] == "ok" else cm
            right = capture(oracle, rm["value"]) if rm["status"] == "ok" else rm
            if name == "cip":
                comparison = compare(left, right, lambda a, b: a["labels"] == b)
                if right["status"] == "ok":
                    native_api = capture(lambda m: cip_labels(m, clear=False), Chem.MolFromSmiles(smiles))
                    if compare(native_api, right) != "match":
                        counts["cip:rdkit_retained_legacy"] += 1
                        result["cip_with_existing_labels"] = native_api
                if left["status"] == "ok" and left["value"]["unresolved"]:
                    counts["cip:unresolved"] += 1
            else:
                comparison = compare(left, right)
            result[name] = {"chematic": left, "rdkit": right, "comparison": comparison}
            if name == "cleanup" and comparison == "mismatch":
                a = normalized_smiles(left["value"])
                b = normalized_smiles(right["value"])
                result[name]["reparse_equal"] = a is not None and a == b
            counts[f"{name}:{comparison}"] += 1
        if cm["status"] == rm["status"] == "ok":
            joint.append((index, cm["value"], rm["value"]))
        out.write(json.dumps(result, sort_keys=True, allow_nan=False) + "\n")
    return joint, counts


def molecule_timings(rows, counts, generator, repeats):
    smiles = [s for _, s in rows]
    c_parse, r_parse = chematic.from_smiles, Chem.MolFromSmiles
    morgan_candidate = lambda mols: [m.rdkit_ecfp4() for m in mols]
    morgan_oracle = lambda mols: generator.GetFingerprints(mols, numThreads=1)
    lanes = {
        "cleanup": (lambda mols: [m.cleanup().rdkit_smiles for m in mols],
                    lambda mols: [Chem.MolToSmiles(rdMolStandardize.Cleanup(m)) for m in mols]),
        "cip": (lambda mols: [m.cip_stereo(mode="accurate") for m in mols],
                lambda mols: [rdCIPLabeler.AssignCIPLabels(m) for m in mols]),
        "morgan": (morgan_candidate, morgan_oracle),
    }
    timing = {}
    for lane, (candidate, oracle) in lanes.items():
        if counts[f"{lane}:match"] != len(rows):
            continue
        if lane == "cip" and (counts["cip:unresolved"] or counts["cip:rdkit_retained_legacy"]):
            continue
        for mode in ("first_use", "parse_inclusive"):
            def prepare(operation, parse):
                if mode == "first_use":
                    mols = list(map(parse, smiles))
                    return lambda: operation(mols)
                return lambda: operation(list(map(parse, smiles)))
            timing[f"{lane}:{mode}"] = time_blocks(
                lambda: prepare(candidate, c_parse), lambda: prepare(oracle, r_parse), repeats
            )
    return timing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--rdkit-version", required=True)
    parser.add_argument("--chematic-wheel", type=Path, required=True)
    parser.add_argument("--rdkit-wheel", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--repeats", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.repeats and (args.repeats < 20 or args.repeats % 2):
        parser.error("use zero or an even number of timing blocks of at least 20")
    if args.repeats and os.environ.get("RAYON_NUM_THREADS") != "1":
        parser.error("set RAYON_NUM_THREADS=1 for single-threaded timing")
    if rdBase.rdkitVersion != args.rdkit_version:
        parser.error(f"loaded RDKit {rdBase.rdkitVersion}; expected {args.rdkit_version}")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    RDLogger.DisableLog("rdApp.*")
    data = args.corpus.read_text(encoding="utf-8")
    rows = [(i, line.split()[0]) for i, line in enumerate(data.splitlines(), 1)
            if line.strip() and not line.lstrip().startswith("#")][:args.limit]
    if not rows:
        parser.error("corpus contains no input rows")
    c_identity = package_identity(args.chematic_wheel, "chematic", metadata.version("chematic"))
    r_identity = package_identity(args.rdkit_wheel, "rdkit", rdBase.rdkitVersion)
    args.output.mkdir(parents=True, exist_ok=False)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    with gzip.open(args.output / "rows.jsonl.gz", "wt", encoding="utf-8") as out:
        joint, counts = score_rows(rows, out, generator)
    report = {
        "corpus": {"path": str(args.corpus), "sha256": digest(args.corpus), "rows": len(rows)},
        "script_sha256": digest(Path(__file__)),
        "rows_sha256": digest(args.output / "rows.jsonl.gz"),
        "runtime": {"python": sys.version, "cpu": cpu_model(), "host": platform.platform()},
        "chematic": {**c_identity, "revision": args.revision},
        "rdkit": r_identity,
        "options": {"morgan_radius": 2, "morgan_bits": 2048, "morgan_chirality": False,
                    "rdkit_library_threads": 1, "rayon_threads": os.environ.get("RAYON_NUM_THREADS"),
                    "repeats": args.repeats},
        "counts": counts, "joint_parsed": len(joint),
        "substructure": list(search_probe(rows, joint, args.repeats)) if joint else [],
        "similarity": similarity_probe(rows, joint, generator, args.repeats) if joint else None,
        "mcs": list(mcs_probe(chematic.from_smiles, Chem.MolFromSmiles)),
        "molecule_timings": molecule_timings(rows, counts, generator, args.repeats)
        if args.repeats else {},
    }
    (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
