#!/usr/bin/env python3
"""Check Cleanup output across deterministic rooted spellings of eight inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "validation/results/rdkit-cleanup-v1.2.2-source-stress-2026-11-21.json"
INDICES = (65, 83, 101, 809, 911, 1664, 2076, 2855)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--all-chemical-differences", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--candidate-commit")
    parser.add_argument("--indices", help="comma-separated stress corpus indices")
    args = parser.parse_args()
    if chematic.__version__ != "1.2.2" or rdBase.rdkitVersion != "2026.09.1":
        raise ValueError("expected chematic 1.2.2 and RDKit 2026.09.1")
    package_root = Path(chematic.__file__).parent.parent
    with zipfile.ZipFile(args.wheel) as archive:
        files = [name for name in archive.namelist() if name.startswith("chematic/") and not name.endswith("/")]
        if not files:
            raise ValueError("wheel contains no chematic package files")
        for name in files:
            installed = package_root / name
            if not installed.is_file() or installed.read_bytes() != archive.read(name):
                raise ValueError(f"installed package differs from wheel: {installed}")
    RDLogger.DisableLog("rdApp.*")
    cases = {row["index"]: row for row in json.loads(INPUT.read_text())["mismatches"]}
    records = []
    indices = INDICES
    if args.indices:
        if args.all_chemical_differences or args.output is None:
            parser.error("explicit indices require --output and cannot be combined with --all-chemical-differences")
        corpus = (ROOT / "validation/cosmolkit_comparison/corpora/stress-v1.smi").read_text().splitlines()
        indices = tuple(int(index) for index in args.indices.split(","))
        cases = {index: {"input": corpus[index].split()[0]} for index in indices}
    elif args.all_chemical_differences:
        audit = ROOT / "validation/results/cleanup-v1.2.2-47-atom-audit-2026-11-21.json"
        indices = tuple(row["index"] for row in json.loads(audit.read_text())["rows"]
                        if (not row["charge_h_connectivity_isomorphic"]
                            or row["bond_isotope_graph_isomorphic"] is False))
    for index in indices:
        mol = Chem.MolFromSmiles(cases[index]["input"])
        forms = sorted({Chem.MolToSmiles(mol, canonical=False, rootedAtAtom=i, isomericSmiles=True)
                        for i in range(mol.GetNumAtoms())})
        variants = []
        excluded = []
        reference_input = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        for form in forms:
            reparsed = Chem.MolFromSmiles(form)
            if reparsed is None:
                raise ValueError(f"rooted spelling did not parse: {index}: {form}")
            try:
                parsed = chematic.from_smiles(form)
                if (Chem.MolToSmiles(reparsed, canonical=True, isomericSmiles=True) != reference_input
                        or parsed.rdkit_smiles != reference_input):
                    excluded.append({"input": form, "reason": "input canonical representations differ"})
                    continue
            except ValueError as error:
                excluded.append({"input": form, "reason": str(error)})
                continue
            try:
                output = parsed.cleanup().rdkit_smiles
            except chematic.ChematicInputError as error:
                if error.category != "ambiguous" or error.code not in ("cleanup_metal_hydrogen", "cleanup_unstable_metal"):
                    raise
                output = None
            variants.append({"input": form,
                             "rdkit": Chem.MolToSmiles(rdMolStandardize.Cleanup(reparsed), isomericSmiles=True),
                             "chematic": output})
        records.append({"index": index, "inputs": len(forms),
                        "excluded": excluded,
                        "rdkit_distinct_outputs": len({row["rdkit"] for row in variants}),
                        "chematic_distinct_outputs": len({row["chematic"] for row in variants if row["chematic"] is not None}),
                        "chematic_refusals": sum(row["chematic"] is None for row in variants),
                        "variants": variants})
    report = {"method": "one deterministic noncanonical SMILES rooted at each atom; duplicate strings removed",
              "rdkit_version": rdBase.rdkitVersion, "chematic_version": chematic.__version__,
              "wheel": args.wheel.name, "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              "installed_package_files_verified": len(files),
              "candidate_commit": args.candidate_commit,
              "source_record": str(INPUT.relative_to(ROOT)), "cases": records}
    filename = ("cleanup-rooted-smiles-31-v1.2.2-2026-11-21.json"
                if args.all_chemical_differences else "cleanup-rooted-smiles-v1.2.2-2026-11-21.json")
    output = args.output or ROOT / "validation/results" / filename
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print([(r["index"], r["inputs"], r["rdkit_distinct_outputs"], r["chematic_distinct_outputs"]) for r in records])


if __name__ == "__main__":
    main()
