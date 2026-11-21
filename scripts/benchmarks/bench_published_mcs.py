#!/usr/bin/env python3
"""Measure prepared MCS search in fresh processes with rotating package order."""

import argparse
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import time
import zipfile

CASES = [
    ["c1ccccc1", "Cc1ccccc1"],
    ["C1CC1", "CCC"],
    ["CC(O)C", "CC(=O)C", "CC(N)C"],
    ["C1CC2CCC1C2", "C1CCC2CCCCC2C1"],
    ["c1ccc2ccccc2c1", "c1ccc2c(c1)cco2"],
]
ARMS = ("v120", "v121", "v122", "rdkit", "cosmolkit")
VERSIONS = {"v120": ("chematic", "1.2.0"), "v121": ("chematic", "1.2.1"),
            "v122": ("chematic", "1.2.2"), "rdkit": ("rdkit", "2026.9.2"),
            "cosmolkit": ("cosmolkit", "0.6.2")}


def verify_wheel(arm, wheel_dir, installed_root):
    package, version = VERSIONS[arm]
    wheels = list(wheel_dir.glob(f"{package}-{version}-*.whl"))
    if len(wheels) != 1:
        raise ValueError(f"expected one {package} {version} wheel, found {len(wheels)}")
    wheel = wheels[0]
    with zipfile.ZipFile(wheel) as archive:
        files = [name for name in archive.namelist() if name.startswith(package + "/") and not name.endswith("/")]
        if not files:
            raise ValueError(f"no package files in {wheel}")
        for name in files:
            installed = installed_root / name
            if not installed.is_file() or installed.read_bytes() != archive.read(name):
                raise ValueError(f"installed file differs from wheel: {installed}")
    return {"filename": wheel.name, "sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
            "installed_files_verified": len(files)}


def worker(arm):
    start = time.perf_counter_ns()
    if arm.startswith("v"):
        lib = importlib.import_module("chematic")
        parse = lib.from_smiles
        def search(mols):
            result, limited = lib.find_mcs_checked(mols, timeout_ms=2000)
            return (result.smiles if result is not None else None), limited, "smiles"
    elif arm == "rdkit":
        lib = importlib.import_module("rdkit")
        from rdkit import Chem
        from rdkit.Chem import rdFMCS
        parse = Chem.MolFromSmiles
        def search(mols):
            result = rdFMCS.FindMCS(mols, timeout=2, maximizeBonds=True)
            return result.smartsString, result.canceled, "smarts"
    else:
        lib = importlib.import_module("cosmolkit")
        parse = lib.Molecule.from_smiles
        def search(mols):
            result = lib.maximum_common_substructure(mols, timeout=2, maximize_bonds=True)
            return result.smarts, not result.completed, "smarts"
    import_ns = time.perf_counter_ns() - start
    start = time.perf_counter_ns()
    molsets = [[parse(smiles) for smiles in case] for case in CASES]
    parse_ns = time.perf_counter_ns() - start
    results, first_ns, repeated, repeated_ns = [], [], [], []
    for mols in molsets:
        start = time.perf_counter_ns()
        results.append(search(mols))
        first_ns.append(time.perf_counter_ns() - start)
    for mols in molsets:
        start = time.perf_counter_ns()
        repeated.append(search(mols))
        repeated_ns.append(time.perf_counter_ns() - start)
    if results != repeated:
        raise ValueError("MCS result differs between repeated searches")
    scale = 1 if sys.platform == "darwin" else 1024
    return {"version": lib.__version__, "module": lib.__file__, "import_ns": import_ns,
            "parse_ns": parse_ns, "first_ns": first_ns, "repeated_ns": repeated_ns,
            "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale,
            "results": results}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=ARMS)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--wheel-dir", type=Path)
    for arm in ARMS[:3]:
        parser.add_argument(f"--{arm}-path", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker)))
        return
    if args.blocks < 20 or not args.output or any(getattr(args, f"{arm}_path") is None for arm in ARMS[:3]):
        parser.error("20 blocks, output and all three installed package paths are required")
    report = {"cases": CASES, "case_sha256": hashlib.sha256(json.dumps(CASES).encode()).hexdigest(),
              "host": platform.platform(), "python": sys.version, "blocks": [], "complete": False}
    if args.validate_only:
        report = json.loads(args.output.read_text())
        if report["cases"] != CASES or len(report["blocks"]) != args.blocks:
            raise ValueError("saved inputs or block count differ")
    else:
        if args.wheel_dir is None:
            parser.error("wheel directory is required for measurement")
        verified = {}
        for arm in ARMS:
            if arm.startswith("v"):
                root = getattr(args, f"{arm}_path")
            else:
                spec = importlib.util.find_spec(arm)
                if spec is None or spec.origin is None:
                    raise ValueError(f"cannot locate {arm}")
                root = Path(spec.origin).parent.parent
            verified[arm] = verify_wheel(arm, args.wheel_dir, root)
        report["wheel_verification"] = verified
    for index in range(0 if args.validate_only else args.blocks):
        order = ARMS[index % len(ARMS):] + ARMS[:index % len(ARMS)]
        block = {"order": order, "load_average": os.getloadavg()}
        for arm in order:
            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            path = getattr(args, f"{arm}_path", None)
            if path:
                env["PYTHONPATH"] = str(path)
            result = subprocess.run([sys.executable, __file__, "--worker", arm], env=env,
                                    text=True, capture_output=True, timeout=120, check=True)
            block[arm] = json.loads(result.stdout)
        report["blocks"].append(block)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"block {index + 1}/{args.blocks}", flush=True)
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.error")
    for block in report["blocks"]:
        for arm in ARMS:
            validations = []
            for case, (query, limited, kind) in zip(CASES, block[arm]["results"], strict=True):
                if not query:
                    validations.append({"error": "no MCS output", "limited": limited})
                    continue
                mol = Chem.MolFromSmiles(query) if kind == "smiles" else Chem.MolFromSmarts(query)
                if mol is None:
                    validations.append({"error": "unparseable output", "limited": limited})
                    continue
                validations.append({"atoms": mol.GetNumAtoms(), "bonds": mol.GetNumBonds(),
                                    "matches": [Chem.MolFromSmiles(s).HasSubstructMatch(mol) for s in case],
                                    "limited": limited})
            block[arm]["validation"] = validations
    report["summary"] = {arm: {key: statistics.median(b[arm][key] for b in report["blocks"])
                                for key in ("import_ns", "parse_ns", "peak_rss_bytes")}
                         for arm in ARMS}
    report["case_summary"] = {}
    for arm in ARMS:
        cases = []
        for index in range(len(CASES)):
            comparable = []
            for block in report["blocks"]:
                actual = block[arm]["validation"][index]
                reference = block["rdkit"]["validation"][index]
                comparable.append(
                    "error" not in actual and "error" not in reference
                    and not actual["limited"] and not reference["limited"]
                    and all(actual["matches"]) and all(reference["matches"])
                    and actual["atoms"] == reference["atoms"]
                    and actual["bonds"] == reference["bonds"]
                )
            cases.append({"case": index, "passes_all_blocks": all(comparable),
                          "first_ns_median": statistics.median(b[arm]["first_ns"][index] for b in report["blocks"]),
                          "repeat_ns_median": statistics.median(b[arm]["repeated_ns"][index] for b in report["blocks"])})
        report["case_summary"][arm] = cases
    report["protocol"] = {
        "order": f"{args.blocks} rotating blocks of five fresh processes",
        "timeout_seconds": 2,
        "search": "prepared molecules; first and second calls timed per case",
        "validation": "RDKit checks output parseability, atom and bond counts, and matches against inputs after timing",
        "memory": "whole process peak RSS",
        "note": "The chematic API returns a concrete molecule; RDKit and COSMolKit return SMARTS. "
                "RDKit parse includes first-use initialization. Matching counts and input matches "
                "do not prove query semantic equivalence. Invalid or differing results cannot support speed ratios.",
    }
    report["complete"] = True
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
