#!/usr/bin/env python3
"""Dump exact R2/R3 Morgan bytes, with and without stereo, for the frontier packet.

Select the source wheel through PYTHONPATH; run RDKit in a separate process.
Exceptions remain per-row errors and never count as agreement.
"""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=("chematic", "rdkit", "cosmolkit"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path('/Users/k_tanabe/.codex/worktrees/rdkit-cosmolkit-speed-20261017/chematic')
    corpora = [root / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
               root / "validation/cosmolkit_comparison/corpora/stress-v1.smi"]
    ops = {}
    if args.engine == "chematic":
        import chematic
        parse = chematic.from_smiles
        for radius in (2, 3):
            for chiral in (False, True):
                ops[f"morgan_r{radius}_chiral{int(chiral)}"] = (
                    lambda mol, r=radius, s=chiral: mol.rdkit_ecfp_config(
                        r, 2048, include_chirality=s).hex()
                )
        version = chematic.__version__
    elif args.engine == "cosmolkit":
        import cosmolkit as ck
        parse = ck.Molecule.from_smiles
        for radius in (2, 3):
            for chiral in (False, True):
                generator = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=radius, fp_size=2048, include_chirality=chiral))
                def evaluate(mol, gen=generator):
                    bits = bytearray(256)
                    for bit in mol.fingerprint_morgan_with_generator(gen).on_bits():
                        bits[bit // 8] |= 1 << (bit % 8)
                    return bits.hex()
                ops[f"morgan_r{radius}_chiral{int(chiral)}"] = evaluate
        version = ck.__version__
    else:
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator
        RDLogger.DisableLog("rdApp.*")
        parse = Chem.MolFromSmiles
        for radius in (2, 3):
            for chiral in (False, True):
                generator = rdFingerprintGenerator.GetMorganGenerator(
                    radius=radius, fpSize=2048, includeChirality=chiral)
                ops[f"morgan_r{radius}_chiral{int(chiral)}"] = (
                    lambda mol, gen=generator: DataStructs.BitVectToBinaryText(
                        gen.GetFingerprint(mol)).hex()
                )
        version = Chem.rdBase.rdkitVersion
    with args.output.open("w") as output:
        for path in corpora:
            data = path.read_bytes()
            rows = [line.split()[0] for line in data.decode().splitlines() if line.strip()]
            output.write(json.dumps({
                "corpus": path.name, "sha256": hashlib.sha256(data).hexdigest(),
                "rows": len(rows), "engine": args.engine, "version": version,
            }) + "\n")
            for index, smiles in enumerate(rows):
                row = {"corpus": path.name, "index": index, "values": {}}
                try:
                    molecule = parse(smiles)
                    if molecule is None:
                        raise ValueError("parse failed")
                    for name, evaluate in ops.items():
                        try:
                            row["values"][name] = {"value": evaluate(molecule)}
                        except Exception as error:
                            row["values"][name] = {"error": type(error).__name__ + ": " + str(error)}
                except Exception as error:
                    row["parse_error"] = type(error).__name__ + ": " + str(error)
                output.write(json.dumps(row, sort_keys=True) + "\n")
            print(path.name, len(rows), "done", flush=True)


if __name__ == "__main__":
    main()
