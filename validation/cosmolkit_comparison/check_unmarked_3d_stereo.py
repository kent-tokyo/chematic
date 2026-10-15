#!/usr/bin/env python3
"""Compare stereo recovered from unmarked RDKit 3D V2000/V3000 SDF records.

Example (Python environment with RDKit and COSMolKit installed):
  cargo build -p chematic-mol --example sdf_v3000_accuracy_dump --offline
  python check_unmarked_3d_stereo.py --chematic ../../target/debug/examples/sdf_v3000_accuracy_dump \\
      --output-prefix /tmp/stereo-a --seed 20261016 --per-category 50
  python check_unmarked_3d_stereo.py --chematic ../../target/debug/examples/sdf_v3000_accuracy_dump \\
      --output-prefix /tmp/stereo-b --seed 20261017 --per-category 100 \\
      --exclude-manifest /tmp/stereo-a-manifest.json

Only molecules whose 3D blocks RDKit itself reads back with the source's
canonical isomeric SMILES are compared. V2000 auto-fallback to V3000 is
rejected. An engine's output is normalized by the same RDKit SMILES parser.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import subprocess
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(smiles: str | None) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None


def unmarked_block(mol: Chem.Mol, v3000: bool) -> str | None:
    lines = Chem.MolToMolBlock(mol, includeStereo=False, forceV3000=v3000).splitlines()
    if not lines[3].endswith("V3000" if v3000 else "V2000"):
        return None
    if v3000:
        lines = [re.sub(r"\sCFG=\d+", "", line) for line in lines]
    else:
        atoms, bonds = int(lines[3][:3]), int(lines[3][3:6])
        for line in range(4, 4 + atoms):
            lines[line] = lines[line][:39] + "  0" + lines[line][42:]
        for line in range(4 + atoms, 4 + atoms + bonds):
            lines[line] = lines[line][:9] + "  0" + lines[line][12:]
    return "\n".join(lines) + "\n"


def prepare(corpus: Path, seed: int, per_category: int, excluded: set[int]):
    sources = [line.split()[0] for line in corpus.read_text().splitlines() if line.strip()]
    indices = [index for index in range(len(sources)) if index not in excluded]
    random.Random(seed).shuffle(indices)
    rows: list[dict] = []
    blocks: dict[str, list[str]] = {"v2000": [], "v3000": []}
    counts = {"tetra": 0, "alkene": 0}
    for index in indices:
        source = sources[index]
        if any(tag in source for tag in ("@SP", "@TB", "@OH")):
            continue
        category = "tetra" if "@" in source else "alkene" if "/" in source or "\\" in source else None
        if category is None or counts[category] >= per_category:
            continue
        mol = Chem.MolFromSmiles(source)
        if mol is None or not 3 <= mol.GetNumAtoms() <= 35:
            continue
        expected = Chem.MolToSmiles(mol, isomericSmiles=True)
        with_h = Chem.AddHs(mol)
        settings = AllChem.ETKDGv3()
        settings.randomSeed = 20261016 + index
        settings.maxIterations = 30
        if AllChem.EmbedMolecule(with_h, settings) != 0:
            continue
        mol = Chem.RemoveHs(with_h)
        candidate = {}
        for fmt, v3000 in (("v2000", False), ("v3000", True)):
            block = unmarked_block(mol, v3000)
            if block is None:
                break
            reread = Chem.MolFromMolBlock(block, removeHs=False)
            if reread is None or Chem.MolToSmiles(reread, isomericSmiles=True) != expected:
                break
            candidate[fmt] = block + "$$$$\n"
        if len(candidate) != 2:
            continue
        rows.append({"source_index": index, "source": source, "expected": expected, "category": category})
        counts[category] += 1
        for fmt in blocks:
            blocks[fmt].append(candidate[fmt])
        if all(count == per_category for count in counts.values()):
            break
    if any(count != per_category for count in counts.values()):
        raise RuntimeError(f"Insufficient source-gated records: {counts}")
    return rows, {fmt: "".join(part) for fmt, part in blocks.items()}, counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=Path(__file__).parent / "corpora/stress-v1.smi")
    parser.add_argument("--chematic", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--per-category", type=int, required=True)
    parser.add_argument("--exclude-manifest", type=Path)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")
    excluded = set()
    if args.exclude_manifest:
        excluded = {row["source_index"] for row in json.loads(args.exclude_manifest.read_text())["rows"]}
    rows, blocks, counts = prepare(args.corpus, args.seed, args.per_category, excluded)
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    binary = args.chematic.resolve()
    manifest = {
        "corpus_sha256": sha256(args.corpus.read_bytes()),
        "seed": args.seed,
        "excluded_source_indices": sorted(excluded),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "counts": counts,
        "rows": rows,
        "sdf_sha256": {fmt: sha256(sdf.encode()) for fmt, sdf in blocks.items()},
        "chematic_binary_sha256": sha256(binary.read_bytes()),
    }
    Path(f"{prefix}-manifest.json").write_text(json.dumps(manifest, indent=2))
    report = {"manifest": str(prefix) + "-manifest.json", "results": {}}
    for fmt, sdf in blocks.items():
        path = Path(f"{prefix}-{fmt}.sdf")
        path.write_text(sdf)
        completed = subprocess.run([str(binary), str(path)], capture_output=True, text=True, check=True)
        chematic_rows = json.loads(completed.stdout)
        cosmolkit_rows = cosmolkit.SdfDataset.open(str(path))
        if len(chematic_rows) != len(rows) or len(cosmolkit_rows) != len(rows):
            raise AssertionError("Reader record count differs from source-gated count")
        hits = {engine: {category: 0 for category in counts} for engine in ("chematic", "cosmolkit")}
        misses = {engine: [] for engine in hits}
        for index, row in enumerate(rows):
            try:
                other = cosmolkit_rows[index].molecule().to_smiles()
                other_error = None
            except Exception as error:
                other, other_error = None, str(error)
            observations = (
                ("chematic", chematic_rows[index].get("smiles"), chematic_rows[index].get("error")),
                ("cosmolkit", other, other_error),
            )
            for engine, raw, error in observations:
                actual = canonical(raw)
                if actual == row["expected"]:
                    hits[engine][row["category"]] += 1
                else:
                    misses[engine].append({"index": index, "source_index": row["source_index"],
                                           "expected": row["expected"], "actual": actual,
                                           "raw": raw, "error": error})
        report["results"][fmt] = {"hits": hits, "misses": misses}
    Path(f"{prefix}-report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({fmt: result["hits"] for fmt, result in report["results"].items()}))


if __name__ == "__main__":
    main()
