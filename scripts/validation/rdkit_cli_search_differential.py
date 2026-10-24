#!/usr/bin/env python3
"""Compare CLI substructure booleans with RDKit on a fixed SMILES corpus.

Run after `cargo build -p chematic-cli`. This measures only the 25 SMARTS
queries below; it does not score arbitrary SMARTS, mapping order, or speed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase


QUERIES = (
    "[OH]", "c", "[#7]", "C=O", "[NX3;H2,H1;!$(NC=O)]", "[r5]", "[r6]",
    "c1ccccc1", "[CX4]", "[!#6;!#1]", "C(=O)[OH]", "[nH]", "[#6]=[#6]",
    "[F,Cl,Br,I]", "[OX2H]", "[#16]", "[c+]", "[C]", "[as]", "[R2]",
    "[x2]", "[D3]", "[H1]", "[N;H2]", "[C;X3]",
)


def score(binary: Path, rows: list[tuple[int, str, str, bool]], parity: bool) -> dict:
    result = {"agree": 0, "mismatch": 0, "rejected": 0, "by_query": {}, "examples": []}
    args = [str(binary), "batch-substructure"]
    if parity:
        args.append("--rdkit-parity")
    for start in range(0, len(rows), 1000):
        batch = rows[start : start + 1000]
        payload = "".join(f"{smiles}\t{query}\n" for _, smiles, query, _ in batch)
        completed = subprocess.run(args, input=payload, text=True, capture_output=True, check=True)
        records = json.loads(completed.stdout)["records"]
        if len(records) != len(batch):
            raise RuntimeError("CLI omitted a search record")
        for (index, smiles, query, expected), record in zip(batch, records, strict=True):
            if record["substructure"] is None:
                status = "rejected"
                actual = None
            else:
                status = "agree" if bool(record["substructure"]["match_count"]) == expected else "mismatch"
                actual = bool(record["substructure"]["match_count"])
            result[status] += 1
            if status != "agree":
                query_counts = result["by_query"].setdefault(query, {"mismatch": 0, "rejected": 0})
                query_counts[status] += 1
                if len(result["examples"]) < 20:
                    result["examples"].append({
                        "index": index, "smiles": smiles, "query": query,
                        "rdkit": expected, "chematic": actual, "error": record["error"],
                    })
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("target/debug/chematic"))
    parser.add_argument(
        "--corpus", type=Path,
        default=Path("validation/cosmolkit_comparison/corpora/stress-v1.smi"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    corpus_bytes = args.corpus.read_bytes()
    queries = [(query, Chem.MolFromSmarts(query)) for query in QUERIES]
    if any(parsed is None for _, parsed in queries):
        raise RuntimeError("an oracle SMARTS failed to parse")
    rows = []
    input_count = 0
    rdkit_rejected = 0
    for index, line in enumerate(corpus_bytes.decode("utf-8").splitlines()):
        if not line.strip() or line.startswith("#"):
            continue
        input_count += 1
        smiles = line.split()[0]
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            rdkit_rejected += 1
            continue
        for query, parsed in queries:
            rows.append((index, smiles, query, bool(molecule.HasSubstructMatch(parsed))))

    report = {
        "operation": "HasSubstructMatch boolean, default chirality off",
        "rdkit_version": rdBase.rdkitVersion,
        "binary": str(args.binary),
        "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
        "corpus": str(args.corpus),
        "corpus_sha256": hashlib.sha256(corpus_bytes).hexdigest(),
        "input_count": input_count,
        "rdkit_rejected": rdkit_rejected,
        "queries": list(QUERIES),
        "comparisons": len(rows),
        "native": score(args.binary, rows, False),
        "rdkit_parity": score(args.binary, rows, True),
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
