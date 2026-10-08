#!/usr/bin/env python3
"""Compare chematic's RDKit-compatible canonical SMILES with RDKit's.

For every SMILES in a corpus (one per line, first whitespace-separated
token), compares ``chematic.from_smiles(s).rdkit_smiles`` with
``Chem.MolToSmiles(Chem.MolFromSmiles(s))`` and prints match / mismatch /
error counts and the first mismatches.

Categories:

* ``match``     both produce the same string;
* ``mismatch``  both produce a string and the strings differ;
* ``error``     RDKit produces a string, chematic raises (unsupported input,
                parse or sanitization failure);
* ``rejected``  RDKit's ``MolFromSmiles`` returns ``None``; counted as a
                match when chematic also raises and as ``mismatch`` when
                chematic returns a string.

Run from the repository root with an environment that has both ``rdkit``
and a ``chematic`` wheel installed::

    python validation/cosmolkit_comparison/check_rdkit_smiles.py \\
        scripts/chembl_accuracy_corpus_4999.smi \\
        validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi

With no corpus arguments those two corpora are used. ``--fail-on-mismatch``
exits with status 1 when any mismatch or error is found.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CORPORA = [
    ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi",
    ROOT / "validation" / "benchmark_corpora" / "rdkit-js-browser-10k-v1.smi",
]


def read_smiles(path: Path) -> list[str]:
    out = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        if tokens:
            out.append(tokens[0])
    return out


def check_corpus(path: Path, max_show: int) -> dict:
    import chematic
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    counts = {"total": 0, "match": 0, "mismatch": 0, "error": 0, "rejected": 0}
    shown = 0
    for smi in read_smiles(path):
        counts["total"] += 1
        rd_mol = Chem.MolFromSmiles(smi)
        expected = Chem.MolToSmiles(rd_mol) if rd_mol is not None else None
        try:
            got = chematic.from_smiles(smi).rdkit_smiles
            failure = None
        except ValueError as exc:
            got = None
            failure = str(exc)
        if expected is None:
            counts["rejected"] += 1
            if got is None:
                counts["match"] += 1
            else:
                counts["mismatch"] += 1
                if shown < max_show:
                    shown += 1
                    print(f"  RDKit rejects {smi}\n    chematic {got}")
        elif got == expected:
            counts["match"] += 1
        elif got is None:
            counts["error"] += 1
            if shown < max_show:
                shown += 1
                print(f"  ERROR    {smi}\n    rdkit    {expected}\n    chematic {failure}")
        else:
            counts["mismatch"] += 1
            if shown < max_show:
                shown += 1
                print(f"  MISMATCH {smi}\n    rdkit    {expected}\n    chematic {got}")
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("corpora", nargs="*", type=Path, help="SMILES files")
    parser.add_argument("--max-show", type=int, default=20,
                        help="mismatches/errors to print per corpus (default 20)")
    parser.add_argument("--fail-on-mismatch", action="store_true",
                        help="exit 1 if any corpus has a mismatch or error")
    args = parser.parse_args()

    import chematic
    from rdkit import rdBase

    print(f"chematic {getattr(chematic, '__version__', '?')}, RDKit {rdBase.rdkitVersion}")
    failed = False
    for path in args.corpora or DEFAULT_CORPORA:
        print(f"== {path}")
        counts = check_corpus(path, args.max_show)
        print(
            f"   total {counts['total']}  match {counts['match']}  "
            f"mismatch {counts['mismatch']}  error {counts['error']}  "
            f"(RDKit rejects {counts['rejected']})"
        )
        failed |= counts["mismatch"] > 0 or counts["error"] > 0
    return 1 if (failed and args.fail_on_mismatch) else 0


if __name__ == "__main__":
    sys.exit(main())
