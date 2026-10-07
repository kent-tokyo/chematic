#!/usr/bin/env python3
"""Write (or check) the BioTransformer corpus's reactant requests.

For each of the corpus's 400 reactants that RDKit 2026.03.6 reads, the
implicit-H and explicit-H SMILES that ``biotransformer_rule_corpus.py`` passes
to chematic, one tab-separated row per reactant. Pinned as
``validation/biotransformer-requests-400.tsv`` so that chematic's side of the
corpus runs without RDKit (``biotransformer_chematic_responses.py``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "validation/biotransformer-requests-400.tsv"


def requests(reactants: Path) -> list[tuple[str, str]]:
    rows = []
    for line in reactants.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        mol = Chem.MolFromSmiles(line.split()[0])
        if mol is None:
            continue
        # As biotransformer_rule_corpus.py: canonical SMILES read back.
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        rows.append((Chem.MolToSmiles(mol), Chem.MolToSmiles(Chem.AddHs(mol))))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reactants", type=Path, required=True)
    ap.add_argument("--output", type=Path, default=DEFAULT)
    ap.add_argument("--check", action="store_true", help="compare with --output instead of writing it")
    args = ap.parse_args()
    if rdBase.rdkitVersion != "2026.03.6":
        ap.error(f"RDKit 2026.03.6 required, got {rdBase.rdkitVersion}")
    text = "# implicit-H\texplicit-H reactant SMILES written by RDKit 2026.03.6\n" + "".join(
        f"{a}\t{b}\n" for a, b in requests(args.reactants))
    if args.check:
        same = args.output.read_text(encoding="utf-8") == text
        print("requests match" if same else "requests differ")
        return 0 if same else 1
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
