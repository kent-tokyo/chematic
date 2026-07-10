#!/usr/bin/env python3
"""
Canonical SMILES structural correctness: does canonical_smiles(parse(x))
represent the SAME molecule as x, for every one of N independently-traversed
input spellings of the same molecule?

This is a stricter, different question than self-stability/idempotency
(does canonical_smiles(parse(x)) == canonical_smiles(parse(canonical_smiles(
parse(x))))). A canonicalization bug can be perfectly idempotent (same wrong
output every time from a given input) while still being WRONG -- it silently
returns a *different stereoisomer* than the one that was parsed. Idempotency
tests cannot detect this class of bug at all; only comparing against the
original molecule's identity (via an independent oracle) can.

Oracle: RDKit's own canonical form (Chem.MolToSmiles), used only to decide
"same molecule or not" -- never to judge whether chematic's specific string
choice is "right," since canonical form is implementation-defined.

Comparison is strict full-isomeric-SMILES equality (Chem.MolToSmiles(m) ==
Chem.MolToSmiles(rd)), not a formula+atomcount fallback -- a fallback that
accepts "same formula, same atom count" as "same molecule" would silently
pass every stereo-inversion case, since inverting stereochemistry never
changes formula or atom count.

Usage:
    python scripts/canonical_structural_correctness.py [SMILES.csv] [-n N_VARIANTS] [--limit N]
"""
import sys
import os
import argparse

import chematic
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def load_smiles(path, limit=None):
    path = os.path.expanduser(path)
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    smis = [
        line.split(",")[0].strip()
        for line in lines
        if line.split(",")[0].strip().lower() != "smiles"
    ]
    return smis[:limit] if limit else smis


def check_molecule(rd, n_variants, max_examples_per_mol=2):
    """Run n_variants independently-traversed spellings of `rd` through
    chematic's canonical_smiles and check each represents the same molecule.
    Returns (is_bad, failure_examples)."""
    orig_canon = Chem.MolToSmiles(rd)
    orig_formula = rdMolDescriptors.CalcMolFormula(rd)
    orig_nostereo = Chem.MolToSmiles(rd, isomericSmiles=False)

    is_bad = False
    examples = []
    for _ in range(n_variants):
        variant = Chem.MolToSmiles(rd, doRandom=True)
        try:
            out = chematic.from_smiles(variant).smiles
        except Exception as e:
            is_bad = True
            if len(examples) < max_examples_per_mol:
                examples.append((variant, None, f"chematic raised: {e}"))
            continue
        m2 = Chem.MolFromSmiles(out)
        if m2 is None:
            is_bad = True
            if len(examples) < max_examples_per_mol:
                examples.append((variant, out, "RDKit could not parse chematic's output"))
            continue
        rt_canon = Chem.MolToSmiles(m2)
        if rt_canon != orig_canon:
            is_bad = True
            if len(examples) < max_examples_per_mol:
                rt_formula = rdMolDescriptors.CalcMolFormula(m2)
                rt_nostereo = Chem.MolToSmiles(m2, isomericSmiles=False)
                if rt_formula != orig_formula:
                    kind = f"FORMULA CHANGED: {orig_formula} -> {rt_formula}"
                elif rt_nostereo != orig_nostereo:
                    kind = "SKELETON CHANGED (not just stereo)"
                else:
                    kind = "STEREO INVERTED (same skeleton/formula)"
                examples.append((variant, out, kind))
    return is_bad, examples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("corpus", nargs="?", default="~/Downloads/SMILES.csv")
    parser.add_argument("-n", type=int, default=10, help="variants per molecule")
    parser.add_argument("--limit", type=int, default=None, help="cap corpus size")
    args = parser.parse_args()

    smis = load_smiles(args.corpus, args.limit)

    n_mol = 0
    mol_bad = 0
    all_examples = []

    for smi in smis:
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        n_mol += 1
        is_bad, examples = check_molecule(rd, args.n)
        if is_bad:
            mol_bad += 1
            if len(all_examples) < 15:
                all_examples.append((smi, examples))

    print(f"corpus: {n_mol}, variants per molecule: {args.n}")
    pct = 100 * mol_bad / n_mol if n_mol else 0.0
    print(f"structural-correctness failures (>=1 of {args.n} variants wrong): {mol_bad}/{n_mol} ({pct:.2f}%)")
    print()
    for smi, examples in all_examples:
        print(f"orig: {smi}")
        for variant, out, kind in examples:
            print(f"  [{kind}]")
            print(f"  variant: {variant}")
            print(f"  chematic canonical: {out}")
        print()

    return 0 if mol_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
