"""Compare chematic's Mol.rdkit_avalon_fp against RDKit's Avalon fingerprint.

For each SMILES in a corpus, the on-bits of
``rdkit.Avalon.pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles(s), nBits=2048)``
are compared with those of ``chematic.from_smiles(s).rdkit_avalon_fp(2048)``.
Prints match / mismatch / error counts and the first mismatches with the bit
differences.

Usage: python check_rdkit_avalon.py CORPUS.smi [--show N] [--nbits N]
"""

import argparse
import sys

import chematic
from rdkit import Chem, RDLogger
from rdkit.Avalon import pyAvalonTools

RDLogger.DisableLog("rdApp.*")


def on_bits(raw):
    return [i * 8 + j for i, b in enumerate(raw) for j in range(8) if b >> j & 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("corpus")
    ap.add_argument("--show", type=int, default=10)
    ap.add_argument("--nbits", type=int, default=2048)
    args = ap.parse_args()

    match = mismatch = rd_none = ch_error = both_fail = 0
    shown = 0
    errors = []
    with open(args.corpus) as f:
        for lineno, line in enumerate(f, 1):
            parts = line.split()
            if not parts:
                continue
            smi = parts[0]
            rm = Chem.MolFromSmiles(smi)
            try:
                cm = chematic.from_smiles(smi)
                got = on_bits(bytes(cm.rdkit_avalon_fp(args.nbits)))
                cerr = None
            except Exception as e:  # noqa: BLE001
                got = None
                cerr = str(e)
            if rm is None:
                if got is None:
                    both_fail += 1
                else:
                    rd_none += 1
                continue
            exp = list(pyAvalonTools.GetAvalonFP(rm, nBits=args.nbits).GetOnBits())
            if got is None:
                ch_error += 1
                errors.append((lineno, smi, cerr))
                continue
            if got == exp:
                match += 1
            else:
                mismatch += 1
                if shown < args.show:
                    shown += 1
                    e, g = set(exp), set(got)
                    print(f"MISMATCH line {lineno}: {smi}")
                    print(f"   only RDKit:    {sorted(e - g)}")
                    print(f"   only chematic: {sorted(g - e)}")
    total = match + mismatch + ch_error
    print(
        f"{args.corpus}: match {match}/{total} mismatch {mismatch} "
        f"chematic-error {ch_error} rdkit-none {rd_none} both-fail {both_fail}"
    )
    for lineno, smi, err in errors[: args.show]:
        print(f"ERROR line {lineno}: {smi}: {err}")
    return 0 if mismatch == 0 and ch_error == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
