#!/usr/bin/env python3
"""Broad-corpus differential validation for P1-S2 (E/Z direction perception
in the MOL V2000 reader).

For each SMILES in the standard corpus:
1. Parse with RDKit, generate 2D coordinates (`AllChem.Compute2DCoords`),
   write a V2000 MOL block (`Chem.MolToMolBlock`) -- this is RDKit's own,
   independently-generated 2D depiction, not anything chematic produced.
2. `rdkit_verdict` = the standard InChI `/b` layer of RDKit's OWN re-parse of
   that MOL block (`Chem.MolFromMolBlock` -> `Chem.MolToInchi`) -- "what E/Z
   does RDKit itself resolve from this 2D depiction".
3. Feed the SAME MOL block through chematic (`chematic.from_mol_block`,
   using the isolated venv's build -- see `--venv-python`).
4. `chematic_verdict` = the InChI `/b` layer of RDKit's re-parse of
   chematic's canonical SMILES output for that molecule.
5. Compare `rdkit_verdict` vs `chematic_verdict` (atom-numbering-independent
   by construction -- InChI's own canonical atom order, not either engine's
   raw index order):
     - both None                          -> no_stereo_both_agree
     - rdkit set, chematic == rdkit       -> semantic_agreement
     - rdkit set, chematic None           -> chematic_abstained
     - rdkit set, chematic != rdkit (set) -> SEMANTIC INVERSION (gate: must be 0)
     - rdkit None, chematic set           -> FALSE POSITIVE (gate: must be 0)

Required gates: semantic inversions == 0, false-positive assignments == 0.
Abstaining (chematic_abstained) is never a gate failure.

Usage:
    /path/to/venv/bin/python scripts/stereo2d_ez_corpus_diagnosis.py \
        ~/Downloads/SMILES.csv --limit 4999
"""

import argparse
import io
import json
import sys
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUMMARY_PATH = ROOT / "validation" / "results" / "stereo2d_ez_corpus_diagnosis_summary.json"
RDKIT_PINNED_COMMIT = "8afba32ec539dcb2369bc84549d802aca3f7eb39"
EXPECTED_RDKIT_VERSION = "2026.03.3"


def inchi_b_layer(rdkit_mol, Chem):
    if rdkit_mol is None:
        return None
    stderr_buf = io.StringIO()
    with redirect_stderr(stderr_buf):
        inchi = Chem.MolToInchi(rdkit_mol)
    if not inchi:
        return None
    for part in inchi.split("/"):
        if part.startswith("b"):
            return part
    return None


def has_cumulene(rdkit_mol, Chem):
    """Cheap structural proxy for allene/cumulene: an sp carbon flanked by
    two double bonds (`*=[#6]=*`)."""
    patt = Chem.MolFromSmarts("*=[#6]=*")
    return rdkit_mol.HasSubstructMatch(patt) if patt is not None else False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("corpus", help="CSV/TXT file with one SMILES per line (first column)")
    parser.add_argument("--limit", type=int, default=4999)
    args = parser.parse_args()

    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem

    installed_version = rdBase.rdkitVersion
    version_mismatch = installed_version != EXPECTED_RDKIT_VERSION
    if version_mismatch:
        print(f"WARNING: installed rdkit=={installed_version} != pinned {EXPECTED_RDKIT_VERSION}", file=sys.stderr)

    try:
        import chematic
    except ImportError:
        print("FAIL: chematic Python module not importable -- build the wheel into this venv first", file=sys.stderr)
        sys.exit(1)

    lines = [
        l.strip().split(",")[0]
        for l in Path(args.corpus).expanduser().read_text().splitlines()
        if l.strip()
    ]
    lines = lines[: args.limit]

    input_count = len(lines)
    parsed_count = 0
    parse_failures = []
    rdkit_resolved_molecules = 0
    resolved_bond_count_proxy = 0  # count of molecules with >=1 resolved b-layer entry (comma-separated)
    chematic_assigned = 0
    chematic_abstained = 0
    no_stereo_both_agree = 0
    semantic_agreements = 0
    semantic_inversions = []
    false_positives = []
    cumulene_cases = 0
    mol_block_write_failures = []
    chematic_parse_failures = []

    for i, smi in enumerate(lines):
        rdkit_mol = Chem.MolFromSmiles(smi)
        if rdkit_mol is None:
            parse_failures.append({"index": i, "smiles": smi})
            continue
        parsed_count += 1

        try:
            AllChem.Compute2DCoords(rdkit_mol)
            mol_block = Chem.MolToMolBlock(rdkit_mol)
        except Exception as e:  # noqa: BLE001 - corpus-scale defensive catch, reported not swallowed
            mol_block_write_failures.append({"index": i, "smiles": smi, "error": str(e)})
            continue

        stderr_buf = io.StringIO()
        with redirect_stderr(stderr_buf):
            rdkit_reparsed = Chem.MolFromMolBlock(mol_block)
        rdkit_verdict = inchi_b_layer(rdkit_reparsed, Chem)
        if rdkit_verdict is not None:
            rdkit_resolved_molecules += 1
            resolved_bond_count_proxy += rdkit_verdict.count(",") + 1

        try:
            chematic_mol = chematic.from_mol_block(mol_block)
            chematic_smiles = chematic_mol.smiles
        except Exception as e:  # noqa: BLE001
            chematic_parse_failures.append({"index": i, "smiles": smi, "error": str(e)})
            continue

        stderr_buf2 = io.StringIO()
        with redirect_stderr(stderr_buf2):
            chematic_reparsed = Chem.MolFromSmiles(chematic_smiles)
        chematic_verdict = inchi_b_layer(chematic_reparsed, Chem)

        if chematic_verdict is not None:
            chematic_assigned += 1
        else:
            chematic_abstained += 1

        if rdkit_verdict is None and chematic_verdict is None:
            no_stereo_both_agree += 1
        elif rdkit_verdict is not None and chematic_verdict == rdkit_verdict:
            semantic_agreements += 1
        elif rdkit_verdict is not None and chematic_verdict is None:
            pass  # chematic_abstained already counted above
        elif rdkit_verdict is not None and chematic_verdict is not None:
            semantic_inversions.append(
                {
                    "index": i,
                    "smiles": smi,
                    "mol_block": mol_block,
                    "rdkit_verdict": rdkit_verdict,
                    "chematic_verdict": chematic_verdict,
                    "chematic_smiles": chematic_smiles,
                }
            )
        elif rdkit_verdict is None and chematic_verdict is not None:
            false_positives.append(
                {
                    "index": i,
                    "smiles": smi,
                    "mol_block": mol_block,
                    "chematic_verdict": chematic_verdict,
                    "chematic_smiles": chematic_smiles,
                }
            )
            if has_cumulene(rdkit_mol, Chem):
                cumulene_cases += 1

    summary = {
        "rdkit_version": installed_version,
        "rdkit_version_mismatch": version_mismatch,
        "rdkit_pinned_commit": RDKIT_PINNED_COMMIT,
        "input_molecules": input_count,
        "parsed_molecules": parsed_count,
        "parse_failures_count": len(parse_failures),
        "parse_failures": parse_failures[:50],
        "mol_block_write_failures_count": len(mol_block_write_failures),
        "chematic_parse_failures_count": len(chematic_parse_failures),
        "chematic_parse_failures": chematic_parse_failures[:50],
        "molecules_with_rdkit_resolved_ez": rdkit_resolved_molecules,
        "resolved_ez_bonds_proxy": resolved_bond_count_proxy,
        "chematic_assigned": chematic_assigned,
        "chematic_abstained": chematic_abstained,
        "no_stereo_both_agree": no_stereo_both_agree,
        "semantic_agreements": semantic_agreements,
        "semantic_inversions_count": len(semantic_inversions),
        "semantic_inversions": semantic_inversions[:50],
        "false_positive_count": len(false_positives),
        "false_positives": false_positives[:50],
        "false_positive_cumulene_subset_count": cumulene_cases,
    }
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2))

    print(f"input={input_count} parsed={parsed_count} parse_failures={len(parse_failures)}")
    print(f"mol_block_write_failures={len(mol_block_write_failures)} chematic_parse_failures={len(chematic_parse_failures)}")
    print(f"rdkit_resolved_molecules={rdkit_resolved_molecules} resolved_bonds_proxy={resolved_bond_count_proxy}")
    print(f"chematic_assigned={chematic_assigned} chematic_abstained={chematic_abstained}")
    print(f"no_stereo_both_agree={no_stereo_both_agree} semantic_agreements={semantic_agreements}")
    print(f"semantic_inversions={len(semantic_inversions)} false_positives={len(false_positives)} (cumulene subset: {cumulene_cases})")

    if semantic_inversions or false_positives:
        print("FAIL: required gates violated (semantic_inversions==0, false_positives==0)", file=sys.stderr)
        sys.exit(1)
    print("OK: both required gates satisfied (0 semantic inversions, 0 false positives).")
    sys.exit(0)


if __name__ == "__main__":
    main()
