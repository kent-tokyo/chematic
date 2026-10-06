#!/usr/bin/env python3
"""Compare a ``chematic_chemistry_dump.py`` dump with RDKit 2026.03.6.

Counts, over the dumped corpus:

* accurate CIP: tetrahedral labels that agree with ``rdCIPLabeler``, typed
  abstentions, and disagreements (rows with ``@``);
* hybridization: atoms whose RDKit-model hybridization equals
  ``GetHybridization`` (SP/SP2/SP3, anything else as 0);
* MMFF94: heavy and hydrogen atoms typed differently from
  ``MMFFGetMoleculeProperties`` (rows RDKit types);
* MOL writer: stereo rows whose block RDKit reads back as the input
  (canonical isomeric SMILES), and ``strict=True`` refusals.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdCIPLabeler

RDLogger.DisableLog("rdApp.*")
HYB = {
    Chem.HybridizationType.SP: 1,
    Chem.HybridizationType.SP2: 2,
    Chem.HybridizationType.SP3: 3,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    lines = args.dump.read_text().splitlines()
    header = json.loads(lines[0])
    cip = Counter()
    hyb = Counter()
    mmff = Counter()
    mol = Counter()
    mol_lost: list[str] = []
    for line in lines[1:]:
        r = json.loads(line)
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None or "error" in r:
            continue
        # CIP
        if "@" in r["smiles"] and "cip" in r:
            rd = Chem.Mol(m)
            rdCIPLabeler.AssignCIPLabels(rd)
            want = {a.GetIdx(): a.GetProp("_CIPCode") for a in rd.GetAtoms() if a.HasProp("_CIPCode")}
            got = {d["atom_idx"]: d["descriptor"] for d in r["cip"] if d["descriptor"] in "RSrs"}
            abstained = {d["atom_idx"] for d in r["cip_unresolved"]}
            for i in set(want) | set(got):
                if i in want and want[i] == got.get(i):
                    cip["agree"] += 1
                elif i in want and i not in got:
                    cip["abstain" if i in abstained else "rdkit_only"] += 1
                elif i not in want:
                    cip["chematic_only"] += 1
                else:
                    cip["differ"] += 1
            cip["kekule_dependent"] += sum(1 for d in r["cip"] if d.get("kekule_dependent"))
        # hybridization
        for a in m.GetAtoms():
            hyb["atoms"] += 1
            hyb["agree"] += HYB.get(a.GetHybridization(), 0) == r["hybridization"][a.GetIdx()]
        # MMFF94
        mh = Chem.AddHs(m)
        props = AllChem.MMFFGetMoleculeProperties(mh, mmffVariant="MMFF94")
        if props is None:
            mmff["rdkit_unsupported"] += 1
        elif "mmff" not in r:
            mmff["chematic_error"] += 1
        elif len(r["mmff"]) != mh.GetNumAtoms():
            mmff["atom_count_mismatch"] += 1
        else:
            mmff["compared_rows"] += 1
            for a in mh.GetAtoms():
                kind = "hydrogen" if a.GetAtomicNum() == 1 else "heavy"
                mmff[kind] += 1
                if props.GetMMFFAtomType(a.GetIdx()) != r["mmff"][a.GetIdx()]:
                    mmff[kind + "_differing"] += 1
        # MOL writer
        want_smiles = Chem.MolToSmiles(m)
        if "molblock" in r and any(c in want_smiles for c in "@/\\"):
            mol["stereo_rows"] += 1
            back = Chem.MolFromMolBlock(r["molblock"])
            if back is not None and Chem.MolToSmiles(back) == want_smiles:
                mol["read_back_as_input"] += 1
            else:
                mol_lost.append(r["smiles"])
            mol["strict_refused"] += bool(r.get("mol_strict_refused"))
    report = {
        "schema": "chematic-chemistry-dump-vs-rdkit/v1",
        "rdkit": rdBase.rdkitVersion,
        "chematic": header,
        "cip_accurate": dict(cip),
        "hybridization": dict(hyb),
        "mmff94": dict(mmff),
        "mol_writer": dict(mol),
        "mol_writer_lost": mol_lost,
    }
    args.output.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: report[k] for k in ("cip_accurate", "hybridization", "mmff94", "mol_writer")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
