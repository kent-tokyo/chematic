"""SDF → descriptor CSV using chematic.rdkit_compat (RDKit-style API).

RDKit equivalent:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    suppl = Chem.SDMolSupplier("in.sdf")
    ... Descriptors.MolWt(mol) ...
"""
import csv
import io
import tempfile
import os

from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors


def sdf_to_rows(sdf_path):
    rows = []
    for mol in Chem.SDMolSupplier(sdf_path):
        if mol is None:
            continue
        rows.append({
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
            "smiles": Chem.MolToSmiles(mol),
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 1),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
        })
    return rows


def main():
    # Build a small demo SDF.
    d = tempfile.mkdtemp()
    sdf = os.path.join(d, "demo.sdf")
    with Chem.SDWriter(sdf) as w:
        for smi, name in [("CCO", "ethanol"), ("c1ccccc1", "benzene"),
                          ("CC(=O)Oc1ccccc1C(=O)O", "aspirin")]:
            m = Chem.MolFromSmiles(smi)
            m.SetProp("_Name", name)
            w.write(m)

    rows = sdf_to_rows(sdf)

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    print(buf.getvalue())

    assert len(rows) == 3
    assert any(r["name"] == "aspirin" and 180.0 < r["mw"] < 180.3 for r in rows)
    print("OK: SDF → CSV (3 molecules)")


if __name__ == "__main__":
    main()
