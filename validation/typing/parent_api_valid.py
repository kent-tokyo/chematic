import chematic as ch

mol = ch.from_smiles("CC(=O)[O-].[Na+]")
parents: list[ch.Mol] = [
    mol.fragment_parent(), mol.charge_parent(), mol.isotope_parent(), mol.stereo_parent(),
]
tautomer, tautomer_status = mol.tautomer_parent(max_transforms=16, timeout_ms=None)
super_mol, super_status = mol.super_parent(max_tautomers=32)
report: ch.SuperParentReport = mol.super_parent_report()
stage: ch.ParentStage = report["stages"][0]
status: str = report["status"]
smiles: str = stage["smiles"]
