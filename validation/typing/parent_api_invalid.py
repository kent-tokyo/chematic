import chematic as ch

mol = ch.from_smiles("CCO")
mol.tautomer_parent(max_transforms="16")
mol.super_parent(timeout_ms="10")
wrong_parent: str = mol.fragment_parent()
report = mol.super_parent_report()
wrong_status: int = report["status"]
wrong_stage: ch.ParentStage = {"name": "cleanup", "smiles": "CCO"}
report["molecule"]
