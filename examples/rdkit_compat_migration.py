"""
RDKit → chematic migration guide (10 common patterns).

Each block shows the original RDKit code (commented) and the chematic
equivalent using chematic.rdkit_compat for minimal diffs.
"""

from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# ── 1. Parse SMILES ─────────────────────────────────────────────────────────
# RDKit:  from rdkit import Chem;  mol = Chem.MolFromSmiles("CCO")
mol = Chem.MolFromSmiles(ASPIRIN)
assert mol is not None

# ── 2. Convert back to SMILES ────────────────────────────────────────────────
# RDKit:  Chem.MolToSmiles(mol)
smi = Chem.MolToSmiles(mol)
print(f"SMILES: {smi}")

# ── 3. Molecular weight / logP / TPSA ───────────────────────────────────────
# RDKit:  from rdkit.Chem import Descriptors
#         Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
mw   = Descriptors.MolWt(mol)
logp = Descriptors.MolLogP(mol)
tpsa = rdMolDescriptors.CalcTPSA(mol)
print(f"MW={mw:.2f}  LogP={logp:.2f}  TPSA={tpsa:.1f}")

# ── 4. H-bond donors / acceptors ─────────────────────────────────────────────
# RDKit:  rdMolDescriptors.CalcNumHBD(mol), rdMolDescriptors.CalcNumHBA(mol)
hbd = rdMolDescriptors.CalcNumHBD(mol)
hba = rdMolDescriptors.CalcNumHBA(mol)
print(f"HBD={hbd}  HBA={hba}")

# ── 5. Heavy atom count ──────────────────────────────────────────────────────
# RDKit:  mol.GetNumHeavyAtoms()
n_heavy = mol.GetNumHeavyAtoms()
print(f"Heavy atoms: {n_heavy}")

# ── 6. Substructure search (SMARTS) ─────────────────────────────────────────
# RDKit:  patt = Chem.MolFromSmarts("[OH]");  mol.HasSubstructMatch(patt)
patt = Chem.MolFromSmarts("[OH]")
has_oh = mol.HasSubstructMatch(patt)
print(f"Has [OH]: {has_oh}")

# ── 7. Morgan fingerprint + Tanimoto similarity ──────────────────────────────
# RDKit:  from rdkit.Chem import rdMolDescriptors, DataStructs
#         fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
#         DataStructs.TanimotoSimilarity(fp1, fp2)
fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2)
fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(
    Chem.MolFromSmiles("CC(=O)O"), 2  # acetic acid
)
sim = DataStructs.TanimotoSimilarity(fp1, fp2)
print(f"Tanimoto(aspirin, acetic acid): {sim:.3f}")

# ── 8. Parse MOL block (e.g., from SDF file) ─────────────────────────────────
# RDKit:  mol2 = Chem.MolFromMolBlock(block)
block = Chem.MolToMolBlock(mol)
mol2 = Chem.MolFromMolBlock(block)
assert mol2 is not None
print(f"MolFromMolBlock: {mol2}")

# ── 9. SDMolSupplier — iterate over SDF file ────────────────────────────────
# RDKit:  suppl = Chem.SDMolSupplier("mols.sdf")
#         for mol in suppl:
#             if mol is not None: process(mol)
#
# (Write a temporary SDF first for the demo)
import tempfile, os
with tempfile.NamedTemporaryFile(suffix=".sdf", mode="w", delete=False) as tmp:
    tmp.write(block + "\n$$$$\n")
    tmp_path = tmp.name

count = 0
with Chem.SDWriter(tmp_path.replace(".sdf", "_out.sdf")) as w:
    for m in Chem.SDMolSupplier(tmp_path):
        if m is not None:
            w.write(m)
            count += 1
os.unlink(tmp_path)
os.unlink(tmp_path.replace(".sdf", "_out.sdf"))
print(f"SDMolSupplier / SDWriter: processed {count} molecule(s)")

# ── 10. AddHs / RemoveHs ────────────────────────────────────────────────────
# RDKit:  mol_h = Chem.AddHs(mol);  mol_noh = Chem.RemoveHs(mol_h)
mol_h   = Chem.AddHs(mol)
mol_noh = Chem.RemoveHs(mol_h)
print(f"AddHs → RemoveHs: {mol_noh.GetNumHeavyAtoms()} heavy atoms (was {n_heavy})")

print("\nAll 10 patterns OK.")
