#!/usr/bin/env python3
"""
Generate RDKit reference values for chematic accuracy benchmarking.

Computes MW, LogP, TPSA, heavy_atom_count for 200 diverse molecules,
and pairwise ECFP4 Tanimoto for a 50-mol subset.

Output files:
  scripts/rdkit_ref_properties.tsv   — per-molecule MW/LogP/TPSA/HAC
  scripts/rdkit_ref_tanimoto.tsv     — pairwise Tanimoto (50×50 upper triangle)
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, AllChem
import sys, json, math

# 200 diverse SMILES covering:
#   - simple organics
#   - heterocycles (pyridine, imidazole, thiophene, indole, furan)
#   - stereocenters, E/Z
#   - multi-ring (naphthalene, quinoline, coumarin, caffeine-family)
#   - drugs (aspirin, ibuprofen, paracetamol, warfarin, morphine, ...)
#   - charged (amines, carboxylates)
#   - heavy atoms (Cl, Br, S, P, N, O)
#   - macrocycles, sugars
SMILES_200 = [
    # --- Simple organics ---
    ("methane",          "C"),
    ("ethane",           "CC"),
    ("propane",          "CCC"),
    ("butane",           "CCCC"),
    ("pentane",          "CCCCC"),
    ("hexane",           "CCCCCC"),
    ("cyclohexane",      "C1CCCCC1"),
    ("benzene",          "c1ccccc1"),
    ("toluene",          "Cc1ccccc1"),
    ("ethylbenzene",     "CCc1ccccc1"),
    ("styrene",          "C=Cc1ccccc1"),
    ("phenol",           "Oc1ccccc1"),
    ("aniline",          "Nc1ccccc1"),
    ("naphthalene",      "c1ccc2ccccc2c1"),
    ("anthracene",       "c1ccc2cc3ccccc3cc2c1"),
    ("acetic_acid",      "CC(=O)O"),
    ("ethanol",          "CCO"),
    ("methanol",         "CO"),
    ("acetone",          "CC(C)=O"),
    ("methylamine",      "CN"),
    # --- Heteroaromatics ---
    ("pyridine",         "c1ccncc1"),
    ("pyrimidine",       "c1ccnc(n1)"),
    ("pyrazine",         "c1cnccn1"),
    ("imidazole",        "c1cnc[nH]1"),
    ("thiophene",        "c1ccsc1"),
    ("furan",            "c1ccco1"),
    ("indole",           "c1ccc2[nH]ccc2c1"),
    ("benzimidazole",    "c1ccc2[nH]cnc2c1"),
    ("purine",           "c1ncc2[nH]cnc2n1"),
    ("quinoline",        "c1ccc2ncccc2c1"),
    ("isoquinoline",     "c1ccc2cnccc2c1"),
    ("pyrrole",          "c1cc[nH]c1"),
    ("oxazole",          "c1cnoc1"),
    ("thiazole",         "c1csc(n1)"),
    ("pyrazole",         "c1cc[nH]n1"),
    ("triazole_1H",      "c1cn[nH]n1"),
    ("tetrazole",        "c1nn[nH]n1"),
    # --- Named drugs ---
    ("aspirin",          "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen",        "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("paracetamol",      "CC(=O)Nc1ccc(O)cc1"),
    ("caffeine",         "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("morphine",         "OC1=CC=C2CC3N(CCC4=C3C2=C1O4)C"),
    ("codeine",          "COC1=CC=C2CC3N(CCC4=C3C2=C1O4)C"),
    ("warfarin",         "OC(=O)CCCC(=O)c1ccc(cc1)Cl"),
    ("lidocaine",        "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    ("diazepam",         "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13"),
    ("metformin",        "CN(C)C(=N)NC(=N)N"),
    ("atorvastatin",     "CC(C)c1c(C(=O)Nc2ccccc2F)c(CC[C@@H](O)C[C@@H](O)CC(=O)O)c(-c2ccc(F)cc2)n1-c1ccccc1"),
    ("methotrexate",     "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)NC(CCC(=O)O)C(=O)O"),
    ("sildenafil",       "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"),
    ("tamoxifen",        "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
    ("propranolol",      "CC(C)NCC(O)COc1cccc2ccccc12"),
    ("naproxen",         "CC(C(=O)O)c1ccc2cc(OC)ccc2c1"),
    ("indomethacin",     "CC1=C(CC(=O)O)c2cc(OC)ccc2N1C(=O)c1ccc(Cl)cc1"),
    ("chlorpromazine",   "CN(C)CCCN1c2ccccc2Sc3ccc(Cl)cc13"),
    ("clonazepam",       "O=C1CN=C(c2ccccc2Cl)c3cc([N+](=O)[O-])ccc13"),
    ("fluoxetine",       "CNCC(COc1ccc(cc1)C(F)(F)F)Oc1ccccc1"),
    # --- Amino acids / peptide-like ---
    ("glycine",          "NCC(=O)O"),
    ("alanine",          "N[C@@H](C)C(=O)O"),
    ("serine",           "N[C@@H](CO)C(=O)O"),
    ("threonine",        "N[C@@H]([C@@H](O)C)C(=O)O"),
    ("valine",           "N[C@@H](C(C)C)C(=O)O"),
    ("leucine",          "N[C@@H](CC(C)C)C(=O)O"),
    ("phenylalanine",    "N[C@@H](Cc1ccccc1)C(=O)O"),
    ("tryptophan",       "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    ("tyrosine",         "N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
    ("histidine",        "N[C@@H](Cc1c[nH]cn1)C(=O)O"),
    ("cysteine",         "N[C@@H](CS)C(=O)O"),
    ("methionine",       "N[C@@H](CCSC)C(=O)O"),
    ("lysine",           "N[C@@H](CCCCN)C(=O)O"),
    ("arginine",         "N[C@@H](CCCNC(=N)N)C(=O)O"),
    ("aspartate",        "N[C@@H](CC(=O)O)C(=O)O"),
    ("glutamate",        "N[C@@H](CCC(=O)O)C(=O)O"),
    # --- Carbohydrates / polyols ---
    ("glucose",          "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("fructose",         "OC[C@@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O"),
    ("ribose",           "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O"),
    ("sucrose",          "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H]1O"),
    ("sorbitol",         "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"),
    # --- Halogens ---
    ("chlorobenzene",    "Clc1ccccc1"),
    ("bromobenzene",     "Brc1ccccc1"),
    ("fluorobenzene",    "Fc1ccccc1"),
    ("iodobenzene",      "Ic1ccccc1"),
    ("dichloromethane",  "ClCCl"),
    ("chloroform",       "ClC(Cl)Cl"),
    ("ccl4",             "ClC(Cl)(Cl)Cl"),
    ("4_chlorophenol",   "Oc1ccc(Cl)cc1"),
    ("4_bromophenol",    "Oc1ccc(Br)cc1"),
    ("trifluoromethylbenzene", "FC(F)(F)c1ccccc1"),
    # --- Sulfur/phosphorus ---
    ("dimethyl_sulfide", "CSC"),
    ("dimethyl_sulfoxide", "CS(C)=O"),
    ("dimethyl_sulfone", "CS(C)(=O)=O"),
    ("methane_sulfonic_acid", "CS(=O)(=O)O"),
    ("thiophenol",       "Sc1ccccc1"),
    ("cystamine",        "NCCSSCCN"),
    ("trimethyl_phosphate", "COP(=O)(OC)OC"),
    # --- Carboxylic acids / esters / amides ---
    ("benzoic_acid",     "OC(=O)c1ccccc1"),
    ("phenylacetic_acid","OC(=O)Cc1ccccc1"),
    ("cinnamic_acid_E",  "OC(=O)/C=C/c1ccccc1"),
    ("methyl_benzoate",  "COC(=O)c1ccccc1"),
    ("ethyl_benzoate",   "CCOC(=O)c1ccccc1"),
    ("benzamide",        "NC(=O)c1ccccc1"),
    ("n_methyl_benzamide","CNC(=O)c1ccccc1"),
    ("urea",             "NC(N)=O"),
    ("dimethyl_urea",    "CN(C(=O)N)C"),
    ("succinic_acid",    "OC(=O)CCC(=O)O"),
    ("maleic_acid",      "OC(=O)/C=C\\C(=O)O"),
    ("fumaric_acid",     "OC(=O)/C=C/C(=O)O"),
    ("phthalic_acid",    "OC(=O)c1ccccc1C(=O)O"),
    ("isophthalic_acid", "OC(=O)c1cccc(c1)C(=O)O"),
    # --- Amines / nitro ---
    ("n_methylaniline",  "CNc1ccccc1"),
    ("diphenylamine",    "c1ccc(cc1)Nc1ccccc1"),
    ("triphenylamine",   "c1ccc(cc1)N(c1ccccc1)c1ccccc1"),
    ("nitrobenzene",     "O=c1ccc([N+](=O)[O-])cc1"),
    ("4_aminophenol",    "Nc1ccc(O)cc1"),
    ("4_nitrophenol",    "Oc1ccc([N+](=O)[O-])cc1"),
    # --- Bicyclic / fused ring systems ---
    ("tetralin",         "C1CCc2ccccc2C1"),
    ("decalin",          "C1CCC2CCCCC2C1"),
    ("coumarin",         "O=c1ccc2ccccc2o1"),
    ("chromone",         "O=c1ccoc2ccccc12"),
    ("xanthene",         "c1ccc2c(c1)Cc1ccccc1O2"),
    ("acridine",         "c1ccc2nc3ccccc3cc2c1"),
    ("carbazole",        "c1ccc2[nH]c3ccccc3c2c1"),
    ("dibenzofuran",     "c1ccc2c(c1)oc1ccccc12"),
    ("dibenzothiophene", "c1ccc2c(c1)sc1ccccc12"),
    ("fluorene",         "C1c2ccccc2-c2ccccc21"),
    ("fluoranthene",     "c1ccc2ccc3cccc4ccc1c2c34"),
    ("pyrene",           "c1cc2ccc3cccc4ccc(c1)c2c34"),
    ("triphenylene",     "c1ccc2ccc3ccccc3c2c1"),
    # --- Heterocyclic drug fragments ---
    ("piperazine",       "C1CNCCN1"),
    ("morpholine",       "C1CNCCO1"),
    ("n_methylpiperazine", "CN1CCNCC1"),
    ("n_boc_piperazine", "O=C(OC(C)(C)C)N1CCNCC1"),
    ("piperidine",       "C1CCNCC1"),
    ("pyrrolidine",      "C1CCNC1"),
    ("azetidine",        "C1CNC1"),
    ("oxetane",          "C1COC1"),
    ("thietane",         "C1CSC1"),
    ("tetrahydrofuran",  "C1CCOC1"),
    ("tetrahydrothiophene", "C1CCSC1"),
    ("tetrahydropyran",  "C1CCOCC1"),
    ("1_4_dioxane",      "C1COCCO1"),
    # --- Natural products / complex ---
    ("nicotine",         "CN1CCC[C@H]1c1cccnc1"),
    ("quinine",          "COc1ccc2nccc(C(O)C3CC4CCN3CC4=C)c2c1"),
    ("capsaicin",        "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O"),
    ("resveratrol",      "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    ("quercetin",        "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"),
    ("curcumin",         "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
    ("penicillin_G",     "O=C(O)[C@@H]1[C@H]2SC(C)(C)[C@@H]2N1C(=O)Cc1ccccc1"),
    ("ampicillin",       "N[C@@H](c1ccccc1)C(=O)N[C@H]2[C@@H]3SC(C)(C)[C@H]3N2C(=O)O"),  # simplified
    # --- Nucleobases ---
    ("adenine",          "Nc1ncnc2[nH]cnc12"),
    ("guanine",          "Nc1nc2[nH]cnc2c(=O)[nH]1"),
    ("cytosine",         "Nc1ccnc(=O)[nH]1"),
    ("thymine",          "Cc1c[nH]c(=O)nc1=O"),
    ("uracil",           "O=c1cc[nH]c(=O)[nH]1"),
    ("hypoxanthine",     "O=c1[nH]cnc2[nH]cnc12"),
    ("xanthine",         "O=c1[nH]c2[nH]cnc2c(=O)[nH]1"),
    # --- Misc / test edges ---
    ("ethylene_glycol",  "OCCO"),
    ("glycerol",         "OCC(O)CO"),
    ("propylene_glycol", "CC(O)CO"),
    ("catechol",         "Oc1ccccc1O"),
    ("resorcinol",       "Oc1cccc(O)c1"),
    ("hydroquinone",     "Oc1ccc(O)cc1"),
    ("salicylic_acid",   "OC(=O)c1ccccc1O"),
    ("vanillin",         "COc1cc(C=O)ccc1O"),
    ("eugenol",          "C=CCc1ccc(O)c(OC)c1"),
    ("safrol",           "C=CCc1ccc2c(c1)OCO2"),
    ("epinephrine",      "CNC[C@@H](O)c1ccc(O)c(O)c1"),
    ("dopamine",         "NCCc1ccc(O)c(O)c1"),
    ("serotonin",        "NCCc1c[nH]c2ccc(O)cc12"),
    ("melatonin",        "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"),
    ("histamine",        "NCCc1c[nH]cn1"),
    ("putrescine",       "NCCCCN"),
    ("spermidine",       "NCCCNCCCCN"),
    ("folic_acid",       "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)nc2c(=O)[nH]1"),
    ("riboflavin",       "Cc1cc2nc3c(=O)[nH]c(=O)nc3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C"),
    ("pyridoxine",       "Cc1ncc(CO)c(CO)c1O"),
    ("thiamine",         "Cc1ncc(C[n+]2csc(CCO)c2C)cn1"),
    ("cobalamin_ring_a", "CC1=C2[C@H](CC(=O)N)..."),  # placeholder — skip complex
]

# Filter to parseable and remove the placeholder
valid = []
for name, smi in SMILES_200:
    if "..." in smi:
        continue
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"[SKIP] {name}: {smi}", file=sys.stderr)
        continue
    valid.append((name, smi, mol))

print(f"# Valid molecules: {len(valid)}", file=sys.stderr)

# ---- Properties ----
rows = []
for name, smi, mol in valid:
    mw    = Descriptors.MolWt(mol)
    logp  = Descriptors.MolLogP(mol)
    tpsa  = rdMolDescriptors.CalcTPSA(mol, includeSandP=True)
    hac   = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    hbd   = rdMolDescriptors.CalcNumHBD(mol)
    hba   = rdMolDescriptors.CalcNumHBA(mol)
    rows.append((name, smi, mw, logp, tpsa, hac, hbd, hba))

with open("scripts/rdkit_ref_properties.tsv", "w") as f:
    f.write("name\tsmiles\tmw\tlogp\ttpsa\thac\thbd\thba\n")
    for name, smi, mw, logp, tpsa, hac, hbd, hba in rows:
        f.write(f"{name}\t{smi}\t{mw:.4f}\t{logp:.4f}\t{tpsa:.4f}\t{hac}\t{hbd}\t{hba}\n")

print(f"Written {len(rows)} rows to scripts/rdkit_ref_properties.tsv", file=sys.stderr)

# ---- ECFP4 Tanimoto (first 50) ----
fp_mols = valid[:50]
fps = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) for _, _, mol in fp_mols]

with open("scripts/rdkit_ref_tanimoto.tsv", "w") as f:
    names_header = "\t".join(n for n, _, _ in fp_mols)
    f.write(f"name\t{names_header}\n")
    for i, (name_i, _, _) in enumerate(fp_mols):
        sims = [DataStructs.TanimotoSimilarity(fps[i], fps[j]) for j in range(len(fp_mols))]
        f.write(name_i + "\t" + "\t".join(f"{s:.6f}" for s in sims) + "\n")

print(f"Written {len(fp_mols)}×{len(fp_mols)} Tanimoto matrix to scripts/rdkit_ref_tanimoto.tsv", file=sys.stderr)

# ---- Summary statistics preview ----
mws   = [r[2] for r in rows]
logps = [r[3] for r in rows]
tpsas = [r[4] for r in rows]
print(f"\nMW   range: {min(mws):.1f} – {max(mws):.1f}  mean={sum(mws)/len(mws):.1f}", file=sys.stderr)
print(f"LogP range: {min(logps):.2f} – {max(logps):.2f}  mean={sum(logps)/len(logps):.2f}", file=sys.stderr)
print(f"TPSA range: {min(tpsas):.1f} – {max(tpsas):.1f}  mean={sum(tpsas)/len(tpsas):.1f}", file=sys.stderr)
