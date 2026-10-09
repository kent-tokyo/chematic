#!/usr/bin/env python3
"""Per-engine operation tables for the corpus-scale chematic / COSMolKit /
RDKit comparison.

Every engine exposes ``parse(smiles)`` and a dict of operations that take the
parsed molecule and return a JSON-serialisable value in a shared normal form.
An operation an engine does not offer is simply absent and is reported as
``unsupported`` -- never as a match. Where chematic has both a native and an
RDKit-compatible API, the RDKit-compatible one is used and named in
``API_NOTES``; nothing is substituted silently.

Normal forms (all indices are input atom/bond order):

* fingerprints: sorted on-bit indices in RDKit's numbering (MACCS keys 1..166)
* CIP: sorted ``["a", atom, label]`` and ``["b", atom, atom, label]`` entries
  (bonds keyed by their two atoms, so engines with different bond numbering
  compare equal);
  chematic abstentions are returned separately by ``cip_abstain``
* SMARTS: per query, sorted list of sorted atom-index tuples (unique matches)
* ``*_rt`` operations: the RDKit isomeric canonical SMILES that RDKit reads back
  from the engine's output (RDKit's own writer included). They are scored
  against the reference's ``input_canonical`` -- the input molecule itself --
  so a match means "the written file still describes the input molecule".
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SMARTS_QUERIES = json.loads(
    (ROOT / "validation/rdkit_rebaseline_smarts_queries.json").read_text(encoding="utf-8")
)["queries"]

# rdMolHash.HashFunction members, in enum order; ``mol_hash_<name>`` ops.
MOL_HASH_FUNCTIONS = [
    "AnonymousGraph", "ElementGraph", "CanonicalSmiles", "MurckoScaffold", "ExtendedMurcko",
    "MolFormula", "AtomBondCounts", "DegreeVector", "Mesomer", "HetAtomTautomer",
    "HetAtomProtomer", "RedoxPair", "Regioisomer", "NetCharge", "SmallWorldIndexBR",
    "SmallWorldIndexBRL", "ArthorSubstructureOrder", "HetAtomTautomerv2", "HetAtomProtomerv2",
]

API_NOTES = {
    "chematic": {
        "canonical_smiles": "Mol.rdkit_smiles (RDKit-compatible writer); canonical_smiles_rt reads back the native Mol.smiles",
        "smarts": "Mol.find_matches_rdkit_parity (opt-in RDKit 2026.03.6 profile; typed refusals)",
        "cip": "Mol.cip_stereo(mode='accurate'); abstentions from cip_stereo_unresolved()",
        "morgan2_2048": "Mol.rdkit_ecfp_config(2, 2048)",
        "mol_wt": "Mol.rdkit_mw", "tpsa": "Mol.rdkit_tpsa",
        "aromatic_rings": "Mol.rdkit_aromatic_ring_count",
        "molblock_rt": "Mol.to_mol_block() (default writer)",
        "murcko_scaffold": "Mol.rdkit_murcko_scaffold (RDKit MurckoDecompose on the RDKit model)",
        "smarts_write": "Mol.rdkit_smarts", "cx_smarts": "Mol.rdkit_cx_smarts",
        "pdb_block": "Mol.rdkit_pdb_block", "stereoisomer_count": "Mol.rdkit_stereoisomer_count",
        "chiral_centers": "Mol.rdkit_chiral_centers(include_unassigned=True)",
        "inchi": "Mol.standard_inchi (needs the native-inchi build feature)",
        "morgan2_bitinfo": "Mol.rdkit_morgan_bit_info(2, 2048)",
        "pdb_read": "chematic.rdkit_pdb_block_to_smiles(Mol.rdkit_pdb_block())",
        "extended_murcko": "Mol.rdkit_mol_hash('ExtendedMurcko')",
        "embed3d": "Mol.add_hydrogens().rdkit_embed(random_seed=42)",
        "mol_hash_*": "Mol.rdkit_mol_hash(<HashFunction name>)",
        "tautomer_*": "Mol.canonical_tautomer() / enumerate_tautomers() (native), read back by RDKit",
        "rxn:*": "chematic.run_smirks_checked(smirks, [mol, *partners], rdkit_compat=True)",
        "distance_matrix": "Mol.rdkit_distance_matrix()",
        "distance_matrix_3d": "Mol.add_hydrogens().rdkit_distance_matrix_3d(RDKit-embedded coords)",
        "cx_smiles": "Mol.rdkit_cx_smiles()",
        "random_smiles5": "Mol.rdkit_random_smiles(5, 42)",
        "fragment_smiles": "Mol.rdkit_fragment_smiles(first half of Mol.rdkit_num_atoms())",
        "chemistry_problems": "Mol.rdkit_chemistry_problems()",
        "chemistry_problems_unsanitized": "chematic.rdkit_detect_chemistry_problems(smiles)",
        "reaction_smarts_roundtrip": "chematic.rdkit_reaction_to_smarts(smirks)",
        "morgan2_sparse_counts": "Mol.rdkit_morgan_sparse_counts(2)",
        "torsion_legacy_counts": "Mol.rdkit_legacy_torsion_counts()",
        "mmff_energy_gradient": "Mol.add_hydrogens()._rdkit_mmff_terms(RDKit-embedded coords)",
        "uff_energy_gradient": "Mol.add_hydrogens().rdkit_uff_energy/rdkit_uff_gradient(RDKit-embedded coords)",
    },
    "cosmolkit": {
        "smarts": "get_substruct_matches(mol, parse_smarts(q))",
        "cip": "Molecule.with_cip_labels()",
        "molblock_rt": "Molecule.to_2d_sdf_string()",
        "morgan2_bitinfo": "fingerprint_morgan_with_output(2, 2048).additional_output().bit_info_map()",
        "pdb_read": "Molecule.from_pdb_block(Molecule.to_pdb_block()).to_smiles()",
        "extended_murcko": "Molecule.net_scaffold().to_smiles()",
        "embed3d": "Molecule.with_hydrogens().with_3d_conformer_result(EmbedParameters.etkdg_v3(), random_seed=42)",
    },
    "cosmolkit-0.5": {
        "smarts": "Molecule.substruct_matches(parse_smarts(q))",
        "cip": "Molecule.with_cip_labels()",
        "molblock_rt": "Molecule.to_sdf_2d()",
        "morgan*": "morgan_fingerprint_with_generator(MorganFingerprintGenerator(params=MorganParams(radius, fp_size=2048, ...)))",
        "morgan2_bitinfo": "morgan_fingerprint_with_generator(..., output=FingerprintAdditionalOutput().bit_info_map())",
        "fp_torsion": "Molecule.legacy_topological_torsion_fingerprint()",
        "smiles_*": "Molecule.to_smiles_with_params(SmilesWriteParams(...))",
        "embed3d": "with_hydrogens().with_3d_conformer_result_with_params(EmbedParams.etkdg_v3() randomSeed=42)",
        "removed": "see COSMOLKIT_05_REMOVED",
        "tautomer_*": "Molecule.canonical_tautomer() / enumerate_tautomers().canonical_smiles()",
        "rxn:*": "parse_smirks(s).run([mol, *partners], ReactionRunParams())",
        "3d ops": "with_hydrogens().with_only_3d_conformer(RDKit-embedded coords)",
        "chemistry_problems": "detect_chemistry_problems(), kind + message mapped to RDKit types",
    },
}


def _rd():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    return Chem


def rdkit_readback(text: str | None, kind: str) -> str | None:
    """RDKit isomeric canonical SMILES of an engine's SMILES or MOL output."""
    Chem = _rd()
    if text is None:
        return None
    mol = Chem.MolFromSmiles(text) if kind == "smiles" else Chem.MolFromMolBlock(text)
    return None if mol is None else Chem.MolToSmiles(mol)


def rdkit_inchi_of(smiles_canonical: str) -> str | None:
    """RDKit's InChI of a molecule, the common input of the InChI-read check."""
    Chem = _rd()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol = Chem.MolFromSmiles(smiles_canonical)
        return Chem.MolToInchi(mol) or None if mol is not None else None


def _bits(raw: bytes, offset: int = 0) -> list[int]:
    return [i + offset for i in range(len(raw) * 8) if raw[i // 8] >> (i % 8) & 1]


def _bit_info(info) -> dict[str, list[list[int]]]:
    """A Morgan bitInfo map as JSON: bit -> sorted [atom, radius] pairs."""
    return {str(bit): sorted([int(a), int(r)] for a, r in env) for bit, env in sorted(info.items())}


def _smarts_sets(matches) -> list[list[int]]:
    return sorted({tuple(sorted(m)) for m in matches})


# --------------------------------------------------------------------- RDKit
def rdkit_engine():
    Chem = _rd()
    from rdkit import rdBase
    from rdkit.Chem import (QED, Crippen, Descriptors, MACCSkeys, rdCIPLabeler,
                            rdFingerprintGenerator, rdMolHash)
    from rdkit.Chem import rdMolDescriptors as D

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    queries = {q: Chem.MolFromSmarts(q) for q in SMARTS_QUERIES}

    def cip(mol):
        m = Chem.Mol(mol)
        rdCIPLabeler.AssignCIPLabels(m)
        out = [["a", a.GetIdx(), a.GetProp("_CIPCode")] for a in m.GetAtoms() if a.HasProp("_CIPCode")]
        out += [["b", *sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), b.GetProp("_CIPCode")]
                for b in m.GetBonds() if b.HasProp("_CIPCode")]
        return sorted(out)

    def inchi(mol):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return Chem.MolToInchi(mol) or None

    ops = {
        "canonical_smiles": lambda m: Chem.MolToSmiles(m),
        "input_canonical": lambda m: Chem.MolToSmiles(m),
        "canonical_smiles_rt": lambda m: rdkit_readback(Chem.MolToSmiles(m), "smiles"),
        "formula": D.CalcMolFormula,
        "mol_wt": Descriptors.MolWt,
        "exact_mw": Descriptors.ExactMolWt,
        "tpsa": D.CalcTPSA,
        "logp": Crippen.MolLogP,
        "mr": Crippen.MolMR,
        "hba": D.CalcNumHBA,
        "hbd": D.CalcNumHBD,
        "rotatable_bonds": Descriptors.NumRotatableBonds,
        "aromatic_rings": D.CalcNumAromaticRings,
        "fsp3": D.CalcFractionCSP3,
        "qed": QED.qed,
        "stereocenters": D.CalcNumAtomStereoCenters,
        "morgan2_2048": lambda m: list(gen.GetFingerprint(m).GetOnBits()),
        "maccs": lambda m: list(MACCSkeys.GenMACCSKeys(m).GetOnBits()),
        "cip": cip,
        "molblock_rt": lambda m: rdkit_readback(Chem.MolToMolBlock(m), "mol"),
        "inchi": inchi,
        "inchikey": lambda m: Chem.InchiToInchiKey(inchi(m)) if inchi(m) else None,
    }
    for q, qm in queries.items():
        ops["smarts:" + q] = (lambda qm: lambda m: _smarts_sets(m.GetSubstructMatches(qm, maxMatches=100000)))(qm)
    from rdkit.Avalon import pyAvalonTools
    from rdkit.Chem import GraphDescriptors, rdMolDescriptors as RD
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, GetStereoisomerCount
    from rdkit.Chem.Scaffolds import MurckoScaffold
    gen3 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    gen_chiral = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    gen_count = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)

    def inchi_rt(m):
        inchi = rdkit_inchi_of(Chem.MolToSmiles(m))
        back = Chem.MolFromInchi(inchi) if inchi else None
        return Chem.MolToSmiles(back) if back is not None else None

    ops.update({
        "smiles_kekule": lambda m: Chem.MolToSmiles(m, kekuleSmiles=True),
        "smiles_noniso": lambda m: Chem.MolToSmiles(m, isomericSmiles=False),
        "smiles_explicit": lambda m: Chem.MolToSmiles(m, allBondsExplicit=True, allHsExplicit=True),
        "smiles_addhs": lambda m: Chem.MolToSmiles(Chem.AddHs(m)),
        "inchi_read": inchi_rt,
        "murcko_scaffold": lambda m: Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)),
        "smarts_write": lambda m: Chem.MolToSmarts(m),
        "cx_smarts": lambda m: Chem.MolToCXSmarts(m),
        "pdb_block": lambda m: Chem.MolToPDBBlock(m),
        "stereoisomer_count": lambda m: GetStereoisomerCount(m),
        "chiral_centers": lambda m: [list(c) for c in
                                     Chem.FindMolChiralCenters(Chem.Mol(m), includeUnassigned=True)],
        "morgan3_2048": lambda m: list(gen3.GetFingerprint(m).GetOnBits()),
        "morgan2_chiral": lambda m: list(gen_chiral.GetFingerprint(m).GetOnBits()),
        "morgan2_countsim": lambda m: list(gen_count.GetFingerprint(m).GetOnBits()),
        "atom_pair_counts": lambda m: sorted([k, v] for k, v in
                                             RD.GetHashedAtomPairFingerprint(m, nBits=2048).GetNonzeroElements().items()),
        "stereoisomers": lambda m: sorted(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)),
        "morgan2_bitinfo": lambda m: morgan_bitinfo(m),
        "pdb_read": lambda m: (lambda r: Chem.MolToSmiles(r) if r is not None else None)(
            Chem.MolFromPDBBlock(Chem.MolToPDBBlock(m))),
        "extended_murcko": lambda m: rdMolHash.MolHash(m, rdMolHash.HashFunction.ExtendedMurcko),
        "embed3d": lambda m: embed3d(m),
    })
    for name in MOL_HASH_FUNCTIONS:
        ops["mol_hash_" + name] = (lambda f: lambda m: rdMolHash.MolHash(m, f))(
            rdMolHash.HashFunction.names[name])
    ops.update({
        "fp_atom_pair": lambda m: list(RD.GetHashedAtomPairFingerprintAsBitVect(m, nBits=2048).GetOnBits()),
        "fp_torsion": lambda m: list(RD.GetHashedTopologicalTorsionFingerprintAsBitVect(m, nBits=2048).GetOnBits()),
        "fp_pattern": lambda m: list(Chem.PatternFingerprint(m, fpSize=2048).GetOnBits()),
        "fp_layered": lambda m: list(Chem.LayeredFingerprint(m, fpSize=2048).GetOnBits()),
        "fp_rdkit": lambda m: list(Chem.RDKFingerprint(m, fpSize=2048).GetOnBits()),
        "fp_avalon": lambda m: list(pyAvalonTools.GetAvalonFP(m, nBits=2048).GetOnBits()),
        "chi0v": RD.CalcChi0v, "chi1v": RD.CalcChi1v, "chi2v": RD.CalcChi2v,
        "chi3v": RD.CalcChi3v, "chi4v": RD.CalcChi4v,
        "chi0": GraphDescriptors.Chi0, "chi1": GraphDescriptors.Chi1,
        "kappa1": RD.CalcKappa1, "kappa2": RD.CalcKappa2, "kappa3": RD.CalcKappa3,
        "hall_kier_alpha": RD.CalcHallKierAlpha,
        "labute_asa": RD.CalcLabuteASA,
        "slogp_vsa": lambda m: list(RD.SlogP_VSA_(m)),
        "smr_vsa": lambda m: list(RD.SMR_VSA_(m)),
        "peoe_vsa": lambda m: list(RD.PEOE_VSA_(m)),
        "mqn": lambda m: list(RD.MQNs_(m)),
        "num_rings": RD.CalcNumRings,
        "num_aliphatic_rings": RD.CalcNumAliphaticRings,
        "num_saturated_rings": RD.CalcNumSaturatedRings,
        "num_aromatic_heterocycles": RD.CalcNumAromaticHeterocycles,
        "num_aliphatic_heterocycles": RD.CalcNumAliphaticHeterocycles,
        "num_saturated_heterocycles": RD.CalcNumSaturatedHeterocycles,
        "num_heteroatoms": RD.CalcNumHeteroatoms,
        "num_amide_bonds": RD.CalcNumAmideBonds,
        "num_bridgehead_atoms": RD.CalcNumBridgeheadAtoms,
        "num_spiro_atoms": RD.CalcNumSpiroAtoms,
        "heavy_atoms": RD.CalcNumHeavyAtoms,
        "unspecified_stereocenters": RD.CalcNumUnspecifiedAtomStereoCenters,
        "bertz_ct": GraphDescriptors.BertzCT,
        "balaban_j": GraphDescriptors.BalabanJ,
        "ipc": GraphDescriptors.Ipc,
    })

    def embed3d(m):
        from rdkit.Chem import AllChem
        mh = Chem.AddHs(m)
        if AllChem.EmbedMolecule(mh, randomSeed=42) != 0:
            return None
        return mh.GetConformer().GetPositions().tolist()

    def morgan_bitinfo(m):
        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        gen.GetFingerprint(m, additionalOutput=ao)
        return _bit_info(ao.GetBitInfoMap())

    def parse(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit could not parse")
        return mol

    ops.update(_rdkit_new_ops(Chem))
    return {"version": rdBase.rdkitVersion, "parse": parse, "ops": ops}


# ------------------------------------------------------------------ chematic
def chematic_engine():
    import chematic as c

    def cip(m):
        out = []
        for d in m.cip_stereo(mode="accurate"):
            out.append(["b", *sorted(d["bond_atoms"]), d["descriptor"]] if "bond_idx" in d
                       else ["a", d["atom_idx"], d["descriptor"]])
        return sorted(out)

    def smarts(q):
        def run(m):
            r = m.find_matches_rdkit_parity(q)
            if r["status"] != "ok":
                raise RefusedError(r.get("reason") or r["status"])
            return _smarts_sets(r["matches"])
        return run

    def molblock(m):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return m.to_mol_block()

    def _native(get):
        try:
            return get()
        except RuntimeError as exc:
            if "native-inchi" in str(exc):  # wheel built without the C library
                raise UnsupportedError(str(exc)) from exc
            raise

    def standard_inchi(m):
        return _native(lambda: m.standard_inchi)

    ops = {
        "canonical_smiles": lambda m: m.rdkit_smiles,
        "canonical_smiles_rt": lambda m: rdkit_readback(m.smiles, "smiles"),
        "formula": lambda m: m.formula,
        "mol_wt": lambda m: m.rdkit_mw,
        "exact_mw": lambda m: m.exact_mass,
        "tpsa": lambda m: m.rdkit_tpsa,
        "logp": lambda m: m.logp,
        "mr": lambda m: m.molar_refractivity,
        "hba": lambda m: m.rdkit_hba,
        "hbd": lambda m: m.hbd,
        "rotatable_bonds": lambda m: m.rotatable_bonds,
        "aromatic_rings": lambda m: m.rdkit_aromatic_ring_count,
        "fsp3": lambda m: m.fsp3,
        "qed": lambda m: m.qed,
        "stereocenters": lambda m: m.num_stereocenters,
        "morgan2_2048": lambda m: _bits(m.rdkit_ecfp_config(2, 2048)),
        "maccs": lambda m: _bits(m.maccs(), 1),
        "cip": cip,
        "cip_abstain": lambda m: sorted(m.cip_stereo_unresolved()),
        "molblock_rt": lambda m: rdkit_readback(molblock(m), "mol"),
        "inchi": standard_inchi,
        "inchikey": lambda m: _native(lambda: m.standard_inchikey),
    }
    for q in SMARTS_QUERIES:
        ops["smarts:" + q] = smarts(q)
    ops.update({
        "fp_atom_pair": lambda m: _bits(m.rdkit_atom_pair_fp()),
        "fp_torsion": lambda m: _bits(m.rdkit_torsion_fp()),
        "fp_pattern": lambda m: _bits(m.rdkit_pattern_fp()),
        "fp_layered": lambda m: _bits(m.rdkit_layered_fp()),
        "fp_rdkit": lambda m: _bits(m.rdkit_rdk_fp()),
        "fp_avalon": lambda m: _bits(m.rdkit_avalon_fp(2048)),
        "chi0v": lambda m: m.chi0v, "chi1v": lambda m: m.chi1v, "chi2v": lambda m: m.chi2v,
        "chi3v": lambda m: m.chi3v, "chi4v": lambda m: m.chi4v,
        "chi0": lambda m: m.chi0, "chi1": lambda m: m.chi1,
        "kappa1": lambda m: m.kappa1, "kappa2": lambda m: m.kappa2, "kappa3": lambda m: m.kappa3,
        "hall_kier_alpha": lambda m: m.hall_kier_alpha,
        "labute_asa": lambda m: m.descriptors()["labute_asa"],
        "slogp_vsa": lambda m: list(m.slogp_vsa()),
        "smr_vsa": lambda m: list(m.smr_vsa()),
        "peoe_vsa": lambda m: list(m.peoe_vsa()),
        "mqn": lambda m: list(m.mqn()),
        "num_rings": lambda m: m.num_rings,
        "num_aliphatic_rings": lambda m: m.num_aliphatic_rings,
        "num_saturated_rings": lambda m: m.num_saturated_rings,
        "num_aromatic_heterocycles": lambda m: m.num_aromatic_heterocycles,
        "num_aliphatic_heterocycles": lambda m: m.num_aliphatic_heterocycles,
        "num_saturated_heterocycles": lambda m: m.num_saturated_heterocycles,
        "num_heteroatoms": lambda m: m.num_heteroatoms,
        "num_amide_bonds": lambda m: m.num_amide_bonds,
        "num_bridgehead_atoms": lambda m: m.num_bridgehead_atoms,
        "num_spiro_atoms": lambda m: m.num_spiro_atoms,
        "heavy_atoms": lambda m: m.heavy_atoms,
        "unspecified_stereocenters": lambda m: m.num_unspecified_stereocenters,
        "bertz_ct": lambda m: m.bertz_ct,
        "balaban_j": lambda m: m.balaban_j,
        "ipc": lambda m: m.ipc,
        "smiles_addhs": lambda m: m.add_hydrogens().rdkit_smiles,
        "inchi_read": lambda m: _inchi_read(m),
        "murcko_scaffold": lambda m: m.rdkit_murcko_scaffold(),
        "smarts_write": lambda m: m.rdkit_smarts(),
        "cx_smarts": lambda m: m.rdkit_cx_smarts(),
        "pdb_block": lambda m: m.rdkit_pdb_block(),
        "stereoisomer_count": lambda m: m.rdkit_stereoisomer_count(),
        "chiral_centers": lambda m: [list(c) for c in m.rdkit_chiral_centers(include_unassigned=True)],
        "morgan3_2048": lambda m: _bits(m.rdkit_ecfp_config(3, 2048)),
        "morgan2_chiral": lambda m: _bits(m.rdkit_ecfp_config(2, 2048, include_chirality=True)),
        "morgan2_countsim": lambda m: _bits(m.rdkit_ecfp_config(2, 2048, count_simulation=True)),
        "atom_pair_counts": lambda m: sorted([k, v] for k, v in m.rdkit_atom_pair_counts(2048)),
        "stereoisomers": lambda m: m.rdkit_stereoisomer_smiles(),
        "morgan2_bitinfo": lambda m: _bit_info(m.rdkit_morgan_bit_info(2, 2048)),
        "pdb_read": lambda m: c.rdkit_pdb_block_to_smiles(m.rdkit_pdb_block()),
        "extended_murcko": lambda m: m.rdkit_mol_hash("ExtendedMurcko"),
        "embed3d": lambda m: _embed3d(m),
    })
    for name in MOL_HASH_FUNCTIONS:
        ops["mol_hash_" + name] = (lambda f: lambda m: m.rdkit_mol_hash(f))(name)
    for name, kwargs in (("smiles_kekule", {"kekule": True}), ("smiles_noniso", {"isomeric": False}),
                         ("smiles_explicit", {"all_bonds_explicit": True, "all_hs_explicit": True})):
        ops[name] = (lambda kw: lambda m: _rdkit_smiles_with(m, kw))(kwargs)

    def _embed3d(m):
        try:
            return m.add_hydrogens().rdkit_embed(random_seed=42)
        except RuntimeError:  # RDKit's EmbedMolecule returns -1 here too
            return None

    def _inchi_read(m):
        inchi = rdkit_inchi_of(rdkit_readback(m.rdkit_smiles, "smiles"))
        return c.from_inchi(inchi).rdkit_smiles if inchi else None
    def _rdkit_smiles_with(m, kwargs):
        writer = getattr(m, "rdkit_smiles_with", None)
        if writer is None:
            raise UnsupportedError("rdkit_smiles options are not exposed")
        return writer(**kwargs)

    raw = {"canonical_smiles": lambda m: m.smiles, "molblock": molblock}
    ops.update(_chematic_new_ops(c))
    return {"version": c.__version__, "parse": _with_input_smiles(c.from_smiles), "ops": ops, "raw": raw}


# ----------------------------------------------------------------- COSMolKit
def _cosmolkit_03_engine():
    def _embed3d(m):
        p = ck.EmbedParameters.etkdg_v3()
        p.random_seed = 42
        try:
            r = m.with_hydrogens().with_3d_conformer_result(p)
        except Exception:  # noqa: BLE001 - embedding failure
            return None
        if not r.ok:
            return None
        mol = r.molecule() if callable(getattr(r, "molecule", None)) else getattr(r, "molecule", None)
        return None if mol is None or mol.num_conformers() == 0 else mol.coordinates_3d().tolist()

    import cosmolkit as ck

    def cip(m):
        lab = m.with_cip_labels()
        out = [["a", a.idx(), a.cip_descriptor()] for a in lab.atoms() if a.cip_descriptor()]
        out += [["b", *sorted((b.begin_atom_idx(), b.end_atom_idx())), b.cip_descriptor()]
                for b in lab.bonds() if b.cip_descriptor()]
        return sorted(out)

    def smarts(q):
        qm = ck.parse_smarts(q)
        return lambda m: _smarts_sets(r.atom_mapping() for r in
                                      ck.get_substruct_matches(m, qm, max_matches=100000))

    ops = {
        "canonical_smiles": lambda m: m.to_smiles(),
        "canonical_smiles_rt": lambda m: rdkit_readback(m.to_smiles(), "smiles"),
        "formula": ck.calc_mol_formula,
        "mol_wt": ck.calc_mol_wt,
        "exact_mw": ck.calc_exact_mol_wt,
        "tpsa": ck.calc_tpsa,
        "logp": lambda m: ck.calc_crippen_descriptors(m)[0],
        "mr": lambda m: ck.calc_crippen_descriptors(m)[1],
        "hba": ck.calc_num_hba,
        "hbd": ck.calc_num_hbd,
        "rotatable_bonds": ck.calc_num_rotatable_bonds,
        "aromatic_rings": ck.calc_num_aromatic_rings,
        "fsp3": ck.calc_fraction_csp3,
        "qed": ck.calc_qed,
        "stereocenters": ck.calc_num_atom_stereo_centers,
        "morgan2_2048": lambda m: sorted(m.fingerprint_morgan(radius=2, n_bits=2048).on_bits()),
        "maccs": lambda m: sorted(i + 1 for i in m.maccs_fingerprint().on_bits()),
        "cip": cip,
        "molblock_rt": lambda m: rdkit_readback(m.to_2d_sdf_string(), "mol"),
        "inchi": lambda m: m.to_inchi(),
        "inchikey": lambda m: m.to_inchi_key(),
    }
    for q in SMARTS_QUERIES:
        ops["smarts:" + q] = smarts(q)
    on = lambda fp: sorted(fp.on_bits())  # noqa: E731
    ops.update({
        "fp_atom_pair": lambda m: on(m.fingerprint_atom_pair(n_bits=2048)),
        "fp_torsion": lambda m: on(ck.get_hashed_topological_torsion_fingerprint_as_bit_vect(m, n_bits=2048)),
        "fp_pattern": lambda m: on(m.pattern_fingerprint(n_bits=2048)),
        "fp_layered": lambda m: on(m.fingerprint_layered(fp_size=2048)),
        "fp_rdkit": lambda m: on(m.topological_fingerprint(fp_size=2048)),
        "fp_avalon": lambda m: on(m.avalon_fingerprint(n_bits=2048)),
        "chi0v": ck.calc_chi_0v, "chi1v": ck.calc_chi_1v, "chi2v": ck.calc_chi_2v,
        "chi3v": ck.calc_chi_3v, "chi4v": ck.calc_chi_4v,
        "chi0": ck.calc_chi_0, "chi1": ck.calc_chi_1,
        "kappa1": ck.calc_kappa_1, "kappa2": ck.calc_kappa_2, "kappa3": ck.calc_kappa_3,
        "hall_kier_alpha": ck.calc_hall_kier_alpha,
        "labute_asa": ck.calc_labute_asa,
        "slogp_vsa": lambda m: list(ck.calc_slogp_vsa(m)),
        "smr_vsa": lambda m: list(ck.calc_smr_vsa(m)),
        "mqn": lambda m: list(ck.calc_mqns(m)),
        "num_rings": ck.calc_num_rings,
        "num_aliphatic_rings": ck.calc_num_aliphatic_rings,
        "num_saturated_rings": ck.calc_num_saturated_rings,
        "num_aromatic_heterocycles": ck.calc_num_aromatic_heterocycles,
        "num_aliphatic_heterocycles": ck.calc_num_aliphatic_heterocycles,
        "num_saturated_heterocycles": ck.calc_num_saturated_heterocycles,
        "num_heteroatoms": ck.calc_num_heteroatoms,
        "num_amide_bonds": ck.calc_num_amide_bonds,
        "num_bridgehead_atoms": ck.calc_num_bridgehead_atoms,
        "num_spiro_atoms": ck.calc_num_spiro_atoms,
        "heavy_atoms": ck.calc_num_heavy_atoms,
        "unspecified_stereocenters": ck.calc_num_unspecified_atom_stereo_centers,
        "smiles_kekule": lambda m: m.to_smiles(kekule=True),
        "smiles_noniso": lambda m: m.to_smiles(isomeric_smiles=False),
        "smiles_explicit": lambda m: m.to_smiles(all_bonds_explicit=True, all_hs_explicit=True),
        "smiles_addhs": lambda m: m.with_hydrogens().to_smiles(),
        "inchi_read": lambda m: (lambda i: ck.Molecule.from_inchi(i).to_smiles() if i else None)(
            rdkit_inchi_of(rdkit_readback(m.to_smiles(), "smiles"))),
        "murcko_scaffold": lambda m: m.murcko_scaffold().to_smiles(),
        "smarts_write": lambda m: m.to_smarts(),
        "cx_smarts": lambda m: m.to_cx_smarts(),
        "pdb_block": lambda m: m.to_pdb_block(),
        "stereoisomer_count": lambda m: m.stereoisomer_count(),
        "chiral_centers": lambda m: [list(c) for c in m.find_chiral_centers(include_unassigned=True)],
        "morgan3_2048": lambda m: on(m.fingerprint_morgan(radius=3, n_bits=2048)),
        "morgan2_chiral": lambda m: on(m.fingerprint_morgan(radius=2, n_bits=2048, include_chirality=True)),
        "morgan2_countsim": lambda m: on(m.fingerprint_morgan(radius=2, n_bits=2048, count_simulation=True)),
        "stereoisomers": lambda m: sorted(x.to_smiles() for x in m.stereoisomers()),
        "morgan2_bitinfo": lambda m: _bit_info(
            m.fingerprint_morgan_with_output(radius=2, n_bits=2048).additional_output().bit_info_map()),
        "pdb_read": lambda m: ck.Molecule.from_pdb_block(m.to_pdb_block()).to_smiles(),
        "extended_murcko": lambda m: m.net_scaffold().to_smiles(),
        "embed3d": _embed3d,
    })
    raw = {"canonical_smiles": lambda m: m.to_smiles(), "molblock": lambda m: m.to_2d_sdf_string()}
    return {"version": getattr(ck, "__version__", "unknown"), "parse": ck.Molecule.from_smiles,
            "ops": ops, "raw": raw}


def cosmolkit_version_tuple(version: str) -> tuple[int, int]:
    """(major, minor) of a COSMolKit version string such as ``0.5.0-rc.15``."""
    parts = version.split(".")
    try:
        return int(parts[0]), int(parts[1])
    except (IndexError, ValueError):
        return (0, 0)


def cosmolkit_engine():
    """COSMolKit adapters for the installed release: the 0.3 API, or the
    reworked 0.5 API (``*_with_params`` names, generator fingerprints)."""
    import cosmolkit as ck

    version = getattr(ck, "__version__", "unknown")
    if cosmolkit_version_tuple(version) >= (0, 5):
        return _cosmolkit_05_engine()
    return _cosmolkit_03_engine()


# Operations COSMolKit 0.5 no longer offers (reported as unsupported).
COSMOLKIT_05_REMOVED = {
    "inchi": "InChI writer removed in 0.5",
    "inchikey": "InChI writer removed in 0.5",
    "inchi_read": "InChI reader removed in 0.5",
    "fp_avalon": "Avalon fingerprint removed in 0.5",
    "murcko_scaffold": "murcko_scaffold removed in 0.5",
    "extended_murcko": "net_scaffold removed in 0.5",
    "pdb_block": "Molecule.to_pdb_block removed in 0.5 (bio PDB writer works on BioStructure)",
    "pdb_read": "Molecule.from_pdb_block removed in 0.5",
    "chiral_centers": "find_chiral_centers now returns chiral-tag names for every atom, "
                      "not FindMolChiralCenters' R/S/? list",
    "peoe_vsa": "not offered", "bertz_ct": "not offered", "balaban_j": "not offered", "ipc": "not offered",
}


def _cosmolkit_05_engine():
    import cosmolkit as ck

    on = lambda fp: sorted(fp.on_bits())  # noqa: E731

    def morgan(radius, **kw):
        gen = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=radius, fp_size=2048, **kw))
        return lambda m: on(m.morgan_fingerprint_with_generator(gen))

    def morgan_bitinfo(m):
        gen = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=2, fp_size=2048))
        out = ck.FingerprintAdditionalOutput()
        out.allocate_bit_info_map()
        m.morgan_fingerprint_with_generator(gen, output=out)
        return _bit_info(out.bit_info_map())

    def smiles_with(**kw):
        p = ck.SmilesWriteParams(**kw)
        return lambda m: m.to_smiles_with_params(p)

    def cip(m):
        lab = m.with_cip_labels()
        out = [["a", a.id(), str(a.cip_descriptor().value)] for a in lab.atoms() if a.cip_descriptor()]
        out += [["b", *sorted((b.begin(), b.end())), str(b.cip_descriptor().value)]
                for b in lab.bonds() if b.cip_descriptor()]
        return sorted(out)

    def smarts(q):
        qm = ck.parse_smarts(q)
        return lambda m: _smarts_sets(_ck05_mapping(r) for r in m.substruct_matches(qm))

    def embed3d(m):
        p = ck.EmbedParams.etkdg_v3()
        p = p.with_json('{"randomSeed": 42}') or p
        try:
            r = m.with_hydrogens().with_3d_conformer_result_with_params(p)
        except Exception:  # noqa: BLE001 - embedding failure
            return None
        if not r.ok:
            return None
        mol = r.molecule
        mol = mol() if callable(mol) else mol
        return None if mol is None or mol.num_3d_conformers() == 0 else _ck05_coords(mol.coordinates_3d())

    ops = {
        "canonical_smiles": lambda m: m.to_smiles(),
        "canonical_smiles_rt": lambda m: rdkit_readback(m.to_smiles(), "smiles"),
        "formula": lambda m: m.molecular_formula(),
        "mol_wt": lambda m: m.molecular_weight(),
        "exact_mw": lambda m: m.exact_molecular_weight(),
        "tpsa": lambda m: m.tpsa(),
        "logp": lambda m: m.crippen_descriptors().logp,
        "mr": lambda m: m.crippen_descriptors().molar_refractivity,
        "hba": lambda m: m.num_hba(),
        "hbd": lambda m: m.num_hbd(),
        "rotatable_bonds": lambda m: m.num_rotatable_bonds(),
        "aromatic_rings": lambda m: m.num_aromatic_rings(),
        "fsp3": lambda m: m.fraction_csp3(),
        "qed": lambda m: m.qed(),
        "stereocenters": lambda m: m.num_atom_stereo_centers(),
        "morgan2_2048": morgan(2),
        "maccs": lambda m: sorted(i + 1 for i in m.maccs_fingerprint().on_bits()),
        "cip": cip,
        "molblock_rt": lambda m: rdkit_readback(m.to_sdf_2d(), "mol"),
    }
    for q in SMARTS_QUERIES:
        ops["smarts:" + q] = smarts(q)
    ops.update({
        "fp_atom_pair": lambda m: on(m.atom_pair_fingerprint()),
        "fp_torsion": lambda m: on(m.legacy_topological_torsion_fingerprint()),
        "fp_pattern": lambda m: on(m.pattern_fingerprint()),
        "fp_layered": lambda m: on(m.layered_fingerprint()),
        "fp_rdkit": lambda m: on(m.topological_fingerprint()),
        "chi0v": lambda m: m.chi_0_v(), "chi1v": lambda m: m.chi_1_v(), "chi2v": lambda m: m.chi_2_v(),
        "chi3v": lambda m: m.chi_3_v(), "chi4v": lambda m: m.chi_4_v(),
        "chi0": lambda m: m.chi_0(), "chi1": lambda m: m.chi_1(),
        "kappa1": lambda m: m.kappa_1(), "kappa2": lambda m: m.kappa_2(), "kappa3": lambda m: m.kappa_3(),
        "hall_kier_alpha": lambda m: m.hall_kier_alpha(),
        "labute_asa": lambda m: m.labute_asa(),
        "slogp_vsa": lambda m: list(m.slogp_vsa()),
        "smr_vsa": lambda m: list(m.smr_vsa()),
        "mqn": lambda m: list(m.mqns()),
        "num_rings": lambda m: m.num_rings(),
        "num_aliphatic_rings": lambda m: m.num_aliphatic_rings(),
        "num_saturated_rings": lambda m: m.num_saturated_rings(),
        "num_aromatic_heterocycles": lambda m: m.num_aromatic_heterocycles(),
        "num_aliphatic_heterocycles": lambda m: m.num_aliphatic_heterocycles(),
        "num_saturated_heterocycles": lambda m: m.num_saturated_heterocycles(),
        "num_heteroatoms": lambda m: m.num_heteroatoms(),
        "num_amide_bonds": lambda m: m.num_amide_bonds(),
        "num_bridgehead_atoms": lambda m: m.num_bridgehead_atoms(),
        "num_spiro_atoms": lambda m: m.num_spiro_atoms(),
        "heavy_atoms": lambda m: m.num_heavy_atoms(),
        "unspecified_stereocenters": lambda m: m.num_unspecified_atom_stereo_centers(),
        "smiles_kekule": smiles_with(do_kekule=True),
        "smiles_noniso": smiles_with(do_isomeric_smiles=False),
        "smiles_explicit": smiles_with(all_bonds_explicit=True, all_hydrogens_explicit=True),
        "smiles_addhs": lambda m: m.with_hydrogens().to_smiles(),
        "smarts_write": lambda m: m.to_smarts(),
        "cx_smarts": lambda m: m.to_cx_smarts(),
        "stereoisomer_count": lambda m: m.stereoisomer_count(),
        "morgan3_2048": morgan(3),
        "morgan2_chiral": morgan(2, include_chirality=True),
        "morgan2_countsim": morgan(2, count_simulation=True),
        "atom_pair_counts": lambda m: sorted([k, v] for k, v in
                                             m.atom_pair_count_fingerprint().nonzero_elements().items()),
        "stereoisomers": lambda m: sorted(x.to_smiles() for x in m.enumerate_stereoisomers()),
        "morgan2_bitinfo": morgan_bitinfo,
        "embed3d": embed3d,
    })
    ops.update(_cosmolkit_05_new_ops(ck))
    raw = {"canonical_smiles": lambda m: m.to_smiles(), "molblock": lambda m: m.to_sdf_2d()}
    return {"version": getattr(ck, "__version__", "unknown"),
            "parse": _with_input_smiles(ck.Molecule.from_smiles),
            "ops": ops, "raw": raw, "unsupported": dict(COSMOLKIT_05_REMOVED)}


def _ck05_mapping(match):
    for name in ("atom_mapping", "atoms", "mapping"):
        v = getattr(match, name, None)
        if v is not None:
            return v() if callable(v) else v
    return match


def _ck05_coords(block):
    v = block() if callable(block) else block
    return v.tolist() if hasattr(v, "tolist") else [list(p) for p in v]


# ------------------------------------------------- COSMolKit 0.5 surfaces
# Fixed reaction set: (name, SMIRKS, partner SMILES for reactant templates
# 2..n). The corpus molecule is reactant 1. Products are compared as sorted,
# de-duplicated product sets, each product written by the engine and read
# back by RDKit (``rdkit_readback``) so that engines with different writers
# compare on the molecules themselves; "<invalid>" marks a product RDKit
# cannot read back (or sanitize, for RDKit's own products).
REACTIONS = [
    ("amide_acid", "[C:1](=[O:2])[OH].[N;!H0;!$(NC=O);!$(N=*);!$(N#*):3]>>[C:1](=[O:2])[N:3]", ["NCc1ccccc1"]),
    ("amide_amine", "[N;!H0;!$(NC=O);!$(N=*);!$(N#*);!$(Nc);!$([N+]):1].[C:2](=[O:3])[OH]>>[N:1][C:2]=[O:3]", ["CC(=O)O"]),
    ("esterification", "[C:1](=[O:2])[OH].[OH:3][CH3:4]>>[C:1](=[O:2])[O:3][C:4]", ["CO"]),
    ("suzuki_like", "[c:1][Br,I].[c:2]B(O)O>>[c:1][c:2]", ["OB(O)c1ccccc1"]),
    ("sulfonamide", "[N;!H0;!$(NC=O);!$(NS(=O)=O);!$([N+]):1].Cl[S:2](=[O:3])(=[O:4])[C:5]>>[N:1][S:2](=[O:3])(=[O:4])[C:5]", ["CS(=O)(=O)Cl"]),
    ("boc_deprotection", "[N:1]C(=O)OC([CH3])([CH3])[CH3]>>[N:1]", []),
    ("ester_hydrolysis", "[C:1](=[O:2])O[CX4]>>[C:1](=[O:2])O", []),
    ("amide_hydrolysis", "[C:1](=[O:2])[NH:3][C:4]>>[C:1](=[O:2])O.[N:3][C:4]", []),
    ("nitro_reduction", "[c:1][N+](=O)[O-]>>[c:1]N", []),
    ("alcohol_oxidation", "[CH2:1][OH:2]>>[CH:1]=[O:2]", []),
    ("ketone_reduction", "[C:1](=[O:2])([#6:3])[#6:4]>>[C:1]([OH:2])([#6:3])[#6:4]", []),
    ("n_methylation", "[NH2:1][c:2]>>C[NH:1][c:2]", []),
    ("halogen_exchange", "[c:1]Cl>>[c:1]F", []),
    ("n_oxidation", "[n;H0;+0:1]>>[n+:1][O-]", []),
    ("lactam_ring_closure", "([C:1](=[O:2])[OH].[NH2:3])>>[C:1](=[O:2])[N:3]", []),
    ("diels_alder", "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1", ["C=CC(=O)OC"]),
    ("stereo_retaining_acylation", "[C@@H:1]([NH2:2])([#6:3])[C:4]=[O:5]>>[C@@H:1]([NH:2]C(C)=O)([#6:3])[C:4]=[O:5]", []),
    ("stereo_inversion", "[C@:1]([N:2])>>[C@@:1]([N:2])", []),
    ("unmapped_methyl_ester", "C(=O)[OH]>>C(=O)OC", []),
    ("partially_mapped_ether", "[c:1][OH]>>[c:1]OC", []),
]

# Deliberately broken SMILES (read with sanitize=False) for the
# DetectChemistryProblems check (check_chemistry_problems.py).
BROKEN_SMILES = [
    "CN(C)(C)(C)C", "C1=CC=CC=C1C(=O)(O)O", "c1ccccc1c", "c1cccc1", "Cc1nccc1",
    "F(C)C", "O(C)(C)C", "CC(C)(C)(C)C", "c1ccc2c(c1)cccn2C", "n1cccc1",
    "C1=CC=CC=C1=C", "Cl(C)C", "[NH4](C)", "c1ccco1C", "B(C)(C)(C)C",
    "C(=O)=O=C", "S(C)(C)(C)(C)(C)(C)C", "C1CC1(C)(C)C", "c1ccccc1.N(C)(C)(C)C", "[O-](C)C",
]


def _round_nested(v, nd=6):
    if isinstance(v, float):
        r = round(v, nd)
        return 0.0 if r == 0 else r
    if isinstance(v, (list, tuple)):
        return [_round_nested(x, nd) for x in v]
    return v


def _distance_matrix_norm(rows):
    """Topological distances as ints (-1 for disconnected atom pairs)."""
    return [[int(round(d)) if d < 1e7 else -1 for d in row] for row in rows]


def rdkit_embedded_h(smiles: str):
    """``Chem.AddHs(Chem.MolFromSmiles(smiles))`` embedded with
    ``EmbedMolecule(randomSeed=42)``: (molecule, coordinates) or None. The
    shared input of every engine's 3D surface ops."""
    Chem = _rd()
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=42) != 0:
        return None
    return mh, mh.GetConformer().GetPositions().tolist()


def _rdkit_reaction(smirks):
    from rdkit.Chem import AllChem
    rxn = AllChem.ReactionFromSmarts(smirks)
    rxn.Initialize()
    return rxn


def _product_sets(sets):
    """Sorted unique product sets: each a '.'-joined list of RDKit-canonical
    product SMILES in product-template order."""
    return sorted({".".join(rdkit_readback(p, "smiles") or "<invalid>" if p is not None else "<invalid>"
                            for p in ps) for ps in sets})


def fragment_atoms(n: int) -> list[int]:
    """The deterministic atom subset of the fragment-SMILES op: the first
    half of the atoms (at least one)."""
    return list(range(max(1, n // 2)))


def _rdkit_new_ops(Chem):
    from rdkit.Chem import AllChem, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.MolStandardize import rdMolStandardize

    enumerator = rdMolStandardize.TautomerEnumerator()
    sparse_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)

    def run_rxn(smirks, partners):
        rxn = _rdkit_reaction(smirks)
        pmols = [Chem.MolFromSmiles(x) for x in partners]

        def run(m):
            sets = []
            for prods in rxn.RunReactants((m, *pmols)):
                smis = []
                for p in prods:
                    try:
                        Chem.SanitizeMol(p)
                        smis.append(Chem.MolToSmiles(p))
                    except Exception:  # noqa: BLE001 - unsanitizable product
                        smis.append(None)
                sets.append(smis)
            return _product_sets(sets)
        return run

    def mmff(m):
        e = _rd_embed_from_mol(m)
        if e is None:
            return None
        mh, _ = e
        props = AllChem.MMFFGetMoleculeProperties(mh)
        if props is None:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(mh, props)
        return _round_nested([ff.CalcEnergy(), list(ff.CalcGrad())])

    def uff(m):
        e = _rd_embed_from_mol(m)
        if e is None:
            return None
        mh, _ = e
        ff = AllChem.UFFGetMoleculeForceField(mh)
        return _round_nested([ff.CalcEnergy(), list(ff.CalcGrad())])

    def dm3d(m):
        e = _rd_embed_from_mol(m)
        if e is None:
            return None
        return _round_nested(Chem.Get3DDistanceMatrix(e[0]).tolist())

    def problems(m):
        return sorted([p.GetType(), _problem_atoms(p)] for p in Chem.DetectChemistryProblems(m))

    ops = {
        "tautomer_canonical": lambda m: Chem.MolToSmiles(enumerator.Canonicalize(m)),
        "tautomer_set": lambda m: sorted({Chem.MolToSmiles(t) for t in enumerator.Enumerate(m)}),
        "cx_smiles": lambda m: Chem.MolToCXSmiles(m),
        "random_smiles5": lambda m: list(Chem.MolToRandomSmilesVect(m, 5, randomSeed=42)),
        "fragment_smiles": lambda m: Chem.MolFragmentToSmiles(m, atomsToUse=fragment_atoms(m.GetNumAtoms())),
        "chemistry_problems": problems,
        "distance_matrix": lambda m: _distance_matrix_norm(Chem.GetDistanceMatrix(m).tolist()),
        "distance_matrix_3d": dm3d,
        "mmff_energy_gradient": mmff,
        "uff_energy_gradient": uff,
        "morgan2_sparse_counts": lambda m: sorted([int(k), v] for k, v in
                                                  sparse_gen.GetSparseCountFingerprint(m).GetNonzeroElements().items()),
        "torsion_legacy_counts": lambda m: sorted([int(k), v] for k, v in
                                                  rdMolDescriptors.GetTopologicalTorsionFingerprint(m).GetNonzeroElements().items()),
    }
    for name, smirks, partners in REACTIONS:
        ops["rxn:" + name] = run_rxn(smirks, partners)
    return ops


_RD_EMBED_MEMO: list = [None, None]


def _rd_embed_from_mol(m):
    """RDKit's embedded explicit-H copy of an RDKit molecule (input order),
    memoized for the molecule the ops are currently running on."""
    if _RD_EMBED_MEMO[0] is m:
        return _RD_EMBED_MEMO[1]
    _RD_EMBED_MEMO[0], _RD_EMBED_MEMO[1] = m, _rd_embed_uncached(m)
    return _RD_EMBED_MEMO[1]


def _rd_embed_uncached(m):
    Chem = _rd()
    from rdkit.Chem import AllChem
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=42) != 0:
        return None
    return mh, mh.GetConformer().GetPositions().tolist()


def _problem_atoms(p):
    if hasattr(p, "GetAtomIndices"):
        return sorted(p.GetAtomIndices())
    return [p.GetAtomIdx()]


class _InputSmiles:
    """The SMILES the engine's ``parse`` last read: 3D surface ops embed the
    input molecule with RDKit (same atom order in every engine)."""
    value: str | None = None


def _with_input_smiles(parse):
    def wrapped(smiles):
        _InputSmiles.value = smiles
        return parse(smiles)
    return wrapped


_EMBED_MEMO: dict = {}


def _embedded_input():
    """``rdkit_embedded_h`` of the current input SMILES (memoized for it)."""
    smi = _InputSmiles.value
    if smi is None:
        return None
    if smi not in _EMBED_MEMO:
        _EMBED_MEMO.clear()
        _EMBED_MEMO[smi] = rdkit_embedded_h(smi)
    return _EMBED_MEMO[smi]


def _chematic_new_ops(c):
    def readback(mol):
        return rdkit_readback(mol.smiles, "smiles")

    def run_rxn(smirks, partners):
        pmols = [c.from_smiles(x) for x in partners]

        def run(m):
            r = c.run_smirks_checked(smirks, [m, *pmols], rdkit_compat=True)
            valence = r["status"] in ("typed_refusal", "partial_products") and r.get("reason") == "product_valence"
            if r["status"] not in ("products", "no_match") and not valence:
                raise RefusedError(f"{r['status']}: {r.get('reason')}")
            # RDKit returns the product sets its sanitization rejects too
            # (``<invalid>`` products); chematic reports them separately.
            sets = [[p.smiles for p in ps] for ps in r.get("products") or []]
            sets += [[p.smiles if ok else None for p, ok in ps] for ps in r.get("rejected_products") or []]
            return _product_sets(sets)
        return run

    def mmff(m):
        e = _embedded_input()
        if e is None:
            return None
        mh = m.add_hydrogens()
        energy, grad = mh._rdkit_mmff_terms(e[1])
        return _round_nested([energy, list(grad)])

    def uff(m):
        e = _embedded_input()
        if e is None:
            return None
        mh = m.add_hydrogens()
        return _round_nested([mh.rdkit_uff_energy(e[1]), list(mh.rdkit_uff_gradient(e[1]))])

    def dm3d(m):
        e = _embedded_input()
        if e is None:
            return None
        return _round_nested(m.add_hydrogens().rdkit_distance_matrix_3d(e[1]))

    def problems(m):
        return sorted([t, sorted(a)] for t, a in m.rdkit_chemistry_problems())

    ops = {
        "tautomer_canonical": lambda m: readback(m.canonical_tautomer()),
        "tautomer_set": lambda m: sorted({readback(t) for t in m.enumerate_tautomers()}),
        "cx_smiles": lambda m: m.rdkit_cx_smiles(),
        "random_smiles5": lambda m: m.rdkit_random_smiles(5, 42),
        "fragment_smiles": lambda m: m.rdkit_fragment_smiles(fragment_atoms(m.rdkit_num_atoms())),
        "chemistry_problems": problems,
        "distance_matrix": lambda m: _distance_matrix_norm(m.rdkit_distance_matrix()),
        "distance_matrix_3d": dm3d,
        "mmff_energy_gradient": mmff,
        "uff_energy_gradient": uff,
        "morgan2_sparse_counts": lambda m: sorted([int(k), v] for k, v in m.rdkit_morgan_sparse_counts(2).items()),
        "torsion_legacy_counts": lambda m: sorted([int(k), v] for k, v in m.rdkit_legacy_torsion_counts().items()),
    }
    for name, smirks, partners in REACTIONS:
        ops["rxn:" + name] = run_rxn(smirks, partners)
    return ops


def _cosmolkit_05_new_ops(ck):
    """Harness operations for surfaces new in COSMolKit 0.5."""
    def run_rxn(smirks, partners):
        rxn = ck.parse_smirks(smirks)
        pmols = [ck.Molecule.from_smiles(x) for x in partners]

        def run(m):
            sets = rxn.run([m, *pmols], ck.ReactionRunParams())
            return _product_sets([[p.to_smiles() for p in ps] for ps in sets])
        return run

    def embedded(m):
        e = _embedded_input()
        if e is None:
            return None
        return m.with_hydrogens().with_only_3d_conformer(e[1])

    def energy_gradient(kind):
        def run(m):
            mh = embedded(m)
            if mh is None:
                return None
            g = getattr(mh, kind + "_energy_gradient")()
            energy = g.energy() if callable(g.energy) else g.energy
            grad = g.gradient() if callable(g.gradient) else g.gradient
            return _round_nested([energy, list(grad)])
        return run

    def dm(values, n):
        return [values[i * n:(i + 1) * n] for i in range(n)]

    def dm3d(m):
        mh = embedded(m)
        if mh is None:
            return None
        d = mh.distance_matrix_3d()
        return _round_nested(dm(d.values(), d.dimension()))

    def topo_dm(m):
        d = m.distance_matrix()
        return _distance_matrix_norm(dm(d.values(), d.dimension()))

    def problems(m):
        found = m.detect_chemistry_problems().problems
        return sorted(_ck05_problem(p) for p in (found() if callable(found) else found))

    sparse_gen = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=2))
    ops = {
        "tautomer_canonical": lambda m: m.canonical_tautomer().to_smiles(),
        "tautomer_set": lambda m: sorted(set(m.enumerate_tautomers().canonical_smiles())),
        "cx_smiles": lambda m: m.to_cx_smiles(),
        "random_smiles5": lambda m: list(m.to_random_smiles(5, 42)),
        "fragment_smiles": lambda m: m.to_fragment_smiles(fragment_atoms(m.num_atoms())),
        "chemistry_problems": problems,
        "distance_matrix": topo_dm,
        "distance_matrix_3d": dm3d,
        "mmff_energy_gradient": energy_gradient("mmff"),
        "uff_energy_gradient": energy_gradient("uff"),
        "morgan2_sparse_counts": lambda m: sorted([int(k), v] for k, v in
                                                  m.morgan_sparse_count_fingerprint_with_generator(sparse_gen)
                                                  .nonzero_elements().items()),
        "torsion_legacy_counts": lambda m: sorted([int(k), v] for k, v in
                                                  m.legacy_topological_torsion_sparse_count_fingerprint()
                                                  .nonzero_elements().items()),
    }
    for name, smirks, partners in REACTIONS:
        ops["rxn:" + name] = run_rxn(smirks, partners)
    return ops


def _ck05_problem(p):
    """A COSMolKit ChemistryProblem as RDKit's ``[GetType(), atoms]`` (the
    0.5 error carries a kind and RDKit's message text)."""
    import re
    err = p.error() if callable(p.error) else p.error
    kind = err.kind() if callable(err.kind) else err.kind
    text = str(err)
    if kind == "Valence":
        m = re.search(r"atom # (\d+)", text)
        return ["AtomValenceException", [int(m.group(1))] if m else []]
    if kind == "Kekulize":
        m = re.search(r"aromatic atom (\d+) is not in a ring", text)
        if m:
            return ["AtomKekulizeException", [int(m.group(1))]]
        return ["KekulizeException", sorted(int(x) for x in re.findall(r"AtomId\((\d+)\)", text))]
    return [str(kind), []]


def chemistry_problems_unsanitized(engine: str, smiles: str):
    """``DetectChemistryProblems`` of ``smiles`` read with sanitize=False, as
    sorted ``[type, atoms]`` pairs; None where the engine cannot read it."""
    if engine == "rdkit":
        Chem = _rd()
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        return None if m is None else sorted(
            [p.GetType(), _problem_atoms(p)] for p in Chem.DetectChemistryProblems(m))
    if engine == "cosmolkit":
        import cosmolkit as ck
        try:
            m = ck.Molecule.from_smiles_with_params(smiles, ck.SmilesParseParams(sanitize=False))
        except Exception:  # noqa: BLE001 - unreadable input
            return None
        problems = m.detect_chemistry_problems().problems
        problems = problems() if callable(problems) else problems
        return sorted(_ck05_problem(p) for p in problems)
    if engine == "chematic":
        import chematic as c
        try:
            found = c.rdkit_detect_chemistry_problems(smiles)
        except ValueError:
            return None
        return sorted([t, sorted(a)] for t, a in found)
    raise UnsupportedError(f"{engine}: no DetectChemistryProblems equivalent")


def reaction_smarts_roundtrip(engine: str, smirks: str) -> str:
    """The engine's reaction-SMARTS writer applied to the parsed reaction
    (RDKit: ``ReactionToSmarts(ReactionFromSmarts(s))``)."""
    if engine == "rdkit":
        from rdkit.Chem import AllChem
        return AllChem.ReactionToSmarts(AllChem.ReactionFromSmarts(smirks))
    if engine == "cosmolkit":
        import cosmolkit as ck
        return ck.parse_smirks(smirks).to_smirks()
    if engine == "chematic":
        import chematic as c
        return c.rdkit_reaction_to_smarts(smirks)
    raise UnsupportedError(f"{engine}: no reaction SMARTS writer")


class RefusedError(Exception):
    """The engine declined this cell with a typed reason (not a wrong answer)."""


class UnsupportedError(Exception):
    """The installed artifact does not offer this operation."""


ENGINES = {"rdkit": rdkit_engine, "chematic": chematic_engine, "cosmolkit": cosmolkit_engine}
