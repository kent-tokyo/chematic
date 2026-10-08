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

API_NOTES = {
    "chematic": {
        "canonical_smiles": "Mol.rdkit_smiles (RDKit-compatible writer); canonical_smiles_rt reads back the native Mol.smiles",
        "smarts": "Mol.find_matches_rdkit_parity (opt-in RDKit 2026.03.6 profile; typed refusals)",
        "cip": "Mol.cip_stereo(mode='accurate'); abstentions from cip_stereo_unresolved()",
        "morgan2_2048": "Mol.rdkit_ecfp_config(2, 2048)",
        "mol_wt": "Mol.rdkit_mw", "tpsa": "Mol.rdkit_tpsa",
        "aromatic_rings": "Mol.rdkit_aromatic_ring_count",
        "molblock_rt": "Mol.to_mol_block() (default writer)",
        "inchi": "Mol.standard_inchi (needs the native-inchi build feature)",
    },
    "cosmolkit": {
        "smarts": "get_substruct_matches(mol, parse_smarts(q))",
        "cip": "Molecule.with_cip_labels()",
        "molblock_rt": "Molecule.to_2d_sdf_string()",
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


def _smarts_sets(matches) -> list[list[int]]:
    return sorted({tuple(sorted(m)) for m in matches})


# --------------------------------------------------------------------- RDKit
def rdkit_engine():
    Chem = _rd()
    from rdkit import rdBase
    from rdkit.Chem import (QED, Crippen, Descriptors, MACCSkeys, rdCIPLabeler,
                            rdFingerprintGenerator)
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
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
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
        "morgan3_2048": lambda m: list(gen3.GetFingerprint(m).GetOnBits()),
        "morgan2_chiral": lambda m: list(gen_chiral.GetFingerprint(m).GetOnBits()),
        "morgan2_countsim": lambda m: list(gen_count.GetFingerprint(m).GetOnBits()),
        "atom_pair_counts": lambda m: sorted([k, v] for k, v in
                                             RD.GetHashedAtomPairFingerprint(m, nBits=2048).GetNonzeroElements().items()),
        "stereoisomers": lambda m: sorted(Chem.MolToSmiles(x) for x in EnumerateStereoisomers(m)),
    })
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

    def parse(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("RDKit could not parse")
        return mol

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
        "murcko_scaffold": lambda m: m.scaffold().rdkit_smiles,
        "morgan3_2048": lambda m: _bits(m.rdkit_ecfp_config(3, 2048)),
        "morgan2_chiral": lambda m: _bits(m.rdkit_ecfp_config(2, 2048, include_chirality=True)),
        "stereoisomers": lambda m: sorted(x.rdkit_smiles for x in m.enumerate_stereoisomers()),
    })
    for name, kwargs in (("smiles_kekule", {"kekule": True}), ("smiles_noniso", {"isomeric": False}),
                         ("smiles_explicit", {"all_bonds_explicit": True, "all_hs_explicit": True})):
        ops[name] = (lambda kw: lambda m: _rdkit_smiles_with(m, kw))(kwargs)

    def _inchi_read(m):
        inchi = rdkit_inchi_of(rdkit_readback(m.rdkit_smiles, "smiles"))
        return c.from_inchi(inchi).rdkit_smiles if inchi else None
    def _rdkit_smiles_with(m, kwargs):
        writer = getattr(m, "rdkit_smiles_with", None)
        if writer is None:
            raise UnsupportedError("rdkit_smiles options are not exposed")
        return writer(**kwargs)

    raw = {"canonical_smiles": lambda m: m.smiles, "molblock": molblock}
    return {"version": c.__version__, "parse": c.from_smiles, "ops": ops, "raw": raw}


# ----------------------------------------------------------------- COSMolKit
def cosmolkit_engine():
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
        "morgan3_2048": lambda m: on(m.fingerprint_morgan(radius=3, n_bits=2048)),
        "morgan2_chiral": lambda m: on(m.fingerprint_morgan(radius=2, n_bits=2048, include_chirality=True)),
        "morgan2_countsim": lambda m: on(m.fingerprint_morgan(radius=2, n_bits=2048, count_simulation=True)),
        "stereoisomers": lambda m: sorted(x.to_smiles() for x in m.stereoisomers()),
    })
    raw = {"canonical_smiles": lambda m: m.to_smiles(), "molblock": lambda m: m.to_2d_sdf_string()}
    return {"version": getattr(ck, "__version__", "unknown"), "parse": ck.Molecule.from_smiles,
            "ops": ops, "raw": raw}


class RefusedError(Exception):
    """The engine declined this cell with a typed reason (not a wrong answer)."""


class UnsupportedError(Exception):
    """The installed artifact does not offer this operation."""


ENGINES = {"rdkit": rdkit_engine, "chematic": chematic_engine, "cosmolkit": cosmolkit_engine}
