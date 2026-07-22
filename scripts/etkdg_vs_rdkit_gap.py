#!/usr/bin/env python3
"""Diagnostic comparison: chematic-3d ETKDG-style conformer generation vs RDKit ETKDGv3.

Research/diagnosis tool for docs/etkdg_3d_gap_rfc.md — NOT a regression gate.
Re-run with: .venv/bin/python scripts/etkdg_vs_rdkit_gap.py [--out validation/results/etkdg_vs_rdkit_summary.json]

For every corpus molecule this measures, per molecule and aggregated:
  - embedding success rate (named failure buckets, no silent drops)
  - chirality retention: declared vs re-perceived-from-3D, split into
    coverage (could the engine even assess it) and match-rate (given coverage)
  - bond-length / bond-angle violations vs an external covalent-radius/
    hybridization-angle reference (same reference applied to both engines)
  - aromatic-ring planarity (RMS deviation from best-fit plane)
  - a *same-yardstick* MMFF94 energy delta: chematic's own
    mmff94_energy_breakdown() applied to (a) chematic's conformer and
    (b) RDKit's conformer projected onto chematic's heavy-atom topology.
    RDKit's own native (all-atom) MMFF/UFF energy is also reported but never
    diffed numerically against chematic's number (different atom sets).
  - Kabsch-aligned heavy-atom RMSD vs RDKit's conformer
  - runtime per engine
  - (flexible-molecule subset) ensemble diversity + duplicate rate for both
    engines' multi-conformer generation

See docs/etkdg_3d_gap_rfc.md for interpretation and the phased roadmap.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

import chematic
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Corpus
# ---------------------------------------------------------------------------
# Curated, not bulk: spans the structural classes that stress a distance-geometry
# embedder (advisor-reviewed scope — tens of molecules, not thousands, since the
# discriminating metrics here are RMSD/energy/chirality, not bond lengths).

CORPUS = [
    # -- rigid rings --------------------------------------------------------
    ("benzene", "c1ccccc1", "rigid_ring"),
    ("naphthalene", "c1ccc2ccccc2c1", "fused_aromatic"),
    ("pyridine", "c1ccncc1", "rigid_ring"),
    ("furan", "c1ccoc1", "rigid_ring"),
    ("thiophene", "c1ccsc1", "rigid_ring"),
    ("adamantane", "C1CC2CC3CC1CC(C2)C3", "rigid_ring"),
    ("cubane", "C1C2C3C1C4C2C3C4", "rigid_ring"),
    ("cyclohexane", "C1CCCCC1", "rigid_ring"),
    ("cyclopentane", "C1CCCC1", "rigid_ring"),
    ("indole", "c1ccc2[nH]ccc2c1", "fused_aromatic"),
    ("purine", "c1ncc2[nH]cnc2n1", "fused_aromatic"),
    ("quinoline", "c1ccc2ncccc2c1", "fused_aromatic"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1", "fused_aromatic"),
    ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34", "fused_aromatic"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1", "fused_aromatic"),
    # -- flexible chains ------------------------------------------------------
    ("butane", "CCCC", "flexible_chain"),
    ("hexane", "CCCCCC", "flexible_chain"),
    ("decane", "CCCCCCCCCC", "flexible_chain"),
    ("triethylene_glycol", "OCCOCCOCCO", "flexible_chain"),
    ("hexanediol", "OCCCCCCO", "flexible_chain"),
    ("hexadecane", "CCCCCCCCCCCCCCCC", "flexible_chain"),
    # -- macrocycles ----------------------------------------------------------
    ("cyclododecane", "C1CCCCCCCCCCC1", "macrocycle"),
    ("crown_12_4", "O1CCOCCOCCOCC1", "macrocycle"),
    ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1", "macrocycle"),
    # -- sp3 stereocenters with an implicit H (the common drug-like case) -----
    ("l_alanine", "N[C@@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("d_alanine", "N[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("l_serine", "N[C@@H](CO)C(=O)O", "stereocenter_implicit_h"),
    ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O", "stereocenter_implicit_h"),
    ("2_butanol_R", "C[C@H](O)CC", "stereocenter_implicit_h"),
    ("2_butanol_S", "C[C@@H](O)CC", "stereocenter_implicit_h"),
    ("2_chlorobutane_R", "C[C@H](Cl)CC", "stereocenter_implicit_h"),
    ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", "stereocenter_implicit_h"),
    ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1", "stereocenter_implicit_h"),
    ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O", "stereocenter_implicit_h"),
    # -- sp3 stereocenters with 4 distinct heavy substituents (assessable) ----
    ("chfclbr_R", "[C@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("chfclbr_S", "[C@@H](F)(Cl)Br", "stereocenter_quaternary"),
    ("quaternary_1_R", "[C@](F)(Cl)(Br)I", "stereocenter_quaternary"),
    ("quaternary_1_S", "[C@@](F)(Cl)(Br)I", "stereocenter_quaternary"),
    ("quaternary_2_R", "[C@](C)(N)(O)F", "stereocenter_quaternary"),
    ("quaternary_2_S", "[C@@](C)(N)(O)F", "stereocenter_quaternary"),
    # -- alkene E/Z (always assessable: no implicit-H ambiguity) -------------
    ("but2ene_E", "C/C=C/C", "alkene_ez"),
    ("but2ene_Z", r"C/C=C\C", "alkene_ez"),
    ("chloropropene_E", "C(/C=C/C)Cl", "alkene_ez"),
    ("chloropropene_Z", r"C(/C=C\C)Cl", "alkene_ez"),
    ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1", "alkene_ez"),
    ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1", "alkene_ez"),
    ("pent2ene_E", "CC/C=C/C", "alkene_ez"),
    ("pent2ene_Z", r"CC/C=C\C", "alkene_ez"),
    # -- drug-like / mixed rigid+flexible -------------------------------------
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", "druglike"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "druglike"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "druglike"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1", "druglike"),
    ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1", "druglike"),
    ("penicillin_core", "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O", "druglike"),
    ("testosterone", "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O", "druglike_rigid"),
    ("cholesterol", "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C", "druglike_stress"),
    ("atorvastatin_fragment", "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", "druglike_stress"),
    ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O", "druglike"),
]

# Flexible-molecule subset used for the ensemble-generation section
# (diversity, duplicate rate, pruning): pick molecules with real rotatable bonds.
ENSEMBLE_SUBSET = [
    "hexane", "decane", "triethylene_glycol", "hexanediol", "hexadecane",
    "ibuprofen", "diphenhydramine", "aspirin",
]

# ---------------------------------------------------------------------------
# External reference tables (independent of chematic's own internal tables —
# used identically for both engines so the comparison is fair).
# ---------------------------------------------------------------------------

_PT = Chem.GetPeriodicTable()
_BOND_ORDER_SCALE = {
    Chem.BondType.SINGLE: 1.00,
    Chem.BondType.DOUBLE: 0.87,
    Chem.BondType.TRIPLE: 0.78,
    Chem.BondType.AROMATIC: 0.93,
}
_HYB_ANGLE_DEG = {
    Chem.HybridizationType.SP: 180.0,
    Chem.HybridizationType.SP2: 120.0,
    Chem.HybridizationType.SP3: 109.47,
}
BOND_LEN_TOL_FRAC = 0.15   # +/-15% of reference length
BOND_ANGLE_TOL_DEG = 12.0  # +/-12 degrees of reference angle
GROSS_CLASH_DIST = 0.5     # any heavy-atom pair closer than this -> degenerate embedding
BOND_BLOWUP_REL_ERROR = 0.5  # >50% off the covalent-radius reference -> torn geometry,
                              # not just an imprecise force-field relaxation


def ref_bond_length(rm: Chem.Mol, i: int, j: int) -> float:
    ai, aj = rm.GetAtomWithIdx(i), rm.GetAtomWithIdx(j)
    r0 = _PT.GetRcovalent(ai.GetAtomicNum()) + _PT.GetRcovalent(aj.GetAtomicNum())
    bond = rm.GetBondBetweenAtoms(i, j)
    scale = _BOND_ORDER_SCALE.get(bond.GetBondType(), 1.00) if bond else 1.00
    return r0 * scale


def ref_angle_deg(rm: Chem.Mol, center: int) -> float:
    hyb = rm.GetAtomWithIdx(center).GetHybridization()
    return _HYB_ANGLE_DEG.get(hyb, 109.47)


# ---------------------------------------------------------------------------
# RDKit side
# ---------------------------------------------------------------------------

def rdkit_embed(smiles: str, seed: int = 0xF00D):
    """Return (rm, rmH, status, t_embed_s, t_opt_s, mmff_pre, mmff_post)."""
    rm = Chem.MolFromSmiles(smiles)
    if rm is None:
        return None, None, "rdkit_parse_failed", None, None, None, None
    rmH = Chem.AddHs(rm)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    t0 = time.perf_counter()
    cid = AllChem.EmbedMolecule(rmH, params)
    if cid == -1:
        params.useRandomCoords = True
        cid = AllChem.EmbedMolecule(rmH, params)
    t_embed = time.perf_counter() - t0
    if cid == -1:
        return rm, rmH, "rdkit_embed_failed", t_embed, None, None, None

    mmff_pre = mmff_post = None
    t_opt = None
    try:
        props = AllChem.MMFFGetMoleculeProperties(rmH)
        ff = AllChem.MMFFGetMoleculeForceField(rmH, props)
        if ff is not None:
            mmff_pre = ff.CalcEnergy()
            t1 = time.perf_counter()
            AllChem.MMFFOptimizeMolecule(rmH)
            t_opt = time.perf_counter() - t1
            ff2 = AllChem.MMFFGetMoleculeForceField(rmH, props)
            mmff_post = ff2.CalcEnergy() if ff2 is not None else None
    except Exception:
        pass  # some elements aren't MMFF-parameterised; energy stays None, not fatal

    return rm, rmH, "ok", t_embed, t_opt, mmff_pre, mmff_post


def rdkit_heavy_coords(rmH: Chem.Mol, n_heavy: int) -> np.ndarray:
    conf = rmH.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(n_heavy)])


def rdkit_declared_stereo(rm: Chem.Mol) -> dict[int, str]:
    out = {}
    for idx, code in Chem.FindMolChiralCenters(
        rm, includeUnassigned=True, useLegacyImplementation=False
    ):
        if code in ("R", "S"):
            out[idx] = code
    Chem.DetectBondStereochemistry(rm)
    Chem.AssignStereochemistry(rm, cleanIt=True, force=True)
    for b in rm.GetBonds():
        if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            code = "E" if b.GetStereo() == Chem.BondStereo.STEREOE else "Z"
            out[b.GetBeginAtomIdx()] = code
    return out


def rdkit_perceived_stereo_from_3d(rmH: Chem.Mol, n_heavy: int) -> dict[int, str]:
    Chem.AssignStereochemistryFrom3D(rmH)
    Chem.AssignStereochemistry(rmH, cleanIt=True, force=True)
    out = {}
    for i in range(n_heavy):
        a = rmH.GetAtomWithIdx(i)
        if a.HasProp("_CIPCode"):
            out[i] = a.GetProp("_CIPCode")
    for b in rmH.GetBonds():
        if b.GetBeginAtomIdx() < n_heavy and b.GetEndAtomIdx() < n_heavy:
            if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
                code = "E" if b.GetStereo() == Chem.BondStereo.STEREOE else "Z"
                out[b.GetBeginAtomIdx()] = code
    return out


# ---------------------------------------------------------------------------
# chematic side
# ---------------------------------------------------------------------------

def chematic_embed(smiles: str):
    """Return (cm, coords(list[[x,y,z]]), status, t_embed_s) via the real ETKDG
    path (etkdg.rs -> minimize_mmff94), matching what the task calls
    'ETKDG-style conformer generation'."""
    try:
        cm = chematic.from_smiles(smiles)
    except Exception as e:
        return None, None, f"chematic_parse_exception:{type(e).__name__}", None

    t0 = time.perf_counter()
    try:
        ens = cm.conformer_ensemble(1, 0.0, "mmff94", 0.0)
    except Exception as e:
        return cm, None, f"chematic_embed_exception:{type(e).__name__}", None
    t_embed = time.perf_counter() - t0

    if not ens or len(ens) != 1:
        return cm, None, "chematic_embed_failed_empty_ensemble", t_embed
    coords = ens[0]
    arr = np.array(coords)
    if not np.all(np.isfinite(arr)):
        return cm, coords, "chematic_nonfinite_coords", t_embed
    n = arr.shape[0]
    if n >= 2:
        d = np.linalg.norm(arr[:, None, :] - arr[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        if d.min() < GROSS_CLASH_DIST:
            return cm, coords, "chematic_gross_clash", t_embed
    return cm, coords, "ok", t_embed


# ---------------------------------------------------------------------------
# Shared geometry metrics (external reference, applied identically to both
# engines' heavy-atom-only coordinate arrays).
# ---------------------------------------------------------------------------

def bond_violations(rm: Chem.Mol, coords: np.ndarray) -> dict:
    n_viol, n_total, max_frac = 0, 0, 0.0
    for b in rm.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        r0 = ref_bond_length(rm, i, j)
        r = float(np.linalg.norm(coords[i] - coords[j]))
        frac = abs(r - r0) / r0
        n_total += 1
        max_frac = max(max_frac, frac)
        if frac > BOND_LEN_TOL_FRAC:
            n_viol += 1
    return {"n_bonds": n_total, "n_violations": n_viol, "max_rel_error": round(max_frac, 4)}


def angle_violations(rm: Chem.Mol, coords: np.ndarray) -> dict:
    n_viol = n_total = 0
    max_dev = 0.0
    for atom in rm.GetAtoms():
        c = atom.GetIdx()
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(nbrs) < 2:
            continue
        theta0 = ref_angle_deg(rm, c)
        pc = coords[c]
        for a in range(len(nbrs)):
            for bb in range(a + 1, len(nbrs)):
                va = coords[nbrs[a]] - pc
                vb = coords[nbrs[bb]] - pc
                na, nb = np.linalg.norm(va), np.linalg.norm(vb)
                if na < 1e-9 or nb < 1e-9:
                    continue
                cos_t = np.clip(np.dot(va, vb) / (na * nb), -1.0, 1.0)
                theta = np.degrees(np.arccos(cos_t))
                dev = abs(theta - theta0)
                n_total += 1
                max_dev = max(max_dev, dev)
                if dev > BOND_ANGLE_TOL_DEG:
                    n_viol += 1
    return {"n_angles": n_total, "n_violations": n_viol, "max_dev_deg": round(max_dev, 2)}


def ring_planarity(rm: Chem.Mol, coords: np.ndarray) -> list[float]:
    """RMS deviation (Å) of each aromatic ring's atoms from its best-fit plane."""
    out = []
    ri = rm.GetRingInfo()
    for ring in ri.AtomRings():
        if not all(rm.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            continue
        pts = coords[list(ring)]
        centroid = pts.mean(axis=0)
        _, _, vt = np.linalg.svd(pts - centroid)
        normal = vt[-1]
        dev = (pts - centroid) @ normal
        out.append(float(np.sqrt(np.mean(dev**2))))
    return out


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    h = p.T @ q
    u, s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    corr = np.diag([1, 1, d])
    r = vt.T @ corr @ u.T
    q_rot = (r @ q.T).T
    return float(np.sqrt(np.mean(np.sum((p - q_rot) ** 2, axis=1))))


def pairwise_rmsd_stats(conf_list: list[np.ndarray]) -> dict:
    n = len(conf_list)
    if n < 2:
        return {"n_conformers": n, "mean_pairwise_rmsd": 0.0, "min_pairwise_rmsd": 0.0, "max_pairwise_rmsd": 0.0}
    vals = [
        kabsch_rmsd(conf_list[i], conf_list[j])
        for i in range(n) for j in range(i + 1, n)
    ]
    return {
        "n_conformers": n,
        "mean_pairwise_rmsd": round(float(np.mean(vals)), 3),
        "min_pairwise_rmsd": round(float(np.min(vals)), 3),
        "max_pairwise_rmsd": round(float(np.max(vals)), 3),
    }


def duplicate_rate(conf_list: list[np.ndarray], threshold: float = 0.2) -> float:
    n = len(conf_list)
    if n < 2:
        return 0.0
    dup = 0
    for i in range(n):
        for j in range(i + 1, n):
            if kabsch_rmsd(conf_list[i], conf_list[j]) < threshold:
                dup += 1
    return round(dup / (n * (n - 1) / 2), 3)


# ---------------------------------------------------------------------------
# Per-molecule pipeline
# ---------------------------------------------------------------------------

def evaluate_one(name: str, smiles: str, category: str) -> dict:
    rec = {"name": name, "smiles": smiles, "category": category, "status": None}

    rm, rmH, rdkit_status, t_embed_rd, t_opt_rd, mmff_pre, mmff_post = rdkit_embed(smiles)
    rec["rdkit_status"] = rdkit_status
    rec["rdkit_embed_time_s"] = round(t_embed_rd, 4) if t_embed_rd is not None else None
    rec["rdkit_mmff_energy_preopt"] = mmff_pre
    rec["rdkit_mmff_energy_postopt"] = mmff_post

    cm, coords_c, chematic_status, t_embed_c = chematic_embed(smiles)
    rec["chematic_status"] = chematic_status
    rec["chematic_embed_time_s"] = round(t_embed_c, 4) if t_embed_c is not None else None

    if rdkit_status != "ok" or chematic_status != "ok":
        rec["status"] = f"rdkit={rdkit_status},chematic={chematic_status}"
        return rec

    n_heavy = rm.GetNumAtoms()
    if cm.heavy_atoms != n_heavy:
        rec["status"] = "atom_count_mismatch"
        return rec

    # Verify atom correspondence (fail loud, per-index element check) before
    # trusting any downstream index-aligned comparison.
    chematic_syms = [row[0] for row in cm.atom_table]
    rdkit_syms = [a.GetSymbol() for a in rm.GetAtoms()]
    if chematic_syms != rdkit_syms:
        rec["status"] = "atom_correspondence_mismatch"
        rec["mismatch_detail"] = {
            "chematic": chematic_syms, "rdkit": rdkit_syms,
        }
        return rec

    coords_c = np.array(coords_c)
    coords_rd = rdkit_heavy_coords(rmH, n_heavy)

    # --- bond length / angle violations, same external reference table -----
    rec["chematic_bond_violations"] = bond_violations(rm, coords_c)
    rec["rdkit_bond_violations"] = bond_violations(rm, coords_rd)
    rec["chematic_angle_violations"] = angle_violations(rm, coords_c)
    rec["rdkit_angle_violations"] = angle_violations(rm, coords_rd)

    # Chirality is computed on whatever geometry exists, blown-up or not: a
    # signed-volume/dihedral-sign read is defined even on a distorted structure,
    # and "does chirality survive" is itself part of what this diagnosis must
    # report (advisor guidance: don't let coverage~0 hide behind other checks).
    # It is tagged with geometry_clean so aggregation can report the clean-only
    # subset separately -- a signed volume computed on a torn molecule conflates
    # "no stereo enforcement" with "garbage geometry" and must not be presented
    # as evidence for the former alone.
    rec["chirality"] = evaluate_chirality(rm, rmH, cm, coords_c, n_heavy)

    # A bond off by >50% of its covalent-radius reference is a torn/degenerate
    # geometry, not force-field imprecision (real relaxation artifacts are a
    # few percent, even strained rings rarely exceed ~20%). Named bucket, not
    # a silent drop: everything computed so far (incl. chirality) is kept in
    # the row, but this record is excluded from the "ok" RMSD/energy/ring
    # aggregates below since a torn molecule makes those numbers meaningless
    # (see docs/etkdg_3d_gap_rfc.md).
    geometry_clean = (
        rec["chematic_bond_violations"]["max_rel_error"] <= BOND_BLOWUP_REL_ERROR
        and rec["rdkit_bond_violations"]["max_rel_error"] <= BOND_BLOWUP_REL_ERROR
    )
    rec["chirality"]["geometry_clean"] = geometry_clean
    if rec["chematic_bond_violations"]["max_rel_error"] > BOND_BLOWUP_REL_ERROR:
        rec["status"] = "chematic_bond_length_blowup"
        return rec
    if rec["rdkit_bond_violations"]["max_rel_error"] > BOND_BLOWUP_REL_ERROR:
        rec["status"] = "rdkit_bond_length_blowup"
        return rec

    # --- ring planarity ------------------------------------------------------
    rec["chematic_ring_planarity_rms"] = ring_planarity(rm, coords_c)
    rec["rdkit_ring_planarity_rms"] = ring_planarity(rm, coords_rd)

    # --- RMSD vs RDKit reference (Kabsch-aligned, heavy atoms only) ---------
    rec["rmsd_vs_rdkit"] = round(kabsch_rmsd(coords_c, coords_rd), 3)

    # --- fair same-yardstick MMFF94 energy delta ----------------------------
    try:
        eb_chematic = cm.mmff94_energy_breakdown(coords_c.tolist())
        eb_rdkit_geom = cm.mmff94_energy_breakdown(coords_rd.tolist())
        rec["chematic_mmff94_own_energy"] = eb_chematic
        rec["chematic_mmff94_on_rdkit_geometry"] = eb_rdkit_geom
        rec["fair_energy_delta"] = round(eb_chematic["total"] - eb_rdkit_geom["total"], 3)
    except Exception as e:
        rec["fair_energy_error"] = f"{type(e).__name__}:{e}"

    rec["status"] = "ok"
    return rec


def evaluate_chirality(rm, rmH, cm, coords_c, n_heavy) -> dict:
    declared_rd = rdkit_declared_stereo(rm)
    perceived_rd = rdkit_perceived_stereo_from_3d(rmH, n_heavy)

    declared_ch = {d["atom_idx"]: d["descriptor"] for d in cm.cip_stereo()}
    perceived_ch = {d["atom_idx"]: d["code"] for d in cm.stereo_from_coords(coords_c.tolist())}

    out = {"centers": []}
    for idx, code in declared_rd.items():
        heavy_degree = rm.GetAtomWithIdx(idx).GetDegree() if code in ("R", "S") else None
        entry = {
            "atom_idx": idx,
            "declared": code,
            "kind": "tetrahedral" if code in ("R", "S") else "alkene",
            "heavy_neighbor_count": heavy_degree,
            "rdkit_covered": idx in perceived_rd,
            "rdkit_match": perceived_rd.get(idx) == code if idx in perceived_rd else None,
            "chematic_declared": declared_ch.get(idx),
            "chematic_covered": idx in perceived_ch,
            "chematic_match": (
                perceived_ch.get(idx) == declared_ch.get(idx))
                if idx in perceived_ch and idx in declared_ch else None,
        }
        out["centers"].append(entry)
    return out


# ---------------------------------------------------------------------------
# Ensemble-generation subset
# ---------------------------------------------------------------------------

def evaluate_ensemble(name: str, smiles: str, n_requested: int = 20) -> dict:
    rec = {"name": name, "n_requested": n_requested}

    # chematic
    try:
        cm = chematic.from_smiles(smiles)
        t0 = time.perf_counter()
        confs = cm.conformer_ensemble(n_requested, 0.3, "mmff94", 30.0)
        rec["chematic_time_s"] = round(time.perf_counter() - t0, 3)
        rec["chematic_n_kept"] = len(confs)
        arrs = [np.array(c) for c in confs]
        rec["chematic_diversity"] = pairwise_rmsd_stats(arrs)
        rec["chematic_duplicate_rate_within_kept"] = duplicate_rate(arrs)
    except Exception as e:
        rec["chematic_error"] = f"{type(e).__name__}:{e}"

    # RDKit
    try:
        rm = Chem.MolFromSmiles(smiles)
        rmH = Chem.AddHs(rm)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        params.pruneRmsThresh = 0.3
        t0 = time.perf_counter()
        cids = AllChem.EmbedMultipleConfs(rmH, numConfs=n_requested, params=params)
        for cid in cids:
            try:
                AllChem.MMFFOptimizeMolecule(rmH, confId=cid)
            except Exception:
                pass
        rec["rdkit_time_s"] = round(time.perf_counter() - t0, 3)
        rec["rdkit_n_kept"] = len(cids)
        n_heavy = rm.GetNumAtoms()
        arrs = [
            np.array([list(rmH.GetConformer(cid).GetAtomPosition(i)) for i in range(n_heavy)])
            for cid in cids
        ]
        rec["rdkit_diversity"] = pairwise_rmsd_stats(arrs)
        rec["rdkit_duplicate_rate_within_kept"] = duplicate_rate(arrs)
    except Exception as e:
        rec["rdkit_error"] = f"{type(e).__name__}:{e}"

    return rec


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(records: list[dict], ensemble_records: list[dict]) -> dict:
    n = len(records)
    ok = [r for r in records if r["status"] == "ok"]
    status_counts: dict[str, int] = {}
    for r in records:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1

    rdkit_embed_ok = sum(1 for r in records if r["rdkit_status"] == "ok")
    chematic_embed_ok = sum(1 for r in records if r["chematic_status"] == "ok")

    rmsds = [r["rmsd_vs_rdkit"] for r in ok if "rmsd_vs_rdkit" in r]
    fair_deltas = [r["fair_energy_delta"] for r in ok if "fair_energy_delta" in r]

    # Runtime is measured on every molecule that reached an embed attempt,
    # independent of the ok/blowup split above (a blown-up geometry still
    # cost real wall-clock time to generate).
    c_times = [r["chematic_embed_time_s"] for r in records if r.get("chematic_embed_time_s") is not None]
    r_times = [r["rdkit_embed_time_s"] for r in records if r.get("rdkit_embed_time_s") is not None]

    chematic_bond_viol = sum(r["chematic_bond_violations"]["n_violations"] for r in ok if "chematic_bond_violations" in r)
    chematic_bond_total = sum(r["chematic_bond_violations"]["n_bonds"] for r in ok if "chematic_bond_violations" in r)
    rdkit_bond_viol = sum(r["rdkit_bond_violations"]["n_violations"] for r in ok if "rdkit_bond_violations" in r)
    rdkit_bond_total = sum(r["rdkit_bond_violations"]["n_bonds"] for r in ok if "rdkit_bond_violations" in r)

    chematic_angle_viol = sum(r["chematic_angle_violations"]["n_violations"] for r in ok if "chematic_angle_violations" in r)
    chematic_angle_total = sum(r["chematic_angle_violations"]["n_angles"] for r in ok if "chematic_angle_violations" in r)
    rdkit_angle_viol = sum(r["rdkit_angle_violations"]["n_violations"] for r in ok if "rdkit_angle_violations" in r)
    rdkit_angle_total = sum(r["rdkit_angle_violations"]["n_angles"] for r in ok if "rdkit_angle_violations" in r)

    # chirality roll-up, split tetrahedral (R/S) vs alkene (E/Z), and split
    # tetrahedral by heavy-neighbor-count (4 = assessable by chematic's current
    # stereo3d code; <4 = has an implicit H, currently invisible to it).
    # Rolled up over every record that reached chirality evaluation (both
    # geometrically clean "ok" rows AND "chematic_bond_length_blowup" rows) --
    # NOT restricted to `ok` -- so coverage figures reflect the whole corpus,
    # not just the subset lucky enough to avoid Phase-0's blow-up bug. A
    # SEPARATE clean-only rollup is also computed, since a signed-volume/
    # dihedral-sign match rate computed on a torn molecule conflates "no
    # stereo enforcement" with "garbage geometry" -- see docs/etkdg_3d_gap_rfc.md.
    def new_bucket():
        return {"declared": 0, "chematic_covered": 0, "chematic_match": 0, "rdkit_covered": 0, "rdkit_match": 0}

    buckets_all = {"tet_4": new_bucket(), "tet_lt4": new_bucket(), "ez": new_bucket()}
    buckets_clean = {"tet_4": new_bucket(), "tet_lt4": new_bucket(), "ez": new_bucket()}
    for r in records:
        chir = r.get("chirality", {})
        clean = chir.get("geometry_clean", False)
        for c in chir.get("centers", []):
            key = "ez" if c["kind"] == "alkene" else ("tet_4" if c["heavy_neighbor_count"] == 4 else "tet_lt4")
            for buckets in ([buckets_all, buckets_clean] if clean else [buckets_all]):
                bucket = buckets[key]
                bucket["declared"] += 1
                if c["chematic_covered"]:
                    bucket["chematic_covered"] += 1
                    if c["chematic_match"]:
                        bucket["chematic_match"] += 1
                if c["rdkit_covered"]:
                    bucket["rdkit_covered"] += 1
                    if c["rdkit_match"]:
                        bucket["rdkit_match"] += 1

    def rate(num, den):
        return round(num / den, 3) if den else None

    def chirality_block(tet_4, tet_lt4, ez):
        return {
            "tetrahedral_4_heavy_neighbors_assessable": {
                "n_declared": tet_4["declared"],
                "chematic_coverage": rate(tet_4["chematic_covered"], tet_4["declared"]),
                "chematic_match_given_covered": rate(tet_4["chematic_match"], tet_4["chematic_covered"]),
                "rdkit_coverage": rate(tet_4["rdkit_covered"], tet_4["declared"]),
                "rdkit_match_given_covered": rate(tet_4["rdkit_match"], tet_4["rdkit_covered"]),
            },
            "tetrahedral_implicit_h_lt4_heavy_neighbors": {
                "n_declared": tet_lt4["declared"],
                "chematic_coverage": rate(tet_lt4["chematic_covered"], tet_lt4["declared"]),
                "chematic_match_given_covered": rate(tet_lt4["chematic_match"], tet_lt4["chematic_covered"]),
                "rdkit_coverage": rate(tet_lt4["rdkit_covered"], tet_lt4["declared"]),
                "rdkit_match_given_covered": rate(tet_lt4["rdkit_match"], tet_lt4["rdkit_covered"]),
            },
            "alkene_ez": {
                "n_declared": ez["declared"],
                "chematic_coverage": rate(ez["chematic_covered"], ez["declared"]),
                "chematic_match_given_covered": rate(ez["chematic_match"], ez["chematic_covered"]),
                "rdkit_coverage": rate(ez["rdkit_covered"], ez["declared"]),
                "rdkit_match_given_covered": rate(ez["rdkit_match"], ez["rdkit_covered"]),
            },
        }

    return {
        "n_molecules": n,
        "status_counts": status_counts,
        "raw_embed_returned_rate": {
            "rdkit": rate(rdkit_embed_ok, n),
            "chematic": rate(chematic_embed_ok, n),
            "note": "fraction where the embedder returned finite, non-clashing "
                    "coordinates at all -- does NOT check geometric validity "
                    "(bond lengths etc). See geometrically_valid_rate for that.",
        },
        "geometrically_valid_rate": {
            "chematic": rate(status_counts.get("ok", 0), n),
            "note": "fraction of the whole corpus that produced a geometry passing "
                    "the >50%-of-covalent-reference bond-length check (status == "
                    "'ok'). This is the metric that should be called 'success' -- "
                    "raw_embed_returned_rate alone is misleadingly high because "
                    "chematic's coordinate generator never refuses to return "
                    "coordinates, it just sometimes returns a torn structure.",
        },
        "n_comparable_ok": len(ok),
        "single_conformer_embed_time_s": {
            "chematic_mean": round(float(np.mean(c_times)), 4) if c_times else None,
            "chematic_median": round(float(np.median(c_times)), 4) if c_times else None,
            "rdkit_mean": round(float(np.mean(r_times)), 4) if r_times else None,
            "rdkit_median": round(float(np.median(r_times)), 4) if r_times else None,
        },
        "rmsd_vs_rdkit": {
            "mean": round(float(np.mean(rmsds)), 3) if rmsds else None,
            "median": round(float(np.median(rmsds)), 3) if rmsds else None,
            "max": round(float(np.max(rmsds)), 3) if rmsds else None,
            "n": len(rmsds),
        },
        "fair_mmff94_energy_delta_chematic_minus_rdkit_geometry": {
            "mean": round(float(np.mean(fair_deltas)), 3) if fair_deltas else None,
            "median": round(float(np.median(fair_deltas)), 3) if fair_deltas else None,
            "n": len(fair_deltas),
            "note": "positive = chematic's own conformer scores WORSE than RDKit's "
                    "geometry under chematic's own full MMFF94 breakdown (same "
                    "topology, same energy function -- an apples-to-apples "
                    "geometry-quality signal, not an absolute energy comparison).",
        },
        "bond_length_violation_rate_ext_reference": {
            "chematic": rate(chematic_bond_viol, chematic_bond_total),
            "rdkit": rate(rdkit_bond_viol, rdkit_bond_total),
        },
        "bond_angle_violation_rate_ext_reference": {
            "chematic": rate(chematic_angle_viol, chematic_angle_total),
            "rdkit": rate(rdkit_angle_viol, rdkit_angle_total),
        },
        "chirality_all_geometry": {
            **chirality_block(buckets_all["tet_4"], buckets_all["tet_lt4"], buckets_all["ez"]),
            "note": "Rolled up over EVERY declared stereocenter reached, including "
                    "torn (chematic_bond_length_blowup) geometry. A signed-volume or "
                    "dihedral-sign read on a torn molecule conflates 'no stereo "
                    "enforcement' with 'garbage geometry' -- use chirality_clean_geometry "
                    "for the unconfounded match-rate signal; this block is provided for "
                    "coverage (which is a structural fact independent of geometry "
                    "quality) and for transparency, not for the match-rate column.",
        },
        "chirality_clean_geometry": {
            **chirality_block(buckets_clean["tet_4"], buckets_clean["tet_lt4"], buckets_clean["ez"]),
            "note": "Same rollup restricted to status=='ok' rows only (no bond >50% "
                    "off the covalent-radius reference). Small n by construction -- "
                    "Phase 0's blow-up bug (see docs/etkdg_3d_gap_rfc.md) removed most "
                    "of the tetrahedral-4-heavy and alkene rows from this subset -- but "
                    "this is the only match-rate signal not confounded by torn geometry.",
        },
        "ensemble_generation": ensemble_records,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="validation/results/etkdg_vs_rdkit_summary.json")
    ap.add_argument("--rows-out", default="validation/results/etkdg_vs_rdkit_rows.jsonl")
    args = ap.parse_args()

    records = []
    for name, smiles, category in CORPUS:
        try:
            rec = evaluate_one(name, smiles, category)
        except Exception as e:  # never silently drop a molecule
            rec = {
                "name": name, "smiles": smiles, "category": category,
                "status": f"harness_exception:{type(e).__name__}:{e}",
            }
        records.append(rec)
        print(f"{rec['status']:<45} {name}")

    ensemble_records = []
    for name in ENSEMBLE_SUBSET:
        smiles = next(s for n, s, _ in CORPUS if n == name)
        ensemble_records.append(evaluate_ensemble(name, smiles))

    summary = aggregate(records, ensemble_records)
    summary["meta"] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python": sys.version,
        "platform": platform.platform(),
        "rdkit_version": __import__("rdkit").__version__,
        "chematic_version": getattr(chematic, "__version__", "unknown"),
        "reproducibility_note": (
            "chematic's Prng (crates/chematic-3d/src/prng.rs) has no public seed "
            "API; it self-seeds from a process-global atomic counter, so ensemble "
            "diversity/duplicate-rate figures are not bit-reproducible run-to-run "
            "(the noise_sigma_deg=0.0 single-conformer path used for all "
            "non-ensemble metrics above IS deterministic). RDKit runs pinned with "
            "randomSeed=0xF00D."
        ),
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(summary, indent=2))

    rows_path = Path(args.rows_out)
    with rows_path.open("w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")

    print("\n--- summary ---")
    print(json.dumps(summary, indent=2)[:4000])
    print(f"\nfull summary written to {out_path}")
    print(f"per-molecule rows written to {rows_path}")


if __name__ == "__main__":
    main()
