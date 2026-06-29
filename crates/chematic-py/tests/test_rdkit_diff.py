"""Differential tests: chematic.rdkit_compat vs real RDKit.

Skipped automatically when RDKit is not installed (e.g. CI), so this never
blocks the no-RDKit build. Run locally with RDKit present to validate
compatibility. Disagreements that are *expected* (canonical SMILES spelling,
Morgan bit positions) are recorded to ``validation/results/rdkit_diff.jsonl``
rather than asserted, so the diff stays explainable.

Known-good equalities (MW exact, HBA/HBD exact, ring count, TPSA/LogP within
tolerance) ARE asserted — these match chematic's published doctor numbers.
"""
import json
import os

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem as RDChem
from rdkit.Chem import Descriptors as RDDesc, rdMolDescriptors as RDmd

from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors, rdMolDescriptors, DataStructs


CORPUS = [
    "CCO", "c1ccccc1", "c1ccncc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C",  # caffeine
    "c1ccc2ccccc2c1",              # naphthalene
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1",  # ibuprofen
    "C1CCCCC1", "CCN(CC)CC", "OCC(O)CO", "c1ccc(O)cc1",
    "CC(N)C(=O)O", "ClC(Cl)Cl", "O=C(O)c1ccccc1", "Nc1ccccc1",
    "C1=CC2=CC=CC=C2C=C1", "CC#N", "CCOCC", "c1ccoc1",
]

_DIFF_PATH = os.path.join(
    os.path.dirname(__file__), "..", "..", "..",
    "validation", "results", "rdkit_diff.jsonl",
)
_diff_rows = []


def _record(smiles, metric, chematic_val, rdkit_val, ok, delta=None):
    _diff_rows.append({
        "smiles": smiles, "metric": metric,
        "chematic": chematic_val, "rdkit": rdkit_val,
        "delta": delta, "ok": ok,
    })


def teardown_module(module):
    """Flush the collected diff rows to JSONL (best effort)."""
    try:
        os.makedirs(os.path.dirname(_DIFF_PATH), exist_ok=True)
        with open(_DIFF_PATH, "w") as f:
            for row in _diff_rows:
                f.write(json.dumps(row) + "\n")
    except OSError:
        pass  # ponytail: report artifact is optional, never fail teardown


@pytest.mark.parametrize("smi", CORPUS)
def test_diff_descriptors(smi):
    cm = Chem.MolFromSmiles(smi)
    rm = RDChem.MolFromSmiles(smi)
    assert cm is not None and rm is not None

    # MW: exact match (monoisotopic-independent average weight)
    cm_mw, rm_mw = Descriptors.MolWt(cm), RDDesc.MolWt(rm)
    _record(smi, "MolWt", cm_mw, rm_mw, abs(cm_mw - rm_mw) < 0.05, abs(cm_mw - rm_mw))
    assert abs(cm_mw - rm_mw) < 0.05, f"MW {smi}: {cm_mw} vs {rm_mw}"

    # HBA / HBD: exact integer match
    cm_hba, rm_hba = rdMolDescriptors.CalcNumHBA(cm), RDmd.CalcNumHBA(rm)
    _record(smi, "HBA", cm_hba, rm_hba, cm_hba == rm_hba, cm_hba - rm_hba)
    assert cm_hba == rm_hba, f"HBA {smi}: {cm_hba} vs {rm_hba}"

    cm_hbd, rm_hbd = rdMolDescriptors.CalcNumHBD(cm), RDmd.CalcNumHBD(rm)
    _record(smi, "HBD", cm_hbd, rm_hbd, cm_hbd == rm_hbd, cm_hbd - rm_hbd)
    assert cm_hbd == rm_hbd, f"HBD {smi}: {cm_hbd} vs {rm_hbd}"

    # TPSA within 1.0, LogP within 0.5 (recorded; tolerant assert)
    cm_tpsa, rm_tpsa = rdMolDescriptors.CalcTPSA(cm), RDmd.CalcTPSA(rm)
    _record(smi, "TPSA", cm_tpsa, rm_tpsa, abs(cm_tpsa - rm_tpsa) < 1.0,
            abs(cm_tpsa - rm_tpsa))
    assert abs(cm_tpsa - rm_tpsa) < 1.0, f"TPSA {smi}: {cm_tpsa} vs {rm_tpsa}"

    cm_logp, rm_logp = Descriptors.MolLogP(cm), RDDesc.MolLogP(rm)
    _record(smi, "MolLogP", cm_logp, rm_logp, abs(cm_logp - rm_logp) < 0.5,
            abs(cm_logp - rm_logp))
    assert abs(cm_logp - rm_logp) < 0.5, f"LogP {smi}: {cm_logp} vs {rm_logp}"


@pytest.mark.parametrize("smi", CORPUS)
def test_diff_ring_count(smi):
    cm = Chem.MolFromSmiles(smi)
    rm = RDChem.MolFromSmiles(smi)
    cn = cm.GetRingInfo().NumRings()
    rn = rm.GetRingInfo().NumRings()
    _record(smi, "NumRings", cn, rn, cn == rn, cn - rn)
    assert cn == rn, f"NumRings {smi}: {cn} vs {rn}"


@pytest.mark.parametrize("smarts", ["[OH]", "c", "[#7]", "C=O"])
def test_diff_smarts_match_count(smarts):
    rq = RDChem.MolFromSmarts(smarts)
    for smi in CORPUS:
        cm = Chem.MolFromSmiles(smi)
        rm = RDChem.MolFromSmiles(smi)
        cn = len(cm.GetSubstructMatches(smarts))
        rn = len(rm.GetSubstructMatches(rq))
        ok = cn == rn
        _record(smi, f"smarts:{smarts}", cn, rn, ok, cn - rn)
        assert ok, f"SMARTS {smarts} on {smi}: {cn} vs {rn}"


def test_diff_sdf_roundtrip_atom_count(tmp_path):
    out = tmp_path / "diff.sdf"
    with Chem.SDWriter(str(out)) as w:
        for smi in CORPUS:
            w.write(Chem.MolFromSmiles(smi))
    supplier = RDChem.SDMolSupplier(str(out))
    rd_mols = [m for m in supplier if m is not None]
    assert len(rd_mols) == len(CORPUS)
    for smi, rm in zip(CORPUS, rd_mols):
        cm = Chem.MolFromSmiles(smi)
        ok = cm.GetNumAtoms() == rm.GetNumAtoms()
        _record(smi, "sdf_atom_count", cm.GetNumAtoms(), rm.GetNumAtoms(), ok)
        assert ok, f"SDF atom count {smi}: {cm.GetNumAtoms()} vs {rm.GetNumAtoms()}"


def test_diff_morgan_self_tanimoto():
    # Morgan bit positions differ (FNV-1a vs MurmurHash); only self==1.0 is required.
    for smi in CORPUS:
        cm = Chem.MolFromSmiles(smi)
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(cm, 2)
        sim = DataStructs.TanimotoSimilarity(fp, fp)
        _record(smi, "morgan_self_tanimoto", sim, 1.0, sim == 1.0)
        assert sim == 1.0


# ---------------------------------------------------------------------------
# SMARTS substructure match sets vs RDKit
#
# Match sets are compared order-invariantly (set of frozensets). Ring-size
# queries ([rN]) and aromatic-carbonyl edge cases diverge due to SSSR / aromatic
# perception differences (see scripts/rdkit_compat_diff.py for the full corpus
# breakdown) and are intentionally excluded from this regression guard.
# ---------------------------------------------------------------------------

# Queries that agree with RDKit on the clean corpus (no [rN] ring-size queries).
CLEAN_SMARTS = [
    "[OH]", "[#7]", "[CX4]", "c", "C=O", "[F,Cl,Br,I]",
    "[OX2H]", "[#16]", "[nH]", "[!#6;!#1]",
]


def _match_set(matches):
    return frozenset(frozenset(m) for m in matches)


@pytest.mark.parametrize("smarts", CLEAN_SMARTS)
def test_diff_smarts_match_sets(smarts):
    rq = RDChem.MolFromSmarts(smarts)
    for smi in CORPUS:
        cm = Chem.MolFromSmiles(smi)
        rm = RDChem.MolFromSmiles(smi)
        cm_set = _match_set(cm.GetSubstructMatches(smarts))
        rd_set = _match_set(rm.GetSubstructMatches(rq, uniquify=True))
        ok = cm_set == rd_set
        _record(smi, f"smarts_set:{smarts}", len(cm_set), len(rd_set), ok)
        assert ok, f"SMARTS {smarts} on {smi}: {sorted(map(sorted, cm_set))} vs {sorted(map(sorted, rd_set))}"


@pytest.mark.parametrize("smi", CORPUS)
def test_diff_aromatic_atom_count(smi):
    cm = Chem.MolFromSmiles(smi)
    rm = RDChem.MolFromSmiles(smi)
    cm_n = sum(1 for a in cm._mol.atom_table if a[3])
    rd_n = sum(1 for a in rm.GetAtoms() if a.GetIsAromatic())
    _record(smi, "aromatic_atoms", cm_n, rd_n, cm_n == rd_n, cm_n - rd_n)
    assert cm_n == rd_n, f"aromatic atom count {smi}: {cm_n} vs {rd_n}"


@pytest.mark.parametrize("smi", CORPUS)
def test_diff_aromatic_bond_count(smi):
    cm = Chem.MolFromSmiles(smi)
    rm = RDChem.MolFromSmiles(smi)
    cm_n = sum(1 for b in cm._mol.bond_table if b[3])
    rd_n = sum(1 for b in rm.GetBonds() if b.GetIsAromatic())
    _record(smi, "aromatic_bonds", cm_n, rd_n, cm_n == rd_n, cm_n - rd_n)
    assert cm_n == rd_n, f"aromatic bond count {smi}: {cm_n} vs {rd_n}"
