"""Canonical SMILES differential validation (chematic vs RDKit).

Canonical SMILES are algorithm-specific, so chematic and RDKit produce
different *strings*. What must hold is **semantic round-trip equivalence**:
chematic's canonical SMILES, re-parsed by RDKit, canonicalizes to the same
molecule as RDKit's native canonicalization of the original input.

Full-corpus measurement lives in ``scripts/canonical_diff.py`` (5,000-mol run;
see docs/rdkit_compat.md for the current agreement rate). This file is a
portable regression guard on a small curated corpus; the RDKit-dependent
parts auto-skip when RDKit is absent.

Formerly-known divergence, now fixed: chematic used to drop directional
(`/`,`\\`) E/Z stereo on **exocyclic C=N** bonds adjacent to an aromatic ring
atom, because the bond's order had to be forced to `Aromatic` (for SMARTS
`:a` matching) with no room left to also record direction. The direction is
now stashed on the side (`bond_directions`, see chematic-core) and consulted
independently by the canonical writer.
"""
import pytest

import chematic
from chematic import rdkit_compat as Chem


# Diverse, non-exocyclic-C=N corpus — every entry must round-trip and be idempotent.
CLEAN_CORPUS = [
    "CCO", "c1ccccc1", "c1ccncc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "c1ccc2ccccc2c1", "OCC(O)CO",
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "C1CCCCC1", "CCN(CC)CC", "c1ccc(O)cc1",
    "C[C@H](N)C(=O)O", "ClC(Cl)Cl", "Nc1ccccc1", "CC#N", "CCOCC",
    "c1ccoc1", "O=C(O)c1ccccc1",
    # E/Z stable skeletons — direction choice is deterministic & idempotent
    "C/C=C/C", "C/C=C\\C", "F/C=C/F", "F/C=C\\F", "CC/C=C/CC", "CC/C=C\\CC",
    "C/C=C/C=C/C", "Cl/C=C/Br", "C/C=C/c1ccccc1",
    "c1ccc2cc3ccccc3cc2c1",        # anthracene
    "c1ccc2ncccc2c1",              # quinoline
    "c1ccc2c(c1)cc[nH]2",          # indole
    "c1cc2ccc3cccc4ccc(c1)c2c34",  # pyrene
    "O=C1CCC(=O)N1", "c1ccsc1", "CN1CCC[C@H]1c1cccnc1",  # nicotine
]

# Exocyclic C=N E/Z adjacent to an aromatic ring atom — formerly a known
# divergence (direction dropped because the ring bond order must stay
# Aromatic); now preserved via the bond_directions side channel (see
# test_exocyclic_cn_stereo_preserved below). Kept out of CLEAN_CORPUS: this
# molecule also happens to hit the separate, pre-existing aromaticity-
# perception round-trip idempotency issue (docs/rdkit_compat.md) — its
# canonical form is not idempotent even with connectivity alone (`/`,`\\`
# stripped), which is unrelated to stereo preservation.
EZ_EXOCYCLIC_FIXED = "CCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O"


# ---------------------------------------------------------------------------
# Idempotency — chematic-internal, no RDKit needed
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smi", CLEAN_CORPUS)
def test_canonical_idempotent(smi):
    c1 = chematic.from_smiles(smi).smiles
    c2 = chematic.from_smiles(c1).smiles
    assert c1 == c2, f"canonical SMILES not idempotent: {c1!r} -> {c2!r}"


# ---------------------------------------------------------------------------
# Aromaticity round-trip consistency (chematic-internal; Sprint 9 diagnostic)
#
# The residual ~1.6% canonical idempotency failures on large fused polycyclics
# are caused by aromaticity perception disagreeing between a molecule and the
# re-parse of its own canonical SMILES (e.g. 16 vs 17 aromatic bonds), which
# shifts Morgan ranks and the emitted atom order. This guards the property on
# fused aromatics that DO round-trip consistently; it can be extended once the
# aromaticity/parser-core fix lands (see docs/rdkit_compat.md).
# ---------------------------------------------------------------------------

AROMATIC_ROUNDTRIP_CORPUS = [
    "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1",
    "c1ccc2ccccc2c1",          # naphthalene
    "c1ccc2ncccc2c1",          # quinoline
    "c1ccc2c(c1)cc[nH]2",      # indole
    "c1ccc2cc3ccccc3cc2c1",    # anthracene
    "c1ccc2[nH]c3ccccc3c2c1",  # carbazole
    "c1ccc2c(c1)oc1ccccc12",   # dibenzofuran
]


def _arom_counts(smi):
    m = chematic.from_smiles(smi)
    return (sum(1 for a in m.atom_table if a[3]),
            sum(1 for b in m.bond_table if b[3]))


@pytest.mark.parametrize("smi", AROMATIC_ROUNDTRIP_CORPUS)
def test_aromaticity_roundtrip_consistent(smi):
    before = _arom_counts(smi)
    canon = chematic.from_smiles(smi).smiles
    after = _arom_counts(canon)
    assert before == after, (
        f"aromatic (atoms, bonds) changed across canonical round-trip for {smi}: "
        f"{before} -> {after}"
    )


# ---------------------------------------------------------------------------
# Round-trip semantic equivalence vs RDKit
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def rdkit_mod():
    rdkit = pytest.importorskip("rdkit")
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from rdkit import Chem as RDChem
    return RDChem


@pytest.mark.parametrize("smi", CLEAN_CORPUS)
def test_chematic_smiles_parseable_by_rdkit(smi, rdkit_mod):
    cm = chematic.from_smiles(smi).smiles
    assert rdkit_mod.MolFromSmiles(cm) is not None, \
        f"chematic SMILES not parseable by RDKit: {cm!r}"


@pytest.mark.parametrize("smi", CLEAN_CORPUS)
def test_roundtrip_equivalence(smi, rdkit_mod):
    cm = chematic.from_smiles(smi).smiles
    rd_native = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(smi))
    rd_of_cm = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(cm))
    assert rd_native == rd_of_cm, \
        f"round-trip mismatch for {smi}: rdkit_native={rd_native!r} rdkit(chematic)={rd_of_cm!r}"


def test_exocyclic_cn_stereo_preserved(rdkit_mod):
    """Exocyclic C=N E/Z direction now survives canonical round-trip.

    Not part of CLEAN_CORPUS (see EZ_EXOCYCLIC_FIXED comment) because this
    molecule separately triggers the pre-existing aromaticity round-trip
    idempotency issue — unrelated to stereo, out of scope here.
    """
    cm = chematic.from_smiles(EZ_EXOCYCLIC_FIXED).smiles
    assert rdkit_mod.MolFromSmiles(cm) is not None, "must be valid SMILES"
    assert "/" in cm or "\\" in cm, "E/Z direction should be written, not dropped"
    rd_native = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(EZ_EXOCYCLIC_FIXED))
    rd_of_cm = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(cm))
    assert rd_native == rd_of_cm, \
        f"round-trip mismatch: rdkit_native={rd_native!r} rdkit(chematic)={rd_of_cm!r}"
