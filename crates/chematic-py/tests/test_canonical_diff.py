"""Canonical SMILES differential validation (chematic vs RDKit).

Canonical SMILES are algorithm-specific, so chematic and RDKit produce
different *strings*. What must hold is **semantic round-trip equivalence**:
chematic's canonical SMILES, re-parsed by RDKit, canonicalizes to the same
molecule as RDKit's native canonicalization of the original input.

Full-corpus measurement lives in ``scripts/canonical_diff.py`` (5,000-mol run:
100% RDKit-parseable, 99.62% round-trip equivalent, the 19 mismatches all from
exocyclic C=N E/Z stereo). This file is a portable regression guard on a small
curated corpus; the RDKit-dependent parts auto-skip when RDKit is absent.

Known divergence (documented, not yet fixed): chematic drops directional
(`/`,`\\`) E/Z stereo on **exocyclic C=N** bonds when writing canonical SMILES.
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
    "c1ccoc1", "O=C(O)c1ccccc1", "C/C=C/C", "F/C=C\\F",
    "c1ccc2cc3ccccc3cc2c1",        # anthracene
    "c1ccc2ncccc2c1",              # quinoline
    "c1ccc2c(c1)cc[nH]2",          # indole
    "c1cc2ccc3cccc4ccc(c1)c2c34",  # pyrene
    "O=C1CCC(=O)N1", "c1ccsc1", "CN1CCC[C@H]1c1cccnc1",  # nicotine
]

# Known E/Z exocyclic-C=N divergence: chematic loses the directional stereo.
EZ_EXOCYCLIC_KNOWN = "CCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O"


# ---------------------------------------------------------------------------
# Idempotency — chematic-internal, no RDKit needed
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smi", CLEAN_CORPUS)
def test_canonical_idempotent(smi):
    c1 = chematic.from_smiles(smi).smiles
    c2 = chematic.from_smiles(c1).smiles
    assert c1 == c2, f"canonical SMILES not idempotent: {c1!r} -> {c2!r}"


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


def test_ez_exocyclic_is_known_divergence(rdkit_mod):
    """Document the one known round-trip divergence class.

    chematic still emits valid, RDKit-parseable SMILES — only the exocyclic
    C=N E/Z stereo is dropped. If chematic ever fixes this, the round-trip
    will match and this test should be promoted into CLEAN_CORPUS.
    """
    cm = chematic.from_smiles(EZ_EXOCYCLIC_KNOWN).smiles
    assert rdkit_mod.MolFromSmiles(cm) is not None, "must still be valid SMILES"
    rd_native = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(EZ_EXOCYCLIC_KNOWN))
    rd_of_cm = rdkit_mod.MolToSmiles(rdkit_mod.MolFromSmiles(cm))
    # Connectivity (ignoring stereo) must still match.
    strip = lambda s: s.replace("/", "").replace("\\", "")
    assert strip(rd_native) == strip(rd_of_cm), \
        "exocyclic-C=N divergence must be stereo-only, not connectivity"
