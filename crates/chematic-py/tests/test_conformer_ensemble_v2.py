"""Regression tests for `Mol.conformer_ensemble_v2()` (A2.1) -- the Python
binding for `chematic_3d::embed_ensemble_v2` (A2, PR #373).

`embed_ensemble_v2` calls `embed_pipeline_v2` under the hood, which is
already proven (see `test_pipeline_v2.py`, issue #172) to generate directly
on the caller's own atom order without reparsing. So this file needs only
a light confirming sweep for atom-index consistency, not the full 6-fixture
regression suite `test_conformer_ensemble.py` needed to catch that bug in
the first place. The focus here is A2.1-specific: determinism, full
attempt provenance (kept / duplicate-pruned / failed), the
conformer<->provenance cross-reference, and the config-validation
error path.
"""
import math

import chematic


MIN_BOND_LEN = 0.6
MAX_BOND_LEN = 2.6
MIN_INTERATOMIC_DIST = 0.4

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
BRANCHED = "CCC(C)C"
DISCONNECTED = "CCO.CC(C)C"
STEREO = "F[C@@](Cl)(Br)I"


def _assert_consistent_indexing(mol, conformer):
    assert len(conformer) == mol.heavy_atoms
    for a1, a2, _btype, _arom in mol.bond_table:
        d = math.dist(conformer[a1], conformer[a2])
        assert MIN_BOND_LEN <= d <= MAX_BOND_LEN, (
            f"bond {a1}-{a2} has length {d:.3f} -- outside sane range "
            f"[{MIN_BOND_LEN}, {MAX_BOND_LEN}]; likely atom-index mismatch"
        )
    n = len(conformer)
    for i in range(n):
        for j in range(i + 1, n):
            d = math.dist(conformer[i], conformer[j])
            assert d >= MIN_INTERATOMIC_DIST, (
                f"atoms {i},{j} are {d:.3f} A apart -- degenerate/clashing geometry"
            )


def _config(count=4, base_seed=1, rmsd_threshold=0.5, max_attempts=1, **overrides):
    per_conformer = chematic.PipelineV2Config.safe(
        force_field="dreiding",
        stereo_policy="verify_only",
        ring_torsion_policy="diagnostic_only",
        max_attempts=max_attempts,
    )
    return chematic.EnsembleV2Config(
        per_conformer,
        count=count,
        base_seed=base_seed,
        rmsd_threshold=rmsd_threshold,
        **overrides,
    )


# ---------------------------------------------------------------------------
# Atom-index consistency (light sweep -- see module doc for why this is
# lighter than test_conformer_ensemble.py's own 6-fixture suite)
# ---------------------------------------------------------------------------


def test_atom_index_consistency_across_structural_classes():
    for smiles in (ASPIRIN, BRANCHED, DISCONNECTED, STEREO):
        mol = chematic.from_smiles(smiles)
        result = mol.conformer_ensemble_v2(_config(count=3, base_seed=7))
        assert len(result["conformers"]) >= 1, smiles
        for conformer in result["conformers"]:
            _assert_consistent_indexing(mol, conformer)


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


def test_same_base_seed_is_deterministic():
    mol = chematic.from_smiles(ASPIRIN)
    config = _config(count=6, base_seed=42)
    r1 = mol.conformer_ensemble_v2(config)
    r2 = mol.conformer_ensemble_v2(config)
    assert r1["conformers"] == r2["conformers"]
    assert r1["conformer_provenance"] == r2["conformer_provenance"]


def test_different_base_seeds_are_not_aliased():
    mol = chematic.from_smiles(ASPIRIN)
    r1 = mol.conformer_ensemble_v2(_config(count=1, base_seed=1))
    r2 = mol.conformer_ensemble_v2(_config(count=1, base_seed=2))
    assert len(r1["conformers"]) == 1
    assert len(r2["conformers"]) == 1
    assert r1["conformers"][0] != r2["conformers"][0]


# ---------------------------------------------------------------------------
# Full provenance shape
# ---------------------------------------------------------------------------


def test_every_attempt_recorded_with_exactly_one_outcome_populated():
    mol = chematic.from_smiles(ASPIRIN)
    result = mol.conformer_ensemble_v2(_config(count=5, base_seed=100))
    assert len(result["attempts"]) == 5
    assert result["requested_count"] == 5
    assert result["termination"] == "completed"
    for i, attempt in enumerate(result["attempts"]):
        assert attempt["attempt_index"] == i
        assert attempt["outcome"] in ("success", "failure")
        if attempt["outcome"] == "success":
            assert attempt["success"] is not None
            assert attempt["failure"] is None
            assert attempt["success"]["disposition"]["kind"] in (
                "kept",
                "pruned_as_duplicate",
            )
        else:
            assert attempt["success"] is None
            assert attempt["failure"] is not None


def test_conformer_provenance_cross_references_kept_attempts():
    mol = chematic.from_smiles(ASPIRIN)
    result = mol.conformer_ensemble_v2(_config(count=8, base_seed=20260823))
    provenance = result["conformer_provenance"]
    assert len(provenance) == len(result["conformers"])

    kept_attempt_indices = {
        a["attempt_index"]
        for a in result["attempts"]
        if a["outcome"] == "success" and a["success"]["disposition"]["kind"] == "kept"
    }
    for entry in provenance:
        assert entry["attempt_index"] in kept_attempt_indices
        matching = result["attempts"][entry["attempt_index"]]
        assert matching["outcome"] == "success"
        assert matching["success"]["disposition"]["kind"] == "kept"
        assert entry["seed"] == matching["seed"]
        assert entry["energy"] == matching["success"]["energy"]


def test_pruned_duplicate_points_back_at_a_kept_representative():
    """Propane's heavy-atom skeleton (C-C-C) has no dihedral degree of
    freedom, so a generous rmsd_threshold reliably collapses every attempt
    to one kept representative -- exercising the real
    PrunedAsDuplicate/Kept round trip end to end, not just its shape."""
    mol = chematic.from_smiles("CCC")
    result = mol.conformer_ensemble_v2(
        _config(count=6, base_seed=99, rmsd_threshold=5.0)
    )
    kept = [
        a
        for a in result["attempts"]
        if a["outcome"] == "success" and a["success"]["disposition"]["kind"] == "kept"
    ]
    pruned = [
        a
        for a in result["attempts"]
        if a["outcome"] == "success"
        and a["success"]["disposition"]["kind"] == "pruned_as_duplicate"
    ]
    assert len(kept) == 1, f"expected exactly one representative, got {len(kept)}"
    assert pruned, "expected at least one pruned duplicate at this generous threshold"
    kept_index = kept[0]["attempt_index"]
    for a in pruned:
        disposition = a["success"]["disposition"]
        assert disposition["representative_attempt_index"] == kept_index
        assert 0.0 <= disposition["rmsd"] < 5.0
        assert disposition["symmetric"] is True


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------


def test_invalid_rmsd_threshold_raises_value_error_at_call_time():
    mol = chematic.from_smiles("CCC")
    per_conformer = chematic.PipelineV2Config.safe(
        force_field="dreiding",
        stereo_policy="verify_only",
        ring_torsion_policy="diagnostic_only",
        max_attempts=1,
    )
    for bad in (-1.0, float("nan"), float("inf")):
        # Construction itself is infallible -- matches the Rust
        # EnsembleV2Config struct, which has no invariant of its own.
        config = chematic.EnsembleV2Config(per_conformer, count=1, base_seed=1, rmsd_threshold=bad)
        try:
            mol.conformer_ensemble_v2(config)
            raise AssertionError(f"rmsd_threshold={bad} should have raised ValueError")
        except ValueError:
            pass


def test_zero_rmsd_threshold_is_accepted_as_pruning_disabled():
    mol = chematic.from_smiles("CCC")
    result = mol.conformer_ensemble_v2(_config(count=2, base_seed=1, rmsd_threshold=0.0))
    assert len(result["attempts"]) == 2


# ---------------------------------------------------------------------------
# Legacy API: deprecation note added, behavior unchanged
# ---------------------------------------------------------------------------


def test_conformer_ensemble_docstring_has_deprecation_note():
    doc = chematic.Mol.conformer_ensemble.__doc__ or ""
    assert "deprecated" in doc.lower()
    assert "conformer_ensemble_v2" in doc


def test_conformer_ensemble_behavior_unchanged():
    """Sanity check that the docstring-only edit didn't touch behavior."""
    mol = chematic.from_smiles(ASPIRIN)
    ensemble = mol.conformer_ensemble(3, 0.5, "dreiding", 30.0)
    assert isinstance(ensemble, list)
    assert len(ensemble) >= 1
    for conformer in ensemble:
        _assert_consistent_indexing(mol, conformer)
