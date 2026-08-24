"""Tests for Mol.embed_pipeline_v2() / PipelineV2Config / PipelineV2Error.

Atom-order-consistency methodology mirrors test_conformer_ensemble.py (issue
#172): embed_pipeline_v2() must be applied directly to Mol.inner, never a
canonicalize-then-reparse copy, or the returned coords desync from the atom/
bond tables the caller already holds.
"""

import ast
import math
from pathlib import Path

import pytest

import chematic

DECANE = "CCCCCCCCCC"
NAPHTHALENE = "c1ccc2ccccc2c1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
BRANCHED = "CCC(C)C"  # 2-methylbutane
# same molecule (aspirin) as ASPIRIN, written starting from the opposite end --
# a genuine atom-order permutation (verified below to actually reorder), not
# a repeat parse of the same traversal.
ASPIRIN_REORDERED = "OC(=O)c1ccccc1OC(C)=O"

MIN_BOND_LEN = 0.6
MAX_BOND_LEN = 2.6
MIN_INTERATOMIC_DIST = 0.4


def _safe_config(**overrides):
    kwargs = dict(
        force_field="none",
        stereo_policy="ignore",
        ring_torsion_policy="fail_closed",
        embed_seed=7,
    )
    kwargs.update(overrides)
    return chematic.PipelineV2Config.safe(**kwargs)


def _assert_consistent_indexing(mol, coords):
    """Every bond in mol.bond_table must have a plausible length in `coords`,
    and no two atoms may coincide -- same property test_conformer_ensemble.py
    checks for issue #172, applied to embed_pipeline_v2()'s own coords."""
    assert len(coords) == mol.heavy_atoms
    for a1, a2, _btype, _arom in mol.bond_table:
        d = math.dist(coords[a1], coords[a2])
        assert MIN_BOND_LEN <= d <= MAX_BOND_LEN, (
            f"bond {a1}-{a2} has length {d:.3f} -- outside sane range "
            f"[{MIN_BOND_LEN}, {MAX_BOND_LEN}]; likely atom-index mismatch"
        )
    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            d = math.dist(coords[i], coords[j])
            assert d >= MIN_INTERATOMIC_DIST, (
                f"atoms {i},{j} are {d:.3f} A apart -- degenerate/clashing geometry"
            )


def _reorders_under_canonicalization(mol) -> bool:
    reparsed = chematic.from_smiles(mol.smiles)
    orig_bonds = {frozenset((a, b)) for a, b, _, _ in mol.bond_table}
    new_bonds = {frozenset((a, b)) for a, b, _, _ in reparsed.bond_table}
    return orig_bonds != new_bonds


# ---------------------------------------------------------------------------
# Atom-order invariant
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("smiles", [DECANE, NAPHTHALENE, ASPIRIN, BRANCHED])
def test_embed_pipeline_v2_atom_index_consistency(smiles):
    """Returned coords must be indexed identically to the Mol the caller
    already holds, for a plain chain, a fused aromatic ring, an
    asymmetrically-substituted ring, and a branched molecule."""
    mol = chematic.from_smiles(smiles)
    result = mol.embed_pipeline_v2(_safe_config())
    _assert_consistent_indexing(mol, result["coords"])
    assert result["final_validation"]["atom_count_unchanged"] is True
    assert result["final_validation"]["all_finite"] is True
    assert len(result["coords"]) == mol.heavy_atoms


def test_naphthalene_aspirin_actually_reorder_under_canonicalization():
    """Sanity-check the test methodology itself (matches issue #172's own
    check): these fixtures must actually reorder under canonicalization, or
    the tests above wouldn't exercise the atom-order bug's trigger condition.
    (DECANE is deliberately excluded here: a plain path graph's edge set is
    frozenset-identical under end-to-end reversal, so this bond-set-based
    check cannot observe decane's reordering either way -- a pre-existing
    limitation of this same methodology in test_conformer_ensemble.py, not
    something introduced or fixed by this PR.)"""
    for smiles in (NAPHTHALENE, ASPIRIN):
        mol = chematic.from_smiles(smiles)
        assert _reorders_under_canonicalization(mol), (
            f"{smiles}: expected canonical_smiles() to reorder atoms"
        )


def test_negative_control_reproduces_old_bug_shape():
    """Reconstruct the OLD conformer_ensemble()-style bug (issue #172) by
    hand: generate on a canonicalize-reparsed copy, then cross-index against
    the ORIGINAL mol's bond_table, as a caller holding only `mol` would.
    Proves _assert_consistent_indexing can actually detect a real atom-order
    mismatch, not just that it happens to pass."""
    mol = chematic.from_smiles(ASPIRIN)
    assert _reorders_under_canonicalization(mol)  # precondition for the bug to bite
    reparsed = chematic.from_smiles(mol.smiles)
    buggy_result = reparsed.embed_pipeline_v2(_safe_config())
    with pytest.raises(AssertionError):
        _assert_consistent_indexing(mol, buggy_result["coords"])


def test_aspirin_reordered_smiles_is_a_genuine_permutation():
    m1 = chematic.from_smiles(ASPIRIN)
    m2 = chematic.from_smiles(ASPIRIN_REORDERED)
    assert m1.smiles == m2.smiles, "expected the same underlying molecule"
    assert m1.bond_table != m2.bond_table, (
        "expected a genuinely different atom numbering, not the same traversal twice"
    )


def test_permuted_input_gives_equivalent_geometry_and_diagnostics():
    """embed_pipeline_v2() is stochastic distance geometry seeded internally
    by atom-traversal order, so a genuine atom permutation is not expected to
    produce bit-identical coordinates. What must hold instead (the actual
    atom-order invariant): each permutation's own coords stay self-consistent
    with THAT Mol's own bond_table (checked independently for each), and the
    permutation-invariant summary diagnostics (declared-stereocenter count,
    atom count, geometric soundness) agree between the two orderings."""
    config = _safe_config()
    m1 = chematic.from_smiles(ASPIRIN)
    m2 = chematic.from_smiles(ASPIRIN_REORDERED)
    r1 = m1.embed_pipeline_v2(config)
    r2 = m2.embed_pipeline_v2(config)

    _assert_consistent_indexing(m1, r1["coords"])
    _assert_consistent_indexing(m2, r2["coords"])

    assert len(r1["coords"]) == len(r2["coords"]) == m1.heavy_atoms == m2.heavy_atoms
    assert r1["stereo_before"]["n_declared"] == r2["stereo_before"]["n_declared"]
    assert r1["final_validation"]["sound"] == r2["final_validation"]["sound"] is True
    assert (
        r1["final_validation"]["atom_count_unchanged"]
        == r2["final_validation"]["atom_count_unchanged"]
        is True
    )


# ---------------------------------------------------------------------------
# Deterministic reproducibility + field-by-field schema check
#
# Grounds every field's expected shape/value in what the Rust source actually
# guarantees for this exact scenario (ForceFieldPolicy::None never fails or
# falls back; decane has no declared stereo; a fixed seed makes the
# stochastic embedder reproducible) -- these are the closest in-process proxy
# for "Rust result and Python result agree field-by-field" available here:
# chematic-py builds only as a cdylib (no rlib target), so a Rust-side
# #[test] cannot itself embed a Python interpreter to call this same binding.
# ---------------------------------------------------------------------------


def test_result_schema_and_values_for_deterministic_case():
    mol = chematic.from_smiles(DECANE)
    config = _safe_config(embed_seed=42)
    result = mol.embed_pipeline_v2(config)

    assert set(result.keys()) == {
        "coords",
        "embed_stats",
        "bound_adjustment_report",
        "torsion_knowledge_report",
        "ring_torsion_evidence",
        "torsion_optimization_report",
        "stereo_before",
        "stereo_repair",
        "stereo_after_repair",
        "force_field",
        "final_stereo",
        "final_validation",
        "elapsed_ms_by_stage",
    }

    # use_macrocycle_14_bounds=False (the .safe() default) -> None, not [].
    assert result["bound_adjustment_report"] is None
    # stereo_policy="ignore" -> never repaired.
    assert result["stereo_repair"] is None
    # no torsion-knowledge flags set -> nothing to optimize.
    assert result["torsion_optimization_report"] is None

    for key in ("stereo_before", "stereo_after_repair", "final_stereo"):
        sv = result[key]
        assert sv == {
            "tetrahedral": [],
            "double_bond": [],
            "n_declared": 0,
            "n_satisfied": 0,
            "n_violations": 0,
            "n_unevaluable": 0,
            "is_fully_satisfied": True,
        }

    ff = result["force_field"]
    assert ff["requested_force_field"] == "none"
    assert ff["actual_force_field_used"] == "none"
    assert ff["fallback_reason"] is None
    assert ff["missing_parameter_classes"] == []
    assert ff["coverage"] is None
    assert ff["energy_before"] == {"kind": "none", "total": 0.0}
    assert ff["energy_after"] == {"kind": "none", "total": 0.0}
    assert ff["converged"] is True
    assert ff["iterations"] == 0
    assert ff["max_residual_force"] == 0.0
    assert ff["starting_geometry"] is None
    assert ff["coords"] == result["coords"]

    stats = result["embed_stats"]
    assert stats["attempts_used"] >= 1
    assert stats["failure_counts"] == {}

    fv = result["final_validation"]
    assert fv["all_finite"] is True
    assert fv["atom_count_unchanged"] is True
    assert fv["sound"] is True
    assert fv["stereo_ok"] is True
    assert isinstance(fv["bounds_conformance"], dict)
    assert set(fv["bounds_conformance"].keys()) == {
        "n_pairs",
        "n_violations",
        "max_rel_violation",
    }

    timings = result["elapsed_ms_by_stage"]
    assert set(timings.keys()) == {
        "torsion_knowledge_ms",
        "bound_adjustment_ms",
        "distance_geometry_ms",
        "torsion_energy_eval_ms",
        "torsion_optimization_ms",
        "stereo_verify_before_ms",
        "stereo_repair_ms",
        "stereo_verify_after_repair_ms",
        "force_field_ms",
        "final_stereo_verify_ms",
        "final_validation_ms",
        "total_ms",
    }


def test_same_seed_is_reproducible():
    mol = chematic.from_smiles(DECANE)
    config = _safe_config(embed_seed=42)
    r1 = mol.embed_pipeline_v2(config)
    r2 = mol.embed_pipeline_v2(config)
    assert r1["coords"] == r2["coords"]
    assert r1["embed_stats"] == r2["embed_stats"]


# ---------------------------------------------------------------------------
# Failure contract
# ---------------------------------------------------------------------------


CYCLOHEXANE_WITH_CHAIN = "C1CCCCC1CCCCCCCCCCCC"  # saturated small ring + acyclic tail


def test_failure_raises_pipeline_v2_error_with_structured_diagnostics():
    # a small-ring torsion potential requested under the fail-closed ring
    # policy on a molecule with a genuine saturated small ring is a typed,
    # reliably-reachable failure (RingTorsionApplicationUnsupported) -- see
    # pipeline_v2.rs's stage-6 gate. (Naphthalene's ring bonds are aromatic,
    # not classified as small-ring torsion candidates at all, so it would
    # not reach this failure -- confirmed empirically before writing this
    # fixture, not assumed.)
    mol = chematic.from_smiles(CYCLOHEXANE_WITH_CHAIN)
    config = _safe_config(use_small_ring_torsions=True, ring_torsion_policy="fail_closed")

    with pytest.raises(chematic.PipelineV2Error) as excinfo:
        mol.embed_pipeline_v2(config)

    err = excinfo.value
    assert isinstance(err, ValueError)
    diag = err.diagnostics
    assert diag["cause"] == {"kind": "ring_torsion_application_unsupported"}
    assert diag["stage"] == "torsion_optimization"
    assert diag["coords_are_diagnostic_only"] is True
    # last_known_coords is present (diagnostic) but distinct from a real
    # "coords" success key -- must never be mistaken for one.
    assert "coords" not in diag
    assert diag["last_known_coords"] is not None
    assert len(diag["last_known_coords"]) == mol.heavy_atoms


def test_diagnostic_only_ring_torsion_policy_avoids_the_failure():
    """Same request as above but under DiagnosticOnly: succeeds, and
    ring_torsion_evidence truthfully reports the potentials as scored-only
    (never silently upgraded to "applied")."""
    mol = chematic.from_smiles(CYCLOHEXANE_WITH_CHAIN)
    config = _safe_config(
        use_small_ring_torsions=True, ring_torsion_policy="diagnostic_only"
    )
    result = mol.embed_pipeline_v2(config)
    evidence = result["ring_torsion_evidence"]
    assert evidence["diagnostic_only"] is True
    assert evidence["n_applied"] == 0
    assert evidence["n_scored_only"] > 0, (
        "expected at least one small-ring potential to actually be matched "
        "and scored -- otherwise this test can't distinguish DiagnosticOnly "
        "from a request that matched nothing at all"
    )


def test_invalid_policy_string_raises_value_error():
    with pytest.raises(ValueError):
        _safe_config(force_field="not_a_real_policy")
    with pytest.raises(ValueError):
        _safe_config(stereo_policy="not_a_real_policy")
    with pytest.raises(ValueError):
        _safe_config(ring_torsion_policy="not_a_real_policy")


# ---------------------------------------------------------------------------
# PipelineV2Config construction
# ---------------------------------------------------------------------------


def test_explicit_constructor_requires_every_field():
    config = chematic.PipelineV2Config(
        embed_seed=1,
        max_attempts=4,
        embed_timeout_ms=None,
        use_exp_torsions=False,
        use_small_ring_torsions=False,
        use_macrocycle_torsions=False,
        use_macrocycle_14_bounds=False,
        include_legacy_torsion_heuristic=False,
        stereo_policy="verify_only",
        fail_on_unevaluable_stereo=True,
        force_field_policy="dreiding",
        force_field_max_iterations=100,
        gate_mmff94_torsion_oop=False,
        gate_mmff94_stretch_bend=False,
        ring_torsion_policy="diagnostic_only",
        total_timeout_ms=5000,
    )
    assert config.embed_seed == 1
    assert config.max_attempts == 4
    assert config.stereo_policy == "verify_only"
    assert config.fail_on_unevaluable_stereo is True
    assert config.force_field_policy == "dreiding"
    assert config.ring_torsion_policy == "diagnostic_only"
    assert config.total_timeout_ms == 5000


def test_safe_convenience_constructor_still_requires_judgment_fields():
    with pytest.raises(TypeError):
        chematic.PipelineV2Config.safe()  # force_field/stereo_policy/ring_torsion_policy missing


def test_safe_uses_conservative_defaults_for_everything_else():
    config = chematic.PipelineV2Config.safe(
        force_field="mmff94_with_uff_fallback",
        stereo_policy="repair_and_verify",
        ring_torsion_policy="fail_closed",
    )
    assert config.force_field_policy == "mmff94_with_uff_fallback"
    assert config.stereo_policy == "repair_and_verify"
    assert config.fail_on_unevaluable_stereo is False
    assert config.use_exp_torsions is False
    assert config.use_small_ring_torsions is False
    assert config.use_macrocycle_torsions is False
    assert config.use_macrocycle_14_bounds is False
    assert config.include_legacy_torsion_heuristic is False
    assert config.max_attempts == 8
    assert config.embed_timeout_ms is None
    assert config.total_timeout_ms is None
    assert config.enforce_chirality is False


# ---------------------------------------------------------------------------
# enforce_chirality (v0.14.0, issue #285's E/Z bound fix) -- Python parity
# ---------------------------------------------------------------------------


def test_enforce_chirality_defaults_false_and_is_backward_compatible():
    """Existing callers that never pass enforce_chirality must keep working
    unchanged -- .safe() and the explicit constructor both default it False."""
    config = _safe_config()
    assert config.enforce_chirality is False

    explicit = chematic.PipelineV2Config(
        embed_seed=1,
        max_attempts=4,
        embed_timeout_ms=None,
        use_exp_torsions=False,
        use_small_ring_torsions=False,
        use_macrocycle_torsions=False,
        use_macrocycle_14_bounds=False,
        include_legacy_torsion_heuristic=False,
        stereo_policy="ignore",
        fail_on_unevaluable_stereo=False,
        force_field_policy="none",
        force_field_max_iterations=100,
        gate_mmff94_torsion_oop=False,
        gate_mmff94_stretch_bend=False,
        ring_torsion_policy="fail_closed",
        total_timeout_ms=None,
    )
    assert explicit.enforce_chirality is False


def test_enforce_chirality_true_fixes_but2ene_z_raw_embed():
    """Direct Python-level confirmation that enforce_chirality reaches
    distance_geometry_v2.rs's apply_declared_ez_bounds (issue #285): but2ene_Z
    (C/C=C\\C) is the exact molecule that fix targets -- raw embedding (no
    force field) must satisfy declared E/Z once enforce_chirality is set,
    across multiple seeds, matching the Rust-level corpus measurement."""
    mol = chematic.from_smiles(r"C/C=C\C")
    for seed in range(5):
        config = _safe_config(
            force_field="none",
            stereo_policy="ignore",
            embed_seed=seed,
            max_attempts=1,
            enforce_chirality=True,
        )
        result = mol.embed_pipeline_v2(config)
        assert result["final_stereo"]["is_fully_satisfied"] is True, (
            f"seed {seed}: raw embed must already satisfy declared E/Z"
        )


def test_enforce_chirality_with_repair_and_verify_is_allowed():
    """Revised 2026-08-24 (issue #291 Step A): enforce_chirality +
    stereo_policy="repair_and_verify" was previously rejected as
    InvalidConfiguration -- now validated (see the Rust-level
    pipeline_v2.rs's revised Stage 1 doc entry and
    crates/chematic-3d/examples/issue291_repair_policy_measurement.rs)."""
    mol = chematic.from_smiles(r"C/C=C\C")
    config = _safe_config(
        force_field="none",
        stereo_policy="repair_and_verify",
        enforce_chirality=True,
    )
    result = mol.embed_pipeline_v2(config)
    assert result["final_stereo"]["is_fully_satisfied"] is True


# ---------------------------------------------------------------------------
# expand_implicit_h_through_pipeline / stereo_safe (issue #291/#383)
# ---------------------------------------------------------------------------

TESTOSTERONE = "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O"


def test_expand_implicit_h_through_pipeline_requires_enforce_chirality():
    mol = chematic.from_smiles(TESTOSTERONE)
    config = _safe_config(
        expand_implicit_h_through_pipeline=True,
        enforce_chirality=False,
    )
    with pytest.raises(chematic.PipelineV2Error) as excinfo:
        mol.embed_pipeline_v2(config)
    diag = excinfo.value.diagnostics
    assert diag["cause"] == {"kind": "invalid_configuration"}
    assert diag["stage"] == "validate_config"


def test_stereo_safe_sets_all_three_flags_together():
    config = chematic.PipelineV2Config.stereo_safe(
        force_field="mmff94_with_uff_fallback",
        ring_torsion_policy="fail_closed",
    )
    assert config.stereo_policy == "repair_and_verify"
    assert config.enforce_chirality is True
    assert config.expand_implicit_h_through_pipeline is True


def test_stereo_safe_fixes_testosterone_via_python_binding():
    """Issue #291/#383: testosterone via the Python binding, same seed/
    configuration already Rust-level tested and cross-checked against an
    independent oracle in pipeline_v2.rs's own
    stereo_safe_matches_the_hand_built_configuration_above test."""
    mol = chematic.from_smiles(TESTOSTERONE)
    config = chematic.PipelineV2Config.stereo_safe(
        force_field="mmff94_with_uff_fallback",
        ring_torsion_policy="diagnostic_only",
        embed_seed=0,
    )
    result = mol.embed_pipeline_v2(config)
    assert result["final_stereo"]["is_fully_satisfied"] is True
    assert result["final_stereo"]["n_violations"] == 0
    assert len(result["coords"]) == mol.heavy_atoms


def test_expand_implicit_h_through_pipeline_is_noop_without_declared_stereo():
    """A molecule with no declared stereo must give an identical result with
    the flag on vs. off -- decane has no stereocenters at all."""
    mol = chematic.from_smiles(DECANE)
    base = mol.embed_pipeline_v2(
        _safe_config(
            force_field="mmff94_with_uff_fallback",
            stereo_policy="repair_and_verify",
            enforce_chirality=True,
            embed_seed=0,
        )
    )
    expanded = mol.embed_pipeline_v2(
        _safe_config(
            force_field="mmff94_with_uff_fallback",
            stereo_policy="repair_and_verify",
            enforce_chirality=True,
            expand_implicit_h_through_pipeline=True,
            embed_seed=0,
        )
    )
    assert expanded["coords"] == base["coords"]


# ---------------------------------------------------------------------------
# .pyi stub existence check ("type stub test")
# ---------------------------------------------------------------------------


def test_pyi_declares_the_new_public_surface():
    """Parse __init__.pyi and confirm the new names it declares actually
    exist at runtime with a matching shape -- catches a stub that drifts
    from the compiled extension (wrong name, missing method, etc.)."""
    pyi_path = (
        Path(__file__).resolve().parents[1] / "python" / "chematic" / "__init__.pyi"
    )
    tree = ast.parse(pyi_path.read_text())
    top_level_classes = {
        node.name: node for node in tree.body if isinstance(node, ast.ClassDef)
    }

    assert "PipelineV2Config" in top_level_classes
    assert "PipelineV2Error" in top_level_classes
    assert hasattr(chematic, "PipelineV2Config")
    assert hasattr(chematic, "PipelineV2Error")
    assert issubclass(chematic.PipelineV2Error, ValueError)

    config_methods = {
        n.name
        for n in top_level_classes["PipelineV2Config"].body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert "safe" in config_methods
    assert hasattr(chematic.PipelineV2Config, "safe")
    assert "stereo_safe" in config_methods
    assert hasattr(chematic.PipelineV2Config, "stereo_safe")

    mol_class = next(
        node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "Mol"
    )
    mol_methods = {
        n.name for n in mol_class.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert "embed_pipeline_v2" in mol_methods
    assert hasattr(chematic.Mol, "embed_pipeline_v2")
