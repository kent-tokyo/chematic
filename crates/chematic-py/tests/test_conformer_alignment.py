"""Batch alignment and atomic RMS output motivated by RDKit #9601."""

import numpy as np
import pytest

from chematic import align_conformers_checked

REF = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
        [2.0, 3.0, 4.0],
    ]
)
ROT = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])


def frames():
    return np.array([REF, REF @ ROT + [3.0, -2.0, 7.0], REF - [5.0, 8.0, 3.0]])


@pytest.mark.parametrize(
    "kind", ["array", "list", "tuple", "strided", "readonly", "endian", "integer"]
)
def test_copy_preserving_rigid_alignment(kind):
    source = frames()
    if kind == "list":
        source = source.tolist()
    elif kind == "tuple":
        source = tuple(tuple(map(tuple, frame)) for frame in source)
    elif kind == "strided":
        backing = np.zeros((3, 5, 6))
        backing[:, :, ::2] = source
        source = backing[:, :, ::2]
    elif kind == "readonly":
        source.flags.writeable = False
    elif kind == "endian":
        source = source.astype(">f8")
    elif kind == "integer":
        source = source.astype(np.int64)
    before = np.array(source, copy=True)
    rms = [12.5]
    aligned = align_conformers_checked(source, rms_list=rms)
    np.testing.assert_allclose(aligned, np.repeat(REF[None], 3, axis=0), atol=1e-13)
    np.testing.assert_array_equal(source, before)
    assert aligned.dtype == np.float64 and not np.shares_memory(aligned, before)
    if isinstance(source, np.ndarray):
        assert not np.shares_memory(aligned, source)
    assert rms[0] == 12.5 and rms[1:] == pytest.approx([0.0, 0.0], abs=1e-13)


def test_selection_order_reference_and_transform_of_unselected_atoms():
    source = frames()
    source[2, 4] += [1.0, -2.0, 3.0]
    before = source.copy()
    rms = []
    result = align_conformers_checked(
        source, atom_ids=[0, 1, 2, 3], conformer_ids=[1, 2], rms_list=rms
    )
    np.testing.assert_array_equal(result[:2], source[:2])
    np.testing.assert_allclose(result[2, :4], source[1, :4], atol=1e-13)
    # The fifth atom receives the fitted rigid transform, not a separate fit.
    np.testing.assert_allclose(
        result[2, 4] - source[1, 4], [1.0, -2.0, 3.0] @ ROT, atol=1e-13
    )
    np.testing.assert_array_equal(source, before)
    assert rms == pytest.approx([0.0], abs=1e-13)


def test_weights_rms_uses_atom_count_denominator():
    source = frames()
    source[1, -1] += [1.0, 2.0, -3.0]
    weights = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    rms = []
    result = align_conformers_checked(source, weights=weights, rms_list=rms)
    direct = np.sqrt(np.sum(weights[:, None] * (result[1] - source[0]) ** 2) / 5)
    assert rms[0] == pytest.approx(direct)
    scaled = []
    other = align_conformers_checked(source, weights=weights * 9, rms_list=scaled)
    np.testing.assert_allclose(other, result, atol=1e-13)
    np.testing.assert_allclose(scaled, np.array(rms) * 3, atol=1e-13)


@pytest.mark.parametrize(
    "kwargs", [{}, {"atom_ids": []}, {"conformer_ids": []}, {"weights": []}]
)
def test_empty_selections_and_weights_mean_all(kwargs):
    np.testing.assert_allclose(
        align_conformers_checked(frames(), **kwargs),
        align_conformers_checked(frames()),
        atol=1e-13,
    )


def test_single_frame_and_single_selected_frame_leave_rms_unchanged():
    rms = [2.0]
    source = frames()
    np.testing.assert_array_equal(
        align_conformers_checked(source[:1], rms_list=rms), source[:1]
    )
    np.testing.assert_array_equal(
        align_conformers_checked(source, conformer_ids=[2], rms_list=rms), source
    )
    assert rms == [2.0]


def test_reflections_are_not_used():
    source = np.array([REF, REF * [-1.0, 1.0, 1.0]])
    rms = []
    align_conformers_checked(source, rms_list=rms)
    assert rms[0] > 0.1


@pytest.mark.parametrize(
    "value,exception",
    [
        (None, TypeError),
        ("abc", TypeError),
        ((x for x in []), TypeError),
        ([], ValueError),
        (np.zeros((0, 2, 3)), ValueError),
        (np.zeros((2, 0, 3)), ValueError),
        (np.zeros((2, 3)), ValueError),
        (np.zeros((2, 2, 4)), ValueError),
        ([REF, REF[:2]], ValueError),
        ([REF, 5], TypeError),
        ([REF, np.array(5)], ValueError),
        (np.full((2, 2, 3), np.nan), ValueError),
        (np.full((2, 2, 3), np.inf), ValueError),
        (np.full((2, 2, 3), 1e101), ValueError),
        (np.zeros((2, 2, 3), dtype=bool), TypeError),
        (np.zeros((2, 2, 3), dtype=complex), TypeError),
        (np.zeros((2, 2, 3), dtype=object), TypeError),
        ([[[True, 0, 0]]], TypeError),
        ([REF, [["1", 2.0, 3.0]] * 5], TypeError),
    ],
)
def test_invalid_coordinates_preserve_rms(value, exception):
    rms = [3.0]
    with pytest.raises(exception):
        align_conformers_checked(value, rms_list=rms)
    assert rms == [3.0]


@pytest.mark.parametrize("name", ["atom_ids", "conformer_ids"])
@pytest.mark.parametrize(
    "value,exception",
    [
        ([0, 0], ValueError),
        ([-1], ValueError),
        ([999], ValueError),
        ([True], TypeError),
        ([1.0], TypeError),
        ("01", TypeError),
        ((i for i in [0]), TypeError),
        (np.array([0]), TypeError),
    ],
)
def test_invalid_indices_preserve_outputs(name, value, exception):
    source = frames()
    before = source.copy()
    rms = [3.0]
    with pytest.raises(exception):
        align_conformers_checked(source, **{name: value}, rms_list=rms)
    assert rms == [3.0]
    np.testing.assert_array_equal(source, before)


@pytest.mark.parametrize(
    "weights,exception",
    [
        ([1.0, 2.0], ValueError),
        ([1.0, 2.0, 3.0, 4.0, 0.0], ValueError),
        ([1.0, 2.0, 3.0, 4.0, -1.0], ValueError),
        ([np.inf] * 5, ValueError),
        ([np.nan] * 5, ValueError),
        ([True] * 5, TypeError),
        (np.ones(5, dtype=complex), TypeError),
        (np.ones((5, 1)), ValueError),
        ((x for x in [1.0] * 5), TypeError),
        ([1e-300, 1e100, 1.0, 1.0, 1.0], ValueError),
    ],
)
def test_invalid_weights_are_atomic(weights, exception):
    rms = [3.0]
    with pytest.raises(exception):
        align_conformers_checked(frames(), weights=weights, rms_list=rms)
    assert rms == [3.0]


@pytest.mark.parametrize(
    "value,exception",
    [
        ((), TypeError),
        (np.array([]), TypeError),
        ([1], TypeError),
        ([True], TypeError),
        ([np.nan], ValueError),
        ([np.inf], ValueError),
    ],
)
def test_invalid_rms_destination(value, exception):
    with pytest.raises(exception):
        align_conformers_checked(frames(), rms_list=value)


def test_rms_list_subclasses_are_rejected():
    class BadList(list):
        def extend(self, values):
            raise AssertionError("must not call arbitrary destination code")

    with pytest.raises(TypeError):
        align_conformers_checked(frames(), rms_list=BadList())


def test_second_fit_failure_does_not_partially_append(monkeypatch):
    from chematic import _conformer_alignment as module

    original = module.rdAlignment.GetAlignmentTransform
    calls = []

    def fit(*args):
        calls.append(1)
        if len(calls) == 2:
            raise ValueError("injected solver failure")
        return original(*args)

    monkeypatch.setattr(module.rdAlignment, "GetAlignmentTransform", fit)
    rms = [5.0]
    source = frames()
    before = source.copy()
    with pytest.raises(ValueError, match="injected"):
        align_conformers_checked(source, rms_list=rms)
    assert rms == [5.0]
    np.testing.assert_array_equal(source, before)


def test_size_limits_before_allocating_large_output():
    for source in [
        np.broadcast_to(0.0, (1025, 1, 3)),
        np.broadcast_to(0.0, (1, 65537, 3)),
        np.broadcast_to(0.0, (1000, 1001, 3)),
    ]:
        with pytest.raises(ValueError):
            align_conformers_checked(source)


def test_single_selected_frame_still_validates_weight_ratios():
    with pytest.raises(ValueError, match="weight ratios"):
        align_conformers_checked(
            frames(), conformer_ids=[0], weights=[1e-300, 1e100, 1.0, 1.0, 1.0]
        )


@pytest.mark.parametrize("n_atoms", [1, 2, 3])
def test_degenerate_fits_preserve_distances_and_residuals(n_atoms):
    source = frames()[:, :n_atoms]
    rms = []
    result = align_conformers_checked(source, rms_list=rms)
    assert rms == pytest.approx([0.0, 0.0], abs=1e-12)
    np.testing.assert_allclose(
        np.linalg.norm(result[:, 1:] - result[:, :1], axis=2),
        np.linalg.norm(source[:, 1:] - source[:, :1], axis=2),
        atol=1e-12,
    )
