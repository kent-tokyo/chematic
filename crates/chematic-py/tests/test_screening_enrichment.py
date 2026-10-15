"""Ranked enrichment factors: RDKit #7981, caller order and error policy."""
import math

import pytest
import chematic


def test_issue_7981_and_full_library():
    scores = sorted(zip([0., 0., 0., 0., 1.], [False, False, False, True, True]), reverse=True)
    labels = [active for _, active in scores]
    assert chematic.enrichment_factors(labels, [.9, 1.]) == [1., 1.]


def test_fraction_order_duplicates_and_rounding():
    assert chematic.enrichment_factors([True, False, True, False, False], [1., .2, .01, .2, .5, .9]) == pytest.approx([1., 2.5, 2.5, 2.5, 5/3, 1.])
    assert chematic.enrichment_factors([True, True, False, False, False], [.21]) == [2.5]


@pytest.mark.parametrize("fraction", [.01, .5, .9, 1.])
def test_single_active_entry(fraction):
    assert chematic.enrichment_factors([True], [fraction]) == [1.]


def test_zero_top_hits():
    assert chematic.enrichment_factors([False, False, True, True], [.5, 1.]) == [0., 1.]


@pytest.mark.parametrize("labels,fractions,message", [
    ([], [1.], "ranking is empty"), ([True], [], "fraction list is empty"),
    ([False, False], [1.], "no actives"),
])
def test_empty_or_undefined_inputs(labels, fractions, message):
    with pytest.raises(ValueError, match=message):
        chematic.enrichment_factors(labels, fractions)


@pytest.mark.parametrize("fraction", [0., -0., -.1, 1.0001, math.nan, math.inf, -math.inf])
def test_invalid_fraction(fraction):
    with pytest.raises(ValueError, match="index 1"):
        chematic.enrichment_factors([True], [1., fraction])


def test_tiny_fraction_selects_one():
    assert chematic.enrichment_factors([True, False], [math.ulp(0.)]) == [2.]


def test_scores_are_not_implicitly_sorted():
    assert chematic.enrichment_factors([True, False], [.5]) == [2.]
    assert chematic.enrichment_factors([False, True], [.5]) == [0.]


@pytest.mark.parametrize("label", [0, 1, .5, "active", None])
def test_nonboolean_labels_rejected(label):
    with pytest.raises(TypeError):
        chematic.enrichment_factors([label], [1.])


def test_numpy_boolean_labels():
    import numpy as np
    assert chematic.enrichment_factors(np.array([True, False]), [.5, 1.]) == [2., 1.]
