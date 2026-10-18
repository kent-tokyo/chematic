"""Ranked-record EF adapter: endpoint correctness and typed staged inputs."""

from fractions import Fraction
import itertools
import math

import numpy as np
import pytest

from chematic.rdkit_compat import Scoring as S


def test_issue_7981_direct_record_api_and_full_library():
    scores = sorted(zip([0.0, 0.0, 0.0, 0.0, 1.0], [0, 0, 0, 1, 1]), reverse=True)
    assert S.CalcEnrichment(scores, col=1, fractions=[0.2, 0.9, 1.0]) == [2.5, 1.0, 1.0]
    assert S.CalcEnrichment([(1.0, 1)], 1, [0.01, 0.9, 1.0]) == [1.0, 1.0, 1.0]


def test_fraction_order_duplicates_selected_counts_and_zero_top_hits():
    scores = [(9, 1), (8, 0), (7, 1), (6, 0), (5, 0)]
    actual = S.CalcEnrichment(scores, 1, [1.0, 0.2, 0.01, 0.2, 0.5, 0.9])
    assert actual == pytest.approx([1.0, 2.5, 2.5, 2.5, 5 / 3, 1.0])
    assert S.CalcEnrichment([(8, 0), (9, 1)], 1, [0.5]) == [0.0]
    assert S.CalcEnrichment([(1, 1), (2, 1), (3, 0), (4, 0), (5, 0)], 1, [0.21]) == [
        2.5
    ]


def test_list_tuple_generator_numpy_and_negative_columns():
    rows = [(0.9, True), (0.8, False), (0.7, True)]
    expected = [1.5, 1.0]
    for records in [rows, tuple(rows), (row for row in rows), np.array(rows)]:
        assert (
            S.CalcEnrichment(records, np.int64(-1), (f for f in [0.1, 1.0])) == expected
        )
    assert S.CalcEnrichment([[True, object()], [False, object()]], 0, [0.5]) == [2.0]


@pytest.mark.parametrize(
    "label",
    [
        False,
        True,
        0,
        1,
        0.0,
        1.0,
        np.bool_(True),
        np.int64(1),
        np.float32(0),
        Fraction(1),
    ],
)
def test_binary_label_values(label):
    assert S.CalcEnrichment([(1, label), (0, True)], 1, [1.0]) == [1.0]


@pytest.mark.parametrize(
    "label,error",
    [
        (2, ValueError),
        (-1, ValueError),
        (0.5, ValueError),
        (math.nan, ValueError),
        (math.inf, ValueError),
        ("1", TypeError),
        (None, TypeError),
        (object(), TypeError),
        (complex(1), TypeError),
    ],
)
def test_nonbinary_or_truthy_labels_are_refused(label, error):
    with pytest.raises(error):
        S.CalcEnrichment([(1, label)], 1, [1.0])


def test_arbitrary_truthiness_is_never_called():
    class TruthTrap:
        def __bool__(self):
            raise AssertionError("truthiness must not run")

    with pytest.raises(TypeError):
        S.CalcEnrichment([(1, TruthTrap())], 1, [1.0])


@pytest.mark.parametrize(
    "col,error",
    [
        (True, TypeError),
        (np.bool_(False), TypeError),
        (1.0, TypeError),
        ("1", TypeError),
        (1024, ValueError),
        (-1025, ValueError),
    ],
)
def test_column_type_and_hard_bounds(col, error):
    with pytest.raises(error):
        S.CalcEnrichment([(1, True)], col, [1.0])


@pytest.mark.parametrize(
    "rows,col",
    [([], 0), ([[]], 0), ([(1, True)], 2), ([(1, True)], -3), ([list(range(1025))], 0)],
)
def test_empty_rows_missing_columns_and_wide_rows(rows, col):
    with pytest.raises(ValueError):
        S.CalcEnrichment(rows, col, [1.0])


@pytest.mark.parametrize("rows", [None, "records", [True], ["10"], [iter([1, 0])]])
def test_records_require_indexed_nontext_rows(rows):
    with pytest.raises(TypeError):
        S.CalcEnrichment(rows, 0, [1.0])


@pytest.mark.parametrize(
    "value,error",
    [
        (0.0, ValueError),
        (-0.1, ValueError),
        (1.0001, ValueError),
        (math.nan, ValueError),
        (math.inf, ValueError),
        (True, TypeError),
        (np.bool_(False), TypeError),
        ("0.5", TypeError),
        (None, TypeError),
    ],
)
def test_fraction_values(value, error):
    with pytest.raises(error):
        S.CalcEnrichment([(1, True)], 1, [value])


def test_empty_fractions_undefined_all_inactive_and_tiny_fraction():
    with pytest.raises(ValueError, match="fraction list is empty"):
        S.CalcEnrichment([(1, True)], 1, [])
    with pytest.raises(ValueError, match="no actives"):
        S.CalcEnrichment([(1, 0), (0, False)], 1, [0.5])
    assert S.CalcEnrichment([(1, 1), (0, 0)], 1, [math.ulp(0.0)]) == [2.0]


def test_input_ownership_and_fresh_results():
    rows = [[1.0, 1], [0.0, 0]]
    fractions = [0.5, 1.0]
    before = [r.copy() for r in rows]
    result = S.CalcEnrichment(rows, 1, fractions)
    result[:] = [99.0]
    assert rows == before and fractions == [0.5, 1.0]
    assert S.CalcEnrichment(rows, 1, fractions) == [2.0, 1.0]


def test_iteration_and_indexing_exceptions_propagate_before_native_call(monkeypatch):
    import chematic._scoring_enrichment as impl

    def forbidden(*args):
        raise AssertionError("native calculation must not start")

    monkeypatch.setattr(impl, "_enrichment_factors", forbidden)
    error = RuntimeError("source failed")

    def rows():
        yield (1.0, 1)
        raise error

    def fractions():
        yield 0.5
        raise error

    class Broken:
        def __len__(self):
            return 2

        def __getitem__(self, index):
            raise error

    for records, cutoffs in [
        (rows(), [0.5]),
        ([(1.0, 1)], fractions()),
        ([Broken()], [0.5]),
    ]:
        with pytest.raises(RuntimeError) as caught:
            S.CalcEnrichment(records, 1, cutoffs)
        assert caught.value is error


def test_invalid_fractions_do_not_consume_records():
    def rows():
        raise AssertionError("rows must not be consumed")
        yield

    with pytest.raises(ValueError):
        S.CalcEnrichment(rows(), 1, [0.0])


def test_actual_hard_record_and_fraction_boundaries():
    assert S.CalcEnrichment(itertools.repeat((1, 1), 100000), 1, [1.0]) == [1.0]
    with pytest.raises(ValueError, match="100000-record"):
        S.CalcEnrichment(itertools.repeat((1, 1)), 1, [1.0])
    assert len(S.CalcEnrichment([(1, 1)], 1, itertools.repeat(1.0, 4096))) == 4096
    with pytest.raises(ValueError, match="4096-entry"):
        S.CalcEnrichment([(1, 1)], 1, itertools.repeat(1.0))
    row = [object()] * 1024
    row[0] = True
    assert S.CalcEnrichment([row], -1024, [1.0]) == [1.0]


def test_exhaustive_small_binary_rankings_against_exact_rational_definition():
    for count in range(1, 9):
        for labels in itertools.product([0, 1], repeat=count):
            active = sum(labels)
            if active == 0:
                continue
            fractions = [0.01, 0.21, 0.5, 0.9, 1.0, 0.21]
            expected = []
            for fraction in fractions:
                selected = max(1, math.ceil(count * fraction))
                expected.append(
                    float(Fraction(sum(labels[:selected]) * count, selected * active))
                )
            rows = [(object(), label) for label in labels]
            assert S.CalcEnrichment(rows, 1, fractions) == pytest.approx(
                expected, rel=1e-14, abs=1e-14
            )


def test_public_namespace_is_exported():
    namespace = {}
    exec("from chematic.rdkit_compat import *", namespace)
    assert namespace["Scoring"] is S
    assert "RankedScoreRow" in namespace
