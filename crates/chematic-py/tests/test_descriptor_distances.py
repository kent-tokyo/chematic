"""Bounded descriptor input and stable Euclidean norms inspired by RDKit #9601."""
import itertools
import math
import pytest
import chematic


def test_lower_triangle_and_duplicates():
    assert chematic.euclidean_distance_matrix_checked([[0,0], [3,4], [3,0]]) == [5,3,4]
    assert chematic.euclidean_distance_matrix_checked([[1], [1], [2]]) == [0,1,1]


@pytest.mark.parametrize('rows,expected', [([],[]), ([[123]],[]), ([[],[],[]],[0,0,0])])
def test_empty_and_zero_dimensions(rows, expected):
    assert chematic.euclidean_distance_matrix_checked(rows) == expected


@pytest.mark.parametrize('value', [math.nan,math.inf,-math.inf])
def test_nonfinite_including_singleton(value):
    with pytest.raises(ValueError, match='row 0, column 0'):
        chematic.euclidean_distance_matrix_checked([[value]])


@pytest.mark.parametrize('rows', [[[0],[1,2]], [[],[0]], [[0],[]]])
def test_ragged(rows):
    with pytest.raises(ValueError, match='descriptor row 1'):
        chematic.euclidean_distance_matrix_checked(rows)


@pytest.mark.parametrize('rows', [None, 1, '12', b'12', bytearray(b'12'), [1], ['12'], [b'12'], [[None]], [['1']], [[1j]]])
def test_wrong_types(rows):
    with pytest.raises(TypeError):
        chematic.euclidean_distance_matrix_checked(rows)


@pytest.mark.parametrize('value', [1e200,1e-200,math.ulp(0.0),float.fromhex('0x1.fffffffffffffp+1023')])
def test_extreme_representable_norm(value):
    assert chematic.euclidean_distance_matrix_checked([[0],[value]]) == [value]


@pytest.mark.parametrize('rows', [[[0,0],[1.7e308,1.7e308]], [[0],[1.7e308],[-1.7e308]]])
def test_unrepresentable_distances(rows):
    with pytest.raises(ValueError, match='exceeds finite float range'):
        chematic.euclidean_distance_matrix_checked(rows)


@pytest.mark.parametrize('rows,message', [([[]]*1025, '1024 rows'), ([[0]*257], '256 features'), ([[0]*256]*1024, 'pair-feature')])
def test_limits(rows, message):
    with pytest.raises(ValueError, match=message):
        chematic.euclidean_distance_matrix_checked(rows)


def test_generators_have_bounded_consumption_and_keep_source_errors():
    consumed = []
    def endless():
        for index in itertools.count():
            consumed.append(index)
            yield []
    with pytest.raises(ValueError, match='1024 rows'):
        chematic.euclidean_distance_matrix_checked(endless())
    assert len(consumed) == 1025
    consumed.clear()
    def features():
        for index in itertools.count():
            consumed.append(index)
            yield 1.0
    with pytest.raises(ValueError, match='256 features'):
        chematic.euclidean_distance_matrix_checked([features()])
    assert len(consumed) == 257
    def broken():
        yield [0]
        raise RuntimeError('descriptor source failed')
    with pytest.raises(RuntimeError, match='descriptor source failed'):
        chematic.euclidean_distance_matrix_checked(broken())
    rows = (iter(row) for row in [[0,0],[3,4],[3,0]])
    assert chematic.euclidean_distance_matrix_checked(rows) == [5,3,4]


def test_numpy_rows_and_strides():
    np = pytest.importorskip('numpy')
    rows = np.array([[0,9,0,9], [3,9,4,9], [3,9,0,9]], dtype=np.float64)[:,::2]
    assert not rows.flags.c_contiguous
    assert chematic.euclidean_distance_matrix_checked(rows) == [5,3,4]
    assert chematic.euclidean_distance_matrix_checked(rows.astype(np.float32)) == [5,3,4]
    assert chematic.euclidean_distance_matrix_checked(np.array(rows, dtype=np.int64)) == [5,3,4]


def test_inputs_are_unchanged_and_keyword_signature():
    rows = [[0,0],[3,4]]
    result = chematic.euclidean_distance_matrix_checked(descriptor_matrix=rows)
    result[0] = -1
    assert rows == [[0,0],[3,4]]
    assert chematic.euclidean_distance_matrix_checked(rows) == [5]
