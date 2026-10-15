import math
import pytest
import chematic
from chematic.rdkit_compat import DataStructs, ExplicitBitVect


def packed(width, on):
    result = bytearray((width + 7) // 8)
    for bit in on:
        result[bit // 8] |= 1 << (bit % 8)
    return bytes(result)


def test_sparse_folds_including_equal_minimum():
    result = chematic.fold_to_target_density(packed(128, range(16)), 0.25, 64)
    assert result == dict(fingerprint=packed(64, range(16)), bit_width=64,
        original_bits=128, on_bit_count=16, density=0.25, folds=1,
        target_reached=True, stop_reason="target_reached")


def test_dense_or_equal_density_is_unchanged():
    for density in [0., 0.5, 1.]:
        raw = packed(256, range(256))
        result = chematic.fold_to_target_density(raw, density)
        assert result["fingerprint"] == raw
        assert result["folds"] == 0


def test_collisions_and_odd_width_padding():
    result = chematic.fold_to_target_density(packed(65, [0, 32, 64]), 1., 32, bit_width=65)
    assert result["bit_width"] == 32
    assert result["on_bit_count"] == 1
    assert result["fingerprint"] == b"\x01\x00\x00\x00"
    assert result["stop_reason"] == "minimum_length"


def test_zero_bits_and_already_short_input_report_unmet_target():
    for width, minimum in [(128, 1), (8, 64)]:
        result = chematic.fold_to_target_density(packed(width, []), 0.3, minimum)
        assert result["bit_width"] == min(width, minimum)
        assert not result["target_reached"]


@pytest.mark.parametrize("density", [math.nan, math.inf, -math.inf, -.1, 1.1])
def test_invalid_density(density):
    with pytest.raises(ValueError, match="finite"):
        chematic.fold_to_target_density(b"\0", density)


@pytest.mark.parametrize("raw,width,message", [
    (b"", None, "at least one bit"), (b"\0", 0, "at least one bit"),
    (b"\0\0", 8, "expected 1 bytes"), (b"\0", 9, "expected 2 bytes"),
    (b"\x80", 7, "padding"), (b"", 1048577, "exceeds limit"),
])
def test_invalid_packed_input(raw, width, message):
    with pytest.raises(ValueError, match=message):
        chematic.fold_to_target_density(raw, bit_width=width)


def test_zero_minimum_and_negative_integer_options():
    with pytest.raises(ValueError, match="positive"):
        chematic.fold_to_target_density(b"\0", min_length=0)
    for kwargs in [dict(min_length=-1), dict(bit_width=-1)]:
        with pytest.raises(OverflowError):
            chematic.fold_to_target_density(b"\0", **kwargs)


def test_borrowed_bytes_type_and_maximum_width():
    with pytest.raises(TypeError):
        chematic.fold_to_target_density([0])
    result = chematic.fold_to_target_density(bytes(131072), 0., 1)
    assert result["bit_width"] == 1048576
    with pytest.raises(ValueError, match="exceeds limit"):
        chematic.fold_to_target_density(bytes(131073))


def test_compat_adapter_folds_without_mutating_input():
    fp = ExplicitBitVect(128)
    for bit in range(16):
        fp.SetBit(bit)
    result = DataStructs.FoldToTargetDensity(fp, 0.25, 64)
    assert len(result) == 64
    assert result.GetNumOnBits() == 16
    assert result.GetOnBits() == list(range(16))
    assert len(fp) == 128
    result.SetBit(63)
    assert fp.GetNumOnBits() == 16


def test_compat_noop_returns_independent_copy_and_checks_type():
    fp = ExplicitBitVect(8)
    fp.SetBit(0)
    result = DataStructs.FoldToTargetDensity(fp, 0., 1)
    assert result is not fp
    assert result.ToBitString() == fp.ToBitString()
    with pytest.raises(TypeError, match="ExplicitBitVect"):
        DataStructs.FoldToTargetDensity(b"\0")


def test_compat_odd_width_and_zero_width():
    fp = ExplicitBitVect(5)
    for bit in [0, 2, 4]: fp.SetBit(bit)
    result = DataStructs.FoldToTargetDensity(fp, 0.8, 2)
    assert len(result) == 2
    assert result.GetOnBits() == [0]
    with pytest.raises(ValueError):
        DataStructs.FoldToTargetDensity(ExplicitBitVect(0))
