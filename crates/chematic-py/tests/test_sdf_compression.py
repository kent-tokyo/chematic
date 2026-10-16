"""Per-file gzip selection motivated by Open Babel #3030."""
import gzip
import os

import pytest
import chematic


def record(smiles="CCO", name="ethanol"):
    block = chematic.from_smiles(smiles).to_mol_block()
    block = name + "\n" + block.split("\n", 1)[1]
    return (block + "> <Activity>\n7.2\n\n$$$$\n").encode("utf-8")


def rows(path, **kwargs):
    return list(chematic.iter_sdf_checked(path, **kwargs))


@pytest.mark.parametrize("mode", ["gzip", "auto"])
@pytest.mark.parametrize("path_kind", ["str", "bytes", "path"])
def test_gzip_records_unicode_and_path_protocols(tmp_path, mode, path_kind):
    path = tmp_path / "not-gzip.dat"
    path.write_bytes(gzip.compress(record(name="日本語") + record("CN", "amine")))
    source = {"str": str(path), "bytes": os.fsencode(path), "path": path}[path_kind]
    result = rows(source, compression=mode)
    assert [index for index, _, _ in result] == [0, 1]
    assert all(error is None for _, _, error in result)
    assert [rec.name for _, rec, _ in result] == ["日本語", "amine"]
    assert [rec.mol.atom_count for _, rec, _ in result] == [3, 2]
    assert result[0][1].get("Activity") == "7.2"


def test_fresh_inputs_after_gzip_full_partial_and_failed_reads(tmp_path):
    plain = tmp_path / "plain.sdf.gz"
    plain.write_bytes(record("CCO"))
    zipped = tmp_path / "compressed.sdf"
    zipped.write_bytes(gzip.compress(record("CN") * 2))
    bad = tmp_path / "bad.sdf"
    bad.write_bytes(b"\x1f\x8bbroken")
    for _ in range(20):
        assert len(rows(zipped, compression="auto")) == 2
        assert rows(plain, compression="auto")[0][1].mol.atom_count == 3
        iterator = chematic.iter_sdf_checked(zipped, compression="auto")
        assert next(iterator)[1].mol.atom_count == 2
        del iterator
        assert rows(bad, compression="auto")[0][2] is not None
        assert rows(plain, compression="auto")[0][1].mol.atom_count == 3


def test_default_plain_behavior_and_explicit_modes(tmp_path):
    path = tmp_path / "input.sdf.gz"
    path.write_bytes(record())
    assert rows(path)[0][1].name == "ethanol"
    assert rows(path, compression="auto")[0][1].name == "ethanol"
    assert rows(path, compression="gzip")[0][2] is not None
    path.write_bytes(gzip.compress(record()))
    assert not any(rec is not None for _, rec, _ in rows(path))
    assert not any(rec is not None for _, rec, _ in rows(path, compression="plain"))
    assert rows(path, compression="auto")[0][1].name == "ethanol"


@pytest.mark.parametrize("mode", ["", "GZIP", "zlib", "zip"])
def test_invalid_compression_refused_before_file_open(tmp_path, mode):
    with pytest.raises(ValueError, match="compression"):
        chematic.iter_sdf_checked(tmp_path / "missing", compression=mode)


@pytest.mark.parametrize("mode", [None, 42, b"gzip"])
def test_nonstring_compression_is_type_error(tmp_path, mode):
    with pytest.raises(TypeError):
        chematic.iter_sdf_checked(tmp_path / "missing", compression=mode)


@pytest.mark.parametrize("split", [1, 13, 180])
def test_concatenated_members_split_inside_record(tmp_path, split):
    data = record(name="日本語")
    path = tmp_path / "members.gz"
    path.write_bytes(gzip.compress(data[:split]) + gzip.compress(data[split:])
                     + gzip.compress(record("CN", "amine")))
    result = rows(path, compression="auto")
    assert len(result) == 2
    assert [rec.name for _, rec, _ in result] == ["日本語", "amine"]
    assert all(error is None for _, _, error in result)


@pytest.mark.parametrize("damage", ["crc", "truncated", "junk", "bad_member"])
def test_integrity_errors_yield_terminal_error(tmp_path, damage):
    data = bytearray(gzip.compress(record()))
    if damage == "crc":
        data[-8] ^= 1
    elif damage == "truncated":
        del data[-4:]
    elif damage == "junk":
        data.extend(b"junk")
    else:
        data.extend(b"\x1f\x8btruncated")
    path = tmp_path / "bad.gz"
    path.write_bytes(data)
    iterator = chematic.iter_sdf_checked(path, compression="auto")
    result = list(iterator)
    assert result[-1][1] is None
    assert result[-1][2] is not None
    assert list(iterator) == []


@pytest.mark.parametrize("limit", ["max_input_bytes", "max_record_bytes", "max_line_bytes", "max_records"])
def test_limits_apply_after_decompression(tmp_path, limit):
    path = tmp_path / "limited.gz"
    path.write_bytes(gzip.compress(record() * 2))
    result = rows(path, compression="auto", **{limit: 1})
    assert result[-1][1] is None
    assert "limit" in result[-1][2].lower()


def test_malformed_record_position_and_later_valid_record(tmp_path):
    path = tmp_path / "mixed.gz"
    path.write_bytes(gzip.compress(record() + b"malformed\n$$$$\n" + record("CN")))
    result = rows(path, compression="auto")
    assert [index for index, _, _ in result] == [0, 1, 2]
    assert result[0][1].mol.atom_count == 3
    assert result[1][1] is None and result[1][2]
    assert result[2][1].mol.atom_count == 2


@pytest.mark.parametrize("data", [b"", gzip.compress(b""), b"\n$$$$\n"])
def test_empty_inputs(tmp_path, data):
    path = tmp_path / "empty"
    path.write_bytes(data)
    assert rows(path, compression="auto") == []


def test_plain_one_byte_is_replayed_for_parsing(tmp_path):
    path = tmp_path / "one"
    path.write_bytes(b"x")
    result = rows(path, compression="auto")
    assert len(result) == 1
    assert result[0][1] is None and result[0][2]


def test_missing_file_is_ioerror(tmp_path):
    with pytest.raises(OSError):
        chematic.iter_sdf_checked(tmp_path / "missing", compression="auto")


def test_v3000_stereo_uses_complete_sdf_reader(tmp_path):
    path = tmp_path / "stereo.sdf"
    with chematic.SDWriter(path) as writer:
        writer.SetForceV3000(True)
        writer.write(chematic.from_smiles("F[C@H](Cl)Br"))
    expected = rows(path)
    path.write_bytes(gzip.compress(path.read_bytes()))
    actual = rows(path, compression="auto")
    assert actual[0][1].smiles == expected[0][1].smiles
    assert "@" in actual[0][1].smiles
