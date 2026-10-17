"""Original source-byte and balanced file-backed SDF partition contracts."""

from dataclasses import FrozenInstanceError, replace
import gzip
import hashlib
import os

import numpy as np
import pytest
import chematic
from chematic import _sdf_split as split


def source(count, newline=b"\n", terminated=True):
    blocks = []
    for i in range(count):
        block = (
            chematic.from_smiles(["CCO", "C[NH3+]", "C[C@H](O)F"][i % 3])
            .to_mol_block()
            .split("\n", 1)[1]
        )
        name = f"original {i}" if i % 2 else ""
        record = (
            name
            + "\n"
            + block
            + "> <sample>\n"
            + str(i)
            + "\n\n> <opaque>\ninline $$$$ text\n未知\n\n"
        ).encode()
        if i < count - 1 or terminated:
            record += b"$$$$\n"
        blocks.append(record.replace(b"\n", newline))
    return b"".join(blocks)


def write(tmp_path, data):
    path = tmp_path / "source.sdf"
    path.write_bytes(data)
    return path


@pytest.mark.parametrize(
    "records,parts",
    [(1, 1), (10, 1), (10, 3), (10, 10), (19, 7), (101, 10), (1003, 400)],
)
def test_balance_and_exact_source(tmp_path, records, parts):
    data = source(records)
    path = write(tmp_path, data)
    plan = chematic.plan_sdf_split(path, parts, prefix="batch")
    assert plan.record_count == records and plan.source_bytes == len(data)
    assert plan.source_sha256 == hashlib.sha256(data).hexdigest()
    expected = [records // parts + int(i < records % parts) for i in range(parts)]
    assert [p.record_count for p in plan.parts] == expected and max(expected) - min(
        expected
    ) <= 1
    assert len(plan.parts) == parts and [p.filename for p in plan.parts] == sorted(
        p.filename for p in plan.parts
    )
    chunks = []
    first = offset = 0
    for i, p in enumerate(plan.parts):
        raw = b"".join(plan.iter_part(i, chunk_bytes=127))
        chunks.append(raw)
        assert raw == data[p.byte_offset : p.byte_offset + p.byte_len]
        assert hashlib.sha256(raw).hexdigest() == p.sha256
        assert p.first_record == first and p.byte_offset == offset
        first += p.record_count
        offset += p.byte_len
    assert b"".join(chunks) == data and offset == len(data)
    assert path.read_bytes() == data and not list(tmp_path.glob("batch*"))


@pytest.mark.parametrize("newline", [b"\n", b"\r\n"])
@pytest.mark.parametrize("terminated", [False, True])
def test_blank_titles_unicode_and_unterminated_eof(tmp_path, newline, terminated):
    data = source(3, newline, terminated) + b" \t\r\n"
    plan = chematic.plan_sdf_split(write(tmp_path, data), 2)
    assert plan.unterminated_final_record is (not terminated)
    assert b"".join(b"".join(plan.iter_part(i)) for i in range(2)) == data


def test_opaque_invalid_record_duplicate_fields_and_eof_delimiter(tmp_path):
    data = b"invalid MOL\n> <x>\nfirst\n\n> <x>\nsecond\n\n$$$$"
    plan = chematic.plan_sdf_split(write(tmp_path, data), 1)
    assert not plan.unterminated_final_record and b"".join(plan.iter_part(0)) == data
    with pytest.raises(FrozenInstanceError):
        plan.parts[0].filename = "changed"


@pytest.mark.parametrize("data", [b"", b" \n\r\n", b"$$$$\n", b"a\n$$$$\n$$$$\n"])
def test_empty_input_or_empty_record(tmp_path, data):
    with pytest.raises(ValueError):
        chematic.plan_sdf_split(write(tmp_path, data), 1)


def test_non_utf8_and_padded_dollars_are_retained(tmp_path):
    data = b"raw\n\xff\n$$$$ \n $$$$\n$$$$\n"
    plan = chematic.plan_sdf_split(write(tmp_path, data), 1)
    assert plan.record_count == 1 and b"".join(plan.iter_part(0)) == data


@pytest.mark.parametrize(
    "prefix", ["", "../x", "/tmp/x", "a/b", "a\\b", ".hidden", "café", "x" * 65]
)
def test_filename_stem(tmp_path, prefix):
    with pytest.raises(ValueError, match="prefix"):
        chematic.plan_sdf_split(write(tmp_path, source(1)), 1, prefix=prefix)


@pytest.mark.parametrize(
    "key,value",
    [
        ("parts", 0),
        ("parts", 10001),
        ("max_input_bytes", 0),
        ("max_input_bytes", 2**30 + 1),
        ("max_records", 0),
        ("max_records", 1_000_001),
        ("max_line_bytes", 0),
        ("max_line_bytes", 2**20 + 1),
        ("max_record_bytes", 0),
        ("max_record_bytes", 64 * 2**20 + 1),
    ],
)
def test_option_ranges(tmp_path, key, value):
    kwargs = {key: value}
    kwargs.setdefault("parts", 1)
    with pytest.raises(ValueError):
        chematic.plan_sdf_split(write(tmp_path, source(1)), **kwargs)


@pytest.mark.parametrize("value", [True, 1.5, "2"])
def test_count_types(tmp_path, value):
    with pytest.raises(TypeError):
        chematic.plan_sdf_split(write(tmp_path, source(2)), value)


def test_exact_limits_and_numpy_indices(tmp_path):
    data = source(3)
    path = write(tmp_path, data)
    plan = chematic.plan_sdf_split(
        path, np.int64(2), max_input_bytes=len(data), max_records=3
    )
    assert b"".join(plan.iter_part(np.int64(0), chunk_bytes=np.int64(1)))
    for kwargs in [
        {"max_input_bytes": len(data) - 1},
        {"max_records": 2},
        {"max_line_bytes": 10},
        {"max_record_bytes": 30},
    ]:
        with pytest.raises(ValueError, match="limit"):
            chematic.plan_sdf_split(path, 1, **kwargs)
    with pytest.raises(ValueError, match="exceed"):
        chematic.plan_sdf_split(path, 4)


def test_iter_bounds_and_invalid_constructed_plan(tmp_path):
    plan = chematic.plan_sdf_split(write(tmp_path, source(2)), 2)
    for index in [-1, 2]:
        with pytest.raises(IndexError):
            plan.iter_part(index)
    with pytest.raises(TypeError):
        plan.iter_part(True)
    for size in [0, 2**20 + 1]:
        with pytest.raises(ValueError):
            plan.iter_part(0, chunk_bytes=size)
    with pytest.raises(TypeError):
        plan.iter_part(0, chunk_bytes=True)
    bad = replace(plan, parts=(replace(plan.parts[0], byte_len=-1),))
    with pytest.raises(ValueError, match="byte range"):
        bad.iter_part(0)


def test_modified_and_replaced_sources(tmp_path):
    path = write(tmp_path, source(3))
    plan = chematic.plan_sdf_split(path, 2)
    stream = plan.iter_part(0, chunk_bytes=5)
    assert len(next(stream)) == 5
    path.write_bytes(b"changed" + source(3))
    with pytest.raises(ValueError, match="changed|truncated"):
        list(stream)
    with pytest.raises(ValueError, match="changed"):
        list(plan.iter_part(1))
    plan = chematic.plan_sdf_split(path, 2)
    other = tmp_path / "replacement.sdf"
    other.write_bytes(path.read_bytes())
    other.replace(path)
    with pytest.raises(ValueError, match="changed"):
        list(plan.iter_part(0))


def test_hash_guard_even_when_stat_check_is_bypassed(tmp_path, monkeypatch):
    path = write(tmp_path, source(1))
    plan = chematic.plan_sdf_split(path, 1)
    monkeypatch.setattr(split, "_identity", lambda f: plan._source_identity)
    data = path.read_bytes()
    path.write_bytes(b"X" + data[1:])
    with pytest.raises(ValueError, match="hash not verified"):
        list(plan.iter_part(0))


def test_changed_source_between_scans(tmp_path, monkeypatch):
    path = write(tmp_path, source(2))
    original = split._parts

    def change(file, ends, count, prefix):
        data = path.read_bytes()
        path.write_bytes(b"X" + data[1:])
        return original(file, ends, count, prefix)

    monkeypatch.setattr(split, "_parts", change)
    with pytest.raises(ValueError, match="changed"):
        chematic.plan_sdf_split(path, 2)


def test_bytes_path_cwd_change_and_symlink(tmp_path, monkeypatch):
    data = source(2)
    path = write(tmp_path, data)
    monkeypatch.chdir(tmp_path)
    plan = chematic.plan_sdf_split(b"source.sdf", 2)
    monkeypatch.chdir(tmp_path.parent)
    assert b"".join(plan.iter_part(0)) == data[: plan.parts[0].byte_len]
    link = tmp_path / "link.sdf"
    try:
        link.symlink_to(path)
    except (OSError, NotImplementedError):
        pytest.skip("symlink creation is unavailable on this platform")
    assert chematic.plan_sdf_split(link, 1).source_sha256 == plan.source_sha256


def test_gzip_directory_and_fifo(tmp_path):
    with pytest.raises(ValueError, match="compressed"):
        chematic.plan_sdf_split(write(tmp_path, gzip.compress(source(2))), 1)
    with pytest.raises((ValueError, OSError)):
        chematic.plan_sdf_split(tmp_path, 1)
    if hasattr(os, "mkfifo"):
        fifo = tmp_path / "fifo.sdf"
        os.mkfifo(fifo)
        with pytest.raises(ValueError, match="regular"):
            chematic.plan_sdf_split(fifo, 1)


def test_maximum_part_count_names(tmp_path):
    plan = chematic.plan_sdf_split(write(tmp_path, b"raw\n$$$$\n" * 10000), 10000)
    assert len(plan.parts) == 10000
    assert (
        plan.parts[0].filename == "split00001.sdf"
        and plan.parts[-1].filename == "split10000.sdf"
    )
    assert all(p.record_count == 1 for p in plan.parts)


def test_record_limit_counts_whitespace_attached_to_last_record(tmp_path):
    path = write(tmp_path, b"raw\n$$$$\n   \n")
    assert chematic.plan_sdf_split(path, 1, max_record_bytes=13).source_bytes == 13
    with pytest.raises(ValueError, match="trailing whitespace"):
        chematic.plan_sdf_split(path, 1, max_record_bytes=12)
