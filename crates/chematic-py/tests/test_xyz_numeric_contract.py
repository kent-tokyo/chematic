"""Finite XYZ/ExtXYZ real reads must not silently turn nonzero text into zero."""
import json
from pathlib import Path
import struct

import chematic
import pytest

PACKET = json.loads((Path(__file__).parents[3] / "validation/xyz_numeric_read/cases.json")
                    .read_text(encoding="utf-8"))


def source(token, field):
    row = ["0", "0", "0"]
    if field < 3:
        row[field] = token
    header = "number"
    extra = ""
    if field == 3:
        header = f'Lattice="1 0 0 0 {token} 0 0 0 1"'
    elif field == 4:
        header = "Properties=species:S:1:pos:R:3:force:R:3"
        extra = f" 0 {token} 0"
    return f"1\n{header}\nC {' '.join(row)}{extra}\n"


def bits(value):
    return struct.pack(">d", value).hex()


@pytest.mark.parametrize("case", PACKET["valid"], ids=lambda c: c["id"])
def test_finite_bits_in_every_binding_field(case):
    for field in range(5):
        text = source(case["token"], field)
        frame = chematic.from_extxyz(text)
        value = (frame["coords"][0][field] if field < 3 else
                 frame["lattice"][4] if field == 3 else frame["properties"]["force"][0][1])
        assert bits(value) == case["bits"]
        if field < 3:
            assert bits(chematic.from_xyz(text)[1][0][field]) == case["bits"]
        if field == 0:
            output = chematic.write_xyz_checked(
                ["C"], frame["coords"], notation="scientific", decimal_places=17
            )
            assert bits(chematic.from_xyz(output)[1][0][0]) == case["bits"]


@pytest.mark.parametrize("case", PACKET["invalid"], ids=lambda c: c["id"])
def test_invalid_scalars_are_refused_in_all_numeric_fields(case):
    for field in range(5):
        text = source(case["token"], field)
        with pytest.raises(ValueError):
            chematic.from_extxyz(text)
        if field < 3:
            with pytest.raises(ValueError):
                chematic.from_xyz(text)


def test_extxyz_all_frames_stops_on_a_bad_real():
    text = source("1", 0) + source("1e-9999", 0) + source("2", 0)
    with pytest.raises(ValueError, match="underflows"):
        chematic.from_extxyz_all(text)


@pytest.mark.parametrize("extended", [False, True])
def test_batch_readers_account_for_bad_frame_and_continue(tmp_path, extended):
    path = tmp_path / "numbers.xyz"
    path.write_text(source("-0e-999999", 0) + source("1e-9999", 0) + source("5e-324", 0),
                    encoding="utf-8")
    reader = chematic.iter_extxyz_batched if extended else chematic.iter_xyz_batched
    stream = reader(str(path), batch_size=2)
    frames = [frame for batch in stream for frame in batch]
    assert [bits(f["coords"][0][0]) for f in frames] == ["8000000000000000", "0000000000000001"]
    manifest = json.loads(stream.manifest_json())
    assert manifest["frames_seen"] == 3
    assert manifest["rejected_frames"] == 1
    assert manifest["status"] == "complete"


def test_non_numeric_string_and_raw_info_are_preserved():
    frame = chematic.from_extxyz(
        "1\nProperties=species:S:1:pos:R:3:label:S:1 energy=1e-9999\nC 0 0 0 1e-9999\n"
    )
    assert frame["properties"]["label"] == [["1e-9999"]]
    assert frame["info"] == {"energy": "1e-9999"}
