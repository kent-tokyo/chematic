"""Preserve declared atom-column kinds through the Python dictionary bridge."""
import json
from pathlib import Path
import struct

import chematic
import pytest

CASES = json.loads(
    (Path(__file__).parents[3] / "validation/extxyz_python_types/cases.json")
    .read_text(encoding="utf-8")
)["valid"]


def exact(value):
    if isinstance(value, float):
        return (float, struct.pack(">d", value))
    if isinstance(value, (list, tuple)):
        return [exact(v) for v in value]
    if isinstance(value, dict):
        return [(k, exact(v)) for k, v in value.items()]
    return (type(value), value)


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_complete_dictionary_roundtrip(case):
    before = chematic.from_extxyz(case["text"])
    after = chematic.from_extxyz(chematic.to_extxyz(**before))
    for name in ("coords", "lattice", "properties", "property_schema", "info"):
        assert exact(after[name]) == exact(before[name])
    assert after["mol"].atom_count == before["mol"].atom_count
    if not (before["lattice"] or before["properties"] or before["info"]):
        assert after["comment"] == before["comment"]


def test_integer_and_empty_frame_oracles():
    frame = chematic.from_extxyz(CASES[1]["text"])
    assert frame["properties"]["id"] == [
        [9007199254740993], [-(1 << 63)], [(1 << 63) - 1]
    ]
    assert frame["property_schema"] == {"id": ("I", 1)}
    frame = chematic.from_extxyz(CASES[5]["text"])
    assert frame["properties"] == {"id": [], "tag": [], "mask": [], "vec": []}
    assert frame["property_schema"] == {
        "id": ("I", 2), "tag": ("S", 1), "mask": ("L", 3), "vec": ("R", 4)
    }


def write(values, schema):
    return chematic.to_extxyz(
        chematic.from_smiles("[C]"), [[0, 0, 0]],
        properties={"data": values}, property_schema=schema,
    )


@pytest.mark.parametrize("schema", [
    {}, {"other": ("R", 1)}, {"data": ("R", 1), "extra": ("I", 1)},
    {"data": "R"}, {"data": ("r", 1)}, {"data": ("X", 1)},
    {"data": ("R", 0)}, {"data": ("R", -1)}, {"data": ("R", True)},
    {"data": ("R", 1.0)}, {"data": ("R", "1")}, {"data": ("R", 1 << 100)},
    {"data": ("R", 2)}, {"data": None},
    {1: ("R", 1)}, {"data": ("R",)}, {"data": ("R", 1, "extra")},
])
def test_inconsistent_or_invalid_schema_is_refused(schema):
    with pytest.raises(ValueError):
        write([[1.0]], schema)


@pytest.mark.parametrize("kind, value", [
    ("I", 1.0), ("I", True), ("I", "1"), ("I", 1 << 63), ("I", -(1 << 63)-1),
    ("L", 0), ("L", "T"), ("L", 1.0),
    ("S", 1), ("S", True), ("S", ""), ("S", "two tokens"), ("S", 'a"b'),
    ("S", "a\\b"), ("S", "a=b"), ("S", "a\nb"),
    ("R", True), ("R", "1.0"), ("R", 9007199254740993),
    ("R", float("nan")), ("R", float("inf")),
])
def test_invalid_or_lossy_scalar_is_refused(kind, value):
    with pytest.raises(ValueError):
        write([[value]], {"data": (kind, 1)})


def test_legacy_real_columns_accept_exact_integers_and_refuse_rounding():
    for value in (0, 1, 1 << 53, 1 << 60):
        frame = chematic.from_extxyz(write([[value]], None))
        assert frame["property_schema"] == {"data": ("R", 1)}
        assert type(frame["properties"]["data"][0][0]) is float
        assert frame["properties"]["data"][0][0] == value
    with pytest.raises(ValueError, match="exactly representable"):
        write([[9007199254740993]], None)


def test_integer_protocol_cannot_bypass_real_rounding_check():
    class WideInteger:
        def __float__(self):
            return float(1 << 53)

        def __index__(self):
            return (1 << 53) + 1

    with pytest.raises(ValueError, match="exactly representable"):
        write([[WideInteger()]], None)


def test_row_counts_and_schema_without_properties_are_checked():
    for rows in ([], [[1], [2]], [[1, 2]], "text", ["text"]):
        with pytest.raises(ValueError):
            write(rows, {"data": ("I", 1)})
    with pytest.raises(ValueError, match="exactly match"):
        chematic.to_extxyz(chematic.from_smiles("C"), [[0, 0, 0]],
                           property_schema={"missing": ("I", 1)})


def test_duplicate_reader_columns_are_refused_before_dict_overwrite():
    for count, rows in ((1, "C 0 0 0 1 2\n"), (0, "")):
        text = (f"{count}\nProperties=species:S:1:pos:R:3:id:I:1:id:R:1\n" + rows)
        with pytest.raises(ValueError, match="duplicate property"):
            chematic.from_extxyz(text)


def test_trajectory_dictionaries_keep_per_frame_schema_and_comments():
    frames = chematic.from_extxyz_all("".join(c["text"] for c in CASES))
    assert len(frames) == len(CASES)
    output = "".join(chematic.to_extxyz(**f) for f in frames)
    after = chematic.from_extxyz_all(output)
    assert [exact(f["property_schema"]) for f in after] == [
        exact(f["property_schema"]) for f in frames
    ]
    assert [exact(f["properties"]) for f in after] == [exact(f["properties"]) for f in frames]


def test_existing_five_positional_arguments_still_work():
    text = chematic.to_extxyz(chematic.from_smiles("C"), [[0, 0, 0]], None,
                             {"weight": [[0.25]]}, {"note": "control"})
    assert chematic.from_extxyz(text)["properties"] == {"weight": [[0.25]]}


def test_streamed_frames_keep_types_and_can_be_written_back(tmp_path):
    path = tmp_path / "typed.extxyz"
    path.write_text("".join(c["text"] for c in CASES), encoding="utf-8")
    stream = chematic.iter_extxyz_batched(str(path), batch_size=3)
    frames = [f for batch in stream for f in batch]
    assert len(frames) == len(CASES)
    for frame in frames:
        result = chematic.from_extxyz(chematic.to_extxyz(**frame))
        assert exact(result["properties"]) == exact(frame["properties"])
        assert result["property_schema"] == frame["property_schema"]
