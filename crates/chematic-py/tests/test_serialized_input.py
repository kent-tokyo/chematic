"""Only serialized format text accepts immutable UTF-8 bytes (RDKit #9601)."""
import ast
import json
from pathlib import Path

import chematic
import pytest

ROOT = Path(__file__).resolve().parents[3]
CASES = json.loads((ROOT / "validation/rdkit_issues/python-serialized-input/cases.json").read_text())["cases"]


def snapshot(value):
    if isinstance(value, chematic.Mol):
        return {"smiles": value.smiles, "formula": value.formula,
                "atoms": value.atom_count, "bonds": value.bond_table}
    if isinstance(value, (list, tuple)):
        return [snapshot(item) for item in value]
    if isinstance(value, dict):
        return {key: snapshot(item) for key, item in value.items()}
    return value


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["function"])
def test_str_bytes_positional_and_existing_keyword_have_identical_results(case):
    reader = getattr(chematic, case["function"])
    expected = snapshot(reader(case["text"]))
    encoded = case["text"].encode("utf-8")
    assert snapshot(reader(encoded)) == expected
    assert snapshot(reader(**{case["keyword"]: encoded})) == expected
    assert snapshot(reader(**{case["keyword"]: case["text"]})) == expected


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["function"])
@pytest.mark.parametrize("payload", [b"\xff", b"C\xe2\x82", b"\xc0\xaf", b"\xed\xa0\x80"])
def test_invalid_utf8_is_a_decode_error_even_for_boolean_validators(case, payload):
    with pytest.raises(UnicodeDecodeError) as expected:
        payload.decode("utf-8")
    with pytest.raises(UnicodeDecodeError) as actual:
        getattr(chematic, case["function"])(payload)
    assert actual.value.encoding == "utf-8"
    assert actual.value.object == payload
    assert (actual.value.start, actual.value.end) == (expected.value.start, expected.value.end)


class NoCoercion:
    def __str__(self):
        raise AssertionError("must not call __str__")

    def __bytes__(self):
        raise AssertionError("must not call __bytes__")

    def __fspath__(self):
        raise AssertionError("must not call __fspath__")


WRONG = [None, True, 1, 1.5, [], {}, (), Path("CO"), bytearray(b"CO"), memoryview(b"CO"), NoCoercion()]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["function"])
@pytest.mark.parametrize("value", WRONG)
def test_wrong_types_are_rejected_without_protocol_conversion(case, value):
    with pytest.raises(TypeError, match="str or bytes"):
        getattr(chematic, case["function"])(value)


class Text(str):
    def __str__(self):
        raise AssertionError("must read native str storage")


class Data(bytes):
    def __bytes__(self):
        raise AssertionError("must read native bytes storage")


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["function"])
def test_subclasses_are_read_without_overridden_conversion_hooks(case):
    reader = getattr(chematic, case["function"])
    expected = snapshot(reader(case["text"]))
    assert snapshot(reader(Text(case["text"]))) == expected
    assert snapshot(reader(Data(case["text"].encode()))) == expected


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["function"])
@pytest.mark.parametrize("text", ["", "not a format", "C\x00("])
def test_valid_utf8_keeps_native_parse_results_and_structured_errors(case, text):
    reader = getattr(chematic, case["function"])
    def result(value):
        try:
            return ("return", snapshot(reader(value)))
        except ValueError as error:
            return ("error", type(error), str(error),
                    tuple(getattr(error, key, None) for key in ["category", "code", "format"]))
    assert result(text) == result(text.encode())


def test_unicode_cxsmiles_labels_are_preserved():
    _, metadata = chematic.from_cxsmiles("CO |$炭素;酸素$|".encode())
    assert metadata["atom_labels"] == ["炭素", "酸素"]


def test_identifiers_paths_and_unselected_apis_keep_their_existing_contract():
    mol = chematic.from_smiles(b"CO")
    with pytest.raises(TypeError):
        mol.SetProp(b"name", "value")
    with pytest.raises(TypeError):
        mol.GetProp(b"name")
    # Filesystem paths have their own public contract and accept native bytes.
    with pytest.raises(OSError):
        chematic.iter_sdf(b"not-a-path")
    with pytest.raises(TypeError):
        chematic.from_condensed(b"CH4")


def test_stub_annotations_cover_exactly_the_selected_thirteen_inputs():
    tree = ast.parse(
        (ROOT / "crates/chematic-py/python/chematic/__init__.pyi").read_text(encoding="utf-8")
    )
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    for case in CASES:
        assert ast.unparse(functions[case["function"]].args.args[0].annotation) == "str | bytes"
    assert len(CASES) == 13
