"""Token preservation and explicit failure boundaries for Open Babel #2668."""
import json
import re
from pathlib import Path

import chematic
import pytest

CASES = json.loads((Path(__file__).parents[3] / "validation/openbabel_issues/mol2-kekule/cases.json").read_text(encoding="utf-8"))["cases"]


def assert_token_only(before, after):
    source_lines = before.splitlines(keepends=True)
    output_lines = after.splitlines(keepends=True)
    assert len(source_lines) == len(output_lines)
    in_bonds = False
    changed = 0
    for source, output in zip(source_lines, output_lines):
        if source.strip().upper().startswith("@<TRIPOS>"):
            in_bonds = source.strip().upper() == "@<TRIPOS>BOND"
        fields = list(re.finditer(r"\S+", source))
        if in_bonds and len(fields) >= 4 and fields[3].group() == "ar":
            start, end = fields[3].span()
            assert output[start] in "12"
            assert output == source[:start] + output[start] + source[end:]
            changed += 1
        else:
            assert source == output
    return changed


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_only_ar_tokens_change_and_repeat_is_an_exact_noop(case, newline):
    source = case["text"].replace("\r\n", "\n").replace("\n", newline)
    result = chematic.kekulize_mol2_record(source)
    assert assert_token_only(source, result) == case["aromatic_bonds"]
    assert chematic.kekulize_mol2_record(result) == result
    assert chematic.roundtrip_mol2_record(source) == source


def ring(n, element="C"):
    atoms = "".join(f"{i} {element}{i} 0 0 0 {element}.ar\n" for i in range(1, n + 1))
    bonds = "".join(f"{i} {i} {i % n + 1} ar\n" for i in range(1, n + 1))
    return f"@<TRIPOS>MOLECULE\nring\n{n} {n}\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n{atoms}@<TRIPOS>BOND\n{bonds}"


@pytest.mark.parametrize("source", [
    ring(5), ring(4), ring(6, "O"), ring(6, "B"),
    ring(6, "N") + "@<TRIPOS>UNITY_ATOM_ATTR\n1 1\ncharge 1\n",
    ring(6).replace("1 1 2 ar", "1 1 2 un"),
    ring(6).replace("1 1 2 ar", "1 1 2 nc"),
    ring(6).replace("1 1 2 ar", "1 1 99 ar"),
    ring(513), ring(6) + ring(6),
    ring(5).replace("1 C1 0 0 0 C.ar", "1 N1 0 0 0 N.pl3"),
])
def test_unsupported_or_invalid_input_raises_without_partial_output(source):
    with pytest.raises(ValueError):
        chematic.kekulize_mol2_record(source)


def test_input_byte_bound_precedes_parse():
    with pytest.raises(ValueError):
        chematic.kekulize_mol2_record("#" + "x" * (8 * 1024 * 1024) + "\n" + ring(6))


def test_sparse_ids_and_interleaved_comments_are_retained():
    source = (Path(__file__).parents[3] / "validation/openbabel_issues/mol2-kekule/rich-source.mol2").read_bytes().decode()
    assert assert_token_only(source, chematic.kekulize_mol2_record(source)) == 6
