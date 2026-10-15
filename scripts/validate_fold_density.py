#!/usr/bin/env python3
"""Corrected #9715 loop around pinned RDKit OR-fold primitives, plus set oracle."""
from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import chematic
from chematic.rdkit_compat import DataStructs as Adapter, ExplicitBitVect as AdapterBits
from rdkit import DataStructs, rdBase

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "validation/rdkit_issues/fold-density/fixtures.json"
OUTPUT = ROOT / "validation/results/rdkit-fold-density-2026-10-15.json"


def normalized(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def packed(width, on):
    data = bytearray((width + 7) // 8)
    for bit in on:
        data[bit // 8] |= 1 << (bit % 8)
    return bytes(data)


def rdkit_reference(width, on, density, minimum, vector_type):
    fp = vector_type(width)
    for bit in on:
        fp.SetBit(bit)
    steps = 0
    # The reported upstream correction is explicit here; the installed
    # FoldToTargetDensity wrapper still has the bug and is not this oracle.
    while fp.GetNumOnBits() / fp.GetNumBits() < density and fp.GetNumBits() // 2 >= minimum:
        old_width, old_on = fp.GetNumBits(), list(fp.GetOnBits())
        try:
            fp = DataStructs.FoldFingerprint(fp, 2)
        except ValueError as error:
            assert old_width == 2 and str(error) == "invalid fold factor", (old_width, str(error))
            return None, steps
        assert fp.GetNumBits() == old_width // 2
        assert list(fp.GetOnBits()) == sorted({bit % (old_width // 2) for bit in old_on})
        steps += 1
    return (fp.GetNumBits(), list(fp.GetOnBits()), steps), steps


def set_reference(width, on, density, minimum):
    bits = set(on)
    steps = 0
    while len(bits) / width < density and width // 2 >= minimum:
        width //= 2
        bits = {bit % width for bit in bits}
        steps += 1
    return width, sorted(bits), steps


def cases(manifest):
    for fixture in manifest["records"]:
        p = fixture["payload"]
        yield fixture["id"], p["width"], p["on_bits"], p["density"], p["min_length"]
    for width in range(1, 9):
        for mask in range(1 << width):
            bits = [bit for bit in range(width) if mask & (1 << bit)]
            for density in [0., .125, .3, .5, .75, 1.]:
                for minimum in [1, 2, 3, 4, 8, 64]:
                    yield f"exhaust-{width}-{mask}-{density}-{minimum}", width, bits, density, minimum
    for width in [9, 31, 63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 1000, 1023, 1024, 2048, 4096]:
        for spacing in [1, 2, 3, 7, 31, 65, width]:
            bits = sorted(set(range(0, width, spacing)) | {width - 1})
            for density in [.3, .5, 1.]:
                for minimum in [1, 32, 64, width // 2, width + 1]:
                    yield f"boundary-{width}-{spacing}-{density}-{minimum}", width, bits, density, minimum


def check_case(width, on, density, minimum):
    expected = set_reference(width, on, density, minimum)
    completed = primitive_steps = boundary_errors = 0
    for vector_type in [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect]:
        reference, steps = rdkit_reference(width, on, density, minimum, vector_type)
        primitive_steps += steps
        if reference is None:
            assert expected[0] == 1 and minimum == 1
            boundary_errors += 1
        else:
            assert reference == expected
            completed += 1
    result = chematic.fold_to_target_density(packed(width, on), density, minimum, bit_width=width)
    final_width, final_on, steps = expected
    assert result["fingerprint"] == packed(final_width, final_on)
    assert result["bit_width"] == final_width
    assert result["original_bits"] == width
    assert result["folds"] == steps
    assert result["on_bit_count"] == len(final_on)
    assert result["density"] == len(final_on) / final_width
    reached = len(final_on) / final_width >= density
    assert result["target_reached"] == reached
    assert result["stop_reason"] == ("target_reached" if reached else "minimum_length")
    adapter = AdapterBits(width)
    for bit in on:
        adapter.SetBit(bit)
    folded = Adapter.FoldToTargetDensity(adapter, density, minimum)
    assert (len(folded), folded.GetOnBits()) == (final_width, final_on)
    assert (len(adapter), adapter.GetOnBits()) == (width, on)
    return result, completed, primitive_steps, boundary_errors


def main():
    assert rdBase.rdkitVersion == "2026.09.1", rdBase.rdkitVersion
    manifest = json.loads(FIXTURES.read_text())
    aggregate = hashlib.sha256()
    count = folds = reached = primitive_cases = primitive_steps = boundary_errors = 0
    for identifier, width, on, density, minimum in cases(manifest):
        result, completed, steps, errors = check_case(width, on, density, minimum)
        primitive_cases += completed
        primitive_steps += steps
        boundary_errors += errors
        aggregate.update(normalized([identifier, width, on, density, minimum,
                                     result["bit_width"], result["fingerprint"].hex()]) + b"\n")
        count += 1
        folds += result["folds"]
        reached += result["target_reached"]
    reproductions = []
    for fixture in manifest["records"]:
        p = fixture["payload"]
        assert hashlib.sha256(normalized(p)).hexdigest() == fixture["sha256"]
        fp = DataStructs.ExplicitBitVect(p["width"])
        for bit in p["on_bits"]:
            fp.SetBit(bit)
        original = DataStructs.FoldToTargetDensity(fp, p["density"], p["min_length"])
        corrected, _, _, _ = check_case(p["width"], p["on_bits"], p["density"], p["min_length"])
        expected = fixture["expected"]
        assert corrected["bit_width"] == expected["width"]
        assert corrected["fingerprint"] == packed(expected["width"], expected["on_bits"])
        reproductions.append(dict(id=fixture["id"], input=p, expected=expected,
            upstream_original=dict(width=original.GetNumBits(), on_bits=list(original.GetOnBits())),
            differs_from_original=(original.GetNumBits(),list(original.GetOnBits())) != (expected["width"],expected["on_bits"])))
    result = dict(schema_version=1, issue=dict(number=9715, state="OPEN", verified_date="2026-10-15",
        updated_at="2026-10-15T08:57:28Z",url="https://github.com/rdkit/rdkit/issues/9715"),
        oracle=dict(rdkit=rdBase.rdkitVersion, wrapper="Explicit corrected loop around FoldFingerprint(fp,2)",
            primitives=["ExplicitBitVect", "SparseBitVect"], independent_oracle="Set union of indices modulo floor(width/2)",
            original_wrapper_sha256=hashlib.sha256(inspect.getsource(DataStructs.FoldToTargetDensity).encode()).hexdigest()),
        fixtures=dict(path=str(FIXTURES.relative_to(ROOT)),sha256=hashlib.sha256(FIXTURES.read_bytes()).hexdigest()),
        results=dict(independent_exact_cases=count, rdkit_completed_exact_comparisons=primitive_cases,
            rdkit_successful_fold_steps=primitive_steps, rdkit_two_to_one_exceptions=boundary_errors,
            adapter_comparisons=count,
            fold_steps=folds, targets_reached=reached, unmet_targets=count-reached,mismatches=0,
            normalized_packet_sha256=aggregate.hexdigest(), bug_reproductions=sum(r["differs_from_original"] for r in reproductions)),
        reproductions=reproductions)
    OUTPUT.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result["results"],indent=2))


if __name__ == "__main__":
    main()
