#!/usr/bin/env python3
"""Cross-read independently authored static XBS scenes with hofmann 0.21.1."""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import random
import struct
import tempfile

import numpy as np
from hofmann import StructureScene
from chematic import read_xbs


def packed(rows):
    return b"".join(struct.pack("!d", value) for row in rows for value in row)


def check_oracle(path, text, labels, points, styles, rules):
    path.write_bytes(text.encode("ascii"))
    reference = StructureScene.from_xbs(path)
    assert tuple(reference.species) == tuple(labels)
    assert packed(np.asarray(reference.frames[0].coords)) == packed(points)
    assert set(reference.atom_styles) == {label for label, _, _ in styles}
    for label, radius, gray in styles:
        style = reference.atom_styles[label]
        assert (style.radius, style.colour) == (radius, gray)
    assert len(reference.bond_specs) == len(rules)
    for got, want in zip(reference.bond_specs, rules):
        assert got.species == want[:2]
        assert (got.min_length, got.max_length, got.radius, got.colour) == want[2:]


def validate():
    assert version("hofmann") == "0.21.1"
    rng = random.Random(25761019)
    styles = [("Site_1", 1.25, 0.75), ("C", 0.5, 0.0), ("O2", 0.75, 1.0)]
    rules = [
        ("C", "O2", 0.0, 2.5, 0.125, 0.5),
        ("Site_1", "Site_1", 0.25, 3.0, 0.25, 0.0),
    ]
    digest = hashlib.sha256()
    scenes = atoms = 0
    with tempfile.TemporaryDirectory(prefix="chematic-xbs-input-") as tmp:
        path = Path(tmp) / "original.bs"
        for size in (1, 2, 3, 7, 64, 257):
            for scale in (1e-200, 1.0, 1e6):
                points = [
                    [rng.uniform(-scale, scale) for _ in range(3)] for _ in range(size)
                ]
                points[0] = [
                    -0.0,
                    float.fromhex("0x0.0000000000001p-1022"),
                    1.2345678901234567,
                ]
                labels = [styles[i % 3][0] for i in range(size)]
                for reverse in (False, True):
                    rows = list(reversed(points)) if reverse else points
                    names = list(reversed(labels)) if reverse else labels
                    for fmt in (".17g", ".16e"):
                        for newline in ("\n", "\r\n"):

                            def numbers(values):
                                return "\t".join(format(v, fmt) for v in values)

                            lines = ["* Original generated static scene", ""]
                            lines += [
                                f"bonds\t{a}\t{b}\t{numbers(values)}"
                                for a, b, *values in rules
                            ]
                            lines += [f"atom\t{names[0]}\t{numbers(rows[0])}"]
                            lines += [
                                f"spec\t{label}\t{numbers((radius, gray))}"
                                for label, radius, gray in styles
                            ]
                            lines += ["* remaining atoms follow"]
                            lines += [
                                f"atom\t{label}\t{numbers(row)}"
                                for label, row in zip(names[1:], rows[1:])
                            ]
                            text = newline.join(lines) + newline
                            r = read_xbs(text.encode("ascii"))
                            assert r.source_text == text == r.to_xbs()
                            assert r.labels == tuple(names)
                            assert packed(r.positions) == packed(rows)
                            assert [
                                (s.label, s.radius, s.gray) for s in r.species
                            ] == styles
                            assert [
                                (
                                    b.species1,
                                    b.species2,
                                    b.min_distance,
                                    b.max_distance,
                                    b.radius,
                                    b.gray,
                                )
                                for b in r.bonds
                            ] == rules
                            assert r.atom_lines == tuple(
                                i
                                for i, line in enumerate(lines, 1)
                                if line.startswith("atom")
                            )
                            for output in (text, r.to_xbs(preserve_source=False)):
                                check_oracle(path, output, names, rows, styles, rules)
                                again = read_xbs(output)
                                assert again.labels == r.labels
                                assert packed(again.positions) == packed(r.positions)
                                assert (
                                    again.species == r.species
                                    and again.bonds == r.bonds
                                )
                                digest.update(output.encode("ascii"))
                            scenes += 1
                            atoms += size
        boundary = "atom C 0 0 0\n" * 10000 + "spec C 1 0\n"
        r = read_xbs(boundary)
        check_oracle(path, boundary, r.labels, r.positions, [("C", 1.0, 0.0)], [])
        digest.update(boundary.encode("ascii"))
    return {
        "scope": "static gray-style XBS input; original and normalized output cross-read, no viewer execution",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "hofmann": version("hofmann"),
        "seed": 25761019,
        "scenes": scenes,
        "atom_rows": atoms,
        "independent_crossreads": scenes * 2 + 1,
        "maximum_atom_scene": 10000,
        "coordinate_check": "binary64 byte equality including negative zero and minimum subnormal",
        "style_checks": scenes * 3 * 2 + 1,
        "rule_checks": scenes * 2 * 2,
        "packet_sha256": digest.hexdigest(),
        "source_sha256": {
            p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in (
                "crates/chematic-py/python/chematic/_xbs_read.py",
                "crates/chematic-py/tests/test_xbs_input.py",
                "validation/validate_xbs_input.py",
            )
        },
        "limitations": [
            "strict subset, not all XBS commands",
            "no chemistry or units inferred",
            "no GUI/render/bond-perception validation",
            "no native rebuild or release claim",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = validate()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps({k: v for k, v in result.items() if k != "source_sha256"}, indent=2)
    )
