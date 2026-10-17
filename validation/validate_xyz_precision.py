#!/usr/bin/env python3
"""Original numeric profiles, independent decimal oracle and chemistry cross-reads."""
from __future__ import annotations

import hashlib
import json
import math
import random
import struct
import subprocess

import chematic
import rdkit
from rdkit import Chem, rdBase


def bits(value):
    return struct.pack(">d", value)


def rows(text):
    return [line.split() for line in text.splitlines()[2:]]


def main():
    cases = 0
    decimals = 0
    binary_checks = 0
    crossreads = 0
    rdkit_scientific_refusals = 0
    digest = hashlib.sha256()
    rng = random.Random(1693)
    ordinary = [[rng.uniform(-100, 100) for _ in range(3)] for _ in range(64)]
    for notation in ["fixed", "scientific"]:
        for places in [0, 1, 5, 6, 8, 12, 16, 17]:
            for width in [0, 1, 20, 64]:
                text = chematic.write_xyz_checked(
                    ["C"] * len(ordinary),
                    ordinary,
                    comment="original numeric corpus",
                    notation=notation,
                    decimal_places=places,
                    min_width=width,
                )
                for row, expected in zip(rows(text), ordinary):
                    assert row[0] == "C" and len(row) == 4
                    for token, value in zip(row[1:], expected):
                        oracle = format(
                            value, f'.{places}{"f" if notation=="fixed" else "e"}'
                        )
                        if notation == "fixed":
                            assert token == oracle
                        else:
                            a, b = token.split("e")
                            c, d = oracle.split("e")
                            assert a == c and int(b) == int(d)
                        decimals += 1
                digest.update(text.encode())
                cases += 1
    values = [
        0.0,
        -0.0,
        5e-324,
        2.2250738585072014e-308,
        1.7976931348623157e308,
        -1.7976931348623157e308,
        1 / 3,
    ]
    while len(values) < 2048:
        value = struct.unpack(">d", rng.getrandbits(64).to_bytes(8, "big"))[0]
        if math.isfinite(value):
            values.append(value)
    for places in [16, 17]:
        text = chematic.write_xyz_checked(
            ["Ne"] * len(values),
            [[v, -v, v] for v in values],
            notation="scientific",
            decimal_places=places,
            min_width=0,
        )
        for row, value in zip(rows(text), values):
            for token, target in zip(row[1:], [value, -value, value]):
                assert bits(float(token)) == bits(target)
                binary_checks += 1
        digest.update(text.encode())
    positions = [
        [0.1234567890123, 0, 0],
        [1.0806913569035, 0, 0],
        [-0.1165432109877, 0.9273456789012, 0],
    ]
    max_crossread_error = 0.0
    obabel = "/opt/homebrew/bin/obabel"
    for notation in ["fixed", "scientific"]:
        for places in [5, 6, 12, 16]:
            text = chematic.write_xyz_checked(
                ["O", "H", "H"],
                positions,
                comment="original water",
                notation=notation,
                decimal_places=places,
            )
            expected = [[float(v) for v in row[1:]] for row in rows(text)]
            with rdBase.BlockLogs():
                reference = Chem.MolFromXYZBlock(text)
            if notation == "scientific":
                assert reference is None
                rdkit_scientific_refusals += 1
            else:
                assert reference is not None
                assert [a.GetSymbol() for a in reference.GetAtoms()] == ["O", "H", "H"]
                conf = reference.GetConformer()
                for i, target in enumerate(expected):
                    p = conf.GetAtomPosition(i)
                    assert [bits(v) for v in [p.x, p.y, p.z]] == [
                        bits(v) for v in target
                    ]
                crossreads += 1
            # EXYZ preserves numeric columns on Open Babel output, unlike its
            # ordinary five-place XYZ writer; no cell/properties are supplied.
            p = subprocess.run(
                [obabel, "-ixyz", "-oexyz", "-ab"],
                input=text,
                capture_output=True,
                text=True,
                check=True,
            )
            actual = [line.split() for line in p.stdout.splitlines()[2:5]]
            assert [r[0] for r in actual] == ["O", "H", "H"]
            for row, target in zip(actual, expected):
                error = max(abs(float(a) - b) for a, b in zip(row[1:4], target))
                assert error <= 1e-5
                max_crossread_error = max(error, max_crossread_error)
            crossreads += 1
    source = "1\noriginal source\nNe 1.2345678901234 -2.3456789012345 0.0000001234567\n"
    p = subprocess.run(
        [obabel, "-ixyz", "-oxyz", "-ab"],
        input=source,
        capture_output=True,
        text=True,
        check=True,
    )
    legacy = chematic.from_smiles("[Ne]").to_xyz(
        [[1.2345678901234, -2.3456789012345, 0.0000001234567]], "legacy"
    )
    result = {
        "rdkit_version": rdkit.__version__,
        "openbabel_version": subprocess.run(
            [obabel, "-V"], capture_output=True, text=True, check=True
        ).stdout.strip(),
        "format_profiles": cases,
        "independent_decimal_checks": decimals,
        "exact_scientific_f64_checks": binary_checks,
        "crossreads": crossreads,
        "rdkit_scientific_refusals": rdkit_scientific_refusals,
        "openbabel_max_coordinate_error": max_crossread_error,
        "corpus_sha256": digest.hexdigest(),
        "baseline": {
            "openbabel_xyz_output": p.stdout,
            "legacy_chematic_to_xyz": legacy,
        },
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
