#!/usr/bin/env python3
"""Cross-read original synthetic XBS scenes with hofmann 0.21.1 (no rendering)."""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import random
import subprocess
import tempfile

import numpy as np
from hofmann import StructureScene
from chematic import XbsSpecies, XbsBondRule, write_xbs


def validate() -> dict:
    if version("hofmann") != "0.21.1":
        raise RuntimeError("requires pinned hofmann 0.21.1")
    rng = random.Random(2576)
    scenes = 0
    atoms = 0
    styles_checked = 0
    rules_checked = 0
    digest = hashlib.sha256()
    labels_pool = ["C", "O2", "Site_1"]
    styles = [
        XbsSpecies(label, 0.4 + i * 0.3, i * 0.5) for i, label in enumerate(labels_pool)
    ]
    rules = [
        XbsBondRule("C", "O2", 0.25, 2.5, 0.15, 0.7),
        XbsBondRule("Site_1", "Site_1", 0, 3.1, 0.2, 0),
    ]
    with tempfile.TemporaryDirectory(prefix="chematic-xbs-crossread-") as tmp:
        path = Path(tmp) / "scene.bs"
        for size in (1, 2, 3, 7, 64, 257):
            for scale in (1e-200, 1.0, 100.0, 1e6):
                points = [
                    [rng.uniform(-scale, scale) for _ in range(3)] for _ in range(size)
                ]
                points[0] = [
                    -0.0,
                    float.fromhex("0x0.0000000000001p-1022"),
                    1.2345678901234567,
                ]
                labels = [labels_pool[i % 3] for i in range(size)]
                for reversed_order in (False, True):
                    q = list(reversed(points)) if reversed_order else points
                    names = list(reversed(labels)) if reversed_order else labels
                    for include_rules in (False, True):
                        selected_rules = rules if include_rules else []
                        text = write_xbs(names, q, species=styles, bonds=selected_rules)
                        digest.update(text.encode("ascii"))
                        path.write_text(text, encoding="ascii")
                        scene = StructureScene.from_xbs(path)
                        assert tuple(scene.species) == tuple(names)
                        actual = np.asarray(scene.frames[0].coords, dtype=np.float64)
                        expected = np.asarray(q, dtype=np.float64)
                        assert actual.tobytes() == expected.tobytes()
                        for style in styles:
                            actual_style = scene.atom_styles[style.label]
                            assert (actual_style.radius, actual_style.colour) == (
                                style.radius,
                                style.gray,
                            )
                        assert len(scene.bond_specs) == len(selected_rules)
                        for got, want in zip(scene.bond_specs, selected_rules):
                            assert got.species == (want.species1, want.species2)
                            assert (
                                got.min_length,
                                got.max_length,
                                got.radius,
                                got.colour,
                            ) == (
                                want.min_distance,
                                want.max_distance,
                                want.radius,
                                want.gray,
                            )
                        scenes += 1
                        atoms += size
                        styles_checked += len(styles)
                        rules_checked += len(selected_rules)
        # One maximum-size scene checks the public bound without bond detection.
        text = write_xbs(["C"] * 10000, [(0.0, 0.0, 0.0)] * 10000, species=[styles[0]])
        path.write_text(text, encoding="ascii")
        assert StructureScene.from_xbs(path).frames[0].coords.shape == (10000, 3)
        digest.update(text.encode("ascii"))
        xyz = "3\nOriginal asymmetric three-site geometry\nO 0 0 0\nH 0.92 0.13 0\nH -0.2 0.88 0.07\n"
        ob = subprocess.run(
            ["/opt/homebrew/bin/obabel", "-ixyz", "-obs"],
            input=xyz,
            text=True,
            capture_output=True,
            check=True,
        )
        listing = subprocess.check_output(
            ["/opt/homebrew/bin/obabel", "-L", "formats"], text=True
        )
        return {
            "schema_version": 1,
            "scope": "static gray-style XBS writer; independent file parser only",
            "python": platform.python_version(),
            "numpy": np.__version__,
            "hofmann": version("hofmann"),
            "seed": 2576,
            "ordinary_scenes": scenes,
            "atom_rows": atoms,
            "style_checks": styles_checked,
            "rule_checks": rules_checked,
            "max_size_scene_atoms": 10000,
            "coordinate_check": "binary64 byte equality, including negative zero and smallest subnormal",
            "output_sha256": digest.hexdigest(),
            "openbabel": {
                "version": subprocess.check_output(
                    ["/opt/homebrew/bin/obabel", "-V"], text=True
                ).strip(),
                "input_xyz": xyz,
                "bs_stdout": ob.stdout,
                "bs_stderr": ob.stderr,
                "has_xbs_format": any(
                    line.split(" -- ")[0].strip() == "xbs"
                    for line in listing.splitlines()
                ),
                "bs_is_xbs_atom_records": ob.stdout.startswith("atom "),
            },
            "limits": [
                "no original XBS GUI execution",
                "no rendering or bond-detection validation",
                "no chemical topology, charge or unit conversion",
                "no native build or performance claim",
            ],
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(validate(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")
