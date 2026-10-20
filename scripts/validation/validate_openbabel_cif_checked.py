#!/usr/bin/env python3
"""Validate checked CIF source bindings and record bounded Open Babel probes."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

import chematic

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "validation/openbabel_issues/poscar/singular-cell.POSCAR"


def probe(args, timeout=10, input_text=None):
    try:
        result = subprocess.run(args, input=input_text, capture_output=True, text=True, timeout=timeout)
        return dict(args=args, returncode=result.returncode, stdout=result.stdout,
                    stderr=result.stderr, timed_out=False)
    except subprocess.TimeoutExpired as exc:
        return dict(args=args, timed_out=True, timeout_seconds=timeout,
                    stdout=(exc.stdout or b"").decode(), stderr=(exc.stderr or b"").decode())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    mol = chematic.from_smiles("[Si]")
    cell = dict(a=4., b=6., c=8., alpha=90., beta=90., gamma=90.)
    rejected = []
    for key in cell:
        values = [0., -1., math.nan, math.inf, -math.inf]
        if key in ("alpha", "beta", "gamma"):
            values += [180., 270.]
        for value in values:
            candidate = dict(cell, **{key: value})
            try:
                chematic.write_cif_checked(mol, [[0., 0., 0.]], candidate)
            except ValueError as exc:
                rejected.append(dict(field=key, value=repr(value), error=str(exc)))
            else:
                raise AssertionError(candidate)
    round_trips = []
    for candidate in [None, cell, dict(cell, alpha=75., beta=83., gamma=64.),
                      dict(cell, a=1e-6, b=2e-6, c=3e-6)]:
        coords = [[.125123456789, -.375987654321, 1.234567891234]]
        text = chematic.write_cif_checked(mol, coords, candidate)
        result = chematic.parse_cif(text)
        error = max(abs(x-y) for x, y in zip(coords[0], result["coords"][0]))
        assert result["cell"] == candidate and error < 1e-12
        round_trips.append(dict(cell=candidate, max_absolute_coordinate_error=error))
    poscar_probe = probe(["obabel", "-iPOSCAR", str(FIXTURE), "-ocif"])
    no_bond_probe = probe(["obabel", "-iPOSCAR", str(FIXTURE), "-ocif", "-ab"])
    text = chematic.write_cif_checked(mol, [[.5, 1.5, 2.5]], cell)
    cross_read = probe(["obabel", "-icif", "-oxyz", "-ab"], input_text=text)
    cross_read["input_sha256"] = hashlib.sha256(text.encode()).hexdigest()
    cross_read_passed = False
    if not cross_read["timed_out"] and cross_read["returncode"] == 0:
        rows = cross_read["stdout"].splitlines()
        assert int(rows[0]) == 1 and rows[2].split()[0] == "Si", rows
        actual = [float(value) for value in rows[2].split()[1:4]]
        assert max(abs(a-b) for a,b in zip(actual, [.5, 1.5, 2.5])) < 1e-6
        cross_read_passed = True
    control = probe(["obabel", "-ixyz", "-ocif"], input_text="1\ncontrol\nSi 0 0 0\n")
    module = Path(chematic.__file__).resolve()
    report = dict(
        issue="https://github.com/openbabel/openbabel/issues/2818",
        source_baseline=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        artifact_lane="local source extension; uncommitted working tree; not a release",
        source_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in [ROOT / "crates/chematic-mol/src/cif.rs", ROOT / "crates/chematic-mol/src/cif/writer.rs",
                                    ROOT / "crates/chematic-py/src/formats/cif.rs"]},
        python=platform.python_version(), module_path=str(module),
        module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),
        fixture_sha256=hashlib.sha256(FIXTURE.read_bytes()).hexdigest(),
        openbabel_version=probe(["obabel", "-V"]),
        singular_poscar=poscar_probe, singular_poscar_bonding_disabled=no_bond_probe,
        invalid_cell_cases_rejected=rejected, round_trips=round_trips,
        valid_cif_openbabel_cross_read=cross_read,
        valid_cif_openbabel_cross_read_passed=cross_read_passed,
        xyz_to_cif_control=control,
        comparison_conclusion=("valid CIF cross-read passed" if cross_read_passed else
                               "Open Babel CIF cross-read unavailable; singular-input timeout does not prove issue-specific failure"),
        scope="plain CIF atomic positions; no POSCAR reader; legacy Rust write_cif remains unchecked",
    )
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(rejected_cells=len(rejected), round_trips=len(round_trips),
                          singular_poscar_timed_out=poscar_probe["timed_out"],
                          cross_read_passed=cross_read_passed)))


if __name__ == "__main__":
    main()
