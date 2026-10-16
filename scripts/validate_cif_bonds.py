#!/usr/bin/env python3
"""Check molecular CIF bond output with Gemmi 0.7.5 and bounded Open Babel probes.

Build the source Python extension, then provide it and Gemmi on PYTHONPATH.
This is an external oracle runner, not a dependency of the hermetic CI tests.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "validation/openbabel_issues/cif-bonds"
SOURCE = ["crates/chematic-mol/src/cif/bonds.rs", "crates/chematic-mol/src/cif.rs",
          "crates/chematic-mol/src/cif/writer.rs",
          "crates/chematic-mol/src/lib.rs", "crates/chematic-py/src/formats/cif.rs",
          "crates/chematic-py/src/formats.rs", "crates/chematic-py/python/chematic/__init__.pyi"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def probe(command, text=None):
    try:
        r = subprocess.run(command, input=text, text=True, capture_output=True, timeout=5)
        return dict(command=command, returncode=r.returncode, stdout=r.stdout, stderr=r.stderr)
    except subprocess.TimeoutExpired:
        return dict(command=command, timeout_seconds=5)


def validate_case(gemmi, chematic, case, cell):
    text = chematic.write_cif_with_bonds_checked(chematic.from_smiles(case["smiles"]), case["coords"], cell)
    block = gemmi.cif.read_string(text).sole_block()
    labels = list(block.find_values("_atom_site_label"))
    expected_labels = [f"{sym}{i+1}" for i, sym in enumerate(case["elements"])]
    assert labels == expected_labels
    assert list(block.find_values("_atom_site_type_symbol")) == case["elements"]
    prefix = "fract" if cell else "Cartn"
    columns = [list(block.find_values(f"_atom_site_{prefix}_{axis}")) for axis in "xyz"]
    checks = 2
    gemmi_cell = gemmi.UnitCell(*(cell[k] for k in ["a","b","c","alpha","beta","gamma"])) if cell else None
    for i, expected in enumerate(case["coords"]):
        actual = [float(col[i]) for col in columns]
        if gemmi_cell:
            p = gemmi_cell.orthogonalize(gemmi.Fractional(*actual))
            actual = [p.x, p.y, p.z]
        for a, b in zip(actual, expected):
            assert math.isclose(a, b, rel_tol=1e-13, abs_tol=1e-12), (a, b)
            checks += 1
    tags = ["_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2", "_geom_bond_distance",
            "_geom_bond_site_symmetry_1", "_geom_bond_site_symmetry_2", "_ccdc_geom_bond_type"]
    rows = list(zip(*(list(block.find_values(tag)) for tag in tags)))
    assert len(rows) == len(case["bonds"])
    expected_edges = {(min(a,b), max(a,b), kind) for a,b,kind in case["bonds"]}
    actual_edges = set()
    for label1, label2, distance, symmetry1, symmetry2, kind in rows:
        a, b = labels.index(label1), labels.index(label2)
        actual_edges.add((min(a,b), max(a,b), kind))
        assert symmetry1 == symmetry2 == "."
        expected = math.dist(case["coords"][a], case["coords"][b])
        assert math.isclose(float(distance), expected, rel_tol=1e-14, abs_tol=0.)
        checks += 3
    assert actual_edges == expected_edges
    checks += 2
    return dict(name=case["name"], cell=cell, atom_count=len(labels), bond_rows=rows,
                checks=checks, cif_sha256=hashlib.sha256(text.encode()).hexdigest()), text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    p.add_argument("--output", type=Path, default=ROOT / "validation/results/openbabel-cif-bonds-2026-10-16.json")
    args = p.parse_args()
    import chematic
    import gemmi
    import importlib
    assert gemmi.__version__ == "0.7.5"
    cells = [None, dict(a=10.,b=10.,c=10.,alpha=90.,beta=90.,gamma=90.),
             dict(a=12.,b=14.,c=15.,alpha=75.,beta=83.,gamma=64.)]
    cases = json.loads((FIXTURES / "cases.json").read_text())["cases"]
    receipts, water = [], None
    for case in cases:
        for cell in cells:
            receipt, text = validate_case(gemmi, chematic, case, cell)
            receipts.append(receipt)
            if case["name"] == "water" and cell is not None:
                water = text
    for length in [1e308, 1e-300, math.ulp(0.)]:
        case = dict(name=f"extreme-{length}",smiles="CC",elements=["C","C"],
                    coords=[[0.,0.,0.],[length,0.,0.]],bonds=[[0,1,"S"]])
        receipt, _ = validate_case(gemmi, chematic, case, None)
        receipts.append(receipt)
    binary = Path(importlib.import_module("chematic.chematic").__file__)
    result = dict(schema="chematic.cif-bonds-validation.v1", utc=datetime.now(timezone.utc).isoformat(),
                  source_base=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                  source_sha256={f:sha(ROOT/f) for f in SOURCE}, extension_sha256=sha(binary),
                  fixture_sha256={f.name:sha(f) for f in FIXTURES.iterdir() if f.is_file()},
                  gemmi_version=gemmi.__version__,
                  gemmi_runtime_sha256={f.name:sha(f) for f in Path(gemmi.__file__).parent.glob("*.so")},
                  cases=receipts,
                  gemmi_checks=sum(r["checks"] for r in receipts), failed_checks=0,
                  openbabel_version=probe([args.obabel,"-V"]), openbabel_binary_sha256=sha(args.obabel),
                  openbabel_writes=[probe([args.obabel,"-ixyz",str(FIXTURES/"water.xyz"),"-ocif",flag]) for flag in ["-ag","-xg"]],
                  openbabel_read=probe([args.obabel,"-icif","-aB","-omolreport"],water),
                  boundary="Source candidate only. #2726 was a CLI option mix-up; Open Babel already supports -xg. Plain chematic CIF reader remains positions-only; no periodic bond perception or full molecular serialization.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(f"Gemmi: {result['gemmi_checks']} checks in {len(receipts)} cases, zero failures")
    print(f"Evidence: {args.output}")


if __name__ == "__main__":
    main()
