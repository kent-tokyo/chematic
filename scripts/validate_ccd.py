#!/usr/bin/env python3
"""Compare the source CCD reader with Gemmi 0.7.5 and diagnose Open Babel.

Install Gemmi only in an oracle environment, build the chematic extension, then
run with --download or supply the eight official files in --source-dir.
The CCD source files are cached locally; their bytes are not copied into the
repository. The original synthetic fixture is the hermetic CI input.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = ["GLY", "HOH", "ALA", "NA", "CL", "BEN", "ATP", "HEM"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("/private/tmp/chematic-ccd-sources"))
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--source-base", default="HEAD", help="Git ref used as the source candidate base")
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/results/openbabel-ccd-2026-10-16.json")
    args = parser.parse_args()
    import chematic
    import gemmi
    assert gemmi.__version__ == "0.7.5", gemmi.__version__
    args.source_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    checks = 0
    version = subprocess.check_output([args.obabel, "-V"], text=True).strip()
    for component_id in COMPONENTS:
        source = args.source_dir / f"{component_id}.cif"
        url = f"https://files.rcsb.org/ligands/download/{component_id}.cif"
        if args.download:
            subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error", "--max-time", "30", url, "-o", str(source)], check=True)
        text = source.read_text(encoding="utf-8")
        block = gemmi.cif.read_string(text).sole_block()
        try:
            record = json.loads(chematic.parse_ccd_json(text))
        except ValueError as error:
            raise AssertionError(f"{component_id}: {error}") from error

        def column(category, field):
            return [gemmi.cif.as_string(v) if v not in (".", "?") else v
                    for v in block.find_values(f"_{category}.{field}")]

        atoms = record["atoms"]
        bonds = record["bonds"]
        assert record["component_id"] == component_id
        assert len(atoms) == len(column("chem_comp_atom", "atom_id"))
        assert len(bonds) == len(column("chem_comp_bond", "atom_id_1"))
        for field, key in [("atom_id", "atom_id"), ("type_symbol", "element"), ("charge", "charge"), ("pdbx_stereo_config", "stereo_config")]:
            reference = column("chem_comp_atom", field)
            expected = [(None if v in (".", "?") else int(v)) if key == "charge" else (gemmi.Element(v).name if key == "element" else (None if v in (".", "?") else v)) for v in reference]
            assert [a[key] for a in atoms] == expected, (component_id, field)
            checks += len(atoms)
        for key, pattern in [("model_coordinates", "model_Cartn_{}"), ("ideal_coordinates", "pdbx_model_Cartn_{}_ideal")]:
            columns = [column("chem_comp_atom", pattern.format(axis)) for axis in "xyz"]
            expected = [None if all(v in (".", "?") for v in values) else list(map(float, values)) for values in zip(*columns)]
            assert [a[key] for a in atoms] == expected, (component_id, key)
            checks += len(atoms) * 3
        for field in ["atom_id_1", "atom_id_2", "value_order"]:
            assert [b[field] for b in bonds] == column("chem_comp_bond", field), (component_id, field)
            checks += len(bonds)
        assert [a["aromatic"] for a in atoms] == [v == "Y" for v in column("chem_comp_atom", "pdbx_aromatic_flag")]
        assert [b["aromatic"] for b in bonds] == [v == "Y" or b["value_order"] == "AROM" for v, b in zip(column("chem_comp_bond", "pdbx_aromatic_flag"), bonds)]
        checks += len(atoms) + len(bonds)
        requires_consent = any(a["stereo_config"] not in (None, "N") for a in atoms) or any(b["stereo_config"] not in (None, "N") for b in bonds)
        if requires_consent:
            try:
                chematic.from_ccd(text)
                raise AssertionError("expected stereo-loss refusal")
            except ValueError as error:
                assert error.code == "ccd_stereo_loss"
        graph = None
        refusal = None
        try:
            graph = chematic.from_ccd(text, allow_stereo_loss=requires_consent)
        except ValueError as error:
            assert any(a["charge"] is None for a in atoms) and error.code == "ccd_unknown_charge"
            refusal = {"code": error.code, "reason": str(error)}
        counts = Counter(gemmi.Element(v).name for v in column("chem_comp_atom", "type_symbol"))
        charge = None if any(a["charge"] is None for a in atoms) else sum(a["charge"] for a in atoms)
        if graph is not None and charge == 0:
            order = (["C", "H"] if "C" in counts else []) + sorted(k for k in counts if "C" not in counts or k not in ("C", "H"))
            formula = "".join(k + (str(counts[k]) if counts[k] != 1 else "") for k in order if counts[k])
            assert graph.formula == formula, (component_id, graph.formula, formula)
        for category, objects in [("chem_comp_atom", atoms), ("chem_comp_bond", bonds)]:
            for index, obj in enumerate(objects):
                for tag, value in obj["extra_fields"].items():
                    raw = block.find_values(tag)[index]
                    expected = raw if raw in (".", "?") else gemmi.cif.as_string(raw)
                    assert value == expected, (component_id, tag, index)
                    checks += 1
            reference = column(category, "pdbx_stereo_config")
            assert [o["stereo_config"] for o in objects] == [None if v in (".", "?") else v for v in reference]
            checks += len(objects)
        diagnostics = {}
        for fmt in ["mmcif", "cif"]:
            command = [args.obabel, f"-i{fmt}", str(source), "-osmi"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=20)
            diagnostics[fmt] = {"returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}
        rows.append({"component_id": component_id, "source_url": url, "input_sha256": sha(source),
                     "atoms": len(atoms), "bonds": len(bonds), "total_formal_charge": charge,
                     "stereo_loss_consent_required": requires_consent, "graph_formula": graph.formula if graph is not None else None, "graph_refusal": refusal,
                     "openbabel": diagnostics, "record": record})
    fixture = ROOT / "validation/openbabel_issues/ccd/water.cif"
    sources = ["crates/chematic-mol/src/ccd.rs", *[str(p.relative_to(ROOT)) for p in sorted((ROOT / "crates/chematic-mol/src/ccd").glob("*.rs"))], "crates/chematic-py/src/formats/ccd.rs", "crates/chematic-wasm/src/ccd.rs"]
    report = {"schema_version": 1, "checked_at_utc": datetime.now(timezone.utc).isoformat(),
              "issue": "https://github.com/openbabel/openbabel/issues/1516", "issue_state_checked": "open",
              "source_base_commit": subprocess.check_output(["git", "rev-parse", args.source_base], cwd=ROOT, text=True).strip(),
              "source_candidate": True, "published_artifact": False, "gemmi_version": gemmi.__version__,
              "openbabel_version": version, "openbabel_binary_sha256": sha(args.obabel),
              "gemmi_extension_sha256": sha(importlib.import_module("gemmi.gemmi_ext").__file__), "extension_sha256": sha(importlib.import_module("chematic.chematic").__file__),
              "source_sha256": {s: sha(ROOT / s) for s in sources}, "field_checks": checks, "failed_cases": 0,
              "fixture_sha256": sha(fixture), "synthetic_record": json.loads(chematic.parse_ccd_json(fixture.read_text(encoding="utf-8"))),
              "official_components": rows}
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"components": len(rows), "field_checks": checks, "failed_cases": 0}))


if __name__ == "__main__":
    main()
