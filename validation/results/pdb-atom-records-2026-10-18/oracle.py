"""Independent Gemmi PDB parsing, exact source retention and local OB control."""

import hashlib
import json
from pathlib import Path
import random
import subprocess

import gemmi
import chematic

HERE = Path(__file__).parent
rng = random.Random(2708)
input_hash = hashlib.sha256()
checked = 0
max_occupancy_error = max_b_error = max_coordinate_error = 0.0

for case in range(96):
    rows = []
    for i in range(6 * (1 + case % 8)):
        serial = i + 1
        kind = "HETATM" if i % 3 == 0 else "ATOM"
        residue = "LIG" if kind == "HETATM" else "VAL"
        alt = "AB"[i % 2]
        chain = ["A", "Z", " "][i % 3]
        insertion = "I" if i % 4 == 0 else " "
        xyz = tuple(rng.uniform(-50, 50) for _ in range(3))
        occupancy = rng.randint(0, 100) / 100
        b_factor = rng.randint(0, 9999) / 100
        row = (
            f"{kind:6}{serial:5d}  C1 {alt}{residue} {chain}{i+1:4d}{insertion}   "
            + "".join(f"{v:8.3f}" for v in xyz)
            + f"{occupancy:6.2f}{b_factor:6.2f}"
            + " " * 10
            + " C  "
        )
        assert len(row) == 80
        rows.append(row)
    rng.shuffle(rows)
    model = None if case % 3 == 0 else (1 if case % 3 == 1 else 7)
    lines = ["HEADER    original synthetic attribute corpus"]
    if model is not None:
        lines.append(f"MODEL     {model:4d}")
    lines.extend(rows)
    if model is not None:
        lines.append("ENDMDL")
    lines.append("END")
    newline = "\r\n" if case % 2 else "\n"
    text = newline.join(lines) + newline
    parsed = gemmi.read_pdb_string(text)
    reference = {cra.atom.serial: cra for cra in parsed[0].all()}
    actual = chematic.read_pdb_atom_records(text)
    assert actual.source_text == text and actual.model_num == model
    assert len(actual.atoms) == len(reference) == len(rows)
    for index, atom in enumerate(actual.atoms):
        assert atom.source_index == index
        assert atom.record_text == rows[index]
        assert text.splitlines()[atom.line_number - 1] == rows[index]
        cra = reference[atom.serial]
        ref = cra.atom
        assert atom.atom_name.strip() == ref.name
        assert atom.alt_loc == ref.altloc
        assert atom.residue_name == cra.residue.name
        assert atom.chain_id.strip() == cra.chain.name
        assert atom.residue_number == cra.residue.seqid.num
        assert atom.insertion_code == cra.residue.seqid.icode
        assert atom.record_type == ("ATOM" if cra.residue.het_flag == "A" else "HETATM")
        xyz = (ref.pos.x, ref.pos.y, ref.pos.z)
        difference = max(abs(a - b) for a, b in zip(atom.position_angstrom, xyz))
        max_coordinate_error = max(max_coordinate_error, difference)
        assert difference == 0
        occupancy_error = abs(atom.occupancy - ref.occ)
        b_error = abs(atom.b_factor - ref.b_iso)
        max_occupancy_error = max(max_occupancy_error, occupancy_error)
        max_b_error = max(max_b_error, b_error)
        # Gemmi stores these attributes as float32; our read keeps Python f64.
        assert occupancy_error < 1e-7 and b_error < 1e-5
        checked += 1
    input_hash.update(text.encode())

control = subprocess.run(
    ["/opt/homebrew/bin/obabel", "-ipdb", str(HERE / "attributes.pdb"), "-opdb", "-b"],
    capture_output=True,
    text=True,
    timeout=30,
)
(HERE / "openbabel.stdout").write_text(control.stdout)
(HERE / "openbabel.stderr").write_text(control.stderr)
source = chematic.read_pdb_atom_records((HERE / "attributes.pdb").read_text())
output = chematic.read_pdb_atom_records(control.stdout)
summary = dict(
    seed=2708,
    gemmi=gemmi.__version__,
    cases=96,
    atoms_checked=checked,
    max_coordinate_difference=max_coordinate_error,
    max_occupancy_difference=max_occupancy_error,
    max_b_factor_difference=max_b_error,
    reference_precision="Gemmi occupancy/B are float32; source-derived Python values are float64",
    inputs_sha256=input_hash.hexdigest(),
    openbabel=dict(
        version=subprocess.check_output(
            ["/opt/homebrew/bin/obabel", "-V"], text=True
        ).strip(),
        returncode=control.returncode,
        input_altloc=[a.alt_loc for a in source.atoms],
        output_altloc=[a.alt_loc for a in output.atoms],
        input_b_factor=[a.b_factor for a in source.atoms],
        output_b_factor=[a.b_factor for a in output.atoms],
        input_occupancy=[a.occupancy for a in source.atoms],
        output_occupancy=[a.occupancy for a in output.atoms],
    ),
)
assert control.returncode == 0
assert summary["openbabel"]["output_b_factor"] == [0.0, 0.0]
assert summary["openbabel"]["output_altloc"] == [" ", " "]
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
