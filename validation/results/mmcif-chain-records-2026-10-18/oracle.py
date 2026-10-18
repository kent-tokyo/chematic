"""Independent CIF parsing and model/label-chain checks using Gemmi 0.7.3."""

from collections import OrderedDict
import hashlib
import json
from pathlib import Path
import random
import subprocess

import gemmi
import chematic

HERE = Path(__file__).parent
SEED = 2458
rng = random.Random(SEED)
base = gemmi.cif.read_string((HERE / "chains.cif").read_text()).sole_block()
table = base.find_mmcif_category("_atom_site.")
columns = [tag.split(".", 1)[1] for tag in table.tags]
base_rows = [
    [row.str(i) if row[i] not in (".", "?") else row[i] for i in range(len(columns))]
    for row in table
]
col = {name: i for i, name in enumerate(columns)}
inputs_hash = hashlib.sha256()
checked_atoms = checked_chains = 0
max_coordinate_error = 0.0


def optional(row, key):
    value = row.str(col[key])
    return value if value else None


def quote(value):
    return "'" + value + "'" if any(c.isspace() for c in value) else value


def compare(text):
    global checked_atoms, checked_chains, max_coordinate_error
    # This parser and row access are independent of the native chematic reader.
    block = gemmi.cif.read_string(text).sole_block()
    rows = block.find("_atom_site.", columns)
    expected = OrderedDict()
    for index, row in enumerate(rows):
        key = (int(row[col["pdbx_PDB_model_num"]]), row.str(col["label_asym_id"]))
        expected.setdefault(key, []).append(index)
    actual = chematic.read_mmcif_chains(text)
    assert actual.source_text == text
    assert [(c.model_num, c.label_asym_id) for c in actual.chains] == list(expected)
    actual_atoms = []
    for chain in actual.chains:
        key = (chain.model_num, chain.label_asym_id)
        assert [a.source_index for a in chain.atoms] == expected[key]
        checked_chains += 1
        actual_atoms.extend(chain.atoms)
    assert sorted(a.source_index for a in actual_atoms) == list(range(len(rows)))
    for atom in actual_atoms:
        row = rows[atom.source_index]
        assert atom.serial == int(row[col["id"]])
        assert atom.group_pdb == row.str(col["group_PDB"])
        assert atom.element == row.str(col["type_symbol"])
        assert atom.label_atom_id == optional(row, "label_atom_id")
        assert atom.auth_atom_id == optional(row, "auth_atom_id")
        assert atom.alt_loc == optional(row, "label_alt_id")
        xyz = tuple(float(row[col[name]]) for name in ["Cartn_x", "Cartn_y", "Cartn_z"])
        assert atom.position_angstrom == xyz
        max_coordinate_error = max(
            max_coordinate_error,
            *[abs(a - b) for a, b in zip(atom.position_angstrom, xyz)]
        )
        assert atom.occupancy == float(row[col["occupancy"]])
        assert atom.b_iso == float(row[col["B_iso_or_equiv"]])
        charge = optional(row, "pdbx_formal_charge")
        assert atom.formal_charge == (int(charge) if charge is not None else None)
        r = atom.residue
        for attr, name in [
            ("auth_asym_id", "auth_asym_id"),
            ("entity_id", "label_entity_id"),
            ("auth_seq_id", "auth_seq_id"),
            ("insertion_code", "pdbx_PDB_ins_code"),
            ("label_comp_id", "label_comp_id"),
            ("auth_comp_id", "auth_comp_id"),
        ]:
            assert getattr(r, attr) == optional(row, name)
        seq = optional(row, "label_seq_id")
        assert r.label_seq_id == (int(seq) if seq is not None else None)
        assert r.model_num == int(row[col["pdbx_PDB_model_num"]])
        assert r.label_asym_id == row.str(col["label_asym_id"])
        checked_atoms += 1
    inputs_hash.update(text.encode())


for case in range(96):
    rows = []
    for replicate in range(1 + case % 12):
        for original in base_rows:
            row = original.copy()
            row[col["id"]] = str(len(rows) + 1)
            row[col["auth_seq_id"]] = str(
                int(row[col["auth_seq_id"]]) + replicate * 100
            )
            seq = row[col["label_seq_id"]]
            if seq not in (".", "?"):
                row[col["label_seq_id"]] = str(int(seq) + replicate * 100)
            for axis in ["Cartn_x", "Cartn_y", "Cartn_z"]:
                row[col[axis]] = format(rng.uniform(-1000, 1000), ".17g")
            if case % 3 == 0:
                row[col["label_asym_id"]] += " chain"
                row[col["auth_atom_id"]] += " atom"
            rows.append(row)
    rng.shuffle(rows)
    permutation = list(range(len(columns)))
    rng.shuffle(permutation)
    header = "data_independent\n# authored randomized corpus\nloop_\n" + "".join(
        "_atom_site." + columns[i] + "\n" for i in permutation
    )
    data = []
    for index, row in enumerate(rows):
        tokens = [quote(row[i]) for i in permutation]
        if index % 5 == 0:
            data.extend([" ".join(tokens[:11]), " ".join(tokens[11:])])
        else:
            data.append(" ".join(tokens))
    compare(header + "\n".join(data) + "\n#\n")

# Local upstream control measures default record granularity, not chemistry.
control = subprocess.run(
    ["/opt/homebrew/bin/obabel", "-immcif", str(HERE / "chains.cif"), "-oxyz", "-b"],
    capture_output=True,
    text=True,
    timeout=30,
)
(HERE / "openbabel.stdout").write_text(control.stdout)
(HERE / "openbabel.stderr").write_text(control.stderr)
version = subprocess.check_output(["/opt/homebrew/bin/obabel", "-V"], text=True).strip()
summary = dict(
    seed=SEED,
    gemmi=gemmi.__version__,
    cases=96,
    checked_atoms=checked_atoms,
    checked_chains=checked_chains,
    max_coordinate_error=max_coordinate_error,
    inputs_sha256=inputs_hash.hexdigest(),
    openbabel=dict(
        version=version,
        returncode=control.returncode,
        reported_one_record="1 molecule converted" in control.stderr,
    ),
    scope="Label-chain/model partitions and typed atom/residue fields, no topology inference or exact source rewriting",
)
assert control.returncode == 0
assert summary["openbabel"]["reported_one_record"]
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
