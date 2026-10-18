"""Independent Gemmi identity/tensor/connectivity checks for PDB block assembly."""

import hashlib
import json
from pathlib import Path
import subprocess

import gemmi

from chematic import merge_pdb_coordinate_blocks

HERE = Path(__file__).parent


def atom(serial, chain, kind, name, residue, number, index, case, source):
    alt = "A" if index == 0 and case % 2 == 0 else " "
    insertion = "I" if case % 3 == 0 else " "
    x, y, z = (
        source * 20 + case * 0.125 + index * 0.5,
        -index * 1.25,
        source * 1.5 + index * 0.125,
    )
    element = "N" if name.strip() == "N" else ("O" if name.strip() == "O" else "C")
    charge = "1+" if index == 0 and case % 4 == 0 else "  "
    return f"{kind:<6}{serial:5} {name:4}{alt}{residue:3} {chain}{number:4}{insertion}   {x:8.3f}{y:8.3f}{z:8.3f}{.25*(index%4+1):6.2f}{10+index*.5:6.2f}          {element:>2}{charge:2}"


def anisou(line):
    return (
        "ANISOU"
        + line[6:28]
        + "".join(f"{value:7}" for value in [100, 200, 300, -4, 5, -6])
        + line[70:]
    )


def snapshot(text):
    structure = gemmi.read_pdb_string(text)
    assert len(structure) == 1
    atoms = []
    for chain in structure[0]:
        for residue in chain:
            for a in residue:
                atoms.append(
                    dict(
                        chain=chain.name,
                        residue=residue.name,
                        number=residue.seqid.num,
                        insertion=residue.seqid.icode,
                        kind=residue.het_flag,
                        name=a.name,
                        altloc=a.altloc,
                        element=a.element.name,
                        charge=a.charge,
                        serial=a.serial,
                        position=(a.pos.x, a.pos.y, a.pos.z),
                        occupancy=a.occ,
                        bfactor=a.b_iso,
                        anisou=tuple(
                            getattr(a.aniso, name)
                            for name in ["u11", "u22", "u33", "u12", "u13", "u23"]
                        ),
                    )
                )
    return atoms, {k: list(v) for k, v in structure.conect_map.items()}


corpus = hashlib.sha256()
assemblies = sources_count = atom_count = anisou_count = connect_count = 0
max_coordinate_error = 0.0
for case in range(96):
    texts = []
    maps = []
    originals = []
    for source in range(2 + case % 5):
        chain = " " if (case + source) % 3 == 0 else "A"
        mapping = {chain: chr(65 + source)}
        kind = "ATOM" if source % 2 == 0 else "HETATM"
        residue = "ALA" if source % 2 == 0 else "LIG"
        count = 2 + (case + source) % 4
        rows = []
        serials = []
        for i in range(count):
            serial = 10 + case * 100 + i * 3
            serials.append(serial)
            line = atom(
                serial,
                chain,
                kind,
                [" CA ", " N  ", " CB ", " O  ", " C2 "][i],
                residue,
                -2 + i // 2,
                i,
                case,
                source,
            )
            rows.append(line)
            if i == 0:
                rows.append(anisou(line))
        if (case + source) % 2 == 0:
            rows.append(f"TER   {serials[-1]+1:5}      {line[17:27]}")
        rows.extend(
            [
                f"CONECT{serials[0]:5}{serials[-1]:5}{serials[-1]:5}",
                f"CONECT{serials[-1]:5}{serials[0]:5}",
                "END",
            ]
        )
        text = "\n".join(rows) + "\n"
        if case % 2:
            text = text.replace("\n", "\r\n")
        texts.append(text)
        maps.append(mapping)
        originals.append(snapshot(text))
    corpus.update(json.dumps([texts, maps]).encode())
    result = merge_pdb_coordinate_blocks(iter(texts), chain_maps=iter(maps))
    assert result.source_texts == tuple(texts) and not result.losses
    actual, connections = snapshot(result.pdb)
    expected = []
    expected_connections = {}
    for source, (before, old_connections) in enumerate(originals):
        correspondences = [a for a in result.atoms if a.source_index == source]
        remap = {a.old_serial: a.new_serial for a in correspondences}
        for a in before:
            a["chain"] = maps[source][a["chain"] or " "]
            a["serial"] = remap[a["serial"]]
            expected.append(a)
        for serial, targets in old_connections.items():
            expected_connections[remap[serial]] = [remap[t] for t in targets]
    assert len(actual) == len(expected) == len(result.atoms)
    for before, after in zip(expected, actual):
        max_coordinate_error = max(
            max_coordinate_error,
            max(abs(a - b) for a, b in zip(before["position"], after["position"])),
        )
        assert before == after, (case, before, after)
    assert connections == expected_connections
    source_rows = [
        line
        for text in texts
        for line in text.splitlines()
        if line.startswith(("ATOM  ", "HETATM", "ANISOU"))
    ]
    output_rows = [
        line
        for line in result.pdb.splitlines()
        if line.startswith(("ATOM  ", "HETATM", "ANISOU"))
    ]
    assert len(source_rows) == len(output_rows)
    for a, b in zip(source_rows, output_rows):
        assert a[:6] + a[11:21] + a[22:] == b[:6] + b[11:21] + b[22:]
    assemblies += 1
    sources_count += len(texts)
    atom_count += len(actual)
    anisou_count += sum(line.startswith("ANISOU") for line in output_rows)
    connect_count += sum(len(v) for v in connections.values())

# Reproduce the reported class/identity/TER problem on an original minimal input.
protein = (
    "\n".join(
        [
            atom(20, "A", "ATOM", " CA ", "ALA", 7, 0, 0, 0).rstrip(),
            atom(35, "A", "ATOM", " N  ", "ALA", 7, 1, 0, 0).rstrip(),
            "END",
        ]
    )
    + "\n"
)
ligand = (
    "\n".join([atom(20, "L", "HETATM", " C1 ", "LIG", 7, 0, 1, 1).rstrip(), "END"])
    + "\n"
)
for name, text in [("protein", protein), ("ligand", ligand)]:
    (HERE / (name + ".pdb")).write_text(text)
merged = merge_pdb_coordinate_blocks([protein, ligand])
(HERE / "merged.pdb").write_text(merged.pdb)
command = [
    "/opt/homebrew/bin/obabel",
    "-ipdb",
    str(HERE / "protein.pdb"),
    str(HERE / "ligand.pdb"),
    "-opdb",
    "-j",
]
process = subprocess.run(command, capture_output=True, text=True)
(HERE / "openbabel.pdb").write_text(process.stdout)
(HERE / "openbabel.stderr").write_text(process.stderr)
assert process.returncode == 0 and "1 molecule converted" in process.stderr
before = snapshot(protein)[0] + snapshot(ligand)[0]
after = snapshot(process.stdout)[0]
assert len(before) == len(after) == 3
control = dict(
    returncode=process.returncode,
    atoms=len(before),
    changed_record_kinds=sum(a["kind"] != b["kind"] for a, b in zip(before, after)),
    changed_residue_names=sum(
        a["residue"] != b["residue"] for a, b in zip(before, after)
    ),
    changed_chains=sum(a["chain"] != b["chain"] for a, b in zip(before, after)),
    openbabel_ter_records=sum(
        line.startswith("TER") for line in process.stdout.splitlines()
    ),
    chematic_ter_records=sum(
        line.startswith("TER") for line in merged.pdb.splitlines()
    ),
    version=subprocess.check_output(
        ["/opt/homebrew/bin/obabel", "-V"], text=True
    ).strip(),
)
assert control["changed_record_kinds"] == 2 and control["changed_residue_names"] == 3
assert control["openbabel_ter_records"] == 0 and control["chematic_ter_records"] == 2
(HERE / "openbabel-summary.json").write_text(json.dumps(control, indent=2) + "\n")
summary = dict(
    gemmi=gemmi.__version__,
    assemblies=assemblies,
    sources=sources_count,
    atoms=atom_count,
    anisou_records=anisou_count,
    conect_target_references=connect_count,
    max_coordinate_error_angstrom=max_coordinate_error,
    inputs_sha256=corpus.hexdigest(),
    scope="Synthetic source identity, numeric spelling, tensor and explicit-connectivity preservation; no chemistry/geometry alignment claim",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
print(json.dumps(control, indent=2))
