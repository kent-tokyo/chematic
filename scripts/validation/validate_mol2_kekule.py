#!/usr/bin/env python3
"""Verify exact token changes, cross-binding outputs and independent chemistry.

Requires the candidate Python extension, RDKit 2026.09.1, Node/WASM, a built
mol2_kekule_dump example and Open Babel. No downloaded upstream fixtures.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import re
import subprocess
import sys
from pathlib import Path

import chematic
from rdkit import Chem, rdBase

ROOT = Path(__file__).resolve().parents[2]


def section(text, name):
    rows = []
    active = False
    for line in text.splitlines(keepends=True):
        if line.strip().upper().startswith("@<TRIPOS>"):
            active = line.strip().upper() == "@<TRIPOS>" + name
        elif active and line.strip() and not line.lstrip().startswith("#"):
            rows.append(line)
    return rows


def variant(text, seed, ending):
    text = text.replace("\r\n", "\n")
    lines = text.splitlines(keepends=True)
    rng = random.Random(seed)
    if seed is not None:
        for name in ["ATOM", "BOND"]:
            rows = section(text, name)
            shuffled = rows.copy()
            rng.shuffle(shuffled)
            # Rows have unique source IDs; keep comments and section positions.
            swaps = dict(zip(rows, shuffled))
            lines = [swaps.get(line, line) for line in lines]
    result = "".join(lines)
    if ending == "none":
        return result.rstrip("\n")
    return result.replace("\n", ending)


def token_proof(source, output):
    before, after = source.splitlines(keepends=True), output.splitlines(keepends=True)
    assert len(before) == len(after)
    active = False
    changes = 0
    for a, b in zip(before, after):
        if a.strip().upper().startswith("@<TRIPOS>"):
            active = a.strip().upper() == "@<TRIPOS>BOND"
        tokens = list(re.finditer(r"\S+", a))
        if active and len(tokens) >= 4 and tokens[3].group() == "ar":
            start, end = tokens[3].span()
            assert b[start] in "12" and b == a[:start] + b[start] + a[end:]
            changes += 1
        else:
            assert a == b
    return changes


def rdkit_check(source, output):
    def parse(text):
        # RDKit's MOL2 reader cannot consume the sparse-ID, comment/status-rich
        # source profile directly. Give it an independently normalized graph
        # projection with sequential IDs, retaining coordinates/types/charges.
        atoms = [row.split() for row in section(text, "ATOM")]
        bonds = [row.split() for row in section(text, "BOND")]
        ids = {int(row[0]): i+1 for i, row in enumerate(atoms)}
        lines = ["@<TRIPOS>MOLECULE", "oracle", f"{len(atoms)} {len(bonds)} 0 0 0", "SMALL", "USER_CHARGES", "", "@<TRIPOS>ATOM"]
        for row in atoms:
            lines.append(" ".join([str(ids[int(row[0])]), *row[1:9]]))
        lines.append("@<TRIPOS>BOND")
        for i, row in enumerate(bonds):
            lines.append(f"{i+1} {ids[int(row[1])]} {ids[int(row[2])]} {row[3]}")
        m = Chem.MolFromMol2Block("\n".join(lines)+"\n", sanitize=False, removeHs=False, cleanupSubstructures=False)
        assert m is not None
        # UNITY charge is explicit source metadata, never inferred from the
        # electrostatic column. Apply it before the reference sanitization.
        rows = iter(section(text, "UNITY_ATOM_ATTR"))
        for row in rows:
            atom, count = map(int, row.split())
            for _ in range(count):
                field, value = next(rows).split()
                if field.lower() == "charge":
                    m.GetAtomWithIdx(ids[atom]-1).SetFormalCharge(int(value))
        Chem.SanitizeMol(m)
        Chem.AssignStereochemistryFrom3D(m)
        return m
    a, b = parse(source), parse(output)
    assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
    assert a.GetNumAtoms() == b.GetNumAtoms() and a.GetNumBonds() == b.GetNumBonds()
    color = lambda atom: (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(), atom.GetTotalNumHs(), atom.GetIsAromatic(), str(atom.GetChiralTag()))
    assert [color(x) for x in a.GetAtoms()] == [color(x) for x in b.GetAtoms()]
    edges = lambda mol: [(x.GetBeginAtomIdx(), x.GetEndAtomIdx(), x.GetBondTypeAsDouble(), str(x.GetStereo())) for x in mol.GetBonds()]
    assert edges(a) == edges(b)
    ids = [int(row.split()[0]) for row in section(source, "ATOM")]
    valences = dict.fromkeys(ids, 0)
    for row in section(output, "BOND"):
        _, x, y, order, *_ = row.split()
        value = 1 if order == "am" else int(order)
        valences[int(x)] += value
        valences[int(y)] += value
    # Independent of re-aromatization: the explicit integer orders must satisfy
    # every reference atom's valence, including explicit and implicit H slots.
    for source_id, atom in zip(ids, a.GetAtoms()):
        assert valences[source_id] + atom.GetTotalNumHs() == atom.GetTotalValence()
    return Chem.MolToSmiles(a)


def obabel(text, executable):
    # Open Babel's stdin reader rejects CRLF here; normalize only line endings
    # in the oracle copy. Native/Python/WASM source-byte checks use original text.
    r = subprocess.run([executable, "-imol2", "-ocan"], input=text.replace("\r\n", "\n"), text=True, capture_output=True, check=True)
    assert r.stdout.strip(), r.stderr
    # The MOL2 reader can keep explicit bonds unflagged as aromatic. Reparse
    # its valid SMILES to compare chemical identity across these encodings.
    normalized = subprocess.run([executable, "-ismi", "-ocan"], input=r.stdout, text=True, capture_output=True, check=True)
    assert normalized.stdout.strip(), normalized.stderr
    return dict(canonical=normalized.stdout.split()[0], mol2_canonical=r.stdout.split()[0], stderr=r.stderr, smiles_reparse_stderr=normalized.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-dump", type=Path, required=True)
    parser.add_argument("--wasm-module", type=Path, required=True)
    parser.add_argument("--obabel", default="/opt/homebrew/bin/obabel")
    parser.add_argument("--output", type=Path, default=ROOT / "validation/results/mol2-kekule-2026-10-16")
    args = parser.parse_args()
    assert rdBase.rdkitVersion == "2026.09.1"
    fixture = ROOT / "validation/openbabel_issues/mol2-kekule/cases.json"
    base = json.loads(fixture.read_text())["cases"]
    cases = [dict(id=f'{c["id"]}-{seed}-{repr(ending)}', text=variant(c["text"], seed, ending), aromatic_bonds=c["aromatic_bonds"]) for c in base for seed in [None, 0, 11, 53] for ending in ["\n", "\r\n", "none"]]
    payload = json.dumps([c["text"] for c in cases])
    native = json.loads(subprocess.run([str(args.native_dump)], input=payload, capture_output=True, text=True, check=True).stdout)
    wasm = json.loads(subprocess.run(["node", str(ROOT / "scripts/validation/validate_mol2_kekule_node.mjs"), str(args.wasm_module.resolve())], input=payload, capture_output=True, text=True, check=True).stdout)
    assert len(native) == len(wasm) == len(cases)
    rows = []
    for case, rust, js in zip(cases, native, wasm):
        output = rust["output"]
        assert output == js["output"] == chematic.kekulize_mol2_record(case["text"])
        assert chematic.kekulize_mol2_record(output) == output
        assert token_proof(case["text"], output) == case["aromatic_bonds"]
        canonical = rdkit_check(case["text"], output)
        rows.append(dict(case, output=output, rdkit_canonical=canonical))
    cross_reads = []
    for case in base:
        output = chematic.kekulize_mol2_record(case["text"])
        a, b = obabel(case["text"], args.obabel), obabel(output, args.obabel)
        assert a["canonical"] == b["canonical"], case["id"]
        cross_reads.append(dict(id=case["id"], before=a, after=b))
    diagnostic = []
    for flags in [[], ["-xk"]]:
        r = subprocess.run([args.obabel, "-:c1ccncc1", "-omol2", *flags], text=True, capture_output=True, check=True)
        count = sum(row.split()[3] == "ar" for row in section(r.stdout, "BOND"))
        diagnostic.append(dict(flags=flags, aromatic_bonds=count, stdout=r.stdout, stderr=r.stderr))
    packet = dict(schema="chematic.mol2-kekule-comparisons.v1", rows=rows, openbabel_cross_reads=cross_reads, openbabel_xk_diagnostic=diagnostic)
    args.output.mkdir(parents=True, exist_ok=True)
    blob = gzip.compress(json.dumps(packet, sort_keys=True).encode(), mtime=0)
    (args.output / "comparisons.json.gz").write_bytes(blob)
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    sources = ["crates/chematic-mol/src/mol2_kekule.rs", "crates/chematic-mol/src/mol2_source.rs", "crates/chematic-core/src/kekulization.rs", "crates/chematic-py/src/formats/interchange.rs", "crates/chematic-wasm/src/mol_io/common_formats.rs", "crates/chematic-mol/examples/mol2_kekule_dump.rs", "scripts/validation/validate_mol2_kekule.py", "scripts/validation/validate_mol2_kekule_node.mjs", str(fixture.relative_to(ROOT))]
    summary = dict(schema="chematic.mol2-kekule-summary.v1", source_base=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), source_sha256={p: sha(ROOT/p) for p in sources}, runtime=dict(python=sys.version.split()[0], rdkit=rdBase.rdkitVersion, openbabel=subprocess.check_output([args.obabel, "-V"], text=True).strip(), python_extension_sha256=sha(chematic.chematic.__file__), native_dump_sha256=sha(args.native_dump), wasm_sha256=sha(args.wasm_module.with_name("chematic_wasm_bg.wasm"))), results=dict(cases=len(base), cross_binding_outputs=len(rows), exact_token_proofs=len(rows), rdkit_identity_and_explicit_valence_checks=len(rows), openbabel_identity_checks=len(cross_reads), mismatches=0), packet_sha256=sha(args.output / "comparisons.json.gz"), rdkit_boundary="Independent normalized MOL2 graph projection with sequential IDs and explicit UNITY formal charges; not a direct raw-record interoperability claim.", scope="Bounded opt-in neutral C/N/O/S MOL2 bond encoding, with source-byte preservation. Not general atom retyping, aromaticity-model parity, 3D generation, charge prediction, performance or release proof.")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["results"], indent=2))


if __name__ == "__main__":
    main()
