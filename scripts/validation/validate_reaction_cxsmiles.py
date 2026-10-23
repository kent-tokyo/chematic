#!/usr/bin/env python3
"""RDKit #9732 fixed-packet writer capture and independent RDKit observation."""
import argparse
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "validation/rdkit-9732-2026-10-24/cases.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix == ".gz" else data)


def fingerprint(text):
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions as reactions
    reaction = reactions.ReactionFromSmarts(text, useSmiles=True)
    if reaction is None:
        raise ValueError("RDKit observer refused reaction")
    grouped = {}
    for role, molecules in enumerate([reaction.GetReactants(), reaction.GetAgents(), reaction.GetProducts()]):
        for mol in molecules:
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            identity = Chem.Mol(mol)
            Chem.RemoveStereochemistry(identity)
            for group in mol.GetStereoGroups():
                for atom in group.GetAtoms():
                    key = (str(group.GetGroupType()), group.GetReadId())
                    grouped.setdefault(key, []).append({
                        "role": role,
                        "map": atom.GetAtomMapNum(),
                        "isotope": atom.GetIsotope(),
                        "tetrahedral": str(atom.GetChiralTag()) in ["CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW"],
                        "cip": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
                        "rooted_graph": Chem.MolToSmiles(mol, rootedAtAtom=atom.GetIdx(), canonical=True, isomericSmiles=True),
                        "rooted_identity": Chem.MolToSmiles(identity, rootedAtAtom=atom.GetIdx(), canonical=True, isomericSmiles=True),
                    })
    # Numeric group names may be changed without changing chemistry. Compare
    # kind and member partitions, using IDs only to assemble each partition.
    result = [{"kind": key[0], "members": sorted(members, key=lambda v: json.dumps(v, sort_keys=True))}
              for key, members in grouped.items()]
    return sorted(result, key=lambda value: json.dumps(value, sort_keys=True))


def member_identity(member):
    return {key: member[key] for key in ["role", "map", "isotope", "tetrahedral", "rooted_identity"]}


def comparison_key(groups):
    result = []
    for group in groups:
        members = sorted(group["members"], key=lambda m: json.dumps(member_identity(m), sort_keys=True))
        codes = tuple(member["cip"] or "?" for member in members)
        if "ABSOLUTE" not in group["kind"]:
            inverse = tuple({"R": "S", "S": "R", "r": "s", "s": "r"}.get(c, c) for c in codes)
            codes = min(codes, inverse)
        result.append(json.dumps([group["kind"], [member_identity(m) for m in members], codes], sort_keys=True))
    return sorted(result)


def flat_members(groups):
    return sorted([json.dumps([group["kind"], member_identity(m)], sort_keys=True)
                   for group in groups for m in group["members"]])


def observe(paths):
    from rdkit import rdBase
    packet = read(FIXTURE)
    original = packet["positive"][0]["input"]
    reordered = packet["positive"][1]["input"]
    assert fingerprint(original) == fingerprint(reordered)
    assert fingerprint(original) == fingerprint(original.replace("&1:", "&7:"))
    inverted = original.replace("[C@H]", "[C@@H]")
    assert comparison_key(fingerprint(original)) == comparison_key(fingerprint(inverted))
    assert comparison_key(fingerprint(original.replace("&1:", "a:"))) != comparison_key(fingerprint(inverted.replace("&1:", "a:")))
    assert fingerprint(original) != fingerprint("CC[C@H](C)O.O>>CCC(C)=O |&1:3|")
    two = "[CH3:10][C@H:11](F)[C@@H:12](O)Cl>>C "
    assert fingerprint(two + "|&1:1|") != fingerprint(two + "|&1:3|")
    reports = []
    for path in paths:
        capture = read(path)
        report = {"capture": str(path), "capture_sha256": sha(path), "engine": capture["engine"], "available": capture["available"], "cases": []}
        for case in capture["cases"]:
            row = {"id": case["id"], "canonical": case["canonical"]}
            try:
                before = fingerprint(case["input"])
                if "error" in case:
                    row.update(status="writer_refused", error=case["error"])
                else:
                    after = fingerprint(case["output"])
                    status = "match" if comparison_key(before) == comparison_key(after) else (
                        "group_partition_difference" if flat_members(before) == flat_members(after)
                        else "atom_membership_difference")
                    row.update(status=status, before=before, after=after)
            except Exception as error:
                row.update(status="observer_refused", error=str(error))
            report["cases"].append(row)
        report["counts"] = {status: sum(c["status"] == status for c in report["cases"]) for status in ["match", "atom_membership_difference", "group_partition_difference", "writer_refused", "observer_refused"]}
        reports.append(report)
    return {"schema": "chematic.reaction-cx-observer.v1", "observer_rdkit": rdBase.rdkitVersion, "self_controls": 6, "reports": reports}


def capture(engine, python_root):
    if python_root:
        sys.path.insert(0, python_root)
    if engine == "chematic":
        import chematic
        from chematic import chematic as native
        available = hasattr(chematic, "normalize_reaction_cxsmiles")
        runtime = {"path": native.__file__, "sha256": sha(native.__file__)}
        normalize = lambda text, canonical: chematic.normalize_reaction_cxsmiles(text, canonical=canonical)
    else:
        from rdkit import rdBase
        from rdkit.Chem import rdChemReactions as reactions
        available = True
        runtime = {"version": rdBase.rdkitVersion, "path": reactions.__file__, "sha256": sha(reactions.__file__)}
        normalize = lambda text, canonical: reactions.ReactionToCXSmiles(reactions.ReactionFromSmarts(text, useSmiles=True), canonical=canonical)
    cases = []
    if available:
        for case in read(FIXTURE)["positive"]:
            for canonical in [False, True]:
                row = {"id": case["id"], "input": case["input"], "canonical": canonical}
                try:
                    row["output"] = normalize(case["input"], canonical)
                except Exception as error:
                    row["error"] = str(error)
                cases.append(row)
    return {"schema": "chematic.reaction-cx-capture.v1", "engine": engine, "available": available, "runtime": runtime, "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--engine", choices=["chematic", "rdkit"])
    parser.add_argument("--python-root")
    mode.add_argument("--observe", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = observe(args.observe) if args.observe else capture(args.engine, args.python_root)
    result.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), platform=platform.platform(), fixture_sha256=sha(FIXTURE), observer_sha256=sha(__file__))
    payload = (json.dumps(result, indent=2, allow_nan=False) + "\n").encode()
    args.output.write_bytes(gzip.compress(payload, mtime=0) if args.output.suffix == ".gz" else payload)
    if args.observe:
        for report in result["reports"]:
            print(report["capture"], report["counts"])
    else:
        print(args.engine, "available:", result["available"], "captured:", len(result["cases"]))


if __name__ == "__main__":
    main()
