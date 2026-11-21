#!/usr/bin/env python3
"""Inspect in-memory Cleanup differences on the recorded 47 rows.

Requires the published chematic 1.2.2 wheel, RDKit 2026.09.1, and NetworkX.
"""

import argparse
import collections
import hashlib
import json
from pathlib import Path
import warnings
import zipfile

import chematic
import networkx as nx
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "validation/results/rdkit-cleanup-v1.2.2-source-stress-2026-11-21.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--wheel", type=Path, required=True)
    args = parser.parse_args()
    if chematic.__version__ != "1.2.2" or rdBase.rdkitVersion != "2026.09.1":
        raise ValueError("expected chematic 1.2.2 and RDKit 2026.09.1")
    package_root = Path(chematic.__file__).parent.parent
    with zipfile.ZipFile(args.wheel) as archive:
        files = [name for name in archive.namelist() if name.startswith("chematic/") and not name.endswith("/")]
        if not files:
            raise ValueError("wheel contains no chematic package files")
        for name in files:
            path = package_root / name
            if not path.is_file() or path.read_bytes() != archive.read(name):
                raise ValueError(f"installed package differs from wheel: {path}")
    RDLogger.DisableLog("rdApp.*")
    rows = []
    for case in json.loads(RECORD.read_text())["mismatches"]:
        source_input = chematic.from_smiles(case["input"])
        oracle_input = Chem.MolFromSmiles(case["input"])
        source_input_edges = {tuple(sorted((a, b))) for a, b, _, _ in source_input.bond_table}
        oracle_input_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in oracle_input.GetBonds()}
        if source_input_edges != oracle_input_edges:
            raise ValueError(f"input atom alignment changed for index {case['index']}")
        source = source_input.cleanup()
        oracle = rdMolStandardize.Cleanup(oracle_input)
        if source.rdkit_smiles != case["chematic_cleanup"]:
            raise ValueError(f"published wheel differs from recorded source output for index {case['index']}")
        source_atoms = source.atom_table
        oracle_atoms = list(oracle.GetAtoms())
        if len(source_atoms) != len(oracle_atoms) or any(
            atom[1] != reference.GetAtomicNum() for atom, reference in zip(source_atoms, oracle_atoms, strict=True)
        ):
            raise ValueError(f"atom alignment changed for index {case['index']}")
        changes = []
        for index, (atom, reference) in enumerate(zip(source_atoms, oracle_atoms, strict=True)):
            source_state = [atom[2], atom[4]]
            oracle_state = [reference.GetFormalCharge(), reference.GetTotalNumHs()]
            if source_state != oracle_state:
                changes.append({"atom": index, "element": atom[0],
                                "chematic_charge_h": source_state, "rdkit_charge_h": oracle_state})
        source_edges = {tuple(sorted((a, b))) for a, b, _, _ in source.bond_table}
        oracle_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in oracle.GetBonds()}
        source_graph = nx.Graph()
        oracle_graph = nx.Graph()
        for index, (atom, reference) in enumerate(zip(source_atoms, oracle_atoms, strict=True)):
            source_graph.add_node(index, state=(atom[1], atom[2], atom[4]))
            oracle_graph.add_node(index, state=(reference.GetAtomicNum(), reference.GetFormalCharge(), reference.GetTotalNumHs()))
        source_graph.add_edges_from(source_edges)
        oracle_graph.add_edges_from(oracle_edges)
        isomorphic = nx.is_isomorphic(source_graph, oracle_graph,
                                      node_match=lambda a, b: a["state"] == b["state"])
        bond_isomorphic = None
        oracle_text = Chem.MolToSmiles(oracle, canonical=True, isomericSmiles=True)
        oracle_reparsed = Chem.MolFromSmiles(oracle_text)
        source_reparsed = Chem.MolFromSmiles(case["chematic_cleanup"])
        source_second = source.cleanup().rdkit_smiles
        oracle_iterations = [oracle_text]
        current_oracle = oracle
        for _ in range(5):
            current_oracle = rdMolStandardize.Cleanup(current_oracle)
            oracle_iterations.append(Chem.MolToSmiles(current_oracle, canonical=True, isomericSmiles=True))
        oracle_second = oracle_iterations[1]
        sanitized_equal = None
        if isomorphic and source_reparsed is not None:
            sanitized_oracle = Chem.Mol(oracle)
            try:
                Chem.SanitizeMol(sanitized_oracle)
                sanitized_equal = (
                    Chem.MolToSmiles(sanitized_oracle, canonical=True, isomericSmiles=True)
                    == Chem.MolToSmiles(source_reparsed, canonical=True, isomericSmiles=True)
                )
            except (ValueError, RuntimeError):
                pass
        if isomorphic:
            # MOL conversion can lose stereo, so this check only compares the
            # molecular graph and leaves stereo and radical review separate.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rendered = Chem.MolFromMolBlock(source.to_mol_block(), sanitize=False, removeHs=False)
            if rendered is not None and rendered.GetNumAtoms() == len(source_atoms):
                try:
                    Chem.SanitizeMol(rendered)
                    for index, (atom, reference) in enumerate(zip(rendered.GetAtoms(), oracle_atoms, strict=True)):
                        if atom.GetAtomicNum() != source_atoms[index][1]:
                            raise ValueError("MOL atom ordering changed")
                        source_graph.nodes[index]["isotope"] = atom.GetIsotope()
                        oracle_graph.nodes[index]["isotope"] = reference.GetIsotope()
                    for bond in rendered.GetBonds():
                        source_graph.edges[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]["order"] = str(bond.GetBondType())
                    for bond in oracle.GetBonds():
                        oracle_graph.edges[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]["order"] = str(bond.GetBondType())
                    bond_isomorphic = nx.is_isomorphic(
                        source_graph, oracle_graph,
                        node_match=lambda a, b: a["state"] == b["state"] and a["isotope"] == b["isotope"],
                        edge_match=lambda a, b: a["order"] == b["order"],
                    )
                except (ValueError, RuntimeError, KeyError):
                    bond_isomorphic = None
        rows.append({"index": case["index"], "recorded_category": case["classification"],
                     "charge_or_h_differences": changes,
                     "chematic_only_bonds": sorted(source_edges - oracle_edges),
                     "rdkit_only_bonds": sorted(oracle_edges - source_edges),
                     "charge_h_connectivity_isomorphic": isomorphic,
                     "bond_isotope_graph_isomorphic": bond_isomorphic,
                     "same_after_rdkit_sanitize": sanitized_equal,
                     "rdkit_serialization_reparseable": oracle_reparsed is not None,
                     "chematic_second_cleanup_matches": source_second == source.rdkit_smiles,
                     "rdkit_second_cleanup_matches": oracle_second == oracle_text,
                     "rdkit_six_cleanup_strings": oracle_iterations,
                     "rdkit_second_cleanup_matches_reparsed_chematic": (
                         None if source_reparsed is None else
                         oracle_second == Chem.MolToSmiles(source_reparsed, canonical=True, isomericSmiles=True)
                     ),
                     "rdkit_in_memory_radical_electrons": sum(atom.GetNumRadicalElectrons() for atom in oracle_atoms),
                     "chematic_serialized_radical_electrons": None if source_reparsed is None else sum(
                         atom.GetNumRadicalElectrons() for atom in source_reparsed.GetAtoms()
                     ),
                     "assessment": ("chemical_difference" if not isomorphic or bond_isomorphic is False
                                    else "representation_difference_after_sanitize" if sanitized_equal
                                    else "requires_stereo_radical_review")})
    counts = collections.Counter(row["assessment"] for row in rows)
    for assessment in ("chemical_difference", "representation_difference_after_sanitize", "refusal", "requires_stereo_radical_review"):
        counts.setdefault(assessment, 0)
    result = {"source_record": str(RECORD.relative_to(ROOT)), "chematic_version": chematic.__version__,
              "wheel": args.wheel.name, "wheel_sha256": hashlib.sha256(args.wheel.read_bytes()).hexdigest(),
              "installed_package_files_verified": len(files),
              "rdkit_version": rdBase.rdkitVersion, "networkx_version": nx.__version__,
              "method": "compare element, formal charge, attached hydrogen and connectivity up to graph isomorphism; for matching graphs compare isotope and aromatic-normalized bond orders after MOL conversion. Representation differences require matching canonical isomeric SMILES after RDKit sanitizes its in-memory Cleanup product and reparses the chematic output; the original unsanitized radical state is recorded separately.",
              "counts": dict(counts), "rows": rows}
    output = ROOT / "validation/results/cleanup-v1.2.2-47-atom-audit-2026-11-21.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(result["counts"])


if __name__ == "__main__":
    main()
