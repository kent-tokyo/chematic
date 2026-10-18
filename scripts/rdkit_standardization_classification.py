"""Shared residual classifiers for RDKit standardization comparisons.

The labels in this module describe compatibility boundaries; they are not an
independent chemical-correctness oracle.  Keeping the source-built and Python
binding lanes on one implementation prevents their residual taxonomies from
drifting apart.
"""

from __future__ import annotations

from rdkit import Chem

METALS = {
    3, 4, 11, 12, 13, 19, 20, *range(21, 32), *range(37, 51),
    *range(55, 84), *range(87, 117),
}

NON_TETRAHEDRAL_TAGS = {
    Chem.ChiralType.CHI_SQUAREPLANAR,
    Chem.ChiralType.CHI_TRIGONALBIPYRAMIDAL,
    Chem.ChiralType.CHI_OCTAHEDRAL,
}


def smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def graph_key(text: str | None, *, keep_charge: bool, keep_bond_order: bool) -> str | None:
    """Return a deterministic RDKit graph key for residual classification."""
    if text is None:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    copy = Chem.RWMol(mol)
    for atom in copy.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        atom.SetIsotope(0)
        if not keep_charge:
            atom.SetFormalCharge(0)
    for bond in copy.GetBonds():
        bond.SetStereo(Chem.BondStereo.STEREONONE)
        if not keep_bond_order:
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for atom in copy.GetAtoms():
        atom.SetIsAromatic(False)
    return Chem.MolToSmiles(copy, canonical=True, isomericSmiles=False)


def classify_cleanup_mismatch(
    source: Chem.Mol,
    normalized: str | None,
    expected: str | None,
    actual: str | None,
    cleanup_error: str | None,
) -> tuple[str, str]:
    """Classify one Cleanup residual without treating it as an exact match."""
    if cleanup_error is not None or actual is None:
        if cleanup_error and "kekulize" in cleanup_error.lower():
            return "aromatic_n_oxide_kekulization", "not_implemented"
        return "operation_error", "unsupported_or_invalid_output"

    has_metal = any(atom.GetAtomicNum() in METALS for atom in source.GetAtoms())
    if has_metal:
        if any(atom.GetChiralTag() in NON_TETRAHEDRAL_TAGS for atom in source.GetAtoms()):
            return "metal_stereo_preservation", "policy_difference"
        exact_graph = graph_key(actual, keep_charge=True, keep_bond_order=True)
        expected_graph = graph_key(expected, keep_charge=True, keep_bond_order=True)
        neutral_graph = graph_key(actual, keep_charge=False, keep_bond_order=True)
        expected_neutral_graph = graph_key(expected, keep_charge=False, keep_bond_order=True)
        if exact_graph is not None and exact_graph == expected_graph:
            return "canonical_representation", "representation_only"
        if neutral_graph is not None and neutral_graph == expected_neutral_graph:
            return "metal_charge_assignment", "not_implemented"
        topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
        expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
        if topology is not None and topology == expected_topology:
            return "metal_bond_order_assignment", "not_implemented"
        return "metal_disconnection_policy", "not_implemented"

    if normalized != expected:
        return "acid_base_reionization", "not_implemented"

    if expected is not None:
        if "[S+]" in actual and "[O-]" in actual and "[S+]" not in expected:
            return "sulfur_oxide_representation", "representation_only"
        if "[H]=C" in smiles(source) or "[CH]" in actual:
            return "explicit_hydrogen_normalization", "not_implemented"

    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "functional_group_normalization", "not_implemented"
    return "unclassified_semantic_difference", "needs_adjudication"


def classify_reionize_mismatch(
    source: str | None, expected: str | None, actual: str | None, error: str | None
) -> tuple[str, str]:
    if error is not None or actual is None:
        return "operation_error", "unsupported_or_invalid_output"
    if expected == source:
        return "independent_ionization_policy", "policy_difference"
    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "acid_base_charge_assignment", "not_implemented"
    return "hydrogen_or_resonance_assignment", "not_implemented"


def classify_normalize_mismatch(
    expected: str | None, actual: str | None, error: str | None
) -> tuple[str, str]:
    if error is not None or actual is None:
        if error and "kekulize" in error.lower():
            return "aromatic_n_oxide_kekulization", "not_implemented"
        return "operation_error", "unsupported_or_invalid_output"
    topology = graph_key(actual, keep_charge=False, keep_bond_order=False)
    expected_topology = graph_key(expected, keep_charge=False, keep_bond_order=False)
    if topology is not None and topology == expected_topology:
        return "functional_group_normalization", "not_implemented"
    return "normalization_graph_assignment", "not_implemented"
