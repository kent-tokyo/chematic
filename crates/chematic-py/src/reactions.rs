//! Reaction (SMIRKS/MMP/BRICS/MCS), enumeration, and reaction-metric bindings.

use crate::Mol;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::sync::Arc;

/// Test whether a reaction SMARTS pattern matches a reaction SMILES.
///
/// ``smarts``: reaction SMARTS (e.g. ``"[OH:1]>>[O-:1]"``).
/// ``reaction_smiles``: reaction SMILES in ``"R>>P"`` or ``"R>A>P"`` format.
///
///     ok = chematic.reaction_smarts_match("[OH]>>[O-]", "CCO>>CC[O-]")
#[pyfunction]
fn reaction_smarts_match(smarts: &str, reaction_smiles: &str) -> PyResult<bool> {
    let query = chematic_rxn::parse_reaction_query(smarts)
        .map_err(|e| PyValueError::new_err(format!("invalid reaction SMARTS: {e}")))?;
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(format!("invalid reaction SMILES: {e}")))?;
    Ok(chematic_rxn::has_reaction_substructure_match(&rxn, &query))
}

/// Find all Matched Molecular Pairs (MMP) in a list of SMILES strings.
///
/// Returns a list of dicts, each with keys:
///   ``mol_a``, ``mol_b`` (canonical SMILES), ``core`` (shared scaffold),
///   ``fragment_a``, ``fragment_b`` (substituent SMILES containing ``[*]``).
///
/// Uses BRICS single-bond cuts. Pairs are deduplicated.
///
///     pairs = chematic.find_mmp(["c1ccccc1", "Cc1ccccc1", "Nc1ccccc1"])
///     for p in pairs:
///         print(f"{p['fragment_a']} → {p['fragment_b']} on {p['core']}")
#[pyfunction]
fn find_mmp<'py>(smiles: Vec<String>, py: Python<'py>) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    let refs: Vec<&chematic_core::Molecule> = mols.iter().collect();
    chematic_chem::find_mmp(&refs)
        .into_iter()
        .map(|pair| {
            let d = PyDict::new(py);
            d.set_item("mol_a", &pair.mol_a)?;
            d.set_item("mol_b", &pair.mol_b)?;
            d.set_item("core", &pair.core)?;
            d.set_item("fragment_a", &pair.fragment_a)?;
            d.set_item("fragment_b", &pair.fragment_b)?;
            Ok(d)
        })
        .collect()
}

/// R-group decomposition — split molecules into a scaffold core and variable R-groups.
///
/// ``scaffold_smarts``: SMARTS pattern defining the common scaffold.
/// ``mols``: list of :class:`Mol` objects to decompose.
///
/// Returns a list of dicts (one per input molecule), or ``None`` when the scaffold
/// does not match a particular molecule.  Each dict contains:
///   ``mol_idx``   (int)  — index in the input list,
///   ``core``      (str)  — scaffold SMILES with ``[*]`` at attachment points,
///   ``R1``, ``R2``, … (str) — SMILES for each R-group (``[*]`` marks attachment).
///
///     mols = [chematic.from_smiles(s) for s in ["CCc1ccccc1", "CCCc1ccccc1"]]
///     results = chematic.rgroup_decompose("c1ccccc1", mols)
///     # [{"mol_idx": 0, "core": "...", "R1": "[*]CC"}, ...]
#[pyfunction]
fn rgroup_decompose<'py>(
    scaffold_smarts: &str,
    mols: Vec<Mol>,
    py: Python<'py>,
) -> PyResult<Vec<Option<Bound<'py, PyDict>>>> {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let results = chematic_chem::rgroup_decompose(scaffold_smarts, &refs)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    results
        .into_iter()
        .map(|opt| match opt {
            None => Ok(None),
            Some(r) => {
                let d = PyDict::new(py);
                d.set_item("mol_idx", r.mol_idx)?;
                d.set_item("core", &r.core_smiles)?;
                for (k, v) in &r.r_groups {
                    d.set_item(format!("R{k}"), v)?;
                }
                Ok(Some(d))
            }
        })
        .collect()
}

/// Detect activity cliffs in a set of molecules with known activity values.
///
/// An activity cliff is a structurally similar pair with a large activity difference —
/// a classic signal of SAR sensitivity. Common in MolScore and mol-eval type analyses.
///
/// ``mols``: list of :class:`Mol` objects.
/// ``activities``: list of floats (one per mol), e.g. pIC50 values.
/// ``sim_threshold``: minimum ECFP4 Tanimoto similarity to consider a pair (default 0.65).
/// ``cliff_delta``: minimum ``|activity_i − activity_j|`` to be a cliff (default 2.0).
///
/// Returns a list of dicts sorted by similarity descending, each containing:
///   ``mol_a_idx`` (int), ``mol_b_idx`` (int), ``similarity`` (float), ``activity_delta`` (float).
///
///     mols = [chematic.from_smiles(s) for s in ["c1ccccc1", "Cc1ccccc1"]]
///     cliffs = chematic.activity_cliffs(mols, [5.0, 8.5], sim_threshold=0.0, cliff_delta=2.0)
///     # [{"mol_a_idx": 0, "mol_b_idx": 1, "similarity": 0.xx, "activity_delta": 3.5}]
#[pyfunction]
#[pyo3(signature = (mols, activities, sim_threshold = 0.65, cliff_delta = 2.0))]
fn activity_cliffs<'py>(
    mols: Vec<Mol>,
    activities: Vec<f64>,
    sim_threshold: f32,
    cliff_delta: f64,
    py: Python<'py>,
) -> PyResult<Vec<Bound<'py, PyDict>>> {
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let cliffs = chematic_chem::activity_cliffs(&refs, &activities, sim_threshold, cliff_delta);
    cliffs
        .into_iter()
        .map(|c| {
            let d = PyDict::new(py);
            d.set_item("mol_a_idx", c.mol_a_idx)?;
            d.set_item("mol_b_idx", c.mol_b_idx)?;
            d.set_item("similarity", c.similarity)?;
            d.set_item("activity_delta", c.activity_delta)?;
            Ok(d)
        })
        .collect()
}

/// Identify the reaction center: bonds broken/formed and atoms changed.
///
/// Returns a dict with keys:
///   ``broken_bonds`` (list of ``[i, j]`` atom index pairs),
///   ``formed_bonds`` (list of ``[i, j]`` atom index pairs),
///   ``changed_atoms`` (list of atom indices).
///
/// Atom indices use reactant-side numbering.
///
///     rc = chematic.find_reaction_center("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
///     print("broken:", rc["broken_bonds"])
///     print("formed:", rc["formed_bonds"])
#[pyfunction]
fn find_reaction_center<'py>(
    reaction_smiles: &str,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let center = chematic_rxn::find_reaction_center(&rxn);
    let d = PyDict::new(py);
    let broken: Vec<[usize; 2]> = center
        .broken_bonds
        .iter()
        .map(|(a, b)| [a.0 as usize, b.0 as usize])
        .collect();
    let formed: Vec<[usize; 2]> = center
        .formed_bonds
        .iter()
        .map(|(a, b)| [a.0 as usize, b.0 as usize])
        .collect();
    let changed: Vec<usize> = center.changed_atoms.iter().map(|a| a.0 as usize).collect();
    d.set_item("broken_bonds", broken)?;
    d.set_item("formed_bonds", formed)?;
    d.set_item("changed_atoms", changed)?;
    Ok(d)
}

/// E-factor (Environmental Factor) — waste-to-product mass ratio.
///
/// E-factor = waste_mass / product_mass.  Lower is greener.
/// Fine chemicals typically E=5–50; pharmaceuticals E=25–100.
///
///     ef = chematic.e_factor(waste_kg=90.0, product_kg=10.0)  # → 9.0
#[pyfunction]
fn e_factor(waste_mass: f64, product_mass: f64) -> f64 {
    chematic_rxn::e_factor(waste_mass, product_mass)
}

/// Process Mass Intensity (PMI) — total mass used per unit product mass.
///
/// PMI = (sum of all input masses) / product_mass. Lower is greener.
///
///     pmi = chematic.pmi_rxn([solvent_kg, reagent1_kg, reagent2_kg], product_kg)
#[pyfunction]
fn pmi_rxn(all_masses: Vec<f64>, product_mass: f64) -> f64 {
    chematic_rxn::pmi_rxn(&all_masses, product_mass)
}

/// Reaction Mass Efficiency (RME) — fraction of reactant mass in the product.
///
/// RME = product_mass / sum(reactant_masses). Range [0, 1].
///
///     rme = chematic.reaction_mass_efficiency([reactant1_g, reactant2_g], product_g)
#[pyfunction]
fn reaction_mass_efficiency(reactant_masses: Vec<f64>, product_mass: f64) -> f64 {
    chematic_rxn::reaction_mass_efficiency(&reactant_masses, product_mass)
}

/// A value of 100% means all atoms in reactants appear in the product.
///
///     ae = chematic.atom_economy("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
#[pyfunction]
fn atom_economy(reaction_smiles: &str) -> PyResult<f64> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_rxn::atom_economy(&rxn))
}

/// Check whether a reaction SMILES is atom-balanced.
///
/// Returns a dict with keys:
///   ``balanced`` (bool), ``diff`` (list of str describing imbalances).
///
///     result = chematic.balance_check("C+O>>CO")
///     print(result["balanced"], result["diff"])
#[pyfunction]
fn balance_check<'py>(reaction_smiles: &str, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = chematic_rxn::balance_check(&rxn);
    let d = PyDict::new(py);
    d.set_item("balanced", result.balanced)?;
    d.set_item("diff", result.diff())?;
    Ok(d)
}

/// Enumerate a combinatorial library from a SMIRKS template and fragment sets.
///
/// Args:
///     smirks: Reaction SMIRKS template (e.g. ``"[C:1]Cl.[N:2]>>[C:1][N:2]"``).
///     fragment_sets: List of SMILES lists — one list per reactant slot.
///                    All combinations across sets are generated.
///     max_size: Maximum library size (default 1_000_000).
///
/// Returns a list of product SMILES strings.
///
///     products = chematic.enumerate_library(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         [["c1ccccc1", "CC"], ["N", "CN"]],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, fragment_sets, max_size = 1_000_000))]
fn enumerate_library(
    smirks: &str,
    fragment_sets: Vec<Vec<String>>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parsed_sets: Vec<Vec<chematic_core::Molecule>> = fragment_sets
        .iter()
        .map(|set| {
            set.iter()
                .filter_map(|s| chematic_smiles::parse(s).ok())
                .collect()
        })
        .collect();
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library(smirks, parsed_sets, &config)
        .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Enumerate a 2-fragment combinatorial library (scaffold × building block).
///
/// Convenience alternative to ``enumerate_library(smirks, [scaffolds, building_blocks])``.
/// The most common combinatorial chemistry pattern: one scaffold set reacted with one
/// building-block set to produce all pairwise products.
///
///     products = chematic.enumerate_library_2way(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         scaffolds=["c1ccccc1C(=O)Cl", "CC(=O)Cl"],
///         building_blocks=["N", "CN"],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, scaffolds, building_blocks, max_size = 1_000_000))]
fn enumerate_library_2way(
    smirks: &str,
    scaffolds: Vec<String>,
    building_blocks: Vec<String>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parse_smiles_set = |set: Vec<String>| -> Vec<chematic_core::Molecule> {
        set.iter()
            .filter_map(|s| chematic_smiles::parse(s).ok())
            .collect()
    };
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library_2way(
        smirks,
        parse_smiles_set(scaffolds),
        parse_smiles_set(building_blocks),
        &config,
    )
    .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
    .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Enumerate a 3-fragment combinatorial library (scaffold × R1 × R2).
///
/// Convenience alternative to ``enumerate_library(smirks, [scaffolds, r1_set, r2_set])``.
/// Covers the common scaffold-decoration pattern with two variable positions.
///
///     products = chematic.enumerate_library_3way(
///         "[C:1]C(=O)Cl.[N:2]>>[C:1]C(=O)[N:2]",
///         scaffolds=["CC(=O)Cl"],
///         r1_set=["N", "CN"],
///         r2_set=["c1ccccc1", "CC"],
///     )
#[pyfunction]
#[pyo3(signature = (smirks, scaffolds, r1_set, r2_set, max_size = 1_000_000))]
fn enumerate_library_3way(
    smirks: &str,
    scaffolds: Vec<String>,
    r1_set: Vec<String>,
    r2_set: Vec<String>,
    max_size: usize,
) -> PyResult<Vec<String>> {
    let parse_smiles_set = |set: Vec<String>| -> Vec<chematic_core::Molecule> {
        set.iter()
            .filter_map(|s| chematic_smiles::parse(s).ok())
            .collect()
    };
    let config = chematic_rxn::LibraryConfig {
        skip_failures: true,
        max_size: Some(max_size),
    };
    chematic_rxn::enumerate_library_3way(
        smirks,
        parse_smiles_set(scaffolds),
        parse_smiles_set(r1_set),
        parse_smiles_set(r2_set),
        &config,
    )
    .map(|mols| mols.iter().map(chematic_smiles::canonical_smiles).collect())
    .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Compute Tanimoto similarity between two reaction SMILES using reaction fingerprints.
///
/// Parses both reaction SMILES, computes reaction fingerprints, and returns
/// the Tanimoto coefficient.
///
///     sim = chematic.tanimoto_reaction_fp("CC>>CO", "c1ccccc1>>c1ccccc1N")
///
/// Raises ``ValueError`` on invalid reaction SMILES.
#[pyfunction]
fn tanimoto_reaction_fp(rxn1: &str, rxn2: &str) -> PyResult<f64> {
    let r1 =
        chematic_rxn::parse_reaction(rxn1).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let r2 =
        chematic_rxn::parse_reaction(rxn2).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_fp::tanimoto_reaction_fp(&r1, &r2))
}

/// Check whether a reaction SMILES matches a reaction SMARTS pattern.
///
/// Returns ``True`` if the reaction matches the query pattern, ``False`` otherwise.
/// Equivalent to ``chematic.reaction_smarts_match`` but returns a simple bool
/// via SMARTS-based pattern matching rather than substructure query.
///
///     matched = chematic.query_reaction("CC>>CO", "[C:1]>>[C:1]O")
///
/// Raises ``ValueError`` on invalid reaction SMILES or SMARTS.
#[pyfunction]
fn query_reaction(reaction_smiles: &str, smarts: &str) -> PyResult<bool> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = chematic_rxn::query_reaction(&rxn, smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(result.is_complete_match)
}

/// Query a list of reaction SMILES against a single SMARTS pattern.
///
/// Returns a dict:
///   - ``total`` (int): total reactions processed
///   - ``matching`` (int): reactions that matched
///   - ``match_pct`` (float): match percentage (0–100)
///   - ``matches`` (list[(int, bool)]): per-reaction results as (original_index, matched)
///
/// Invalid SMILES are silently skipped (their indices will not appear in ``matches``).
///
/// Raises ``ValueError`` on invalid SMARTS.
///
///     rxns = ["CC>>CO", "CCCC>>CCCCO", "c1ccccc1>>c1ccccc1N"]
///     r = chematic.batch_query_reactions(rxns, "[C:1]>>[C:1]O")
///     print(r["matching"], "/", r["total"])
#[pyfunction]
fn batch_query_reactions<'py>(
    reactions: Vec<String>,
    smarts: &str,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut valid: Vec<chematic_rxn::Reaction> = Vec::new();
    let mut original_indices: Vec<usize> = Vec::new();
    for (i, s) in reactions.iter().enumerate() {
        if let Ok(rxn) = chematic_rxn::parse_reaction(s) {
            valid.push(rxn);
            original_indices.push(i);
        }
    }
    let result = chematic_rxn::batch_query_reactions(&valid, smarts)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Remap indices from filtered slice back to original input indices
    let matches: Vec<(usize, bool)> = result
        .matches
        .iter()
        .map(|&(idx, matched)| (original_indices[idx], matched))
        .collect();
    let d = PyDict::new(py);
    d.set_item("total", result.total_reactions)?;
    d.set_item("matching", result.matching_reactions)?;
    d.set_item("match_pct", result.match_percentage)?;
    d.set_item("matches", matches)?;
    Ok(d)
}

/// Render a reaction SMILES as an SVG diagram.
///
/// Returns an SVG string showing reactants → products with an arrow.
/// Equivalent to RDKit's ``Draw.ReactionToImage(rxn)``.
///
///     svg = chematic.reaction_svg("CC(=O)Cl.[NH3]>>CC(=O)N.HCl")
///     with open("reaction.svg", "w") as f:
///         f.write(svg)
///
/// Raises ``ValueError`` on invalid reaction SMILES.
#[pyfunction]
fn reaction_svg(reaction_smiles: &str) -> PyResult<String> {
    let rxn = chematic_rxn::parse_reaction(reaction_smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chematic_depict::depict_reaction_svg(&rxn))
}

/// Compute scaffold network statistics across a molecule library.
///
/// Returns a dict with three parallel lists:
///   - ``scaffolds``: canonical SMILES of each unique scaffold
///   - ``counts``: how many input molecules contain each scaffold
///   - ``parents``: index of the parent (simpler) scaffold, or ``None`` for root
///
/// Invalid SMILES are silently skipped.
///
///     result = chematic.scaffold_network_counts(smiles_list)
///     for smi, n in zip(result["scaffolds"], result["counts"]):
///         print(smi, n)
#[pyfunction]
fn scaffold_network_counts<'py>(
    smiles: Vec<String>,
    py: Python<'py>,
) -> PyResult<Bound<'py, PyDict>> {
    let mols: Vec<chematic_core::Molecule> = smiles
        .iter()
        .filter_map(|s| chematic_smiles::parse(s).ok())
        .collect();
    let net = chematic_chem::scaffold_network_with_counts(&mols);
    let scaffolds: Vec<String> = net
        .scaffolds
        .iter()
        .map(chematic_smiles::canonical_smiles)
        .collect();
    let d = PyDict::new(py);
    d.set_item("scaffolds", scaffolds)?;
    d.set_item("counts", net.counts)?;
    d.set_item("parents", net.parents)?;
    Ok(d)
}

/// Apply a SMIRKS reaction template to a list of reactant molecules.
///
/// Returns a list of product sets; each set is a list of Mol.
/// Raises ``ValueError`` on SMIRKS parse failure or reactant count mismatch.
///
/// **Stereochemistry**: when the reactant template contains ``@``/``@@`` stereo
/// descriptors, only reactant molecules whose chiral centres match the template
/// configuration are accepted (parity-aware comparison, SMILES write-order
/// independent). Templates without stereo descriptors match both enantiomers.
///
///     products = chematic.run_smirks("[OH:1]>>[O-:1]", [mol])
///     # → [[product_mol], ...]
///
///     # Stereo-selective: only L-amino acids match this template
///     l_products = chematic.run_smirks("[N:1][C@@H:2](C)C(=O)O>>[N:1].[C@@H:2](C)C(=O)O", [mol])
#[pyfunction]
fn run_smirks(smirks: &str, reactants: Vec<Mol>) -> PyResult<Vec<Vec<Mol>>> {
    for mol in &reactants {
        if mol.inner.atom_count() > 300 {
            return Err(PyValueError::new_err(
                "reactant too large for run_smirks (max 300 heavy atoms)",
            ));
        }
    }
    let refs: Vec<&chematic_core::Molecule> = reactants.iter().map(|m| m.inner.as_ref()).collect();
    chematic_rxn::run_reactants(smirks, &refs)
        .map(|sets| {
            sets.into_iter()
                .map(|set| set.into_iter().map(Mol::bare).collect())
                .collect()
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Apply a SMIRKS reaction template to a list of reactant molecules (strict mode).
///
/// Like :func:`run_smirks` but **does not carry substituents** into products.
/// Only atoms explicitly mapped in the product template are included.
/// Stereo filtering behaviour is identical to :func:`run_smirks`.
///
///     products = chematic.run_smirks_strict("[N:1][C:2]>>[N:1].[C:2]", [mol])
///     # → only the mapped N and C atoms; no R-groups attached
#[pyfunction]
fn run_smirks_strict(smirks: &str, reactants: Vec<Mol>) -> PyResult<Vec<Vec<Mol>>> {
    for mol in &reactants {
        if mol.inner.atom_count() > 300 {
            return Err(PyValueError::new_err(
                "reactant too large for run_smirks_strict (max 300 heavy atoms)",
            ));
        }
    }
    let refs: Vec<&chematic_core::Molecule> = reactants.iter().map(|m| m.inner.as_ref()).collect();
    chematic_rxn::run_reactants_strict(smirks, &refs)
        .map(|sets| {
            sets.into_iter()
                .map(|set| set.into_iter().map(Mol::bare).collect())
                .collect()
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Build a full `McsConfig` from the individual keyword arguments shared by
/// [`find_mcs`] and [`find_mcs_checked`] -- kept as one function so the two
/// entry points can never silently drift apart on how a given argument set
/// maps to `McsConfig`.
#[allow(clippy::too_many_arguments)]
fn build_mcs_config(
    match_bonds: bool,
    min_atoms: usize,
    timeout_ms: Option<u64>,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    atom_compare: &str,
    bond_compare: &str,
    match_chiral_tag: bool,
    match_charge: bool,
    match_isotope: bool,
    maximize_bonds: bool,
) -> PyResult<chematic_smarts::McsConfig> {
    use chematic_smarts::{AtomCompare, BondCompare, McsConfig};

    let atom_compare = match atom_compare {
        "elements" => AtomCompare::Elements,
        "any_heavy_atom" => AtomCompare::AnyHeavyAtom,
        "any" => AtomCompare::Any,
        other => {
            return Err(PyValueError::new_err(format!(
                "invalid atom_compare {other:?}: expected \"elements\", \"any_heavy_atom\", or \"any\""
            )));
        }
    };
    let bond_compare = match bond_compare {
        "order_or_aromatic" => BondCompare::OrderOrAromatic,
        "any" => BondCompare::Any,
        other => {
            return Err(PyValueError::new_err(format!(
                "invalid bond_compare {other:?}: expected \"order_or_aromatic\" or \"any\""
            )));
        }
    };
    Ok(McsConfig {
        match_bonds,
        min_atoms,
        timeout_ms,
        ring_matches_ring_only,
        complete_rings_only,
        atom_compare,
        bond_compare,
        match_chiral_tag,
        match_charge,
        match_isotope,
        maximize_bonds,
    })
}

/// Reconstruct a concrete `Mol` from a `QueryMolecule` produced by MCS search
/// (atom queries are always `AtomicNum` primitives; bond queries are typed
/// primitives) -- `None` when the query has no atoms (no common substructure).
fn qmol_to_mol(qmol: &chematic_smarts::QueryMolecule) -> Option<Mol> {
    use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
    use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery};

    if qmol.atoms.is_empty() {
        return None;
    }

    fn extract_atomic_num(q: &AtomQuery) -> Option<u8> {
        match q {
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
            AtomQuery::And(lhs, rhs) => extract_atomic_num(lhs).or_else(|| extract_atomic_num(rhs)),
            _ => None,
        }
    }

    // `build_query`/`molecule_to_query` never encode aromaticity as a per-atom
    // constraint (matches RDKit's own `CompareElements` representation) -- it's
    // carried entirely by the aromatic bond queries, so an atom is aromatic here
    // iff at least one of its query bonds is `BondPrimitive::Aromatic`.
    let mut aromatic_atoms = vec![false; qmol.atoms.len()];
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if matches!(
                qmol.bonds[*bond_idx].query,
                BondQuery::Primitive(BondPrimitive::Aromatic)
            ) {
                aromatic_atoms[atom_idx] = true;
                aromatic_atoms[*neighbor_idx] = true;
            }
        }
    }

    let mut builder = MoleculeBuilder::new();
    for (idx, qa) in qmol.atoms.iter().enumerate() {
        let elem = extract_atomic_num(&qa.query)
            .and_then(Element::from_atomic_number)
            .unwrap_or(Element::C);
        let mut atom = Atom::new(elem);
        atom.aromatic = aromatic_atoms[idx];
        builder.add_atom(atom);
    }
    for (atom_idx, neighbors) in qmol.adj.iter().enumerate() {
        for (bond_idx, neighbor_idx) in neighbors {
            if atom_idx < *neighbor_idx {
                let order = match &qmol.bonds[*bond_idx].query {
                    BondQuery::Primitive(BondPrimitive::Double) => BondOrder::Double,
                    BondQuery::Primitive(BondPrimitive::Triple) => BondOrder::Triple,
                    BondQuery::Primitive(BondPrimitive::Aromatic) => BondOrder::Aromatic,
                    _ => BondOrder::Single,
                };
                let _ = builder.add_bond(
                    AtomIdx(atom_idx as u32),
                    AtomIdx(*neighbor_idx as u32),
                    order,
                );
            }
        }
    }
    Some(Mol {
        inner: Arc::new(builder.build()),
        props: Default::default(),
    })
}

/// Find the Maximum Common Substructure (MCS) of a list of molecules.
///
/// Returns the MCS as a Mol, or ``None`` when there is no common substructure.
/// If ``timeout_ms`` is reached before the search finishes, returns the best
/// result found so far -- indistinguishable here from an exhaustive result;
/// use [`find_mcs_checked`] when that distinction matters.
///
///     mcs = chematic.find_mcs([mol1, mol2])
///     if mcs: print(mcs.smiles)
///
///     # Ring-aware scaffold extraction
///     scaffold = chematic.find_mcs(mols, ring_matches_ring_only=True, complete_rings_only=True)
///
///     # Ignore element identity, match any heavy atom (scaffold hopping)
///     core = chematic.find_mcs(mols, atom_compare="any_heavy_atom")
#[pyfunction]
#[pyo3(signature = (
    mols,
    match_bonds=true,
    min_atoms=1,
    timeout_ms=None,
    ring_matches_ring_only=false,
    complete_rings_only=false,
    atom_compare="elements",
    bond_compare="order_or_aromatic",
    match_chiral_tag=false,
    match_charge=false,
    match_isotope=false,
    maximize_bonds=true,
))]
#[allow(clippy::too_many_arguments)]
fn find_mcs(
    mols: Vec<Mol>,
    match_bonds: bool,
    min_atoms: usize,
    timeout_ms: Option<u64>,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    atom_compare: &str,
    bond_compare: &str,
    match_chiral_tag: bool,
    match_charge: bool,
    match_isotope: bool,
    maximize_bonds: bool,
) -> PyResult<Option<Mol>> {
    let config = build_mcs_config(
        match_bonds,
        min_atoms,
        timeout_ms,
        ring_matches_ring_only,
        complete_rings_only,
        atom_compare,
        bond_compare,
        match_chiral_tag,
        match_charge,
        match_isotope,
        maximize_bonds,
    )?;
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let qmol = chematic_smarts::find_mcs_with_config(&refs, &config);
    Ok(qmol_to_mol(&qmol))
}

/// Like [`find_mcs`], but also reports whether ``timeout_ms`` was reached
/// before the search finished exhaustively.
///
/// Returns a tuple ``(mcs, was_timed_out)``: ``mcs`` is the MCS as a Mol (or
/// ``None`` if there is no common substructure), and ``was_timed_out`` is
/// ``True`` if the search was cut off before proving ``mcs`` optimal.
///
///     mcs, timed_out = chematic.find_mcs_checked(mols, timeout_ms=500)
///     if timed_out:
///         print("warning: MCS may not be optimal")
#[pyfunction]
#[pyo3(signature = (
    mols,
    match_bonds=true,
    min_atoms=1,
    timeout_ms=None,
    ring_matches_ring_only=false,
    complete_rings_only=false,
    atom_compare="elements",
    bond_compare="order_or_aromatic",
    match_chiral_tag=false,
    match_charge=false,
    match_isotope=false,
    maximize_bonds=true,
))]
#[allow(clippy::too_many_arguments)]
fn find_mcs_checked(
    mols: Vec<Mol>,
    match_bonds: bool,
    min_atoms: usize,
    timeout_ms: Option<u64>,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    atom_compare: &str,
    bond_compare: &str,
    match_chiral_tag: bool,
    match_charge: bool,
    match_isotope: bool,
    maximize_bonds: bool,
) -> PyResult<(Option<Mol>, bool)> {
    let config = build_mcs_config(
        match_bonds,
        min_atoms,
        timeout_ms,
        ring_matches_ring_only,
        complete_rings_only,
        atom_compare,
        bond_compare,
        match_chiral_tag,
        match_charge,
        match_isotope,
        maximize_bonds,
    )?;
    let refs: Vec<&chematic_core::Molecule> = mols.iter().map(|m| m.inner.as_ref()).collect();
    let outcome = chematic_smarts::find_mcs_with_config_checked(&refs, &config);
    let was_timed_out = outcome.was_timed_out();
    let qmol = outcome.into_query();
    Ok((qmol_to_mol(&qmol), was_timed_out))
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reaction_smarts_match, m)?)?;
    m.add_function(wrap_pyfunction!(find_mmp, m)?)?;
    m.add_function(wrap_pyfunction!(rgroup_decompose, m)?)?;
    m.add_function(wrap_pyfunction!(activity_cliffs, m)?)?;
    m.add_function(wrap_pyfunction!(find_reaction_center, m)?)?;
    m.add_function(wrap_pyfunction!(e_factor, m)?)?;
    m.add_function(wrap_pyfunction!(pmi_rxn, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_mass_efficiency, m)?)?;
    m.add_function(wrap_pyfunction!(atom_economy, m)?)?;
    m.add_function(wrap_pyfunction!(balance_check, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library_2way, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_library_3way, m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto_reaction_fp, m)?)?;
    m.add_function(wrap_pyfunction!(query_reaction, m)?)?;
    m.add_function(wrap_pyfunction!(batch_query_reactions, m)?)?;
    m.add_function(wrap_pyfunction!(reaction_svg, m)?)?;
    m.add_function(wrap_pyfunction!(scaffold_network_counts, m)?)?;
    m.add_function(wrap_pyfunction!(run_smirks, m)?)?;
    m.add_function(wrap_pyfunction!(run_smirks_strict, m)?)?;
    m.add_function(wrap_pyfunction!(find_mcs, m)?)?;
    m.add_function(wrap_pyfunction!(find_mcs_checked, m)?)?;
    Ok(())
}
