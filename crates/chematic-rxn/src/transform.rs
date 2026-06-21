use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, STEREO_H_SENTINEL,
    validate_valence,
};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule, find_matches,
};

use crate::reaction::{RxnError, parse_reaction};

/// Error type for SMIRKS transformation.
#[derive(Debug)]
pub enum TransformError {
    SmirksParse(RxnError),
    ReactantCountMismatch { expected: usize, got: usize },
}

impl core::fmt::Display for TransformError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::SmirksParse(e) => write!(f, "SMIRKS parse error: {e}"),
            Self::ReactantCountMismatch { expected, got } => {
                write!(f, "reactant count mismatch: expected {expected}, got {got}")
            }
        }
    }
}

impl std::error::Error for TransformError {}

impl From<RxnError> for TransformError {
    fn from(e: RxnError) -> Self {
        Self::SmirksParse(e)
    }
}

/// Apply a SMIRKS template to input reactant molecules.
///
/// Returns all combinations of product sets — one per unique match across all
/// reactant templates.  Each inner `Vec<Molecule>` contains one product per
/// product component in the SMIRKS right-hand side.
///
/// Returns `Ok(vec![])` when no match is found.
///
/// Unmapped atoms attached to a mapped core atom (substituents) are
/// automatically carried through to the matching product template.
/// Use [`run_reactants_strict`] to return only mapped atoms.
pub fn run_reactants(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_impl(smirks, reactants, true)
}

/// Like [`run_reactants`] but **does not carry through substituents**.
///
/// Only atoms that appear explicitly in the product template (via atom maps or
/// new template atoms) are included in each product.  Unmapped neighbors of
/// core atoms are **not** collected via BFS.
///
/// Useful when the SMIRKS describes a complete molecule transformation and
/// you do not want R-group carry-through behaviour.
pub fn run_reactants_strict(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    run_reactants_impl(smirks, reactants, false)
}

fn run_reactants_impl(
    smirks: &str,
    reactants: &[&Molecule],
    carry_substituents: bool,
) -> Result<Vec<Vec<Molecule>>, TransformError> {
    let rxn = parse_reaction(smirks)?;

    let n_templates = rxn.reactants.len();
    if reactants.len() != n_templates {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: reactants.len(),
        });
    }

    // Build a QueryMolecule from each reactant template, and record the
    // atom-map number for each query atom index.
    let queries: Vec<QueryMolecule> = rxn.reactants.iter().map(mol_to_query).collect();
    let template_atom_maps: Vec<Vec<Option<u16>>> = rxn
        .reactants
        .iter()
        .map(|tmpl| {
            (0..tmpl.atom_count())
                .map(|i| tmpl.atom(AtomIdx(i as u32)).atom_map)
                .collect()
        })
        .collect();

    // Detect whether any reactant template carries @/@@ stereo, so we can apply
    // the parity-aware post-check after VF2 completes.  Chirality is NOT encoded
    // into the VF2 query because the raw flag comparison in eval_chirality is
    // SMILES-write-order-dependent; the correct check requires the full mapping
    // (see smirks_chirality_ok below).
    let has_stereo = rxn.reactants.iter().any(|r| {
        r.atoms().any(|(_, a)| a.chirality != Chirality::None)
    });

    // VF2 match: for each (template_query, input_mol) pair.
    let all_match_sets: Vec<Vec<HashMap<usize, AtomIdx>>> = queries
        .iter()
        .zip(reactants.iter())
        .map(|(q, mol)| find_matches(q, mol))
        .collect();

    // No products when any template has no match.
    if all_match_sets.iter().any(|ms| ms.is_empty()) {
        return Ok(vec![]);
    }

    let mut results: Vec<Vec<Molecule>> = Vec::new();

    for combo in cartesian_product(&all_match_sets) {
        // global_map: atom_map_number → (reactant_mol_idx, matched_AtomIdx)
        let mut global_map: HashMap<u16, (usize, AtomIdx)> = HashMap::new();
        for (ri, match_map) in combo.iter().enumerate() {
            for (&qi, &t_idx) in match_map {
                if let Some(am) = template_atom_maps[ri][qi] {
                    global_map.insert(am, (ri, t_idx));
                }
            }
        }

        // all_template_atoms: every (mol_idx, AtomIdx) matched by any reactant template atom.
        // Used as BFS walls to prevent substituent collection from crossing into the
        // template region, and to identify bonds that the product template replaces.
        let mut all_template_atoms: HashSet<(usize, AtomIdx)> = HashSet::new();
        for (ri, match_map) in combo.iter().enumerate() {
            for &t_idx in match_map.values() {
                all_template_atoms.insert((ri, t_idx));
            }
        }

        // Parity-aware chirality post-check.  Runs only when the SMIRKS has @/@@.
        // This must happen after the complete VF2 mapping is known, because
        // correct chirality comparison requires the full neighbor permutation.
        if has_stereo {
            let ok = (0..rxn.reactants.len()).all(|ri| {
                smirks_chirality_ok(&rxn.reactants[ri], reactants[ri], &combo[ri])
            });
            if !ok {
                continue;
            }
        }

        let products: Vec<Molecule> = rxn
            .products
            .iter()
            .map(|pt| {
                build_product(pt, &global_map, reactants, &all_template_atoms, carry_substituents)
            })
            .collect();

        // Skip product sets that contain any over-valenced atom.
        if products.iter().all(|p| validate_valence(p).is_empty()) {
            results.push(products);
        }
    }

    Ok(results)
}

// ---------------------------------------------------------------------------
// SMIRKS chirality post-check (parity-aware)
// ---------------------------------------------------------------------------

/// Returns the parity of the permutation P that maps `from_seq` to `to_seq`:
/// `Some(true)` = even (same chirality sense), `Some(false)` = odd (inverted).
/// Returns `None` when the sequences differ in length or contain elements that
/// cannot be aligned (e.g. an unmapped neighbour).
fn permutation_parity(from_seq: &[u32], to_seq: &[u32]) -> Option<bool> {
    let n = from_seq.len();
    if n != to_seq.len() {
        return None;
    }
    // Build perm[j] = index i in from_seq where from_seq[i] == to_seq[j].
    let mut perm = Vec::with_capacity(n);
    for &t in to_seq {
        let pos = from_seq.iter().position(|&f| f == t)?;
        perm.push(pos);
    }
    // Count inversions to determine parity.
    let mut inv = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            if perm[i] > perm[j] {
                inv += 1;
            }
        }
    }
    Some(inv.is_multiple_of(2)) // true = even = chirality flags must agree for same config
}

/// Parity-aware stereo check for one (template, reactant, mapping) triple.
///
/// For each chiral atom in `tmpl`, maps the template's recorded SMILES stereo
/// neighbor order through the VF2 `match_map` into the reactant atom-index space,
/// then computes the parity of the permutation relative to the reactant's recorded
/// SMILES stereo neighbor order.
///
/// - Even parity → chirality flags must agree for same absolute configuration.
/// - Odd parity  → chirality flags must differ for same absolute configuration.
///
/// Returns `true` if all chiral centres are consistent, `false` if any mismatch.
fn smirks_chirality_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &HashMap<usize, AtomIdx>,
) -> bool {
    for i in 0..tmpl.atom_count() {
        let tmpl_atom = tmpl.atom(AtomIdx(i as u32));
        if tmpl_atom.chirality == Chirality::None {
            continue;
        }

        // Template atom's SMILES stereo neighbour order (template atom indices).
        let Some(tmpl_order) = tmpl.stereo_neighbor_order(AtomIdx(i as u32)) else {
            continue; // No recorded order — skip this centre.
        };

        // Corresponding matched reactant atom.
        let Some(&react_idx) = match_map.get(&i) else {
            continue; // Template atom not in mapping (shouldn't happen for complete match).
        };

        let react_atom = reactant.atom(react_idx);
        if react_atom.chirality == Chirality::None {
            return false; // Template requires stereo; reactant atom has none.
        }

        // Map each template stereo-neighbour index to the corresponding reactant atom index.
        let mut mapped: Vec<u32> = Vec::with_capacity(tmpl_order.len());
        let mut all_mapped = true;
        for &t in tmpl_order {
            if t == STEREO_H_SENTINEL {
                mapped.push(STEREO_H_SENTINEL);
            } else {
                match match_map.get(&(t as usize)) {
                    Some(ri) => mapped.push(ri.0),
                    None => {
                        all_mapped = false;
                        break;
                    }
                }
            }
        }
        if !all_mapped {
            continue; // Partial substructure match — cannot verify chirality.
        }

        // Reactant atom's SMILES stereo neighbour order (reactant atom indices).
        let Some(react_order) = reactant.stereo_neighbor_order(react_idx) else {
            // No recorded order in reactant — fall back to raw flag comparison.
            if react_atom.chirality != tmpl_atom.chirality {
                return false;
            }
            continue;
        };

        let Some(even_parity) = permutation_parity(&mapped, react_order) else {
            continue; // Alignment failed — skip.
        };

        // Check: even parity → flags must agree; odd parity → flags must differ.
        let same_flag = tmpl_atom.chirality == react_atom.chirality;
        if same_flag != even_parity {
            return false;
        }
    }
    true
}

/// Convert a SMIRKS reactant-template `Molecule` to a `QueryMolecule` for VF2.
///
/// Constraints included:
/// - `AtomicNum` and `Aromatic` (always)
/// - `Charge` when non-zero
/// - `HCount` when a bracket atom specifies H > 0 (e.g. `[NH2:1]`)
///   Zero-H bracket atoms (`[N:1]`) are treated as "any H count" because
///   the parser returns 0 for both unspecified and explicit-zero H.
///
/// Chirality (`@`/`@@`) is NOT encoded into the query here; it is checked after
/// VF2 matching via `smirks_chirality_ok`, which uses a permutation-parity
/// comparison so that the same absolute configuration is recognised regardless
/// of how the reactant molecule was written as SMILES.
fn mol_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let mut q = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(
                atom.element.atomic_number(),
            ))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(atom.aromatic))),
        );

        if atom.charge != 0 {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge))),
            );
        }

        if let Some(h) = atom.hydrogen_count
            && h > 0
        {
            q = AtomQuery::And(
                Box::new(q),
                Box::new(AtomQuery::Primitive(AtomPrimitive::HCount(h))),
            );
        }

        qmol.add_atom_with_map(q, atom.atom_map);
    }

    for (_bidx, bond) in mol.bonds() {
        let bq = match bond.order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => {
                BondQuery::Primitive(BondPrimitive::Single)
            }
            BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
            BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
            BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
            BondOrder::QuerySingleOrDouble => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
            ),
            BondOrder::QuerySingleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Single)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::QueryDoubleOrAromatic => BondQuery::Or(
                Box::new(BondQuery::Primitive(BondPrimitive::Double)),
                Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
            ),
            BondOrder::Quadruple | BondOrder::Zero | BondOrder::QueryAny => {
                BondQuery::Primitive(BondPrimitive::Any)
            }
        };
        qmol.add_bond(bond.atom1.0 as usize, bond.atom2.0 as usize, bq);
    }

    qmol
}

/// Clear Up/Down stereo markers from bonds that have no adjacent double bond.
///
/// After a SMIRKS reaction, C=C → C=O conversions can leave stale Up/Down
/// markers (E/Z direction indicators) on single bonds that are no longer
/// adjacent to any double bond.  Such orphaned markers produce invalid SMILES
/// (`/C=O` is nonsensical) and must be demoted to plain Single bonds (RDKit #9339).
fn clear_orphaned_stereo_bonds(mol: Molecule) -> Molecule {
    let orphaned: Vec<BondIdx> = mol
        .bonds()
        .filter_map(|(bidx, bond)| {
            if bond.order != BondOrder::Up && bond.order != BondOrder::Down {
                return None;
            }
            // Up/Down is valid only when at least one endpoint has an adjacent
            // double bond (the one that the Up/Down bond specifies direction for).
            let has_double = [bond.atom1, bond.atom2].iter().any(|&a| {
                mol.neighbors(a)
                    .any(|(_, nb_bidx)| mol.bond(nb_bidx).order == BondOrder::Double)
            });
            if has_double { None } else { Some(bidx) }
        })
        .collect();

    if orphaned.is_empty() {
        return mol;
    }

    let mut builder = chematic_core::MoleculeBuilder::new();
    for (_, atom) in mol.atoms() {
        builder.add_atom(atom.clone());
    }
    for (bidx, bond) in mol.bonds() {
        let order = if orphaned.contains(&bidx) {
            BondOrder::Single
        } else {
            bond.order
        };
        let _ = builder.add_bond(bond.atom1, bond.atom2, order);
    }
    builder.copy_stereo_from(&mol);
    builder.build()
}

/// Build one product molecule applying full SMIRKS semantics.
///
/// 1. Atom-mapped product atoms: copy source atom + override aromatic/charge/H from template.
/// 2. New product atoms (no map): clone from template.
/// 3. BFS from core (mapped) atoms through input molecules, collecting substituents
///    (non-template atoms reachable without crossing template-atom walls).
/// 4. Add product-template bonds (new/changed bonds).
/// 5. Carry through bonds from source molecules where at least one endpoint is a substituent.
fn build_product(
    product_template: &Molecule,
    global_map: &HashMap<u16, (usize, AtomIdx)>,
    input_mols: &[&Molecule],
    all_template_atoms: &HashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    // template_idx_to_new[i]: new AtomIdx for product template atom i.
    let mut template_idx_to_new: Vec<Option<AtomIdx>> = vec![None; product_template.atom_count()];
    // src_to_new: (mol_idx, src_AtomIdx) → new AtomIdx in the product.
    let mut src_to_new: HashMap<(usize, AtomIdx), AtomIdx> = HashMap::new();

    // --- Step 1: add product template atoms ---
    // core_keys: only source atoms that are mapped by THIS product template.
    // Using global_map.values() (all matched atoms across all templates) would
    // seed the BFS in Step 2 from atoms belonging to *other* product templates,
    // causing their substituents to leak into this product (issue #13).
    let product_maps: HashSet<u16> = (0..product_template.atom_count())
        .filter_map(|i| product_template.atom(AtomIdx(i as u32)).atom_map)
        .collect();
    let core_keys: HashSet<(usize, AtomIdx)> = global_map
        .iter()
        .filter(|(am, _)| product_maps.contains(am))
        .map(|(_, &src)| src)
        .collect();

    for (i, slot) in template_idx_to_new.iter_mut().enumerate() {
        let tmpl_atom = product_template.atom(AtomIdx(i as u32));
        let new_idx = if let Some(am) = tmpl_atom.atom_map {
            if let Some(&(mol_idx, src_idx)) = global_map.get(&am) {
                // Core atom: copy source, then override electronic state from template.
                let src_atom = input_mols[mol_idx].atom(src_idx);
                let mut new_atom = src_atom.clone();
                new_atom.aromatic = tmpl_atom.aromatic;
                new_atom.charge = tmpl_atom.charge;
                // Copy H count only when template specifies > 0 (e.g. [NH2:1]).
                // A bare bracket atom ([O:1], [N:1]) has hydrogen_count=Some(0) which
                // means "unspecified" in a product context — clear it so implicit
                // valence rules determine H count (fixes issue #18).
                new_atom.hydrogen_count = tmpl_atom.hydrogen_count.filter(|&h| h > 0);
                // Apply product-template chirality when explicitly specified (@/@@).
                // When the template has Chirality::None:
                //   • If the non-mapped substituents written into the product template
                //     have a DIFFERENT element composition from what the reactant
                //     template matched (substituent replacement, e.g. Br→I), the
                //     neighbour topology changes and the source stereo is invalid.
                //   • If the substituent element sets are identical (pass-through,
                //     e.g. [C:1](F)(Cl)Br >> [C:1](F)(Cl)Br)) preserve the source
                //     chirality — the common case for remote-reaction SMIRKS.
                if tmpl_atom.chirality != Chirality::None {
                    new_atom.chirality = tmpl_atom.chirality;
                } else {
                    // Element multiset of non-mapped atoms in the product template
                    // adjacent to this core atom.
                    let mut prod_elems: HashMap<u8, usize> = HashMap::new();
                    for (nb, _) in product_template.neighbors(AtomIdx(i as u32)) {
                        if product_template.atom(nb).atom_map.is_none() {
                            *prod_elems
                                .entry(product_template.atom(nb).element.atomic_number())
                                .or_insert(0) += 1;
                        }
                    }
                    if !prod_elems.is_empty() {
                        // Element multiset of this core atom's neighbours that were
                        // matched by the reactant template (heavy atoms only).
                        let mut rxn_elems: HashMap<u8, usize> = HashMap::new();
                        for (nb, _) in input_mols[mol_idx].neighbors(src_idx) {
                            if all_template_atoms.contains(&(mol_idx, nb)) {
                                *rxn_elems
                                    .entry(
                                        input_mols[mol_idx].atom(nb).element.atomic_number(),
                                    )
                                    .or_insert(0) += 1;
                            }
                        }
                        if prod_elems != rxn_elems {
                            new_atom.chirality = Chirality::None;
                        }
                    }
                }
                new_atom.atom_map = None;
                let idx = builder.add_atom(new_atom);
                src_to_new.insert((mol_idx, src_idx), idx);
                idx
            } else {
                // Map number not in reactants — new atom from template.
                let mut new_atom = tmpl_atom.clone();
                new_atom.atom_map = None;
                builder.add_atom(new_atom)
            }
        } else {
            // No atom_map — entirely new atom from template.
            let mut new_atom = tmpl_atom.clone();
            new_atom.atom_map = None;
            builder.add_atom(new_atom)
        };
        *slot = Some(new_idx);
    }

    // --- Step 2: BFS from core atoms to collect substituents ---
    // Skipped when carry_substituents = false (run_reactants_strict mode).
    // Seed visited with all template atoms so BFS cannot cross into the template region.
    let mut visited: HashSet<(usize, AtomIdx)> = all_template_atoms.clone();
    if carry_substituents {
        let mut queue: VecDeque<(usize, AtomIdx)> = core_keys.iter().cloned().collect();

        while let Some((mol_idx, cur_idx)) = queue.pop_front() {
            for (nb_idx, _bond_idx) in input_mols[mol_idx].neighbors(cur_idx) {
                let key = (mol_idx, nb_idx);
                if visited.contains(&key) {
                    continue;
                }
                visited.insert(key);
                let src_atom = input_mols[mol_idx].atom(nb_idx);
                let mut new_atom = src_atom.clone();
                new_atom.atom_map = None;
                let new_idx = builder.add_atom(new_atom);
                src_to_new.insert(key, new_idx);
                queue.push_back(key);
            }
        }
    }

    // --- Step 3: add product template bonds ---
    let mut added_bond_pairs: HashSet<(AtomIdx, AtomIdx)> = HashSet::new();

    for (_bidx, bond) in product_template.bonds() {
        let a_new = template_idx_to_new[bond.atom1.0 as usize].unwrap();
        let b_new = template_idx_to_new[bond.atom2.0 as usize].unwrap();
        let _ = builder.add_bond(a_new, b_new, bond.order);
        added_bond_pairs.insert((a_new.min(b_new), a_new.max(b_new)));
    }

    // --- Step 4: carry-through bonds from source molecules ---
    // Bonds where both endpoints are template atoms are replaced or broken by the template;
    // bonds where at least one endpoint is a substituent are carried through.
    for (&(mol_idx, src_idx), &a_new) in &src_to_new {
        for (nb_idx, bond_idx) in input_mols[mol_idx].neighbors(src_idx) {
            let nb_key = (mol_idx, nb_idx);
            let Some(&b_new) = src_to_new.get(&nb_key) else {
                continue;
            };
            if all_template_atoms.contains(&(mol_idx, src_idx))
                && all_template_atoms.contains(&nb_key)
            {
                continue;
            }
            let pair = (a_new.min(b_new), a_new.max(b_new));
            if added_bond_pairs.contains(&pair) {
                continue;
            }
            added_bond_pairs.insert(pair);
            let bond_order = input_mols[mol_idx].bond(bond_idx).order;
            let _ = builder.add_bond(a_new, b_new, bond_order);
        }
    }

    // Clear any Up/Down stereo markers left on bonds that are no longer adjacent
    // to a double bond (e.g. after C=C → C=O conversion via SMIRKS).
    clear_orphaned_stereo_bonds(builder.build())
}

/// Standard Cartesian product: given `sets[0], sets[1], …`, return all
/// ordered selections of one element from each set.
fn cartesian_product<T: Clone>(sets: &[Vec<T>]) -> Vec<Vec<T>> {
    let mut result: Vec<Vec<T>> = vec![vec![]];
    for set in sets {
        result = result
            .into_iter()
            .flat_map(|combo| {
                set.iter().map(move |item| {
                    let mut new_combo = combo.clone();
                    new_combo.push(item.clone());
                    new_combo
                })
            })
            .collect();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn identity_single_atom() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[C:1]>>[C:1]", &[&mol]).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 1);
        assert_eq!(results[0][0].atom_count(), 1);
    }

    #[test]
    fn no_match_returns_empty() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert!(
            results.is_empty(),
            "nitrogen template must not match methane"
        );
    }

    #[test]
    fn multiple_matches_in_single_mol() {
        let mol = parse("NCCN").unwrap();
        let results = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert_eq!(results.len(), 2, "two N atoms in NCCN → two product sets");
    }

    #[test]
    fn bond_formation_two_mols() {
        let n_mol = parse("N").unwrap();
        let c_mol = parse("C").unwrap();
        let results = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&n_mol, &c_mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2, "product must have 2 atoms");
        assert_eq!(prod.bonds().count(), 1, "product must have 1 bond");
    }

    #[test]
    fn bond_cleavage_two_products() {
        let mol = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let products = &results[0];
        assert_eq!(products.len(), 2, "two product templates → two products");
        assert_eq!(products[0].atom_count(), 1);
        assert_eq!(products[1].atom_count(), 1);
    }

    #[test]
    fn reactant_count_mismatch_error() {
        let mol = parse("C").unwrap();
        let err = run_reactants("[N:1].[C:2]>>[N:1][C:2]", &[&mol]);
        assert!(
            matches!(
                err,
                Err(TransformError::ReactantCountMismatch {
                    expected: 2,
                    got: 1
                })
            ),
            "two-template SMIRKS with one reactant must error"
        );
    }

    #[test]
    fn invalid_smirks_error() {
        let mol = parse("C").unwrap();
        let err = run_reactants("[X]>>[X]", &[&mol]);
        assert!(
            matches!(err, Err(TransformError::SmirksParse(_))),
            "unknown element must yield SmirksParse error"
        );
    }

    #[test]
    fn overvalent_product_filtered_oxygen() {
        // O normally has max valence 2.
        // SMIRKS adds two carbons to an oxygen that already has one bond → 3 bonds on O → invalid.
        // CCO: the O is bonded to 1 C (bond_sum=1). Template [O:1]>>[O:1](C)C adds 2 more.
        let ethanol = parse("CCO").unwrap();
        let results = run_reactants("[O:1]>>[O:1](C)C", &[&ethanol]).unwrap();
        // The O that already had 1 bond would get 3 → over-valenced → filtered out.
        // The only match is the terminal O (1 bond → +2 = 3 bonds, invalid).
        assert!(
            results.is_empty(),
            "product with O having 3 bonds must be filtered out, got {} sets",
            results.len()
        );
    }

    #[test]
    fn valid_charged_product_kept() {
        // N with charge +1 can have up to 4 bonds (normal valences [3,5], +1 allows 4).
        // trimethylamine N(C)(C)C has N with bond_sum=3, charge=0.
        // Template [N:1]>>[N+:1] just changes charge, keeps 3 bonds → valid.
        let tma = parse("N(C)(C)C").unwrap();
        let results = run_reactants("[N:1]>>[N+:1]", &[&tma]).unwrap();
        assert!(
            !results.is_empty(),
            "N+ with 3 bonds must be valid and kept"
        );
    }

    #[test]
    fn new_atom_in_product() {
        let mol = parse("C").unwrap();
        let results = run_reactants("[C:1]>>[C:1]=O", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2, "C + new O = 2 atoms");
    }

    #[test]
    fn amide_bond_formation() {
        // NH3 + H-C(=O)-Cl → H-C(=O)-NH2 (formamide)
        let nh3 = parse("N").unwrap();
        let hcocl = parse("C(=O)Cl").unwrap();
        let results = run_reactants("[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]", &[&nh3, &hcocl]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 3, "C + O(new) + N = 3 atoms");
    }

    #[test]
    fn double_bond_product() {
        let mol = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        let prod = &results[0][0];
        assert_eq!(prod.atom_count(), 2);
        let bond_orders: Vec<BondOrder> = prod.bonds().map(|(_, b)| b.order).collect();
        assert!(
            bond_orders.contains(&BondOrder::Double),
            "product must contain a double bond"
        );
    }

    #[test]
    fn substituent_carry_through() {
        // Methylamine + acetyl chloride → N-methylacetamide (5 heavy atoms)
        // CH3-NH2 + CH3-C(=O)-Cl → CH3-C(=O)-NH-CH3
        let methylamine = parse("NC").unwrap();
        let acetyl_cl = parse("CC(=O)Cl").unwrap();
        let results = run_reactants(
            "[N:1].[C:2](=O)Cl>>[C:2](=O)[N:1]",
            &[&methylamine, &acetyl_cl],
        )
        .unwrap();
        assert!(!results.is_empty(), "must produce at least one product set");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom_count(),
            5,
            "N-methylacetamide has 5 heavy atoms, got {}",
            prod.atom_count()
        );
    }

    #[test]
    fn bfs_no_leakage_into_other_product_template_atoms() {
        // Issue #13: in diethylamine (CCNCC), the SMIRKS [N:1][C:2]>>[N:1].[C:2]
        // should cleave the N-C bond and produce:
        //   product1 [N:1] = N + right ethyl chain  (3 atoms: N, C, C)
        //   product2 [C:2] = left ethyl fragment     (2 atoms: C, C)
        //
        // Before the #13 fix, the BFS for product2 was seeded from BOTH N and C:2,
        // causing the right ethyl chain (atoms beyond N) to leak into product2
        // → product2 would have 4 atoms instead of 2.
        let diethylamine = parse("CCNCC").unwrap(); // C-C-N-C-C, 5 heavy atoms
        let results = run_reactants("[N:1][C:2]>>[N:1].[C:2]", &[&diethylamine]).unwrap();
        assert!(
            !results.is_empty(),
            "should find at least one N-C bond match"
        );

        // Find a result where product2 ([C:2]) has exactly 2 atoms (ethyl fragment)
        // — this is only possible when BFS does NOT leak the other ethyl chain.
        let clean_cleavage = results.iter().find(|ps| {
            ps.len() == 2
                && ((ps[0].atom_count() == 3 && ps[1].atom_count() == 2)
                    || (ps[0].atom_count() == 2 && ps[1].atom_count() == 3))
        });
        assert!(
            clean_cleavage.is_some(),
            "expected at least one product set with sizes {{3, 2}} (N+ethyl, ethyl); \
             all sets: {:?}",
            results
                .iter()
                .map(|ps| ps.iter().map(|p| p.atom_count()).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn single_product_no_leakage_from_other_template_core() {
        // Ethane cleavage: each product should be a single carbon atom.
        let ethane = parse("CC").unwrap();
        let results = run_reactants("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
        assert!(!results.is_empty());
        for ps in &results {
            assert_eq!(ps.len(), 2, "two product templates → two products");
            assert_eq!(ps[0].atom_count(), 1, "each product is a single carbon");
            assert_eq!(ps[1].atom_count(), 1, "each product is a single carbon");
        }
    }

    // ── Stereo SMIRKS tests ───────────────────────────────────────────────────

    #[test]
    fn stereo_preserved_when_template_has_no_spec() {
        // Product template [C:1] has no chirality → source @@ is preserved via clone.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)Br", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        // The core C atom is first in the builder (index 0).
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        // Template has None → source Clockwise (@@ in SMILES = Clockwise) is preserved.
        assert_eq!(
            core_chirality,
            Chirality::Clockwise,
            "source @@ chirality must be preserved when template has no stereo spec"
        );
    }

    #[test]
    fn stereo_inverted_by_template() {
        // Product template [C@H:1] has @ (CounterClockwise) → overrides source @@ (Clockwise).
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::CounterClockwise,
            "product template @ must override source @@ → CounterClockwise"
        );
    }

    // ── run_reactants_strict tests ────────────────────────────────────────────

    #[test]
    fn strict_mode_excludes_substituents() {
        // Methylamine (NC): in normal mode [N:1]>>[N:1] carries C through as substituent.
        // In strict mode only N is returned (no C).
        let mol = parse("NC").unwrap();
        let normal = run_reactants("[N:1]>>[N:1]", &[&mol]).unwrap();
        let strict = run_reactants_strict("[N:1]>>[N:1]", &[&mol]).unwrap();
        assert!(!normal.is_empty());
        assert!(!strict.is_empty());
        let normal_atoms = normal[0][0].atom_count();
        let strict_atoms = strict[0][0].atom_count();
        assert!(
            normal_atoms > strict_atoms,
            "normal mode carries substituent C (got {normal_atoms}), \
             strict mode only mapped N (got {strict_atoms})"
        );
        assert_eq!(strict_atoms, 1, "strict mode: only the mapped N atom");
    }

    #[test]
    fn strict_mode_bond_cleavage() {
        // Ethane cleavage: strict mode gives 1-atom products, same as normal here
        // (no unmapped substituents on either C).
        let ethane = parse("CC").unwrap();
        let results = run_reactants_strict("[C:1][C:2]>>[C:1].[C:2]", &[&ethane]).unwrap();
        assert!(!results.is_empty());
        for ps in &results {
            assert_eq!(ps[0].atom_count(), 1);
            assert_eq!(ps[1].atom_count(), 1);
        }
    }

    // ── Issue #18: product bracket notation cleanup ───────────────────────────

    #[test]
    fn product_removes_bracket_from_bare_bracket_atoms() {
        // Issue #18: [O:1] in product template (hydrogen_count=Some(0)) must produce
        // clean `O` SMILES, not `[O]`.
        use chematic_smiles::canonical_smiles;
        let mol = parse("OCC").unwrap();
        let results = run_reactants("[OH:1]>>[O:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match hydroxyl");
        let prod_smi = canonical_smiles(&results[0][0]);
        assert!(
            !prod_smi.contains("[O]"),
            "bare [O:1] product must write as O, not [O], got: {prod_smi}"
        );
    }

    #[test]
    fn product_preserves_explicit_h_from_template() {
        // [NH2:1] in product template specifies 2H explicitly — must be kept.
        use chematic_smiles::canonical_smiles;
        let mol = parse("NC").unwrap();
        let results = run_reactants("[N:1]>>[NH2:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match amine N");
        let smi = canonical_smiles(&results[0][0]);
        // With explicit H2, the N must be in brackets.
        assert!(
            smi.contains("[NH2]"),
            "explicit [NH2:1] in product must keep [NH2], got: {smi}"
        );
    }

    // ── Issue #20: SMIRKS stereo filtering ───────────────────────────────────

    #[test]
    fn stereo_filter_rejects_wrong_enantiomer() {
        // Issue #20: SMIRKS with @@ reactant template must match @@ but NOT @ reactant.
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine (@@)
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (@)

        let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";
        let results_l = run_reactants(smirks, &[&l_ala]).unwrap();
        let results_d = run_reactants(smirks, &[&d_ala]).unwrap();

        assert!(!results_l.is_empty(), "L-alanine (@@) must match @@ template");
        assert!(
            results_d.is_empty(),
            "D-alanine (@) must NOT match @@ template (stereo filter, issue #20)"
        );
    }

    #[test]
    fn stereo_neutral_smirks_matches_both_enantiomers() {
        // SMIRKS without @/@@ must still match both enantiomers (backward compat).
        let l_ala = parse("N[C@@H](C)C(=O)O").unwrap();
        let d_ala = parse("N[C@H](C)C(=O)O").unwrap();
        let smirks = "[N:1][CH:2](C)C(=O)O>>[N:1][CH:2](C)C(=O)O";
        let r_l = run_reactants(smirks, &[&l_ala]).unwrap();
        let r_d = run_reactants(smirks, &[&d_ala]).unwrap();
        assert!(!r_l.is_empty(), "L-alanine must match non-stereo template");
        assert!(!r_d.is_empty(), "D-alanine must match non-stereo template");
    }

    #[test]
    fn stereo_filter_same_config_different_write_order() {
        // Parity-aware matching must accept the same absolute configuration
        // regardless of SMILES atom write order (the confirmed bug in raw flag
        // comparison). Both molecules ARE L-alanine.
        //   Form A: N[C@@H](C)C(=O)O  — N written first, stored as Clockwise
        //   Form B: C[C@H](N)C(=O)O   — C_methyl first, stored as CounterClockwise
        // A raw flag comparison would reject Form B against an @@ template.
        let l_form_a = parse("N[C@@H](C)C(=O)O").unwrap();
        let l_form_b = parse("C[C@H](N)C(=O)O").unwrap(); // same absolute config, diff write order
        let d_form   = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (opposite config)

        let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";

        let r_a = run_reactants(smirks, &[&l_form_a]).unwrap();
        let r_b = run_reactants(smirks, &[&l_form_b]).unwrap();
        let r_d = run_reactants(smirks, &[&d_form]).unwrap();

        assert!(!r_a.is_empty(), "L-alanine form A (N-first @@) must match");
        assert!(!r_b.is_empty(),
            "L-alanine form B (C-first @, same absolute config) must also match \
             — parity-aware comparison required");
        assert!(r_d.is_empty(), "D-alanine must still be rejected");
    }

    // ── RDKit #9339: orphaned stereo bonds cleared from products ─────────────

    #[test]
    fn smirks_reaction_clears_orphaned_stereo_bonds() {
        // (E)-2-butene C/C=C/C has Up/Down bonds adjacent to the C=C double bond.
        // After SMIRKS [C:1]=[C:2]>>[C:1][C:2] the double bond is reduced to a
        // single bond.  The Up/Down single bonds on C0-C1 and C2-C3 are no longer
        // adjacent to ANY double bond → orphaned → must be cleared (RDKit PR #9339).
        let mol = parse("C/C=C/C").unwrap(); // (E)-2-butene
        let results = run_reactants("[C:1]=[C:2]>>[C:1][C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should produce at least one product");
        for prod_set in &results {
            for prod in prod_set {
                for (_, bond) in prod.bonds() {
                    assert_ne!(
                        bond.order, BondOrder::Up,
                        "stray Up bond in product after C=C→C-C (RDKit #9339)"
                    );
                    assert_ne!(
                        bond.order, BondOrder::Down,
                        "stray Down bond in product after C=C→C-C (RDKit #9339)"
                    );
                }
            }
        }
    }

    #[test]
    fn smirks_preserves_stereo_bonds_adjacent_to_remaining_double() {
        // If the double bond is kept unchanged, Up/Down must be preserved.
        // [C:1]=[C:2]>>[C:1]=[C:2] is the identity transformation for alkenes.
        let mol = parse("C/C=C/C").unwrap(); // (E)-2-butene
        let results = run_reactants("[C:1]=[C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty());
        // At least one Up or Down bond must remain (E/Z stereo preserved).
        let has_stereo_bond = results[0]
            .iter()
            .flat_map(|p| p.bonds())
            .any(|(_, b)| b.order == BondOrder::Up || b.order == BondOrder::Down);
        assert!(has_stereo_bond, "stereo bonds adjacent to retained double bond must be preserved");
    }
}
