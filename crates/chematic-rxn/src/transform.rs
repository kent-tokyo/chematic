use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, validate_valence};
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

/// Convert a SMIRKS reactant-template `Molecule` to a `QueryMolecule` for VF2.
///
/// Constraints included:
/// - `AtomicNum` and `Aromatic` (always)
/// - `Charge` when non-zero
/// - `HCount` when a bracket atom specifies H > 0 (e.g. `[NH2:1]`)
///   Zero-H bracket atoms (`[N:1]`) are treated as "any H count" because
///   the parser returns 0 for both unspecified and explicit-zero H.
fn mol_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));

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

        qmol.add_atom(q);
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
                if tmpl_atom.hydrogen_count.is_some() {
                    new_atom.hydrogen_count = tmpl_atom.hydrogen_count;
                }
                // Apply product-template chirality when explicitly specified (@/@@).
                // When the template has Chirality::None, the source chirality is
                // preserved (inherited from the clone above) — this is the common
                // case for reactions that don't change the stereocentre.
                if tmpl_atom.chirality != Chirality::None {
                    new_atom.chirality = tmpl_atom.chirality;
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

    builder.build()
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
}
