use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

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
    crate::perf_counters::record_run_reactants_call();
    let prepared = prepare_reaction(smirks)?;
    let matches = find_matches_impl(&prepared, reactants)?;
    Ok(matches
        .iter()
        .filter_map(|m| apply_match_impl(&prepared, reactants, m, carry_substituents))
        .collect())
}

/// One accepted match of `smirks`'s reactant-side pattern(s) against a set of
/// input molecules — one map per reactant-template slot, keyed by the query
/// atom index (position within that reactant template) to the matched
/// [`AtomIdx`] in the corresponding input molecule. This is exactly the
/// per-combination mapping [`run_reactants_impl`] already builds internally,
/// given a name and a public seam between "enumerate matches" and
/// "apply one match" (issue #225).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionMatch {
    pub per_reactant: Vec<FxHashMap<usize, AtomIdx>>,
}

impl ReactionMatch {
    /// `atom_map` number → (reactant slot index, matched `AtomIdx`),
    /// resolved against `smirks`'s own atom-map annotations — so a caller
    /// can look up "where did mapped atom N end up in my input molecules"
    /// without re-deriving the reactant templates' atom maps itself.
    pub fn atom_map_positions(
        &self,
        smirks: &str,
    ) -> Result<FxHashMap<u16, (usize, AtomIdx)>, TransformError> {
        let rxn = parse_reaction(smirks)?;
        let n_templates = rxn.reactants.len();
        if self.per_reactant.len() != n_templates {
            return Err(TransformError::ReactantCountMismatch {
                expected: n_templates,
                got: self.per_reactant.len(),
            });
        }
        Ok(global_map_of(
            &self.per_reactant,
            &template_atom_maps_of(&rxn),
        ))
    }
}

/// Enumerate every match of `smirks`'s reactant-side pattern(s) against
/// `reactants`, without applying the transformation or building any
/// product. One entry per combination that passes the existing
/// chirality/E-Z stereo post-checks — i.e. exactly the matches
/// [`run_reactants`] would go on to build products for (issue #225).
pub fn find_reaction_matches(
    smirks: &str,
    reactants: &[&Molecule],
) -> Result<Vec<ReactionMatch>, TransformError> {
    let prepared = prepare_reaction(smirks)?;
    find_matches_impl(&prepared, reactants)
}

/// Apply the reaction for exactly one match (as returned by
/// [`find_reaction_matches`]), producing the product set for that match
/// alone. `Ok(None)` means this match's product set failed the existing
/// valence filter — the same case [`run_reactants`] silently drops today
/// (issue #225).
///
/// Does not re-run the chirality/E-Z stereo post-checks that
/// [`find_reaction_matches`] already applied — `m` is expected to be one
/// of the matches it returned (or otherwise already known to satisfy
/// them). Re-parses `smirks`, matching [`run_reactants`]'s own existing
/// per-call behavior.
pub fn apply_reaction_match(
    smirks: &str,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Result<Option<Vec<Molecule>>, TransformError> {
    let prepared = prepare_reaction(smirks)?;
    let n_templates = prepared.rxn.reactants.len();
    if reactants.len() != n_templates {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: reactants.len(),
        });
    }
    if m.per_reactant.len() != n_templates {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: m.per_reactant.len(),
        });
    }
    Ok(apply_match_impl(
        &prepared,
        reactants,
        m,
        carry_substituents,
    ))
}

/// Parsed SMIRKS plus everything derived from it that matching needs —
/// shared by [`run_reactants_impl`], [`find_reaction_matches`], and
/// [`apply_reaction_match`] so the three can never compute it
/// inconsistently.
struct PreparedReaction {
    rxn: crate::reaction::Reaction,
    queries: Vec<QueryMolecule>,
    template_atom_maps: Vec<Vec<Option<u16>>>,
    has_stereo: bool,
    has_ez_stereo: bool,
}

fn prepare_reaction(smirks: &str) -> Result<PreparedReaction, TransformError> {
    crate::perf_counters::record_reaction_parse_call();
    let rxn = parse_reaction(smirks)?;

    // Build a QueryMolecule from each reactant template, and record the
    // atom-map number for each query atom index.
    let queries: Vec<QueryMolecule> = rxn.reactants.iter().map(mol_to_query).collect();
    let template_atom_maps = template_atom_maps_of(&rxn);

    // Detect whether any reactant template carries @/@@ stereo, so we can apply
    // the parity-aware post-check after VF2 completes.  Chirality is NOT encoded
    // into the VF2 query because the raw flag comparison in eval_chirality is
    // SMILES-write-order-dependent; the correct check requires the full mapping
    // (see smirks_chirality_ok below).
    let has_stereo = rxn
        .reactants
        .iter()
        .any(|r| r.atoms().any(|(_, a)| a.chirality != Chirality::None));
    // Similarly, E/Z double-bond stereo (/ and \) is NOT encoded into the VF2
    // query; it is checked post-VF2 via smirks_ez_stereo_ok.
    let has_ez_stereo = rxn.reactants.iter().any(|r| {
        r.bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
    });

    Ok(PreparedReaction {
        rxn,
        queries,
        template_atom_maps,
        has_stereo,
        has_ez_stereo,
    })
}

fn template_atom_maps_of(rxn: &crate::reaction::Reaction) -> Vec<Vec<Option<u16>>> {
    rxn.reactants
        .iter()
        .map(|tmpl| {
            (0..tmpl.atom_count())
                .map(|i| tmpl.atom(AtomIdx(i as u32)).atom_map)
                .collect()
        })
        .collect()
}

/// `atom_map` number → (reactant slot index, matched `AtomIdx`), built from
/// one match's per-reactant maps and the reactant templates' own atom-map
/// annotations. Shared by [`ReactionMatch::atom_map_positions`] and
/// [`apply_match_impl`].
fn global_map_of(
    per_reactant: &[FxHashMap<usize, AtomIdx>],
    template_atom_maps: &[Vec<Option<u16>>],
) -> FxHashMap<u16, (usize, AtomIdx)> {
    let mut global_map: FxHashMap<u16, (usize, AtomIdx)> = FxHashMap::default();
    for (ri, match_map) in per_reactant.iter().enumerate() {
        for (&qi, &t_idx) in match_map {
            if let Some(am) = template_atom_maps[ri][qi] {
                global_map.insert(am, (ri, t_idx));
            }
        }
    }
    global_map
}

/// Steps 1–2 of the original `run_reactants_impl`: VF2-match every reactant
/// template against its input molecule, take the cartesian product across
/// template slots, and keep only combinations that survive the
/// chirality/E-Z stereo post-checks. Does not build any product.
fn find_matches_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
) -> Result<Vec<ReactionMatch>, TransformError> {
    let n_templates = prepared.rxn.reactants.len();
    if reactants.len() != n_templates {
        return Err(TransformError::ReactantCountMismatch {
            expected: n_templates,
            got: reactants.len(),
        });
    }

    // VF2 match: for each (template_query, input_mol) pair.
    let all_match_sets: Vec<Vec<FxHashMap<usize, AtomIdx>>> = prepared
        .queries
        .iter()
        .zip(reactants.iter())
        .map(|(q, mol)| {
            let matches = find_matches(q, mol);
            crate::perf_counters::record_reactant_query_match_call(matches.len());
            matches
        })
        .collect();

    // No matches when any template has no match.
    if all_match_sets.iter().any(|ms| ms.is_empty()) {
        return Ok(vec![]);
    }

    let mut matches: Vec<ReactionMatch> = Vec::new();

    for combo in cartesian_product(&all_match_sets) {
        crate::perf_counters::record_match_combination();

        // Parity-aware chirality post-check.  Runs only when the SMIRKS has @/@@.
        // This must happen after the complete VF2 mapping is known, because
        // correct chirality comparison requires the full neighbor permutation.
        if prepared.has_stereo {
            let ok = (0..prepared.rxn.reactants.len()).all(|ri| {
                smirks_chirality_ok(&prepared.rxn.reactants[ri], reactants[ri], &combo[ri])
            });
            if !ok {
                continue;
            }
        }
        // E/Z double-bond stereo post-check.  Runs only when the SMIRKS has /\.
        if prepared.has_ez_stereo {
            let ok = (0..prepared.rxn.reactants.len()).all(|ri| {
                smirks_ez_stereo_ok(&prepared.rxn.reactants[ri], reactants[ri], &combo[ri])
            });
            if !ok {
                continue;
            }
        }

        matches.push(ReactionMatch {
            per_reactant: combo,
        });
    }

    Ok(matches)
}

/// Step 3 of the original `run_reactants_impl`: build the product set for
/// one already-accepted match and apply the valence filter. `None` means
/// the product set contained an over-valenced atom.
fn apply_match_impl(
    prepared: &PreparedReaction,
    reactants: &[&Molecule],
    m: &ReactionMatch,
    carry_substituents: bool,
) -> Option<Vec<Molecule>> {
    // global_map: atom_map_number → (reactant_mol_idx, matched_AtomIdx)
    let global_map = global_map_of(&m.per_reactant, &prepared.template_atom_maps);

    // all_template_atoms: every (mol_idx, AtomIdx) matched by any reactant template atom.
    // Used as BFS walls to prevent substituent collection from crossing into the
    // template region, and to identify bonds that the product template replaces.
    let mut all_template_atoms: FxHashSet<(usize, AtomIdx)> = FxHashSet::default();
    for (ri, match_map) in m.per_reactant.iter().enumerate() {
        for &t_idx in match_map.values() {
            all_template_atoms.insert((ri, t_idx));
        }
    }

    let products: Vec<Molecule> = prepared
        .rxn
        .products
        .iter()
        .map(|pt| {
            let product = build_product(
                pt,
                &global_map,
                reactants,
                &all_template_atoms,
                carry_substituents,
            );
            crate::perf_counters::record_build_product_call(
                product.atom_count(),
                product.bond_count(),
            );
            product
        })
        .collect();

    // Skip product sets that contain any over-valenced atom.
    if products.iter().all(|p| validate_valence(p).is_empty()) {
        crate::perf_counters::record_product_set();
        Some(products)
    } else {
        None
    }
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
    match_map: &FxHashMap<usize, AtomIdx>,
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

// ---------------------------------------------------------------------------
// SMIRKS E/Z double-bond stereo post-check
// ---------------------------------------------------------------------------

/// Returns the "outward" stereo-bond direction from `atom` (one endpoint of a
/// double bond) toward its Up/Down-annotated substituent, or `None` if no
/// such bond exists.
///
/// E/Z stereo is encoded on the *substituent* bonds adjacent to a C=C, not on
/// the double bond itself.  Whether the stored bond goes *into* or *out of*
/// `atom` determines how to read its direction:
/// - Outgoing (`atom` is `bond.atom1`): direction is as stored.
/// - Incoming (`atom` is `bond.atom2`): direction is flipped (Up ↔ Down).
///
/// The `other` parameter is the opposite endpoint of the double bond; that
/// bond is skipped so we look only at substituents.
fn ez_stereo_outward(mol: &Molecule, atom: AtomIdx, other: AtomIdx) -> Option<BondOrder> {
    for (nb, bidx) in mol.neighbors(atom) {
        if nb == other {
            continue; // skip the double bond itself
        }
        let bond = mol.bond(bidx);
        match bond.order {
            BondOrder::Up | BondOrder::Down => {
                let outward = if bond.atom1 == atom {
                    // bond goes FROM atom outward → direction as stored
                    bond.order
                } else {
                    // bond comes INTO atom → flip to get outward direction
                    match bond.order {
                        BondOrder::Up => BondOrder::Down,
                        _ => BondOrder::Up,
                    }
                };
                return Some(outward);
            }
            _ => {}
        }
    }
    None
}

/// E/Z stereo post-check for one (template, reactant, mapping) triple.
///
/// For each double bond in `tmpl` that has Up/Down substituent bonds on **both**
/// sides, verify that the corresponding double bond in `reactant` (found via
/// the VF2 `match_map`) encodes the same E/Z parity.
///
/// Parity is determined by comparing the "outward direction" from each end of
/// the double bond (see [`ez_stereo_outward`]).  Two outward directions that are
/// equal → same-side (Z/cis); directions that differ → opposite-side (E/trans).
///
/// If only one side of a template double bond has a stereo bond (or the
/// reactant doesn't annotate stereo at all), the constraint is skipped — this
/// matches the behaviour of SMIRKS templates extracted from rdchiral where a
/// single-sided annotation is common.
///
/// Returns `true` if all constrained double bonds are consistent.
fn smirks_ez_stereo_ok(
    tmpl: &Molecule,
    reactant: &Molecule,
    match_map: &FxHashMap<usize, AtomIdx>,
) -> bool {
    for (_, bond) in tmpl.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let ta = bond.atom1;
        let tb = bond.atom2;

        // Outward directions from each end of the template double bond.
        let sa = ez_stereo_outward(tmpl, ta, tb);
        let sb = ez_stereo_outward(tmpl, tb, ta);

        // Both sides must be specified to establish an E/Z constraint.
        let (sa, sb) = match (sa, sb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue, // no constraint on this double bond
        };

        // Map template atoms to reactant atoms.
        let Some(&ra) = match_map.get(&(ta.0 as usize)) else {
            continue;
        };
        let Some(&rb) = match_map.get(&(tb.0 as usize)) else {
            continue;
        };

        // Outward directions from the corresponding reactant double bond ends.
        let ma = ez_stereo_outward(reactant, ra, rb);
        let mb = ez_stereo_outward(reactant, rb, ra);

        // If the reactant doesn't encode stereo on either end, skip (don't reject).
        let (ma, mb) = match (ma, mb) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        // Compare E/Z parity: same outward direction = Z, different = E.
        // If template and reactant disagree, reject this mapping.
        if (sa == sb) != (ma == mb) {
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
    // copy_stereo_from copies stereo_neighbor_order but NOT stereo_groups.
    // Preserve both by applying each separately.
    builder.copy_stereo_from(&mol);
    let mut result = builder.build();
    // Restore enhanced stereo groups (ABS/OR/AND) that copy_stereo_from omits.
    result.set_stereo_groups(mol.stereo_groups().to_vec());
    result
}

// ---------------------------------------------------------------------------
// Product-side chirality correction (parity-aware)
// ---------------------------------------------------------------------------
//
// `build_product`'s Step 1 leaves every mapped core atom's chirality as
// whatever `src_atom.clone()` produced -- the reactant's own flag, if any.
// Two things can make that flag wrong or meaningless by the time the
// product molecule's bonds are fully built (Steps 3-4):
//
//   - An EXPLICIT product-template flag (`[C@:1]`/`[C@@:1]`) describes an
//     absolute configuration relative to the TEMPLATE's own neighbour
//     write-order, not the reactant's. Naively copying the symbol without
//     accounting for a reordered template ignores exactly the kind of
//     reorder that inverts configuration (mirroring the parity math
//     `smirks_chirality_ok`/`permutation_parity` already does for
//     REACTANT-side match validation, just never applied to the product).
//   - An INHERITED flag (no explicit template chirality) can survive on an
//     atom whose real bonding topology changed within the reaction's own
//     mapped/core region (e.g. a ring bond broken by the reaction) --
//     `build_product`'s old invalidation heuristic only ever compared
//     *unmapped* substituent element sets, never mapped-neighbour identity.
//
// Both are handled here, in a single post-build pass (run after all of
// `build_product`'s bonds exist, since validating either case requires the
// atom's REAL final neighbour set): an explicit template flag is
// transcribed together with its own neighbour order (validated against the
// atom's real adjacency, not re-derived); an inherited flag is kept only if
// both an order exists to trust it by and the core neighbourhood provably
// didn't change, and cleared otherwise. Two mechanisms, not one unified
// permutation-parity re-expression, because an unmapped/template-literal
// substituent has no reactant identity to map an order back to at all.

/// Map a template-space stereo order into product-index space via
/// `template_idx_to_new`. [`STEREO_H_SENTINEL`] passes through unchanged.
/// `None` if any real (non-sentinel) template atom didn't make it into the
/// product (shouldn't happen for a well-formed template, but fails closed).
fn map_stereo_order(order: &[u32], template_idx_to_new: &[Option<AtomIdx>]) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                template_idx_to_new
                    .get(t as usize)
                    .copied()
                    .flatten()
                    .map(|a| a.0)
            }
        })
        .collect()
}

/// True when `order` (already in product-index space) matches `new_idx`'s
/// REAL final neighbour set in `product`: at most one H-sentinel, no
/// duplicate real entries, right length, and the real entries are exactly
/// `new_idx`'s actual adjacency. Deliberately a set/degree check, not a
/// permutation-parity one -- `order` is stored verbatim once validated;
/// `corrected_chirality` (the SMILES writer) and `chematic-cip` already do
/// the write-order-independent parity math whenever they consume it.
fn order_matches_final_topology(product: &Molecule, new_idx: AtomIdx, order: &[u32]) -> bool {
    let sentinels = order.iter().filter(|&&x| x == STEREO_H_SENTINEL).count();
    if sentinels > 1 {
        return false;
    }
    let real: Vec<u32> = order
        .iter()
        .copied()
        .filter(|&x| x != STEREO_H_SENTINEL)
        .collect();
    let real_set: FxHashSet<u32> = real.iter().copied().collect();
    if real.len() != real_set.len() {
        return false; // Duplicate entry -- malformed order.
    }
    if order.len() != product.degree(new_idx) + sentinels {
        return false;
    }
    let actual: FxHashSet<u32> = product.neighbors(new_idx).map(|(nb, _)| nb.0).collect();
    real_set == actual
}

/// Bug-B validity gate for an atom that inherited its chirality flag from
/// the reactant (no explicit product-template chirality): the original
/// unmapped-substituent element-multiset comparison, PLUS a symmetric
/// mapped-neighbour atom-map-number-set comparison -- closing the gap where
/// a core/mapped neighbour's topology changes (e.g. a ring bond broken by
/// the reaction) invisibly to the element-multiset-only check.
fn bug_b_topology_unchanged(
    product_template: &Molecule,
    ti: AtomIdx,
    reactant: &Molecule,
    src_idx: AtomIdx,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    mol_idx: usize,
    atom_to_map: &FxHashMap<(usize, AtomIdx), u16>,
) -> bool {
    let mut prod_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut prod_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in product_template.neighbors(ti) {
        match product_template.atom(nb).atom_map {
            None => {
                *prod_elems
                    .entry(product_template.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
            Some(am) => {
                prod_mapped.insert(am);
            }
        }
    }

    let mut rxn_elems: FxHashMap<u8, usize> = FxHashMap::default();
    let mut rxn_mapped: FxHashSet<u16> = FxHashSet::default();
    for (nb, _) in reactant.neighbors(src_idx) {
        if !all_template_atoms.contains(&(mol_idx, nb)) {
            continue;
        }
        match atom_to_map.get(&(mol_idx, nb)) {
            Some(&am) => {
                rxn_mapped.insert(am);
            }
            None => {
                *rxn_elems
                    .entry(reactant.atom(nb).element.atomic_number())
                    .or_insert(0) += 1;
            }
        }
    }

    // NOTE: `prod_elems.is_empty()` is deliberately NOT special-cased here.
    // An empty `prod_elems` legitimately means "product-template atom has no
    // unmapped neighbours" -- which is also true when the reactant's
    // template-matched unmapped substituents (F/Cl/Br etc., matched by the
    // reactant SMARTS itself, not carried-through remote substituents --
    // those never enter `rxn_elems` at all, since `rxn_elems` only counts
    // neighbours inside `all_template_atoms`) were silently DELETED by the
    // product template. `prod_elems == rxn_elems` alone correctly requires
    // both to be empty together, or both to hold the same multiset --
    // exactly what "unchanged" means; special-casing empty-on-one-side
    // would let a genuine deletion through undetected.
    prod_elems == rxn_elems && prod_mapped == rxn_mapped
}

/// Map a REACTANT-space stereo order (atom indices into
/// `input_mols[mol_idx]`) into product-index space via `src_to_new`.
/// [`STEREO_H_SENTINEL`] passes through unchanged. `None` if any real
/// (non-sentinel) reactant neighbour never made it into the product (e.g. it
/// was part of a leaving group the reaction consumed) -- the correspondence
/// is not derivable, so the caller must fail closed rather than guess.
fn remap_reactant_stereo_order(
    order: &[u32],
    mol_idx: usize,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Option<Vec<u32>> {
    order
        .iter()
        .map(|&t| {
            if t == STEREO_H_SENTINEL {
                Some(STEREO_H_SENTINEL)
            } else {
                src_to_new.get(&(mol_idx, AtomIdx(t))).map(|a| a.0)
            }
        })
        .collect()
}

/// Post-build chirality correction. Must run after `build_product`'s Steps
/// 3-4 (all bonds added), since validation needs each atom's real final
/// degree/neighbour set. See the module-level doc above for why this is two
/// mechanisms (transcribe-and-validate for an explicit template flag,
/// gate-and-preserve-or-clear for an inherited one), not one.
fn correct_product_stereo(
    mut product: Molecule,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> Molecule {
    let atom_to_map: FxHashMap<(usize, AtomIdx), u16> =
        global_map.iter().map(|(&am, &k)| (k, am)).collect();

    for i in 0..product_template.atom_count() {
        let ti = AtomIdx(i as u32);
        let Some(new_idx) = template_idx_to_new[i] else {
            continue;
        };
        let tmpl_atom = product_template.atom(ti);

        if tmpl_atom.chirality != Chirality::None {
            // Explicit product-template chirality: the template is the sole
            // source of truth here regardless of whether this atom is a
            // matched core atom, an unmatched map number, or fully new --
            // no global_map lookup needed on this branch.
            let resolved = product_template
                .stereo_neighbor_order(ti)
                .and_then(|order| map_stereo_order(order, template_idx_to_new))
                .filter(|order| order_matches_final_topology(&product, new_idx, order));
            match resolved {
                Some(order) => {
                    product.set_chirality(new_idx, tmpl_atom.chirality);
                    product.set_stereo_neighbor_order(new_idx, order);
                }
                None => product.set_chirality(new_idx, Chirality::None),
            }
        } else if let Some(am) = tmpl_atom.atom_map
            && let Some(&(mol_idx, src_idx)) = global_map.get(&am)
        {
            // Matched core atom, no explicit template chirality: keep the
            // flag `src_atom.clone()` gave it in Step 1 only if it's
            // actually usable and provably still valid. Three-way rule,
            // all conditions required:
            //   1. bug_b_topology_unchanged: the TEMPLATE-level neighbour
            //      composition (unmapped element multiset + mapped atom-map
            //      set) is unchanged reactant-template -> product-template.
            //   2. remap_reactant_stereo_order succeeds: every atom the
            //      reactant's own recorded order references is uniquely
            //      identifiable in the product (via src_to_new) -- fails
            //      closed (None) for a leaving-group neighbour the reaction
            //      consumed, where no correspondence is derivable at all.
            //   3. order_matches_final_topology: the remapped order's real
            //      entries are EXACTLY this atom's real final adjacency in
            //      the built product molecule -- catches a neighbour that
            //      survived into the product (so (2) succeeds) but whose
            //      specific BOND to this atom did not (e.g. reattached
            //      elsewhere by the reaction), which a template-level
            //      heuristic alone cannot see.
            // Never keep a raw flag without a validated order alongside it.
            let src_atom = input_mols[mol_idx].atom(src_idx);
            if src_atom.chirality != Chirality::None {
                let topology_unchanged = bug_b_topology_unchanged(
                    product_template,
                    ti,
                    input_mols[mol_idx],
                    src_idx,
                    all_template_atoms,
                    mol_idx,
                    &atom_to_map,
                );
                let remapped_order = topology_unchanged
                    .then(|| input_mols[mol_idx].stereo_neighbor_order(src_idx))
                    .flatten()
                    .and_then(|order| remap_reactant_stereo_order(order, mol_idx, src_to_new))
                    .filter(|order| order_matches_final_topology(&product, new_idx, order));
                match remapped_order {
                    Some(order) => {
                        product.set_chirality(new_idx, src_atom.chirality);
                        product.set_stereo_neighbor_order(new_idx, order);
                    }
                    None => product.set_chirality(new_idx, Chirality::None),
                }
            }
        }
        // else: no template chirality, and not a matched core atom --
        // new_atom.chirality is already Chirality::None from Step 1's clone
        // of a template-literal atom.
    }

    product
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
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    // template_idx_to_new[i]: new AtomIdx for product template atom i.
    let mut template_idx_to_new: Vec<Option<AtomIdx>> = vec![None; product_template.atom_count()];
    // src_to_new: (mol_idx, src_AtomIdx) → new AtomIdx in the product.
    let mut src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx> = FxHashMap::default();

    // --- Step 1: add product template atoms ---
    // core_keys: only source atoms that are mapped by THIS product template.
    // Using global_map.values() (all matched atoms across all templates) would
    // seed the BFS in Step 2 from atoms belonging to *other* product templates,
    // causing their substituents to leak into this product (issue #13).
    let product_maps: FxHashSet<u16> = (0..product_template.atom_count())
        .filter_map(|i| product_template.atom(AtomIdx(i as u32)).atom_map)
        .collect();
    let core_keys: FxHashSet<(usize, AtomIdx)> = global_map
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
                // Chirality is intentionally left as whatever src_atom.clone()
                // produced above (i.e. the reactant's own flag, if any) --
                // correct_product_stereo, run after all bonds are added
                // below, is the sole authority on the final chirality/order
                // for every atom, whether it inherits from the reactant or
                // carries an explicit product-template @/@@.
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
    let mut visited: FxHashSet<(usize, AtomIdx)> = all_template_atoms.clone();
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
    let mut added_bond_pairs: FxHashSet<(AtomIdx, AtomIdx)> = FxHashSet::default();

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
            let ob = input_mols[mol_idx].bond(bond_idx);
            // Preserve the original atom1→atom2 orientation. Up/Down (E/Z)
            // bond semantics are direction-dependent, so adding the bond with
            // endpoints swapped relative to the source would flip the geometry.
            let (a, b) = if ob.atom1 == src_idx {
                (a_new, b_new)
            } else {
                (b_new, a_new)
            };
            let _ = builder.add_bond(a, b, ob.order);
        }
    }

    // Parity-aware atom chirality correction (see the doc above
    // correct_product_stereo) -- must run after all bonds above are added,
    // since it validates against each atom's real final neighbour set.
    let product = correct_product_stereo(
        builder.build(),
        product_template,
        &template_idx_to_new,
        global_map,
        all_template_atoms,
        input_mols,
        &src_to_new,
    );

    // Clear any Up/Down stereo markers left on bonds that are no longer adjacent
    // to a double bond (e.g. after C=C → C=O conversion via SMIRKS).
    clear_orphaned_stereo_bonds(product)
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

    fn canonical(mol: &Molecule) -> String {
        chematic_smiles::canonical_smiles(mol)
    }

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
    fn stereo_cleared_when_all_neighbors_are_unmapped_template_literals() {
        // Product template [C:1] has no chirality, and all three substituents
        // (F, Cl, Br) are unmapped literal atoms in BOTH templates -- SMIRKS
        // gives no per-atom reactant<->product correspondence for an unmapped
        // atom (only `:n` atom-map numbers establish identity), so there is no
        // derivable stereo_neighbor_order to carry the inherited flag with.
        // Fail-closed: clear rather than keep a flag with an undefined order
        // (chematic-cip already treats order-less chirality as unassigned, so
        // keeping the raw flag here was never actually meaningful downstream).
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)Br", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        // The core C atom is first in the builder (index 0).
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::None,
            "no derivable stereo_neighbor_order for an all-unmapped-substituent \
             inherited flag must clear chirality, not keep an order-less flag"
        );
    }

    #[test]
    fn stereo_preserved_when_all_neighbors_are_mapped_and_remappable() {
        // Every substituent around the stereocenter is atom-mapped, so each
        // has a real, unique reactant->product correspondence via src_to_new.
        // This is the case the remap mechanism *can* and must handle: keep
        // both the flag and a validated, remapped stereo_neighbor_order.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants(
            "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
            &[&mol],
        )
        .unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::Clockwise,
            "source @@ chirality must be preserved when every neighbor is \
             mapped and remappable, with an identity-order template"
        );
        assert!(
            prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
            "a kept inherited flag must always carry a validated stereo_neighbor_order"
        );
    }

    #[test]
    fn stereo_unmapped_leaving_groups_fully_removed_clears_chirality() {
        // Blocker #1 from review: the reactant's literal unmapped
        // substituents (F, Cl, Br) are removed entirely by the product
        // template (bare [C:1], zero neighbours) -- bug_b_topology_unchanged
        // must not special-case an empty product-side unmapped set as
        // "unchanged" when the reactant side was non-empty. In the current
        // architecture this is enforced twice over: the topology gate itself
        // (fixed here) AND independently by remap_reactant_stereo_order,
        // since a literal/unmapped reactant-template neighbour never has a
        // src_to_new entry regardless of what the product does with it. The
        // observable behavior this test pins is the one the reviewer asked
        // for either way: chirality must clear, not survive.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom(AtomIdx(0)).chirality,
            Chirality::None,
            "deleting all of the stereocenter's unmapped substituents must clear chirality"
        );
    }

    #[test]
    fn stereo_remote_reaction_preserves_stereocenter() {
        // The bond change (Br leaving, N arriving) happens two bonds away
        // from the stereocenter, through a mapped, unchanged carrier chain
        // (F, Cl, C4 all keep the same atom-map numbers and the same
        // relative template order on both sides) -- the stereocenter's own
        // immediate neighbor set and order are completely untouched by the
        // remote reaction. Configuration must survive.
        let mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
        let amine = parse("N").unwrap();
        let results = run_reactants(
            "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
             >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
            &[&mol, &amine],
        )
        .unwrap();
        assert!(!results.is_empty(), "should match and produce a product");
        let prod = &results[0][0];
        assert_eq!(
            prod.atom(AtomIdx(0)).chirality,
            Chirality::Clockwise,
            "a remote bond change 2 bonds away must not disturb the \
             stereocenter's own unchanged, fully-mapped neighbor set"
        );
        assert!(
            prod.stereo_neighbor_order(AtomIdx(0)).is_some(),
            "preserved chirality must carry a validated stereo_neighbor_order"
        );
    }

    #[test]
    fn stereo_survives_16_plus_atom_map_relabelings_smiles_and_cip_invariant() {
        // The chemistry must not depend on which integers a SMIRKS author
        // picks for `:n` map numbers -- only the correspondence they encode.
        // Re-run the same reaction with 20 distinct map-number assignments
        // (same template text shape and order each time, only the four
        // integers change) and assert canonical SMILES and accurate CIP
        // agree with the first run every time.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        // Central atom is fixed at map `:1`; offset the other three ranges
        // so none of them ever collides with it or each other.
        let labelings: Vec<[u32; 3]> = (2..=21).map(|i| [i, i + 100, i + 200]).collect();
        assert!(labelings.len() >= 16, "need at least 16 relabelings");

        let mut baseline_smiles: Option<String> = None;
        let mut baseline_cip: Option<Option<chematic_core::CipCode>> = None;
        for [b, c, d] in labelings {
            let smirks =
                format!("[C@@H:1]([F:{b}])([Cl:{c}])[Br:{d}]>>[C:1]([F:{b}])([Cl:{c}])[Br:{d}]");
            let results = run_reactants(&smirks, &[&mol]).unwrap();
            assert!(!results.is_empty(), "labeling {b},{c},{d} must still match");
            let prod = &results[0][0];
            let smi = canonical(prod);
            let cip = chematic_chem::assign_cip_with_mode(prod, chematic_chem::CipMode::Accurate)
                .unwrap()
                .get(AtomIdx(0));

            match (&baseline_smiles, &baseline_cip) {
                (None, None) => {
                    baseline_smiles = Some(smi);
                    baseline_cip = Some(cip);
                }
                (Some(base_smi), Some(base_cip)) => {
                    assert_eq!(
                        &smi, base_smi,
                        "canonical SMILES must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                    );
                    assert_eq!(
                        &cip, base_cip,
                        "accurate CIP code must be invariant to atom-map relabeling \
                         (labeling {b},{c},{d})"
                    );
                }
                _ => unreachable!(),
            }
        }
    }

    #[test]
    fn stereo_chirality_implies_valid_stereo_neighbor_order_invariant() {
        // Cross-cutting invariant the reviewer required: whenever an atom
        // carries a non-None chirality flag, it must always have a valid
        // stereo_neighbor_order alongside it -- never a raw flag with no
        // order (which chematic-cip and the SMILES writer cannot interpret
        // meaningfully). Checked across every atom of a representative
        // spread of already-covered product-generating scenarios.
        let assert_invariant_holds = |mol: &Molecule, label: &str| {
            for i in 0..mol.atom_count() {
                let idx = AtomIdx(i as u32);
                if mol.atom(idx).chirality != Chirality::None {
                    assert!(
                        mol.stereo_neighbor_order(idx).is_some(),
                        "{label}: atom {i} has chirality {:?} but no stereo_neighbor_order",
                        mol.atom(idx).chirality
                    );
                }
            }
        };

        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let chain_mol = parse("[C@@H](F)(Cl)CCBr").unwrap();
        let amine = parse("N").unwrap();

        let scenarios: Vec<(&str, Vec<Molecule>)> = vec![
            (
                "identity template order",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "reordered template inverts",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "all neighbors mapped and remappable",
                vec![
                    run_reactants(
                        "[C@@H:1]([F:2])([Cl:3])[Br:4]>>[C:1]([F:2])([Cl:3])[Br:4]",
                        &[&mol],
                    )
                    .unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "all unmapped substituents deleted",
                vec![run_reactants("[C@@H:1](F)(Cl)Br>>[C:1]", &[&mol]).unwrap()[0][0].clone()],
            ),
            (
                "substituent replacement",
                vec![
                    run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap()[0][0]
                        .clone(),
                ],
            ),
            (
                "remote reaction",
                vec![
                    run_reactants(
                        "[C@@H:1]([F:2])([Cl:3])[C:4][C:5][Br:6].[N:7]\
                         >>[C:1]([F:2])([Cl:3])[C:4][C:5][N:7]",
                        &[&chain_mol, &amine],
                    )
                    .unwrap()[0][0]
                        .clone(),
                ],
            ),
        ];

        for (label, products) in &scenarios {
            for prod in products {
                assert_invariant_holds(prod, label);
            }
        }
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

    #[test]
    fn stereo_identity_template_order_preserves_configuration() {
        // Product template repeats the exact same neighbour order as the
        // reactant template (F, Cl, Br) with an explicit @@ -- the simplest
        // "nothing changed" case for the parity-aware correction.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap();
        let prod = &results[0][0];
        let smi = canonical(prod);
        let input_smi = canonical(&mol);
        assert_eq!(
            smi, input_smi,
            "identity template order should reproduce the input unchanged"
        );
    }

    #[test]
    fn stereo_reordered_template_inverts_absolute_configuration() {
        // Issue found while surveying RDKit's open issues (analogous to
        // RDKit #9257): reordering two substituents in the product template
        // while keeping the SAME @@ symbol must invert the absolute
        // configuration -- the symbol describes an order-relative sense,
        // not an absolute one. Cross-checked against a live RDKit oracle
        // (`rdkit.Chem.rdChemReactions`): RDKit's `reordered` output equals
        // its own `inverted` output, both differing from `identity`.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let identity = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0],
        );
        let reordered = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@@H:1](Cl)(F)Br", &[&mol]).unwrap()[0][0],
        );
        let explicit_invert = canonical(
            &run_reactants("[C@@H:1](F)(Cl)Br>>[C@H:1](F)(Cl)Br", &[&mol]).unwrap()[0][0],
        );
        assert_ne!(
            reordered, identity,
            "reordering two substituents under the same @@ symbol must change \
             the absolute configuration"
        );
        assert_eq!(
            reordered, explicit_invert,
            "reordering two substituents under @@ must match the explicit-@ \
             (same order, opposite symbol) result -- both express the same \
             inverted configuration"
        );
    }

    #[test]
    fn stereo_substituent_replacement_clears_chirality() {
        // Neighbour SET changes (Br -> I): the old flag can no longer mean
        // anything -- must clear to Chirality::None, not carry a stale @@
        // onto a differently-substituted center. Pre-existing behavior,
        // now driven by correct_product_stereo's Bug-B gate instead of the
        // deleted ad hoc heuristic; must not regress.
        let mol = parse("[C@@H](F)(Cl)Br").unwrap();
        let results = run_reactants("[C@@H:1](F)(Cl)Br>>[C:1](F)(Cl)I", &[&mol]).unwrap();
        let prod = &results[0][0];
        let core_chirality = prod.atom(AtomIdx(0)).chirality;
        assert_eq!(
            core_chirality,
            Chirality::None,
            "substituent replacement (Br->I) must clear chirality, not preserve a stale flag"
        );
    }

    #[test]
    fn stereo_lost_mapped_neighbor_clears_spurious_chirality() {
        // Issue found while surveying RDKit's open issues: a mapped/core
        // neighbour (N) is dropped by the product template, while the
        // template's own C:1 atom carries no explicit chirality. The old
        // invalidation heuristic only ever compared *unmapped* substituent
        // element sets (both empty here, since every neighbour of C:1 is
        // mapped) and missed this case entirely, producing a spurious
        // chirality tag on a carbon with 3 identical implicit hydrogens.
        let mol = parse("[C@H](N)(O)Cl").unwrap();
        let results = run_reactants("[C@H:1]([N:2])([O:3])Cl>>[C:1][O:3]", &[&mol]).unwrap();
        let prod = &results[0][0];
        let smi = canonical(prod);
        assert!(
            !smi.contains('@'),
            "carbon losing its N neighbour is no longer a stereocenter (3 identical \
             implicit H's) -- product must carry no chirality symbol, got {smi}"
        );
    }

    #[test]
    fn stereo_polyol_ring_contraction_no_spurious_ch2oh_chirality() {
        // The originally reported repro (analogous to RDKit #9257), minimized
        // versions of which are the two tests directly above. A pyranose ->
        // furanose ring contraction where one ring carbon (mapped :2) loses
        // its own explicit chirality in the product template and becomes a
        // plain CH2OH (degree 2, bonded to a degree-1 O) -- pre-fix this
        // atom kept a spurious inherited chirality tag ([C@H2], a carbon
        // with 3 identical implicit hydrogens). Does NOT assert a specific
        // sign at the two atoms whose absolute configuration this fix does
        // not by itself resolve either way (mapped :11, unchanged reactant->
        // product template and therefore architecturally untouched by this
        // fix; and :15, the debated center in the original RDKit report) --
        // only that the reaction produces exactly 4 real stereocenters, not
        // 5, and that the correct atom (the one that structurally lost its
        // stereocenter) is the one that lost its tag.
        let smirks = "[O:1][C@H:2]1[O:3][C@H:4]([C:5][C:6])[C@@H:11]([O:12])[C@H:13]([O:14])\
                       [C@H:15]1[O:16]>>[O:1][C:2][C@:15]([O:16])1[O:3][C@H:4]([C:5][C:6])\
                       [C@@H:11]([O:12])[C@H:13]1[O:14]";
        let mol = parse("CC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O").unwrap();
        let results = run_reactants(smirks, &[&mol]).unwrap();
        let prod = &results[0][0];

        let n = prod.atom_count();
        let chiral_atoms: Vec<AtomIdx> = (0..n)
            .map(|i| AtomIdx(i as u32))
            .filter(|&a| prod.atom(a).chirality != Chirality::None)
            .collect();
        assert_eq!(
            chiral_atoms.len(),
            4,
            "expected exactly 4 real stereocenters in the product, got {}: {:?}",
            chiral_atoms.len(),
            chiral_atoms
        );

        // Structurally identify the CH2OH carbon: degree 2, bonded to a
        // degree-1 oxygen -- not by a hardcoded index, since atom order
        // depends on build_product's own construction order.
        let ch2oh = (0..n).map(|i| AtomIdx(i as u32)).find(|&a| {
            let atom = prod.atom(a);
            atom.element == chematic_core::Element::C
                && prod.degree(a) == 2
                && prod.neighbors(a).any(|(nb, _)| {
                    prod.atom(nb).element == chematic_core::Element::O && prod.degree(nb) == 1
                })
        });
        let ch2oh = ch2oh.expect("product must contain a CH2OH carbon");
        assert_eq!(
            prod.atom(ch2oh).chirality,
            Chirality::None,
            "the CH2OH carbon (lost its own chirality in the product template, degree 2, \
             no longer has 4 distinct substituents) must not carry a stereo tag"
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
        // [NH2:1] in product template specifies 2H explicitly — the built
        // atom's `hydrogen_count` must be `Some(2)`, not silently dropped or
        // recomputed.
        //
        // Checks the data-level field directly rather than the canonical
        // SMILES string's bracket notation (issue #205): after
        // `initial_invariant`/`emit_atom`'s explicit/implicit-H-count
        // unification fix, `canonical_smiles` correctly stops forcing
        // brackets for an atom whose explicit H count merely repeats what
        // organic-subset valence inference would already give -- this atom
        // (N bonded to one carbon) infers to 2 implicit H anyway, so its
        // canonical form is now the fully standard, bracket-free "CN"
        // (methylamine), not "C[NH2]". The template's explicit
        // specification is still honored -- it is what the "2" came from --
        // just no longer forced into visible bracket notation once it's
        // redundant with inference.
        let mol = parse("NC").unwrap();
        let results = run_reactants("[N:1]>>[NH2:1]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "should match amine N");
        let product = &results[0][0];
        let n_atom = product
            .atoms()
            .find(|(_, a)| a.element == chematic_core::Element::N)
            .map(|(_, a)| a)
            .expect("product must contain the templated N atom");
        assert_eq!(
            n_atom.hydrogen_count,
            Some(2),
            "explicit [NH2:1] in product must set hydrogen_count = Some(2)"
        );
    }

    #[test]
    fn reaction_derived_matches_direct_parse_chlorobenzene() {
        // The exact fixture from kent-tokyo/renkin PR #65: a reaction-
        // derived molecule (built via `run_reactants`, whose product-
        // template atom comes from bracket notation `[Cl]` in the SMIRKS)
        // must canonicalize identically to a directly-parsed organic-subset
        // "Clc1ccccc1" of the same compound (issue #205).
        use chematic_smiles::canonical_smiles;
        let fwd = "[c:1][Br]>>[c:1][Cl]";
        let known = parse("Brc1ccccc1").unwrap();
        let results = run_reactants(fwd, &[&known]).unwrap();
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].len(), 1);
        let reaction_derived = &results[0][0];
        let direct = parse("Clc1ccccc1").unwrap();

        // Structural sanity first: same atom count and element multiset,
        // so canonical-string equality below is a real invariance proof,
        // not an assumption that the two are the same molecule.
        assert_eq!(reaction_derived.atom_count(), direct.atom_count());
        let mut els_a: Vec<_> = reaction_derived.atoms().map(|(_, a)| a.element).collect();
        let mut els_b: Vec<_> = direct.atoms().map(|(_, a)| a.element).collect();
        els_a.sort();
        els_b.sort();
        assert_eq!(els_a, els_b, "same element multiset required");

        assert_eq!(
            canonical_smiles(reaction_derived),
            canonical_smiles(&direct),
            "reaction-derived and directly-parsed chlorobenzene must canonicalize identically"
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

        assert!(
            !results_l.is_empty(),
            "L-alanine (@@) must match @@ template"
        );
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
        let d_form = parse("N[C@H](C)C(=O)O").unwrap(); // D-alanine (opposite config)

        let smirks = "[N:1][C@@H:2](C)C(=O)O>>[N:1][C@@H:2](C)C(=O)O";

        let r_a = run_reactants(smirks, &[&l_form_a]).unwrap();
        let r_b = run_reactants(smirks, &[&l_form_b]).unwrap();
        let r_d = run_reactants(smirks, &[&d_form]).unwrap();

        assert!(!r_a.is_empty(), "L-alanine form A (N-first @@) must match");
        assert!(
            !r_b.is_empty(),
            "L-alanine form B (C-first @, same absolute config) must also match \
             — parity-aware comparison required"
        );
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
                        bond.order,
                        BondOrder::Up,
                        "stray Up bond in product after C=C→C-C (RDKit #9339)"
                    );
                    assert_ne!(
                        bond.order,
                        BondOrder::Down,
                        "stray Down bond in product after C=C→C-C (RDKit #9339)"
                    );
                }
            }
        }
    }

    #[test]
    fn smirks_preserves_stereo_bonds_adjacent_to_remaining_double() {
        // If the double bond is kept unchanged, the *exact* E/Z geometry must be
        // preserved — not merely "some Up/Down bond survives" (which passed even
        // while the geometry was flipping E↔Z, issue #50). Verify by comparing the
        // product's canonical SMILES to the canonical SMILES of the known input.
        use chematic_smiles::canonical_smiles;
        for input in ["C/C=C/C", "C/C=C\\C"] {
            let mol = parse(input).unwrap();
            let results = run_reactants("[C:1]=[C:2]>>[C:1]=[C:2]", &[&mol]).unwrap();
            assert!(!results.is_empty());
            let expected = canonical_smiles(&mol);
            let got = canonical_smiles(&results[0][0]);
            assert_eq!(
                got, expected,
                "identity SMIRKS must preserve exact E/Z geometry for {input}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // E/Z double-bond stereo filtering (issue #21)
    // -----------------------------------------------------------------------

    #[test]
    fn ez_stereo_e_template_matches_e_alkene() {
        // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
        // E-2-butene (C/C=C/C) should produce 1 result.
        let e_alkene = parse("C/C=C/C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
        let results = run_reactants(smirks, &[&e_alkene]).unwrap();
        assert!(!results.is_empty(), "E-template must match E-alkene");
    }

    #[test]
    fn ez_stereo_e_template_rejects_z_alkene() {
        // Template specifies E: [C:1]/[C:2]=[C:3]/[C:4]
        // Z-2-butene (C/C=C\C) should produce 0 results.
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]/[C:4]>>[C:1][C:2][C:3][C:4]";
        let results = run_reactants(smirks, &[&z_alkene]).unwrap();
        assert!(results.is_empty(), "E-template must reject Z-alkene");
    }

    #[test]
    fn ez_stereo_neutral_template_matches_both_geometries() {
        // Template without stereo: [C:1][C:2]=[C:3][C:4]>>[C:1]
        // Both E and Z alkenes should match.
        let e_alkene = parse("C/C=C/C").unwrap();
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1][C:2]=[C:3][C:4]>>[C:1]";
        assert!(
            !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "neutral template must match E-alkene"
        );
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "neutral template must match Z-alkene"
        );
    }

    #[test]
    fn ez_stereo_one_sided_template_matches_both_geometries() {
        // Single-sided stereo bond in template: [C:1]/[C:2]=[C:3][C:4]
        // Without both sides specified, E/Z is ambiguous → no filtering.
        let e_alkene = parse("C/C=C/C").unwrap();
        let z_alkene = parse("C/C=C\\C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3][C:4]>>[C:1]";
        assert!(
            !run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "one-sided template must match E-alkene"
        );
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "one-sided template must match Z-alkene"
        );
    }

    #[test]
    fn ez_stereo_retro_wittig_z_matches_z_hexene() {
        // Retro-Wittig (Z-alkene → two carbonyls).
        // SMIRKS: [C:1]/[C:2]=[C:3]\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]
        //   reads: C:2 and C:3 on OPPOSITE sides (E/trans for those two)
        //   but the substituents C:1 and C:4 are on the SAME side (Z-selectivity)
        //
        // Z-3-hexene (CC/C=C\CC) should match; E-3-hexene (CC/C=C/CC) should not.
        let z_hexene = parse("CC/C=C\\CC").unwrap();
        let e_hexene = parse("CC/C=C/CC").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2]=O.[O:3]=[C:4]";
        assert!(
            !run_reactants(smirks, &[&z_hexene]).unwrap().is_empty(),
            "Z-template must match Z-3-hexene"
        );
        assert!(
            run_reactants(smirks, &[&e_hexene]).unwrap().is_empty(),
            "Z-template must reject E-3-hexene"
        );
    }

    #[test]
    fn ez_stereo_z_template_matches_z_alkene() {
        // Template specifies Z: [C:1]/[C:2]=[C:3]\[C:4]
        // Z-2-butene (C/C=C\C) should match.
        let z_alkene = parse("C/C=C\\C").unwrap();
        let e_alkene = parse("C/C=C/C").unwrap();
        let smirks = "[C:1]/[C:2]=[C:3]\\[C:4]>>[C:1][C:2][C:3][C:4]";
        assert!(
            !run_reactants(smirks, &[&z_alkene]).unwrap().is_empty(),
            "Z-template must match Z-alkene"
        );
        assert!(
            run_reactants(smirks, &[&e_alkene]).unwrap().is_empty(),
            "Z-template must reject E-alkene"
        );
    }

    // -----------------------------------------------------------------------
    // E/Z stereo transfer & creation in products (issue #50)
    //
    // Geometry is verified by comparing the product's canonical SMILES to the
    // canonical SMILES of a reference molecule of known E/Z — exact, not "some
    // Up/Down survives". (CipCode-based verification lives in the Python tests;
    // chematic-rxn cannot depend on chematic-chem without a dependency cycle.)
    // -----------------------------------------------------------------------

    /// Canonical SMILES of the single product of `smirks` applied to `inputs`.
    fn product_canon(smirks: &str, inputs: &[&str]) -> String {
        use chematic_smiles::canonical_smiles;
        let mols: Vec<Molecule> = inputs.iter().map(|s| parse(s).unwrap()).collect();
        let refs: Vec<&Molecule> = mols.iter().collect();
        let results = run_reactants(smirks, &refs).unwrap();
        assert!(!results.is_empty(), "no product for {smirks} on {inputs:?}");
        canonical_smiles(&results[0][0])
    }

    fn canon(smiles: &str) -> String {
        chematic_smiles::canonical_smiles(&parse(smiles).unwrap())
    }

    /// Writer-invariant E/Z of the first C=C double bond: `Some(true)` = E,
    /// `Some(false)` = Z, `None` = no specified geometry. Reuses the crate's
    /// own `ez_stereo_outward` (same convention as the #21 filter): equal
    /// outward directions → Z, opposite → E.
    fn double_bond_is_e(smiles: &str) -> Option<bool> {
        let mol = parse(smiles).unwrap();
        let (a1, a2) = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(_, b)| (b.atom1, b.atom2))?;
        let sa = ez_stereo_outward(&mol, a1, a2)?;
        let sb = ez_stereo_outward(&mol, a2, a1)?;
        Some(sa != sb)
    }

    #[test]
    fn issue50_transfer_identity_preserves_e() {
        // Identity SMIRKS on E-2-butene must yield an E product (was Z before Fix A).
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C/C"]),
            canon("C/C=C/C"),
        );
    }

    #[test]
    fn issue50_transfer_identity_preserves_z() {
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["C/C=C\\C"]),
            canon("C/C=C\\C"),
        );
    }

    #[test]
    fn issue50_create_e_from_template() {
        // Product template introduces an E double bond from a saturated chain.
        assert_eq!(
            product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]/[C:4]", &["CCCC"]),
            canon("C/C=C/C"),
        );
    }

    #[test]
    fn issue50_create_z_from_template() {
        assert_eq!(
            product_canon("[C:1][C:2][C:3][C:4]>>[C:1]/[C:2]=[C:3]\\[C:4]", &["CCCC"]),
            canon("C/C=C\\C"),
        );
    }

    #[test]
    fn issue50_transfer_remote_reaction_keeps_e() {
        // Reaction at a remote site (aldehyde→alcohol) must not disturb the
        // E geometry of a carried-through alkene. The canonical writer may pick
        // `/C=C/` or `\C=C\` (both E) depending on traversal, so assert geometry
        // directly rather than the exact string.
        let got = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C/CC=O"]);
        assert_eq!(
            double_bond_is_e(&got),
            Some(true),
            "E geometry must survive a remote edit"
        );
        // And a Z input stays Z.
        let got_z = product_canon("[CH:1]=O>>[C:1]O", &["CC/C=C\\CC=O"]);
        assert_eq!(
            double_bond_is_e(&got_z),
            Some(false),
            "Z geometry must survive a remote edit"
        );
    }

    #[test]
    fn issue50_geometry_is_deterministic() {
        // The pre-fix bug was nondeterministic (FxHashMap iteration order).
        // The same transform must give the same geometry on every run.
        let first = product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]);
        for _ in 0..6 {
            assert_eq!(
                product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC"]),
                first,
                "product geometry must be deterministic across runs"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Reaction-transform performance regression witnesses (see
    // docs/rfcs/reaction_transform_perf.md). Root cause: `chematic-smiles`'s
    // `canonical_smiles()` wrote the winning individualize-refine branch's
    // string, threw it away, and had `winning_individualized_ranks`'s caller
    // write it a *second* time -- one fully redundant DFS-and-format pass on
    // every single call, tied or not. Fixed by returning the already-written
    // string instead of recomputing it. These three cases mirror the ones
    // used to characterize and fix the regression:
    // (a) a highly symmetric molecule (many individualize-refine branches --
    //     this is the case that actually reproduces a large, measured slowdown
    //     between chematic 0.4.25 and 0.4.30, NOT run_reactants match/product
    //     volume, which stayed flat across versions);
    // (b) an E/Z stereo control, since the fix touches the same
    //     canonical-writer code path `resolve_ez_markers` depends on;
    // (c) a negative control (asymmetric, no ties) that should show only the
    //     universal (small, single-redundant-write) improvement, not the
    //     symmetric-molecule-specific one.
    // -----------------------------------------------------------------------

    #[test]
    fn perf_witness_a_symmetric_molecule_product_is_correct() {
        // Positive witness: adamantane (Td cage symmetry, 24
        // individualize-refine branches at time of writing) run through a
        // simple ring-opening SMIRKS. The fix must not change *which* string
        // wins -- only how many times it gets written -- so the product's
        // canonical SMILES must still round-trip to the same structure.
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let results = run_reactants("[C:1][C:2]>>[C:1][C:2]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "expected at least one C-C bond match");
        let canon = chematic_smiles::canonical_smiles(&results[0][0]);
        // Adamantane's own canonical form is a fixed point of this
        // identity-shaped SMIRKS: it must reparse to the exact same molecule
        // (same atom/bond count -- the transform doesn't add/remove atoms).
        let reparsed = parse(&canon).unwrap();
        assert_eq!(reparsed.atom_count(), mol.atom_count());
        assert_eq!(reparsed.bond_count(), mol.bond_count());
    }

    #[test]
    fn perf_witness_b_ez_stereo_control_survives_symmetric_fix() {
        // Stereo control: identity/remote transforms on the E/Z pair used in
        // the perf investigation's own fixture
        // (crates/chematic-rxn/fixtures/witness_molecules.smi) must still
        // preserve exact geometry -- this is the non-negotiable issue #50
        // gate, re-run here against the specific molecules this perf fix
        // touched.
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C/CC(=O)O"]),
            canon("CC/C=C/CC(=O)O"),
            "(E)-hex-3-enoic acid must keep its E geometry"
        );
        assert_eq!(
            product_canon("[C:1]=[C:2]>>[C:1]=[C:2]", &["CC/C=C\\CC(=O)O"]),
            canon("CC/C=C\\CC(=O)O"),
            "(Z)-hex-3-enoic acid must keep its Z geometry"
        );
    }

    #[test]
    fn perf_witness_c_negative_control_asymmetric_molecule() {
        // Negative control: aspirin has no automorphism ties (every ring
        // carbon has a distinct substitution environment), so
        // `winning_individualized_ranks` takes exactly one branch. This
        // exercises the same code path as (a) but should show only the
        // universal single-redundant-write saving, not a
        // symmetric-molecule-specific one -- included so a future reader can
        // tell the two effects apart empirically, not just by reasoning.
        let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap(); // aspirin
        let results = run_reactants("[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]", &[&mol]).unwrap();
        assert!(!results.is_empty(), "expected the carboxylic acid to match");
        let canon = chematic_smiles::canonical_smiles(&results[0][0]);
        let reparsed = parse(&canon).unwrap();
        assert_eq!(reparsed.atom_count(), mol.atom_count() + 1); // +1 methyl carbon
    }

    // -------------------------------------------------------------------
    // Match-level reaction application (issue #225)
    // -------------------------------------------------------------------

    /// `run_reactants(smirks, reactants)` must equal
    /// `find_reaction_matches(...).filter_map(|m| apply_reaction_match(...))`
    /// -- the exact equivalence issue #225's proposed API is built on.
    /// Checked against several existing SMIRKS/reactant pairs already used
    /// elsewhere in this file, not just one.
    #[test]
    fn find_and_apply_match_equals_run_reactants() {
        let cases: Vec<(&str, Vec<chematic_core::Molecule>)> = vec![
            ("[N:1]>>[N:1]", vec![parse("NCCN").unwrap()]),
            (
                "[N:1].[C:2]>>[N:1][C:2]",
                vec![parse("N").unwrap(), parse("C").unwrap()],
            ),
            ("[C:1][C:2]>>[C:1].[C:2]", vec![parse("CC").unwrap()]),
            (
                "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]",
                vec![parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap()],
            ),
        ];

        for (smirks, mols) in &cases {
            let reactants: Vec<&Molecule> = mols.iter().collect();
            let direct = run_reactants(smirks, &reactants).unwrap();

            let matches = find_reaction_matches(smirks, &reactants).unwrap();
            let via_matches: Vec<Vec<Molecule>> = matches
                .iter()
                .filter_map(|m| apply_reaction_match(smirks, &reactants, m, true).unwrap())
                .collect();

            assert_eq!(
                direct.len(),
                via_matches.len(),
                "{smirks}: product-set count must match"
            );
            for (d, v) in direct.iter().zip(via_matches.iter()) {
                assert_eq!(
                    d.len(),
                    v.len(),
                    "{smirks}: product count per set must match"
                );
                for (dp, vp) in d.iter().zip(v.iter()) {
                    assert_eq!(
                        chematic_smiles::canonical_smiles(dp),
                        chematic_smiles::canonical_smiles(vp),
                        "{smirks}: product molecule must match run_reactants exactly"
                    );
                }
            }
        }
    }

    /// `run_reactants_strict` (carry_substituents=false) must also compose
    /// the same way as the `carry_substituents=true` case above.
    #[test]
    fn find_and_apply_match_equals_run_reactants_strict() {
        let mol = parse("CC(=O)OC1=CC=CC=C1C(=O)O").unwrap();
        let reactants = [&mol];
        let smirks = "[OH:1]-[C:2]=[O:3]>>C-[O:1]-[C:2]=[O:3]";

        let direct = run_reactants_strict(smirks, &reactants).unwrap();
        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        let via_matches: Vec<Vec<Molecule>> = matches
            .iter()
            .filter_map(|m| apply_reaction_match(smirks, &reactants, m, false).unwrap())
            .collect();

        assert_eq!(direct.len(), via_matches.len());
        for (d, v) in direct.iter().zip(via_matches.iter()) {
            for (dp, vp) in d.iter().zip(v.iter()) {
                assert_eq!(
                    chematic_smiles::canonical_smiles(dp),
                    chematic_smiles::canonical_smiles(vp)
                );
            }
        }
    }

    /// The core motivating use case from issue #225: enumerate matches
    /// independently of applying them, reject some based on a property of
    /// the match itself, and apply only the accepted ones -- without
    /// discarding the legitimate matches along with the rejected one.
    #[test]
    fn selective_match_application() {
        let mol = parse("NCCN").unwrap();
        let reactants = [&mol];
        let smirks = "[N:1]>>[N:1]";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        assert_eq!(matches.len(), 2, "two N atoms in NCCN → two matches");

        // Reject the match touching the higher-numbered atom (arbitrary
        // match-specific policy standing in for RENKIN's real ring-bond
        // rejection rule), keep the other.
        let positions: Vec<_> = matches
            .iter()
            .map(|m| m.atom_map_positions(smirks).unwrap()[&1])
            .collect();
        let keep_idx = if positions[0].1.0 < positions[1].1.0 {
            0
        } else {
            1
        };

        let applied = apply_reaction_match(smirks, &reactants, &matches[keep_idx], true).unwrap();
        assert!(applied.is_some(), "the accepted match must still apply");

        // Applying only one match must yield exactly one of the two
        // products `run_reactants` would have returned for the whole call,
        // not both and not neither.
        let full = run_reactants(smirks, &reactants).unwrap();
        assert_eq!(full.len(), 2);
        let applied_canon = chematic_smiles::canonical_smiles(&applied.unwrap()[0]);
        assert!(
            full.iter()
                .any(|ps| chematic_smiles::canonical_smiles(&ps[0]) == applied_canon),
            "the selectively-applied product must be one of run_reactants's own outputs"
        );
    }

    /// `ReactionMatch::atom_map_positions` must resolve atom_map:1 to the
    /// actual matched N atom, for each of the two NCCN matches separately.
    #[test]
    fn atom_map_positions_resolves_matched_atom() {
        let mol = parse("NCCN").unwrap();
        let reactants = [&mol];
        let smirks = "[N:1]>>[N:1]";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        let mut matched_atoms: Vec<AtomIdx> = matches
            .iter()
            .map(|m| {
                let positions = m.atom_map_positions(smirks).unwrap();
                let (reactant_slot, atom_idx) = positions[&1];
                assert_eq!(reactant_slot, 0, "single-reactant SMIRKS: slot must be 0");
                assert_eq!(mol.atom(atom_idx).element.symbol(), "N");
                atom_idx
            })
            .collect();
        matched_atoms.sort_by_key(|a| a.0);
        assert_eq!(
            matched_atoms,
            vec![AtomIdx(0), AtomIdx(3)],
            "NCCN's two N atoms are at index 0 and 3"
        );
    }

    /// `apply_reaction_match` must return `Ok(None)` -- not an error and not
    /// a product -- for a match whose product set fails the existing
    /// valence filter, matching the case [`run_reactants`] silently drops
    /// (`overvalent_product_filtered_oxygen` above).
    #[test]
    fn apply_reaction_match_none_on_valence_violation() {
        let ethanol = parse("CCO").unwrap();
        let reactants = [&ethanol];
        let smirks = "[O:1]>>[O:1](C)C";

        let matches = find_reaction_matches(smirks, &reactants).unwrap();
        assert_eq!(matches.len(), 1, "exactly one O in ethanol");

        let applied = apply_reaction_match(smirks, &reactants, &matches[0], true).unwrap();
        assert!(
            applied.is_none(),
            "over-valenced product must come back as Ok(None), not Some(..) or Err(..)"
        );
    }

    /// `find_reaction_matches` and `apply_reaction_match` must propagate
    /// the same `ReactantCountMismatch` error `run_reactants` does when the
    /// number of input molecules doesn't match the SMIRKS's reactant-slot
    /// count.
    #[test]
    fn find_and_apply_match_reactant_count_mismatch_errors() {
        let mol = parse("C").unwrap();
        let smirks = "[N:1].[C:2]>>[N:1][C:2]";

        let find_err = find_reaction_matches(smirks, &[&mol]);
        assert!(matches!(
            find_err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));

        let dummy_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let apply_err = apply_reaction_match(smirks, &[&mol], &dummy_match, true);
        assert!(matches!(
            apply_err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }

    /// `apply_reaction_match` must also reject a `ReactionMatch` whose own
    /// `per_reactant` shape doesn't match the SMIRKS being applied against
    /// (e.g. a match obtained from a different SMIRKS), rather than
    /// panicking on an out-of-bounds index into `template_atom_maps`.
    #[test]
    fn apply_reaction_match_rejects_mismatched_match_shape() {
        let n_mol = parse("N").unwrap();
        let c_mol = parse("C").unwrap();
        let reactants = [&n_mol, &c_mol];
        let smirks = "[N:1].[C:2]>>[N:1][C:2]";

        // A match shaped for a single-reactant SMIRKS, applied against a
        // two-reactant one.
        let mismatched_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let err = apply_reaction_match(smirks, &reactants, &mismatched_match, true);
        assert!(matches!(
            err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }

    /// `ReactionMatch::atom_map_positions` must likewise reject a
    /// shape mismatch rather than panicking.
    #[test]
    fn atom_map_positions_rejects_mismatched_match_shape() {
        let mismatched_match = ReactionMatch {
            per_reactant: vec![FxHashMap::default()],
        };
        let err = mismatched_match.atom_map_positions("[N:1].[C:2]>>[N:1][C:2]");
        assert!(matches!(
            err,
            Err(TransformError::ReactantCountMismatch {
                expected: 2,
                got: 1
            })
        ));
    }
}
