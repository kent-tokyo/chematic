use super::*;
use chematic_core::ez_markers::{EzFact, EzFactEnd, raw_direction, read_ez_facts, write_ez_facts};

/// RDKit's legacy CIP ranks (`assignStereochemistry`'s atom invariants,
/// which pick a double bond's stereo atoms): atomic number and isotope,
/// refined by the sorted ranks of the neighbours, each counted by bond order
/// (twice the order, so a double bond's partner counts twice), and implicit
/// H as the lowest rank. Explicit H atoms are ordinary neighbours, so an
/// explicit-H molecule ranks as its implicit-H form.
pub(super) fn rdkit_legacy_cip_ranks(mol: &Molecule) -> Vec<u32> {
    let n = mol.atom_count();
    let key0 = |i: usize| {
        let a = mol.atom(AtomIdx(i as u32));
        (
            u32::from(a.element.atomic_number()),
            u32::from(a.isotope.unwrap_or(0)),
        )
    };
    let mut keys0: Vec<(u32, u32)> = (0..n).map(key0).collect();
    keys0.sort_unstable();
    keys0.dedup();
    let mut ranks: Vec<u32> = (0..n)
        .map(|i| keys0.binary_search(&key0(i)).unwrap_or(0) as u32 + 1)
        .collect();
    let mut distinct = keys0.len();
    for _ in 0..n {
        let keys: Vec<(u32, Vec<u32>)> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut nbrs: Vec<u32> = Vec::new();
                for (nb, b) in mol.neighbors(idx) {
                    let count = match mol.bond(b).order {
                        BondOrder::Double => 4,
                        BondOrder::Triple => 6,
                        BondOrder::Aromatic => 3,
                        BondOrder::Zero => 0,
                        _ => 2,
                    };
                    nbrs.extend(std::iter::repeat_n(ranks[nb.0 as usize], count));
                }
                // Implicit H: rank 0, below every atom.
                let h = chematic_core::implicit_hcount(mol, idx) as usize;
                nbrs.extend(std::iter::repeat_n(0, 2 * h));
                nbrs.sort_unstable_by(|a, b| b.cmp(a));
                (ranks[i], nbrs)
            })
            .collect();
        let mut sorted: Vec<&(u32, Vec<u32>)> = keys.iter().collect();
        sorted.sort_unstable();
        sorted.dedup();
        let next: Vec<u32> = keys
            .iter()
            .map(|k| sorted.binary_search(&k).unwrap_or(0) as u32 + 1)
            .collect();
        let now = sorted.len();
        ranks = next;
        if now == distinct {
            break;
        }
        distinct = now;
    }
    ranks
}

/// RDKit's reading of a reactant double bond's E/Z at one end: its stereo
/// atom (the substituent of higher legacy CIP rank, which RDKit records
/// whatever the markers), that atom's side seen from the end, and the end's
/// other substituent. `None` when the end has no stereo atom.
pub(super) fn rdkit_reactant_ez_end(
    mol: &Molecule,
    ranks: &[u32],
    double: BondIdx,
    end: AtomIdx,
) -> Option<(AtomIdx, bool, Option<AtomIdx>)> {
    use chematic_core::ez_markers::marker_side;
    let subs: Vec<(AtomIdx, BondIdx)> = mol
        .neighbors(end)
        .filter(|&(_, b)| b != double && mol.bond(b).order != BondOrder::Double)
        .collect();
    let (stereo, other) = match subs.as_slice() {
        [one] => (*one, None),
        [a, b] => match ranks[a.0.0 as usize].cmp(&ranks[b.0.0 as usize]) {
            std::cmp::Ordering::Less => (*b, Some(*a)),
            std::cmp::Ordering::Greater => (*a, Some(*b)),
            std::cmp::Ordering::Equal => return None,
        },
        _ => return None,
    };
    let side = marker_side(mol, stereo.1, end)
        .or_else(|| other.and_then(|o| marker_side(mol, o.1, end)).map(|up| !up))?;
    Some((stereo.0, side, other.map(|o| o.0)))
}

/// E/Z of product double bonds copied from a reactant double bond, as
/// RDKit's reaction runner carries it: the reactant bond's stereo atoms are
/// followed into the product. A stereo atom left out of the product is
/// replaced by the end's other substituent, on the other side (the E/Z
/// survives losing the bond that carried its marker); a stereo atom still in
/// the product but no longer bonded to its end drops the E/Z (a ring
/// expansion inserting an atom on that side). A template double bond with
/// template neighbours at both ends takes the template's stereo: the
/// template's own when it spells one, else none.
pub(super) fn rdkit_carried_ez(
    product: &mut TracedProduct,
    reactants: &[&Molecule],
    template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
) {
    if !reactants
        .iter()
        .any(|r| r.bonds().any(|(b, _)| raw_direction(r, b).is_some()))
    {
        return;
    }
    let mol = &product.molecule;
    let template_atom_of: FxHashMap<AtomIdx, AtomIdx> = template_idx_to_new
        .iter()
        .enumerate()
        .filter_map(|(t, p)| p.map(|p| (p, AtomIdx(t as u32))))
        .collect();
    let mut dropped: Vec<BondIdx> = Vec::new();
    let mut derived: Vec<EzFact> = Vec::new();
    let mut rank_cache: FxHashMap<usize, Vec<u32>> = FxHashMap::default();
    for (pb, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let (pa, pz) = (bond.atom1, bond.atom2);
        let (Some(sa), Some(sz)) = (
            product.atom_sources[pa.0 as usize],
            product.atom_sources[pz.0 as usize],
        ) else {
            continue;
        };
        if sa.reactant != sz.reactant {
            continue;
        }
        let r = reactants[sa.reactant];
        let Some((rb, rbond)) = r.bond_between(sa.atom, sz.atom) else {
            continue;
        };
        if rbond.order != BondOrder::Double || r.atom(sa.atom).aromatic && r.atom(sz.atom).aromatic
        {
            continue;
        }
        let ranks = rank_cache
            .entry(sa.reactant)
            .or_insert_with(|| rdkit_legacy_cip_ranks(r));
        let (Some(ea), Some(ez)) = (
            rdkit_reactant_ez_end(r, ranks, rb, sa.atom),
            rdkit_reactant_ez_end(r, ranks, rb, sz.atom),
        ) else {
            continue;
        };
        match template_double_stereo(template, &template_atom_of, pa, pz) {
            // The template spells this double bond's stereo: keep it.
            Some((_, true)) => continue,
            // `[C:4][C:1]=[C:2][C:3]>>[C:4][C:1]=[C:2][C:3]` on an alkene,
            // BTMR0072 on a bis-allylic CH2: no E/Z.
            Some((true, false)) => {
                dropped.push(pb);
                continue;
            }
            _ => {}
        }
        let mut ends: Vec<EzFactEnd> = Vec::with_capacity(2);
        let mut drop = false;
        for (p_end, (stereo, side, other)) in [(pa, ea), (pz, ez)] {
            let reference = match product_atom_of(product, sa.reactant, stereo) {
                Some(ps) => match mol.bond_between(p_end, ps) {
                    Some((b, _)) => Some((b, side)),
                    None => {
                        drop = true;
                        None
                    }
                },
                // Only the end's other reactant substituent replaces it; an
                // atom the template adds does not (a vinyl iodide reduced to
                // `C=C[H]` has no E/Z, BTMR0671).
                None => other
                    .and_then(|o| product_atom_of(product, sa.reactant, o))
                    .and_then(|po| mol.bond_between(p_end, po))
                    .map(|(b, _)| (b, !side)),
            };
            match reference {
                Some((reference, up)) => ends.push(EzFactEnd {
                    atom: p_end,
                    reference,
                    up,
                }),
                None => {
                    drop = true;
                    break;
                }
            }
        }
        if drop {
            dropped.push(pb);
        } else if let [a, b] = ends.as_slice() {
            derived.push(EzFact {
                double: pb,
                ends: [*a, *b],
            });
        }
    }
    apply_carried_ez_changes(product, &dropped, derived);
}

/// Drop carried E/Z when the product template explicitly provides a
/// substituent at both ends of a double bond but no direction markers.
///
/// This is a property of the complete product template, so it must be applied
/// once after all reactants have been assembled. Applying it during each
/// reactant pass makes the result depend on reactant order (RDKit #9691).
pub(super) fn clear_product_template_unspecified_ez(
    product: &mut TracedProduct,
    template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
) {
    let dropped: Vec<BondIdx> = template
        .bonds()
        .filter(|(_, bond)| bond.order == BondOrder::Double)
        .filter(|&(template_bond, bond)| {
            [bond.atom1, bond.atom2].iter().all(|&end| {
                template
                    .neighbors(end)
                    .any(|(_, neighbour_bond)| neighbour_bond != template_bond)
            }) && ![bond.atom1, bond.atom2].iter().all(|&end| {
                template.neighbors(end).any(|(_, neighbour_bond)| {
                    neighbour_bond != template_bond
                        && raw_direction(template, neighbour_bond).is_some()
                })
            })
        })
        .filter_map(|(_, bond)| {
            let a = template_idx_to_new
                .get(bond.atom1.0 as usize)
                .copied()
                .flatten()?;
            let b = template_idx_to_new
                .get(bond.atom2.0 as usize)
                .copied()
                .flatten()?;
            product
                .molecule
                .bond_between(a, b)
                .filter(|(_, product_bond)| product_bond.order == BondOrder::Double)
                .map(|(product_bond, _)| product_bond)
        })
        .collect();
    apply_carried_ez_changes(product, &dropped, Vec::new());
}

/// The template double bond between two product atoms, with whether each
/// template end has another template neighbour and a marked one.
fn template_double_stereo(
    template: &Molecule,
    template_atom_of: &FxHashMap<AtomIdx, AtomIdx>,
    a: AtomIdx,
    b: AtomIdx,
) -> Option<(bool, bool)> {
    let (&ta, &tb) = (template_atom_of.get(&a)?, template_atom_of.get(&b)?);
    let (tbond, entry) = template.bond_between(ta, tb)?;
    (entry.order == BondOrder::Double).then(|| {
        let ends = [ta, tb];
        let neighbours = ends
            .iter()
            .all(|&end| template.neighbors(end).any(|(_, nb)| nb != tbond));
        let marked = ends.iter().all(|&end| {
            template
                .neighbors(end)
                .any(|(_, nb)| nb != tbond && raw_direction(template, nb).is_some())
        });
        (neighbours, marked)
    })
}

fn product_atom_of(product: &TracedProduct, reactant: usize, atom: AtomIdx) -> Option<AtomIdx> {
    product
        .atom_sources
        .iter()
        .position(|source| source.is_some_and(|s| s.reactant == reactant && s.atom == atom))
        .map(|index| AtomIdx(index as u32))
}

fn apply_carried_ez_changes(
    product: &mut TracedProduct,
    dropped: &[BondIdx],
    derived: Vec<EzFact>,
) {
    if dropped.is_empty() && derived.is_empty() {
        return;
    }
    // Markers of dropped double bonds go first, so they cannot be read as
    // another double bond's (or conflict with its).
    let mut trial = product.molecule.clone();
    let keep: FxHashSet<BondIdx> = derived
        .iter()
        .flat_map(|f| f.ends.iter().map(|e| e.reference))
        .collect();
    for &pb in dropped {
        let bond = trial.bond(pb);
        for end in [bond.atom1, bond.atom2] {
            let subs: Vec<BondIdx> = trial
                .neighbors(end)
                .filter(|&(_, b)| b != pb && !keep.contains(&b))
                .map(|(_, b)| b)
                .collect();
            for b in subs {
                if matches!(trial.bond(b).order, BondOrder::Up | BondOrder::Down) {
                    trial.set_bond_order(b, BondOrder::Single);
                }
            }
        }
    }
    let mut facts = read_ez_facts(&trial)
        .map(|(facts, _)| facts)
        .unwrap_or_default();
    facts.retain(|f| !dropped.contains(&f.double) && !derived.iter().any(|d| d.double == f.double));
    facts.extend(derived);
    if write_ez_facts(&mut trial, &facts) {
        product.molecule = trial;
    }
}

/// Atoms flagged aromatic with a plain double bond between two of them or
/// with fewer than two aromatic bonds, or an aromatic bond on an atom not
/// flagged aromatic: a template rewrote part of an aromatic system.
pub(super) fn mixed_aromatic_state(mol: &Molecule) -> bool {
    mol.bonds().any(|(_, b)| {
        b.order == BondOrder::Double && mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic
    }) || mol.atoms().any(|(idx, a)| {
        let aromatic_bonds = mol
            .neighbors(idx)
            .filter(|&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
            .count();
        (a.aromatic && aromatic_bonds < 2) || (!a.aromatic && aromatic_bonds > 0)
    })
}

/// Per atom, whether it is in a mixed aromatic state (see
/// [`mixed_aromatic_state`]).
pub(super) fn mixed_aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    let mut out = vec![false; mol.atom_count()];
    for (_, b) in mol.bonds() {
        if b.order == BondOrder::Double && mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic
        {
            out[b.atom1.0 as usize] = true;
            out[b.atom2.0 as usize] = true;
        }
    }
    for (idx, a) in mol.atoms() {
        let aromatic_bonds = mol
            .neighbors(idx)
            .filter(|&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
            .count();
        if (a.aromatic && aromatic_bonds < 2) || (!a.aromatic && aromatic_bonds > 0) {
            out[idx.0 as usize] = true;
        }
    }
    out
}

/// Whether the template made a product's mixed aromatic state: some mixed
/// atom is new or was not mixed in its reactant. A state carried over from
/// the reactant unchanged (an aromatic ring written with `#`) stays as
/// RDKit's product keeps it.
pub(super) fn mixed_state_introduced(product: &TracedProduct, reactants: &[&Molecule]) -> bool {
    let mixed = mixed_aromatic_atoms(&product.molecule);
    let mut reactant_mixed: Vec<Option<Vec<bool>>> = vec![None; reactants.len()];
    mixed.iter().enumerate().any(|(i, &m)| {
        if !m {
            return false;
        }
        let Some(src) = product.atom_sources.get(i).copied().flatten() else {
            return true;
        };
        let rm = reactant_mixed[src.reactant]
            .get_or_insert_with(|| mixed_aromatic_atoms(reactants[src.reactant]));
        !rm[src.atom.0 as usize]
    })
}

/// `mol` after RDKit's sanitize steps for a mixed aromatic state: the
/// remaining aromatic bonds kekulized, aromatic flags cleared, aromaticity
/// perceived again with RDKit's model. `None` when the aromatic bonds do
/// not kekulize (RDKit's sanitize fails too).
pub(super) fn rdkit_sanitized_view(mol: &Molecule) -> Option<Molecule> {
    // RDKit's kekulization refuses an aromatic atom outside a ring.
    let in_ring = atoms_in_rings(mol);
    if mol
        .atoms()
        .any(|(idx, a)| a.aromatic && !in_ring[idx.0 as usize])
    {
        return None;
    }
    // RDKit kekulizes the aromatic bonds; see [`rdkit_kekule_candidate`].
    let mut kekule = if mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        let k = chematic_core::kekulize_with(mol, |idx| rdkit_kekule_candidate(mol, idx)).ok()?;
        chematic_core::apply_kekule(mol, &k)
    } else {
        mol.clone()
    };
    let flagged: Vec<AtomIdx> = kekule
        .atoms()
        .filter(|(_, a)| a.aromatic)
        .map(|(i, _)| i)
        .collect();
    for idx in flagged {
        kekule.set_atom_aromatic(idx, false);
    }
    let mut view = match &*chematic_perception::apply_aromaticity_rdkit_parity_shared(&kekule) {
        Ok(view) => view.clone(),
        Err(_) => kekule,
    };
    // RDKit recomputes the H count of a carbon read without brackets: one
    // left non-aromatic takes the hydrogens its Kekulé bonds leave room for
    // (an explicit-H ring CH next to a carbon the template made sp3 becomes
    // CH2, BTMR0655).
    for (idx, atom) in mol.atoms() {
        let now = view.atom(idx);
        if atom.aromatic
            && atom.element == chematic_core::Element::C
            && atom.charge == 0
            && !now.aromatic
            && now.hydrogen_count.is_some()
        {
            view.set_hydrogen_count(idx, None);
        }
    }
    Some(view)
}

/// Whether RDKit's kekulization would give atom `idx` of `mol` a double
/// bond. RDKit makes an atom a candidate when its bond orders (aromatic
/// bonds as 1) and explicit H come to one less than its default valence,
/// `sbo == dv - 1`. That rule decides for an atom that already has a
/// non-aromatic double or triple bond (a benzene carbon whose ring bond a
/// template rewrote as `=` still is a candidate, RDKit's `C1=C=CC=CC=1`; a
/// ring carbon with a substituent or a second ring bond rewritten is not)
/// for a chalcogen cation (`[o+]`/`[s+]` with two ring bonds, as in
/// pyrylium, is one; a thiophene S-oxide `[O-][s+]1cccc1` is not) and for a
/// neutral N or P (pyridine's N is one; a ring-fusion N with three aromatic
/// bonds is not, also when its third ring neighbour is a `c(=O)` carbon, as
/// in imidazo[1,2-a]pyrimidin-5-one, where chematic's bridgehead fallback
/// does not apply). Other atoms follow chematic's ordinary rule.
pub(super) fn rdkit_kekule_candidate(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    let pi = mol
        .neighbors(idx)
        .any(|(_, b)| matches!(mol.bond(b).order, BondOrder::Double | BondOrder::Triple));
    let chalcogen_cation =
        matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52) && atom.charge == 1;
    let neutral_pnictogen = matches!(atom.element.atomic_number(), 7 | 15)
        && atom.charge == 0
        && !mol
            .neighbors(idx)
            .any(|(nb, _)| mol.atom(nb).element.atomic_number() == 1);
    // A zero-order bond (RDKit's unspecified bond from a product `~` or
    // ring-led `@`) is not counted for a neutral N or P: an [nH] of an
    // explicit-H indole left with one aromatic bond and a zero-order bond
    // to its fused ring is then a candidate, as in RDKit, and the product
    // fails to kekulize (BTMR1031). Carbon keeps the ordinary rule there
    // (RDKit kekulizes `c1ccc~cc1`).
    let pnictogen_beside_zero = matches!(atom.element.atomic_number(), 7 | 15)
        && atom.charge == 0
        && mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Zero);
    if !pi && !chalcogen_cation && !neutral_pnictogen && !pnictogen_beside_zero {
        // A three-connected aromatic P, As, S, Se or Te is judged by RDKit's
        // own rule; chematic's never makes it a candidate.
        if matches!(atom.element.atomic_number(), 15 | 16 | 33 | 34 | 52)
            && atom.charge == 0
            && mol.degree(idx) >= 3
            && let Some(candidate) = rdkit_mark_dbond_cand(mol, idx)
        {
            return candidate;
        }
        return chematic_core::atom_must_be_matched(mol, idx);
    }
    let sbo: i16 = mol
        .neighbors(idx)
        .map(|(_, b)| match mol.bond(b).order {
            BondOrder::Aromatic => 1,
            order => i16::from(order.order_int()),
        })
        .sum::<i16>()
        + i16::from(atom.hydrogen_count.unwrap_or(0));
    crate::rdkit_valence::default_valence(atom.element.atomic_number(), atom.charge)
        .is_some_and(|dv| sbo == dv - 1)
}

/// RDKit's `markDbondCands` decision for an aromatic atom with no double
/// or triple bond (`Kekulize.cpp`, RDKit 2026.03): aromatic bonds count 1
/// toward `sbo`, other bonds their order, plus all hydrogens; the default
/// valence, shifted by the charge (negated for a carbocation and for the
/// early elements), is raised through the element's valence list while the
/// total valence (aromatic bonds as 1.5) exceeds it; an atom whose degree
/// plus implicit hydrogens reaches that valence is no candidate; otherwise
/// it is one when the valence is `sbo + 1`, or `sbo + 2` for a bracket
/// atom. `None` for an element RDKit leaves unrestricted or a non-aromatic
/// atom. [`rdkit_kekule_candidate`] uses it for three-connected aromatic
/// P, As, S, Se and Te: a ring carbon a template rewrites to `[S:1]` keeps
/// the reactant's aromatic flag in RDKit, and 2-chloropyridine gives the
/// hypervalent `ClS1=CC=CC=N1` (#754 fuzz lane). Applied to every aromatic
/// atom it accepted products RDKit's sanitize rejects, so the other atoms
/// keep chematic's rule.
pub(super) fn rdkit_mark_dbond_cand(mol: &Molecule, idx: AtomIdx) -> Option<bool> {
    let atom = mol.atom(idx);
    if !atom.aromatic || atom.wildcard {
        return None;
    }
    let z = atom.element.atomic_number();
    let list = crate::rdkit_valence::valence_list(z)?;
    let early = matches!(z, 1..=5 | 11..=13 | 19..=31 | 37..=49 | 55..=81 | 87..=112);
    let mut chrg = i16::from(atom.charge);
    if early || (z == 6 && chrg > 0) {
        chrg = -chrg;
    }
    let hydrogens = i16::from(chematic_core::implicit_hcount(mol, idx));
    let mut sbo = hydrogens;
    let mut total2 = 2 * hydrogens;
    let mut ignored = 0i16;
    for (_, b) in mol.neighbors(idx) {
        let bond = mol.bond(b);
        match bond.order {
            BondOrder::Aromatic => {
                sbo += 1;
                total2 += 3;
            }
            BondOrder::Zero | BondOrder::QueryAny => ignored += 1,
            BondOrder::Dative if bond.atom1 == idx => ignored += 1,
            order => {
                let o = i16::from(order.order_int());
                sbo += o;
                total2 += 2 * o;
            }
        }
    }
    let tbo = (total2 + 1) / 2;
    let mut dv = i16::from(list[0]) + chrg;
    for &v in &list[1..] {
        if tbo <= dv {
            break;
        }
        dv = i16::from(v) + chrg;
    }
    let implicit = if atom.hydrogen_count.is_none() {
        hydrogens
    } else {
        0
    };
    let total_degree = mol.degree(idx) as i16 + implicit - ignored;
    if total_degree >= dv {
        return Some(false);
    }
    Some(dv == sbo + 1 || (atom.hydrogen_count.is_some() && dv == sbo + 2))
}

/// Whether RDKit's sanitize step would accept `mol` (issue #734).
///
/// Follows RDKit's `SanitizeMol`: aromatic atoms must lie in rings and
/// aromatic systems must kekulize; RDKit's clean-up spellings are accepted
/// (a neutral five-valent N with `=O` or `#N` is read as `[N+][O-]` /
/// `[N+]=[N-]`, a neutral Cl/Br/I of valence 3, 5 or 7 bonded only to O as
/// `[X+k]([O-])…`); and each atom's explicit valence — Kekulé bond orders
/// (a dative bond counts for its acceptor only) plus its explicit H count —
/// must not exceed RDKit's largest allowed valence for its element and
/// charge (see [`crate::rdkit_valence`]).
pub(super) fn product_sanitization_diagnostic(
    product: &TracedProduct,
    product_index: usize,
) -> ProductSanitizationDiagnostic {
    if let Some(rejection) = aromaticity_rejection(product, product_index) {
        return rejection;
    }
    let orders = match product_kekule_orders(product, product_index) {
        Ok(orders) => orders,
        Err(rejection) => return rejection,
    };
    let invalid = invalid_valence_atoms(product, &orders);
    if !invalid.is_empty() {
        return rejected_product(product_index, ProductRejectionReason::Valence, invalid);
    }
    ProductSanitizationDiagnostic {
        product_index,
        accepted: true,
        reason: None,
        atoms: Vec::new(),
    }
}

const MAX_DIAGNOSTIC_ATOMS: usize = 32;

fn atom_rejection_context(
    product: &TracedProduct,
    atom: AtomIdx,
    observed_valence: Option<i16>,
    max_allowed_valence: Option<i16>,
) -> ProductAtomRejectionDiagnostic {
    let value = product.molecule.atom(atom);
    ProductAtomRejectionDiagnostic {
        atom,
        template_map: product
            .template_maps
            .get(atom.0 as usize)
            .copied()
            .flatten(),
        element: value.element,
        explicit_hydrogen_count: value.hydrogen_count,
        observed_valence,
        max_allowed_valence,
    }
}

fn rejected_product(
    product_index: usize,
    reason: ProductRejectionReason,
    mut atoms: Vec<ProductAtomRejectionDiagnostic>,
) -> ProductSanitizationDiagnostic {
    atoms.sort_unstable_by_key(|atom| atom.atom.0);
    atoms.dedup_by_key(|atom| atom.atom.0);
    atoms.truncate(MAX_DIAGNOSTIC_ATOMS);
    ProductSanitizationDiagnostic {
        product_index,
        accepted: false,
        reason: Some(reason),
        atoms,
    }
}

fn aromaticity_rejection(
    product: &TracedProduct,
    product_index: usize,
) -> Option<ProductSanitizationDiagnostic> {
    let mol = &product.molecule;
    if !mol.atoms().any(|(_, atom)| atom.aromatic) {
        return None;
    }
    let in_ring = atoms_in_rings(mol);
    let outside_ring: Vec<_> = mol
        .atoms()
        .filter(|(idx, atom)| atom.aromatic && !in_ring[idx.0 as usize])
        .take(MAX_DIAGNOSTIC_ATOMS)
        .map(|(idx, _)| atom_rejection_context(product, idx, None, None))
        .collect();
    if !outside_ring.is_empty() {
        return Some(rejected_product(
            product_index,
            ProductRejectionReason::Aromaticity,
            outside_ring,
        ));
    }

    // A plain double bond between aromatic ring atoms is not a form RDKit's
    // kekulization accepts; the written SMILES would not parse.
    mol.bonds()
        .find(|(_, bond)| {
            bond.order == BondOrder::Double
                && mol.atom(bond.atom1).aromatic
                && mol.atom(bond.atom2).aromatic
                && in_ring[bond.atom1.0 as usize]
                && in_ring[bond.atom2.0 as usize]
        })
        .map(|(_, bond)| {
            let atoms = [bond.atom1, bond.atom2]
                .into_iter()
                .map(|idx| atom_rejection_context(product, idx, None, None))
                .collect();
            rejected_product(product_index, ProductRejectionReason::Aromaticity, atoms)
        })
}

fn product_kekule_orders(
    product: &TracedProduct,
    product_index: usize,
) -> Result<Vec<BondOrder>, ProductSanitizationDiagnostic> {
    let mol = &product.molecule;
    let mut orders: Vec<BondOrder> = mol.bonds().map(|(_, bond)| bond.order).collect();
    if !orders.contains(&BondOrder::Aromatic) {
        return Ok(orders);
    }
    let kekule = chematic_core::kekulize_with(mol, |idx| rdkit_kekule_candidate(mol, idx))
        .map_err(|_| {
            let atoms = mol
                .atoms()
                .filter(|(_, atom)| atom.aromatic)
                .take(MAX_DIAGNOSTIC_ATOMS)
                .map(|(idx, _)| atom_rejection_context(product, idx, None, None))
                .collect();
            rejected_product(product_index, ProductRejectionReason::Kekulization, atoms)
        })?;
    for (bond, order) in kekule {
        orders[bond.0 as usize] = order;
    }
    Ok(orders)
}

fn explicit_valence(mol: &Molecule, orders: &[BondOrder], idx: AtomIdx) -> i16 {
    let bonds: i16 = mol
        .neighbors(idx)
        .map(|(_, bond_idx)| {
            let bond = mol.bond(bond_idx);
            match orders[bond_idx.0 as usize] {
                BondOrder::Dative if bond.atom1 == idx => 0,
                order => i16::from(order.order_int()),
            }
        })
        .sum();
    bonds + i16::from(mol.atom(idx).hydrogen_count.unwrap_or(0))
}

fn cleaned_charge_and_valence(mol: &Molecule, orders: &[BondOrder], idx: AtomIdx) -> (i8, i16) {
    let atom = mol.atom(idx);
    let mut charge = atom.charge;
    let mut used = explicit_valence(mol, orders, idx);
    if charge != 0 {
        return (charge, used);
    }
    let bond_to = |order: BondOrder, atomic_number: u8| {
        mol.neighbors(idx).any(|(neighbor, bond)| {
            orders[bond.0 as usize] == order
                && mol.atom(neighbor).element.atomic_number() == atomic_number
                && mol.atom(neighbor).charge == 0
        })
    };
    let atomic_number = atom.element.atomic_number();
    if atomic_number == 7
        && used == 5
        && (bond_to(BondOrder::Double, 8) || bond_to(BondOrder::Triple, 7))
    {
        return (1, 4);
    }
    if matches!(atomic_number, 17 | 35 | 53)
        && matches!(used, 3 | 5 | 7)
        && mol
            .neighbors(idx)
            .all(|(neighbor, _)| mol.atom(neighbor).element.atomic_number() == 8)
    {
        let doubles = mol
            .neighbors(idx)
            .filter(|&(_, bond)| orders[bond.0 as usize] == BondOrder::Double)
            .count() as i16;
        charge = doubles as i8;
        used -= doubles;
    }
    (charge, used)
}

fn invalid_valence_atoms(
    product: &TracedProduct,
    orders: &[BondOrder],
) -> Vec<ProductAtomRejectionDiagnostic> {
    use crate::rdkit_valence::max_valence;

    let mol = &product.molecule;
    let mut invalid = Vec::new();
    for (idx, atom) in mol.atoms().filter(|(_, atom)| !atom.wildcard) {
        let (charge, used) = cleaned_charge_and_valence(mol, orders, idx);
        let max = max_valence(atom.element.atomic_number(), charge);
        if max.is_some_and(|max| used > max) {
            invalid.push(atom_rejection_context(product, idx, Some(used), max));
            if invalid.len() == MAX_DIAGNOSTIC_ATOMS {
                break;
            }
        }
    }
    invalid
}
