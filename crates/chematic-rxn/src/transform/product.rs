use super::*;

/// RDKit's bond order around a carried, non-ring atom of the product: the
/// bond to `discoverer` (the neighbour the walk reached it from) first, then
/// its other bonds in reactant order. `None` for an atom in a ring (other
/// neighbours may add their bonds earlier) or with an implicit-H slot.
fn rdkit_carried_bond_order(
    mol: &Molecule,
    atom: AtomIdx,
    discoverer: Option<AtomIdx>,
) -> Option<Vec<u32>> {
    let d = discoverer?;
    if atoms_in_rings(mol)[atom.0 as usize] {
        return None;
    }
    let reactant: Vec<u32> = match mol.stereo_neighbor_order(atom) {
        Some(o) => o.to_vec(),
        None => mol.neighbors(atom).map(|(n, _)| n.0).collect(),
    };
    if reactant.contains(&chematic_core::STEREO_H_SENTINEL) || !reactant.contains(&d.0) {
        return None;
    }
    let mut order = vec![d.0];
    order.extend(reactant.into_iter().filter(|&n| n != d.0));
    Some(order)
}

fn add_product_template_bonds(
    builder: &mut MoleculeBuilder,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    source_of: &[Option<(usize, AtomIdx)>],
    input_mols: &[&Molecule],
) -> Vec<AtomIdx> {
    // Atoms of a product `:` bond are flagged aromatic, as RDKit does.
    let mut colon_atoms = Vec::new();
    for (_bidx, bond) in product_template.bonds() {
        let a_new = template_idx_to_new[bond.atom1.0 as usize].unwrap();
        let b_new = template_idx_to_new[bond.atom2.0 as usize].unwrap();
        if bond.order == BondOrder::Aromatic {
            colon_atoms.extend([a_new, b_new]);
        }
        let mut order = bond.order;
        let (mut from, mut to) = (a_new, b_new);
        // A template single bond over a reactant `/`/`\` bond keeps the
        // direction: it carries the E/Z of an adjacent double bond the
        // template leaves alone.
        if order == BondOrder::Single
            && let (Some((ma, sa)), Some((mb, sb))) =
                (source_of[a_new.0 as usize], source_of[b_new.0 as usize])
            && ma == mb
            && let Some((_, ob)) = input_mols[ma].bond_between(sa, sb)
            && matches!(ob.order, BondOrder::Up | BondOrder::Down)
        {
            order = ob.order;
            if ob.atom1 != sa {
                (from, to) = (b_new, a_new);
            }
        }
        // A product `~` keeps the reactant bond between mapped atoms and is
        // a zero-order bond otherwise.
        if order == BondOrder::QueryAny {
            order = BondOrder::Zero;
            if let (Some((ma, sa)), Some((mb, sb))) =
                (source_of[a_new.0 as usize], source_of[b_new.0 as usize])
                && ma == mb
                && let Some((_, ob)) = input_mols[ma].bond_between(sa, sb)
            {
                order = ob.order;
                if ob.atom1 != sa {
                    (from, to) = (b_new, a_new);
                }
            }
        }
        let _ = builder.add_bond(from, to, order);
    }
    colon_atoms
}

#[allow(clippy::too_many_arguments)]
fn carry_source_bonds(
    builder: &mut MoleculeBuilder,
    source_of: &[Option<(usize, AtomIdx)>],
    new_of: &[Vec<Option<AtomIdx>>],
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    reactant_bond_maps: &FxHashSet<(u16, u16)>,
) {
    let map_of = |key: (usize, AtomIdx)| {
        global_map
            .iter()
            .find(|&(_, &mapped)| mapped == key)
            .map(|(&atom_map, _)| atom_map)
    };
    let carried = source_of
        .iter()
        .enumerate()
        .filter_map(|(new, src)| src.map(|src| (src, AtomIdx(new as u32))));
    for ((mol_idx, src_idx), a_new) in carried {
        for (nb_idx, bond_idx) in input_mols[mol_idx].neighbors(src_idx) {
            let nb_key = (mol_idx, nb_idx);
            let Some(b_new) = new_of[mol_idx][nb_idx.0 as usize] else {
                continue;
            };
            if all_template_atoms.contains(&(mol_idx, src_idx))
                && all_template_atoms.contains(&nb_key)
            {
                // Both ends matched: templates decide this bond, except for
                // a reactant bond between mapped atoms the template does not
                // itself connect (RDKit, github #1387).
                let kept = match (map_of((mol_idx, src_idx)), map_of(nb_key)) {
                    (Some(a), Some(b)) => !reactant_bond_maps.contains(&(a.min(b), a.max(b))),
                    _ => false,
                };
                if !kept {
                    continue;
                }
            }
            let source_bond = input_mols[mol_idx].bond(bond_idx);
            let (a, b) = if source_bond.atom1 == src_idx {
                (a_new, b_new)
            } else {
                (b_new, a_new)
            };
            let _ = builder.add_bond(a, b, source_bond.order);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn restore_carried_atom_stereo(
    product: &mut Molecule,
    product_template: &Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
    input_mols: &[&Molecule],
    source_of: &[Option<(usize, AtomIdx)>],
    discoverer: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    keep_carried_maps: bool,
) -> Vec<((usize, AtomIdx), AtomIdx)> {
    let sourced_atoms: Vec<_> = source_of
        .iter()
        .enumerate()
        .filter_map(|(new, src)| src.map(|src| (src, AtomIdx(new as u32))))
        .collect();
    for &((mol_idx, src_idx), new_idx) in &sourced_atoms {
        let src = input_mols[mol_idx];
        if keep_carried_maps
            && let Some((class, _)) = src.atom(src_idx).chirality.nontetrahedral()
            && let Some(template_idx) = template_idx_to_new
                .iter()
                .position(|idx| *idx == Some(new_idx))
            && product_template
                .atom(AtomIdx(template_idx as u32))
                .chirality
                == Chirality::None
        {
            product.set_chirality(new_idx, Chirality::from_nontetrahedral(class, 0).unwrap());
            continue;
        }
        if keep_carried_maps
            && src.atom(src_idx).chirality.nontetrahedral().is_some()
            && !template_idx_to_new.contains(&Some(new_idx))
        {
            let order = rdkit_carried_bond_order(
                src,
                src_idx,
                discoverer.get(&(mol_idx, src_idx)).copied(),
            )
            .and_then(|order| remap_reactant_stereo_order(&order, mol_idx, src_to_new))
            .filter(|order| order_matches_final_topology(product, new_idx, order));
            match order {
                Some(order) => {
                    product.set_chirality(new_idx, src.atom(src_idx).chirality);
                    product.set_stereo_neighbor_order(new_idx, order);
                }
                None => product.set_chirality(new_idx, Chirality::None),
            }
            continue;
        }
        if !src.atom(src_idx).chirality.is_tetrahedral()
            || template_idx_to_new.contains(&Some(new_idx))
        {
            continue;
        }
        match src
            .stereo_neighbor_order(src_idx)
            .and_then(|order| remap_reactant_stereo_order(order, mol_idx, src_to_new))
            .filter(|order| order_matches_final_topology(product, new_idx, order))
        {
            Some(order) => product.set_stereo_neighbor_order(new_idx, order),
            None => product.set_chirality(new_idx, Chirality::None),
        }
    }
    sourced_atoms
}

fn transfer_folded_h_stereo(
    product: &mut Molecule,
    folded_h: &FxHashSet<(usize, AtomIdx)>,
    input_mols: &[&Molecule],
    new_of: &[Vec<Option<AtomIdx>>],
) {
    let from_core = |order: BondOrder, atom1: AtomIdx, at: AtomIdx| match (order, atom1 == at) {
        (o, true) => o,
        (BondOrder::Up, false) => BondOrder::Down,
        (BondOrder::Down, false) => BondOrder::Up,
        (o, false) => o,
    };
    for &(mol_idx, h) in folded_h {
        let src = input_mols[mol_idx];
        let Some((core, h_bond_idx)) = src.neighbors(h).next() else {
            continue;
        };
        let h_bond = src.bond(h_bond_idx);
        let h_dir = from_core(h_bond.order, h_bond.atom1, core);
        if !matches!(h_dir, BondOrder::Up | BondOrder::Down) {
            continue;
        }
        let Some(core_new) = new_of[mol_idx][core.0 as usize] else {
            continue;
        };
        let on_double = product
            .neighbors(core_new)
            .any(|(_, bond)| product.bond(bond).order == BondOrder::Double);
        let marked = product
            .neighbors(core_new)
            .any(|(_, bond)| matches!(product.bond(bond).order, BondOrder::Up | BondOrder::Down));
        if !on_double || marked {
            continue;
        }
        let other = if h_dir == BondOrder::Up {
            BondOrder::Down
        } else {
            BondOrder::Up
        };
        let carrier = product
            .neighbors(core_new)
            .find(|&(_, bond)| product.bond(bond).order == BondOrder::Single)
            .map(|(_, bond)| bond);
        if let Some(bond) = carrier {
            let atom1 = product.bond(bond).atom1;
            product.set_bond_order(bond, from_core(other, atom1, core_new));
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn carry_substituent_atoms(
    builder: &mut MoleculeBuilder,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    folded_h: &FxHashSet<(usize, AtomIdx)>,
    core_keys: Vec<(usize, AtomIdx)>,
    carry_substituents: bool,
    keep_carried_maps: bool,
    src_to_new: &mut FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> FxHashMap<(usize, AtomIdx), AtomIdx> {
    let mut visited: Vec<Vec<bool>> = input_mols
        .iter()
        .map(|molecule| vec![false; molecule.atom_count()])
        .collect();
    for &(mol_idx, atom) in all_template_atoms.iter().chain(folded_h.iter()) {
        visited[mol_idx][atom.0 as usize] = true;
    }

    // Each carried atom's neighbour from which the walk reached it. RDKit's
    // product bond order around that atom starts with this bond.
    let mut discoverer = FxHashMap::default();
    if !carry_substituents {
        return discoverer;
    }

    // The seeds are in reactant atom order so product atom order never
    // depends on hash-table iteration order.
    let mut queue: VecDeque<(usize, AtomIdx)> = core_keys.into();
    while let Some((mol_idx, current)) = queue.pop_front() {
        for (neighbor, _bond) in input_mols[mol_idx].neighbors(current) {
            let key = (mol_idx, neighbor);
            if std::mem::replace(&mut visited[mol_idx][neighbor.0 as usize], true) {
                continue;
            }
            discoverer.insert(key, current);
            let mut atom = input_mols[mol_idx].atom(neighbor).clone();
            if !keep_carried_maps {
                atom.atom_map = None;
            }
            let product_atom = builder.add_atom(atom);
            src_to_new.insert(key, product_atom);
            queue.push_back(key);
        }
    }
    discoverer
}

struct DenseAtomMaps {
    source_of: Vec<Option<(usize, AtomIdx)>>,
    new_of: Vec<Vec<Option<AtomIdx>>>,
}

fn dense_atom_maps(
    builder_atom_count: usize,
    input_mols: &[&Molecule],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
) -> DenseAtomMaps {
    let mut source_of = vec![None; builder_atom_count];
    let mut new_of: Vec<Vec<Option<AtomIdx>>> = input_mols
        .iter()
        .map(|molecule| vec![None; molecule.atom_count()])
        .collect();
    for (&source, &product) in src_to_new {
        source_of[product.0 as usize] = Some(source);
        new_of[source.0][source.1.0 as usize] = Some(product);
    }
    DenseAtomMaps { source_of, new_of }
}

fn restore_product_aromaticity(
    molecule: &mut Molecule,
    colon_atoms: &[AtomIdx],
    dearomatized: &[AtomIdx],
) {
    for &atom in colon_atoms {
        if !molecule.atom(atom).aromatic {
            molecule.set_atom_aromatic(atom, true);
        }
    }
    for &atom in dearomatized {
        let aromatic_bonds = molecule
            .neighbors(atom)
            .filter(|&(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
            .count();
        if aromatic_bonds >= 2 {
            molecule.set_atom_aromatic(atom, true);
        }
    }
}

fn reactant_has_explicit_h(
    input_mols: &[&Molecule],
    cache: &[std::cell::OnceCell<bool>],
    mol_idx: usize,
) -> bool {
    *cache[mol_idx].get_or_init(|| {
        let molecule = input_mols[mol_idx];
        molecule.atoms().any(|(atom_idx, atom)| {
            atom.element == chematic_core::Element::H && molecule.degree(atom_idx) == 1
        })
    })
}

#[allow(clippy::too_many_arguments)]
fn refill_product_hydrogens(
    molecule: &mut Molecule,
    sourced_atoms: &[((usize, AtomIdx), AtomIdx)],
    input_mols: &[&Molecule],
    template_idx_to_new: &[Option<AtomIdx>],
    explicit_h_cache: &[std::cell::OnceCell<bool>],
    mut folded_cores: Vec<AtomIdx>,
    source_of: &[Option<(usize, AtomIdx)>],
    folded_h: &FxHashSet<(usize, AtomIdx)>,
) {
    // A carried atom of an explicit-H reactant that lost a bond gets its
    // hydrogens back, matching RDKit's refill after an edit.
    for &((mol_idx, source_idx), product_idx) in sourced_atoms {
        let source = input_mols[mol_idx];
        let atom = source.atom(source_idx);
        if atom.hydrogen_count != Some(0)
            || atom.charge != 0
            || atom.aromatic
            || !atom.element.is_organic_subset()
            || atom.element == chematic_core::Element::H
            || molecule.degree(product_idx) >= source.degree(source_idx)
            || template_idx_to_new.contains(&Some(product_idx))
        {
            continue;
        }
        if reactant_has_explicit_h(input_mols, explicit_h_cache, mol_idx) {
            molecule.set_hydrogen_count(product_idx, None);
            folded_cores.push(product_idx);
        }
    }

    for core in folded_cores {
        let atom = molecule.atom(core);
        let mixed_aromatic = !atom.aromatic
            && molecule
                .neighbors(core)
                .any(|(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic);
        let element_changed = source_of[core.0 as usize].is_some_and(|(mol_idx, source_idx)| {
            input_mols[mol_idx].atom(source_idx).element != atom.element
        });
        if mixed_aromatic && element_changed {
            molecule.set_hydrogen_count(core, None);
            continue;
        }
        let count = if mixed_aromatic {
            source_of[core.0 as usize].map_or(0, |(mol_idx, source_idx)| {
                input_mols[mol_idx]
                    .neighbors(source_idx)
                    .filter(|&(neighbor, _)| folded_h.contains(&(mol_idx, neighbor)))
                    .count() as u8
            })
        } else if atom.charge != 0 && !atom.aromatic {
            let bond_order_sum: i16 = molecule
                .neighbors(core)
                .map(|(_, bond)| i16::from(molecule.bond(bond).order.order_int()))
                .sum();
            crate::rdkit_valence::implicit_hydrogens(
                atom.element.atomic_number(),
                atom.charge,
                bond_order_sum,
            )
        } else {
            chematic_core::implicit_hcount(molecule, core)
        };
        molecule.set_hydrogen_count(core, Some(0));
        for _ in 0..count {
            let mut hydrogen = chematic_core::Atom::new(chematic_core::Element::H);
            hydrogen.hydrogen_count = Some(0);
            let hydrogen_idx = molecule.add_atom(hydrogen);
            let _ = molecule.add_bond(core, hydrogen_idx, BondOrder::Single);
        }
    }
}

fn apply_rdkit_template_hydrogens(
    molecule: &mut Molecule,
    template_idx_to_new: &[Option<AtomIdx>],
) {
    for &idx in template_idx_to_new.iter().flatten() {
        let atom = molecule.atom(idx);
        let aromatic_bonds = molecule
            .neighbors(idx)
            .filter(|&(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
            .count();
        if atom.aromatic
            && atom.charge != 0
            && atom.hydrogen_count.is_none()
            && !atom.wildcard
            && aromatic_bonds == 2
        {
            let other_bonds: i16 = molecule
                .neighbors(idx)
                .filter(|&(_, bond)| molecule.bond(bond).order != BondOrder::Aromatic)
                .map(|(_, bond)| i16::from(molecule.bond(bond).order.order_int()))
                .sum();
            let rdkit_h = crate::rdkit_valence::implicit_hydrogens(
                atom.element.atomic_number(),
                atom.charge,
                3 + other_bonds,
            );
            if rdkit_h != chematic_core::implicit_hcount(molecule, idx) {
                molecule.set_hydrogen_count(idx, Some(rdkit_h));
            }
            continue;
        }
        if (atom.charge == 0 && atom.element.is_organic_subset())
            || atom.hydrogen_count.is_some()
            || atom.wildcard
            || atom.aromatic
            || molecule
                .neighbors(idx)
                .any(|(_, bond)| molecule.bond(bond).order == BondOrder::Aromatic)
        {
            continue;
        }
        let bond_order_sum: i16 = molecule
            .neighbors(idx)
            .map(|(_, bond)| {
                let bond = molecule.bond(bond);
                match bond.order {
                    BondOrder::Dative if bond.atom1 == idx => 0,
                    order => i16::from(order.order_int()),
                }
            })
            .sum();
        let rdkit_h = crate::rdkit_valence::implicit_hydrogens(
            atom.element.atomic_number(),
            atom.charge,
            bond_order_sum,
        );
        if rdkit_h != chematic_core::implicit_hcount(molecule, idx) {
            molecule.set_hydrogen_count(idx, Some(rdkit_h));
        }
    }
}

fn traced_product(
    molecule: Molecule,
    product_template: &Molecule,
    template_idx_to_new: Vec<Option<AtomIdx>>,
    sourced_atoms: &[((usize, AtomIdx), AtomIdx)],
) -> (TracedProduct, Vec<Option<AtomIdx>>) {
    let mut atom_sources = vec![None; molecule.atom_count()];
    for &((reactant, atom), product_atom) in sourced_atoms {
        atom_sources[product_atom.0 as usize] = Some(ReactantAtom { reactant, atom });
    }
    let mut template_maps = vec![None; molecule.atom_count()];
    for (template_atom, product_atom) in template_idx_to_new.iter().enumerate() {
        if let Some(product_atom) = product_atom {
            template_maps[product_atom.0 as usize] = product_template
                .atom(AtomIdx(template_atom as u32))
                .atom_map;
        }
    }
    (
        TracedProduct {
            molecule,
            atom_sources,
            template_maps,
        },
        template_idx_to_new,
    )
}

/// Build one product molecule applying full SMIRKS semantics.
///
/// Atom-mapped and new template atoms are assembled first. Substituents and
/// source bonds are then carried into the product before stereo, explicit-H,
/// aromaticity, and trace metadata are finalized.
#[allow(clippy::too_many_arguments)]
pub(super) fn build_product(
    product_template: &Molecule,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
    reactant_template_atoms: &FxHashMap<u16, ReactantMapAtom>,
    reactant_bond_maps: &FxHashSet<(u16, u16)>,
    input_mols: &[&Molecule],
    all_template_atoms: &FxHashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
    specs: &[ProductAtomSpec],
    keep_carried_maps: bool,
) -> (TracedProduct, Vec<Option<AtomIdx>>) {
    // A product holds at most the reactants' atoms and bonds plus the
    // template's: reserve them once instead of growing per atom and bond.
    let mut builder = MoleculeBuilder::with_capacity(
        input_mols.iter().map(|m| m.atom_count()).sum::<usize>() + product_template.atom_count(),
        input_mols.iter().map(|m| m.bond_count()).sum::<usize>() + product_template.bond_count(),
    );
    // What each template atom's spelling specifies; without spelling
    // information a charge always applies and only `H<n>` with n > 0 pins
    // the H count (a bare bracket atom's `Some(0)` is "unspecified", #18).
    let specs = (specs.len() == product_template.atom_count()).then_some(specs);
    let spec_of = |i: usize| -> ProductAtomSpec {
        let atom = product_template.atom(AtomIdx(i as u32));
        specs.map_or(
            ProductAtomSpec {
                charge: true,
                hcount: atom.hydrogen_count.is_some_and(|h| h > 0),
            },
            |s| s[i],
        )
    };
    let template_h = |i: usize| -> Option<u8> {
        let atom = product_template.atom(AtomIdx(i as u32));
        spec_of(i).hcount.then(|| atom.hydrogen_count.unwrap_or(0))
    };

    // template_idx_to_new[i]: new AtomIdx for product template atom i.
    let mut template_idx_to_new: Vec<Option<AtomIdx>> = vec![None; product_template.atom_count()];
    // src_to_new: (mol_idx, src_AtomIdx) → new AtomIdx in the product.
    let mut src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx> = FxHashMap::with_capacity_and_hasher(
        input_mols.iter().map(|m| m.atom_count()).sum(),
        Default::default(),
    );

    // --- Step 1: add product template atoms ---
    // core_keys: only source atoms that are mapped by THIS product template.
    // Using global_map.values() (all matched atoms across all templates) would
    // seed the BFS in Step 2 from atoms belonging to *other* product templates,
    // causing their substituents to leak into this product (issue #13).
    let mut core_keys: Vec<(usize, AtomIdx)> = product_template
        .atoms()
        .filter_map(|(_, atom)| atom.atom_map.and_then(|am| global_map.get(&am).copied()))
        .collect();
    core_keys.sort_unstable();
    core_keys.dedup();

    // Explicit hydrogen atoms (issue #734): a mapped, non-stereo atom's plain
    // H-atom neighbours that the template did not match are not carried as
    // substituents. The atom's hydrogens are re-derived from its product
    // valence and added back as explicit H atoms after assembly, so an
    // `add_hydrogens` reactant gives the same products as its implicit-H
    // form instead of keeping every H and failing the valence check
    // (`[C:1]-[C:2]>>[C:1]=[C:2]`, `[O:1]>>[O:1]C`, `[CH3:1]>>[CH3:1]`).
    let mut folded_h: FxHashSet<(usize, AtomIdx)> = FxHashSet::default();
    // Mapped aromatic atoms the template spells aliphatic (see below).
    let mut dearomatized: Vec<AtomIdx> = Vec::new();
    let mut folded_cores: Vec<AtomIdx> = Vec::new();
    let plain_h = |mol: &Molecule, idx: AtomIdx| {
        let a = mol.atom(idx);
        a.element.atomic_number() == 1
            && !a.wildcard
            && a.isotope.is_none()
            && a.charge == 0
            && a.atom_map.is_none()
            && mol.degree(idx) == 1
    };

    // Reactants given with explicit H atoms (`add_hydrogens`).
    let explicit_h_cache: Vec<std::cell::OnceCell<bool>> = input_mols
        .iter()
        .map(|_| std::cell::OnceCell::new())
        .collect();
    for (i, slot) in template_idx_to_new.iter_mut().enumerate() {
        let tmpl_atom = product_template.atom(AtomIdx(i as u32));
        let new_idx = if let Some(am) = tmpl_atom.atom_map {
            if let Some(&(mol_idx, src_idx)) = global_map.get(&am) {
                // Core atom: copy source, then override electronic state from template.
                let src_atom = input_mols[mol_idx].atom(src_idx);
                let mut new_atom = src_atom.clone();
                // Element change (issue #734): as in RDKit, a mapped atom takes
                // the product template's element when that element differs
                // from the reactant template's (`[C:1]>>[N:1]`). A wildcard
                // product atom (`[*:1]`) or an unchanged element keeps the
                // matched atom's element; a reactant query naming no single
                // element (`[C,N:1]`, `[*:1]`) takes any concrete product one.
                let reactant_tmpl = reactant_template_atoms.get(&am);
                if !tmpl_atom.wildcard {
                    let changes =
                        reactant_tmpl.is_none_or(|r| r.wildcard || r.element != tmpl_atom.element);
                    if changes && new_atom.element != tmpl_atom.element {
                        new_atom.element = tmpl_atom.element;
                        new_atom.isotope = tmpl_atom.isotope;
                    }
                }
                // A wildcard product atom (`[*:1]`) keeps the matched atom's
                // aromaticity along with its element.
                if !tmpl_atom.wildcard {
                    new_atom.aromatic = tmpl_atom.aromatic;
                }
                // Unspecified charge keeps the reactant atom's (RDKit:
                // `[O-:1]>>[O:1]` leaves the alkoxide charged), unless the
                // element changes: RDKit then takes the template atom's
                // charge too (`[N:1]>>[S:1]` on a quaternary ammonium gives
                // a neutral sulfur, not `[SH+]`; found by the SMIRKS fuzz
                // lane, #754).
                if spec_of(i).charge || new_atom.element != src_atom.element {
                    new_atom.charge = tmpl_atom.charge;
                }
                // An isotope in the template applies; none keeps the
                // reactant's (`[C:1]>>[13C:1]` labels the atom).
                if tmpl_atom.isotope.is_some() {
                    new_atom.isotope = tmpl_atom.isotope;
                }
                // H count: a template `H<n>` with n > 0 pins it (`[NH2:1]`); a
                // bare bracket atom's `Some(0)` means "unspecified" (issue #18).
                // Unspecified, RDKit keeps the reactant atom's explicit H count
                // when the atom's template degree is unchanged (`[n:1]>>[n:1]`
                // keeps pyrrole's `[nH]`) and re-derives it from valence when
                // bonds were added or removed (`[n:1]>>[n:1]C`).
                let degree_unchanged = reactant_tmpl
                    .is_some_and(|r| r.degree == product_template.degree(AtomIdx(i as u32)));
                // An element change re-derives H as well (`[C:1]>>[N:1]` on
                // `[13CH3]C` gives `CN`, not a five-valent N).
                let element_unchanged = new_atom.element == src_atom.element;
                new_atom.hydrogen_count = match template_h(i) {
                    Some(h) => Some(h),
                    None if degree_unchanged && element_unchanged => src_atom.hydrogen_count,
                    None => None,
                };
                // Chirality is intentionally left as whatever src_atom.clone()
                // produced above (i.e. the reactant's own flag, if any) --
                // correct_product_stereo, run after all bonds are added
                // below, is the sole authority on the final chirality/order
                // for every atom, whether it inherits from the reactant or
                // carries an explicit product-template @/@@.
                new_atom.atom_map = None;
                // Only organic-subset atoms: their H count can be re-derived
                // from valence.
                let explicit_h: Vec<AtomIdx> = if carry_substituents
                    && src_atom.chirality == Chirality::None
                    && src_atom.element.is_organic_subset()
                    && new_atom.element.is_organic_subset()
                {
                    input_mols[mol_idx]
                        .neighbors(src_idx)
                        .map(|(nb, _)| nb)
                        .filter(|&nb| {
                            plain_h(input_mols[mol_idx], nb)
                                && !all_template_atoms.contains(&(mol_idx, nb))
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                // An H-free mapped atom of an explicit-H reactant has its
                // count pinned to 0 by `add_hydrogens`; RDKit refills it when
                // the edit lowers its valence (`[C:1]=[O:2]>>[C:1]-[O:2]` on
                // explicit-H acetic acid gives `CC(O)O`, #754).
                let refill = explicit_h.is_empty()
                    && carry_substituents
                    && src_atom.chirality == Chirality::None
                    && src_atom.hydrogen_count == Some(0)
                    && src_atom.charge == 0
                    // An aromatic carbon the template spells aliphatic
                    // (`[#6;a:1]` in RDKit's reading) is refilled too.
                    && (!src_atom.aromatic || src_atom.element == chematic_core::Element::C)
                    && !new_atom.aromatic
                    && src_atom.element.is_organic_subset()
                    && new_atom.element.is_organic_subset()
                    && template_h(i).is_none()
                    && reactant_has_explicit_h(input_mols, &explicit_h_cache, mol_idx);
                if !explicit_h.is_empty() || refill {
                    // The folded H atoms are re-derived from valence below.
                    new_atom.hydrogen_count = template_h(i);
                    // An aromatic atom keeps its H count when the template
                    // leaves its degree, element and charge alone, as the
                    // implicit form keeps `[nH]`: valence inference cannot
                    // tell pyrrole's N from pyridine's.
                    if new_atom.hydrogen_count.is_none()
                        && new_atom.aromatic
                        && src_atom.aromatic
                        && degree_unchanged
                        && element_unchanged
                        && new_atom.charge == src_atom.charge
                    {
                        new_atom.hydrogen_count = Some(explicit_h.len() as u8);
                    }
                }
                let dearomatize = src_atom.aromatic && !new_atom.aromatic;
                let idx = builder.add_atom(new_atom);
                if dearomatize {
                    dearomatized.push(idx);
                }
                if refill {
                    folded_cores.push(idx);
                }
                if !explicit_h.is_empty() {
                    folded_h.extend(explicit_h.into_iter().map(|h| (mol_idx, h)));
                    folded_cores.push(idx);
                }
                src_to_new.insert((mol_idx, src_idx), idx);
                idx
            } else {
                // Map number not in reactants — new atom from template.
                let mut new_atom = tmpl_atom.clone();
                new_atom.atom_map = None;
                // RDKit does not pin H0 on a product-only atom: it uses
                // ordinary valence inference, as for a bare `[C]`/`[O]`.
                // A mapped reactant atom above still honours explicit H0.
                new_atom.hydrogen_count = template_h(i).filter(|&h| h != 0);
                builder.add_atom(new_atom)
            }
        } else {
            // No atom_map — entirely new atom from template. As for core atoms,
            // a bare bracket atom's `Some(0)` means "unspecified" in a product
            // template (`[*:1][C](=[O])[C]` adds an acetyl group, not three
            // radicals — issue #679); only an explicit `H<n>` with n > 0 pins
            // the count, as in SMIRKS. RDKit also leaves an explicit H0 on
            // a new atom unspecified; it is not a radical marker.
            let mut new_atom = tmpl_atom.clone();
            new_atom.atom_map = None;
            new_atom.hydrogen_count = template_h(i).filter(|&h| h != 0);
            builder.add_atom(new_atom)
        };
        *slot = Some(new_idx);
    }

    // --- Step 2: collect substituents and index source/product atoms ---
    let discoverer = carry_substituent_atoms(
        &mut builder,
        input_mols,
        all_template_atoms,
        &folded_h,
        core_keys,
        carry_substituents,
        keep_carried_maps,
        &mut src_to_new,
    );
    let DenseAtomMaps { source_of, new_of } =
        dense_atom_maps(builder.atom_count(), input_mols, &src_to_new);

    // --- Steps 3 and 4: template bonds, then carried source bonds ---
    let colon_atoms = add_product_template_bonds(
        &mut builder,
        product_template,
        &template_idx_to_new,
        &source_of,
        input_mols,
    );
    carry_source_bonds(
        &mut builder,
        &source_of,
        &new_of,
        input_mols,
        all_template_atoms,
        global_map,
        reactant_bond_maps,
    );

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

    // Carried stereocentres keep their geometry: re-express the reactant's
    // stereo record in product atoms, since the product's bond order around
    // them can differ from the reactant's (a ring re-entered from another
    // side).
    let mut product = product;
    let sourced_atoms = restore_carried_atom_stereo(
        &mut product,
        product_template,
        &template_idx_to_new,
        input_mols,
        &source_of,
        &discoverer,
        &src_to_new,
        keep_carried_maps,
    );

    // A folded H atom (explicit-H reactant, see above) can be the only
    // `/`/`\` carrier of its double-bond end (`[H]/C(=C(/[H])C)C`). Move
    // its direction to another substituent of that end, on the other side,
    // so the E/Z survives as it does in RDKit.
    transfer_folded_h_stereo(&mut product, &folded_h, input_mols, &new_of);

    // A `/`/`\` read on an aromatic ring bond is stashed beside the bond
    // (`Molecule::bond_direction`): it carries the E/Z of an exocyclic
    // double bond (`C/N=c1\c(O)c(O)c1`). Carry it onto the product bond
    // built from that ring bond while the double bond is still there.
    carry_direction_stash(&mut product, input_mols, &src_to_new);

    // Template markers on a double bond the template creates can sit next
    // to a carried marker of a conjugated double bond, asking the shared
    // end for the same side twice; readers then drop or misread one E/Z
    // (BioTransformer BTMR0709 on an enyne gave both E and Z). Rewrite the
    // markers consistently, each double bond keeping its own geometry.
    chematic_core::ez_markers::reconcile_ez_markers(&mut product);

    // Clear any Up/Down stereo markers left on bonds that are no longer adjacent
    // to a double bond (e.g. after C=C → C=O conversion via SMIRKS).
    let mut molecule = clear_orphaned_stereo_bonds(product);

    refill_product_hydrogens(
        &mut molecule,
        &sourced_atoms,
        input_mols,
        &template_idx_to_new,
        &explicit_h_cache,
        folded_cores,
        &source_of,
        &folded_h,
    );

    restore_product_aromaticity(&mut molecule, &colon_atoms, &dearomatized);
    apply_rdkit_template_hydrogens(&mut molecule, &template_idx_to_new);
    traced_product(
        molecule,
        product_template,
        template_idx_to_new,
        &sourced_atoms,
    )
}
