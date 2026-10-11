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
    let inputs = ProductBuildInputs {
        product_template,
        global_map,
        reactant_template_atoms,
        reactant_bond_maps,
        input_mols,
        all_template_atoms,
        carry_substituents,
        specs: (specs.len() == product_template.atom_count()).then_some(specs),
        keep_carried_maps,
    };
    let assembly = add_template_atoms(&inputs);
    finish_product(&inputs, assembly)
}

struct ProductBuildInputs<'a> {
    product_template: &'a Molecule,
    global_map: &'a FxHashMap<u16, (usize, AtomIdx)>,
    reactant_template_atoms: &'a FxHashMap<u16, ReactantMapAtom>,
    reactant_bond_maps: &'a FxHashSet<(u16, u16)>,
    input_mols: &'a [&'a Molecule],
    all_template_atoms: &'a FxHashSet<(usize, AtomIdx)>,
    carry_substituents: bool,
    specs: Option<&'a [ProductAtomSpec]>,
    keep_carried_maps: bool,
}

impl ProductBuildInputs<'_> {
    fn spec_of(&self, index: usize) -> ProductAtomSpec {
        let atom = self.product_template.atom(AtomIdx(index as u32));
        self.specs.map_or(
            ProductAtomSpec {
                charge: true,
                hcount: atom.hydrogen_count.is_some_and(|count| count > 0),
            },
            |specs| specs[index],
        )
    }

    fn template_h(&self, index: usize) -> Option<u8> {
        let atom = self.product_template.atom(AtomIdx(index as u32));
        self.spec_of(index)
            .hcount
            .then(|| atom.hydrogen_count.unwrap_or(0))
    }
}

struct TemplateAssembly {
    builder: MoleculeBuilder,
    template_idx_to_new: Vec<Option<AtomIdx>>,
    src_to_new: FxHashMap<(usize, AtomIdx), AtomIdx>,
    core_keys: Vec<(usize, AtomIdx)>,
    folded_h: FxHashSet<(usize, AtomIdx)>,
    dearomatized: Vec<AtomIdx>,
    folded_cores: Vec<AtomIdx>,
    explicit_h_cache: Vec<std::cell::OnceCell<bool>>,
}

fn add_template_atoms(inputs: &ProductBuildInputs<'_>) -> TemplateAssembly {
    let atom_capacity = inputs
        .input_mols
        .iter()
        .map(|m| m.atom_count())
        .sum::<usize>()
        + inputs.product_template.atom_count();
    let bond_capacity = inputs
        .input_mols
        .iter()
        .map(|m| m.bond_count())
        .sum::<usize>()
        + inputs.product_template.bond_count();
    let mut state = TemplateAssembly {
        builder: MoleculeBuilder::with_capacity(atom_capacity, bond_capacity),
        template_idx_to_new: vec![None; inputs.product_template.atom_count()],
        src_to_new: FxHashMap::with_capacity_and_hasher(
            inputs.input_mols.iter().map(|m| m.atom_count()).sum(),
            Default::default(),
        ),
        core_keys: mapped_core_keys(inputs.product_template, inputs.global_map),
        folded_h: FxHashSet::default(),
        dearomatized: Vec::new(),
        folded_cores: Vec::new(),
        explicit_h_cache: inputs
            .input_mols
            .iter()
            .map(|_| std::cell::OnceCell::new())
            .collect(),
    };
    for index in 0..inputs.product_template.atom_count() {
        let new_index = add_template_atom(inputs, &mut state, index);
        state.template_idx_to_new[index] = Some(new_index);
    }
    state
}

fn mapped_core_keys(
    template: &Molecule,
    global_map: &FxHashMap<u16, (usize, AtomIdx)>,
) -> Vec<(usize, AtomIdx)> {
    let mut keys: Vec<_> = template
        .atoms()
        .filter_map(|(_, atom)| atom.atom_map.and_then(|map| global_map.get(&map).copied()))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    keys
}

fn add_template_atom(
    inputs: &ProductBuildInputs<'_>,
    state: &mut TemplateAssembly,
    index: usize,
) -> AtomIdx {
    let template_atom = inputs.product_template.atom(AtomIdx(index as u32));
    let source = template_atom
        .atom_map
        .and_then(|map| inputs.global_map.get(&map).map(|&source| (map, source)));
    match source {
        Some((map, (molecule, atom))) => {
            add_mapped_template_atom(inputs, state, index, map, molecule, atom)
        }
        None => add_new_template_atom(inputs, &mut state.builder, index),
    }
}

fn add_new_template_atom(
    inputs: &ProductBuildInputs<'_>,
    builder: &mut MoleculeBuilder,
    index: usize,
) -> AtomIdx {
    let mut atom = inputs.product_template.atom(AtomIdx(index as u32)).clone();
    atom.atom_map = None;
    // RDKit leaves H0 on a product-only atom unspecified.
    atom.hydrogen_count = inputs.template_h(index).filter(|&count| count != 0);
    builder.add_atom(atom)
}

fn add_mapped_template_atom(
    inputs: &ProductBuildInputs<'_>,
    state: &mut TemplateAssembly,
    index: usize,
    atom_map: u16,
    molecule_index: usize,
    source_index: AtomIdx,
) -> AtomIdx {
    let source = inputs.input_mols[molecule_index].atom(source_index);
    let template = inputs.product_template.atom(AtomIdx(index as u32));
    let reactant_template = inputs.reactant_template_atoms.get(&atom_map);
    let (mut atom, degree_unchanged, element_unchanged) =
        mapped_atom_electronics(inputs, index, source, template, reactant_template);
    let explicit_h = foldable_hydrogens(inputs, molecule_index, source_index, source, &atom);
    let refill = should_refill_mapped_hydrogens(
        inputs,
        state,
        index,
        molecule_index,
        source,
        &atom,
        &explicit_h,
    );
    apply_folded_hydrogen_policy(
        inputs,
        index,
        source,
        &mut atom,
        degree_unchanged,
        element_unchanged,
        &explicit_h,
        refill,
    );

    let dearomatize = source.aromatic && !atom.aromatic;
    let new_index = state.builder.add_atom(atom);
    if dearomatize {
        state.dearomatized.push(new_index);
    }
    if refill || !explicit_h.is_empty() {
        state.folded_cores.push(new_index);
    }
    state
        .folded_h
        .extend(explicit_h.into_iter().map(|atom| (molecule_index, atom)));
    state
        .src_to_new
        .insert((molecule_index, source_index), new_index);
    new_index
}

fn mapped_atom_electronics(
    inputs: &ProductBuildInputs<'_>,
    index: usize,
    source: &chematic_core::Atom,
    template: &chematic_core::Atom,
    reactant_template: Option<&ReactantMapAtom>,
) -> (chematic_core::Atom, bool, bool) {
    let mut atom = source.clone();
    if !template.wildcard {
        let changes = reactant_template
            .is_none_or(|value| value.wildcard || value.element != template.element);
        if changes && atom.element != template.element {
            atom.element = template.element;
            atom.isotope = template.isotope;
        }
        atom.aromatic = template.aromatic;
    }
    if inputs.spec_of(index).charge || atom.element != source.element {
        atom.charge = template.charge;
    }
    if template.isotope.is_some() {
        atom.isotope = template.isotope;
    }
    let degree_unchanged = reactant_template
        .is_some_and(|value| value.degree == inputs.product_template.degree(AtomIdx(index as u32)));
    let element_unchanged = atom.element == source.element;
    atom.hydrogen_count = match inputs.template_h(index) {
        Some(count) => Some(count),
        None if degree_unchanged && element_unchanged => source.hydrogen_count,
        None => None,
    };
    atom.atom_map = None;
    (atom, degree_unchanged, element_unchanged)
}

fn plain_hydrogen(mol: &Molecule, index: AtomIdx) -> bool {
    let atom = mol.atom(index);
    atom.element.atomic_number() == 1
        && !atom.wildcard
        && atom.isotope.is_none()
        && atom.charge == 0
        && atom.atom_map.is_none()
        && mol.degree(index) == 1
}

fn foldable_hydrogens(
    inputs: &ProductBuildInputs<'_>,
    molecule_index: usize,
    source_index: AtomIdx,
    source: &chematic_core::Atom,
    product: &chematic_core::Atom,
) -> Vec<AtomIdx> {
    if !inputs.carry_substituents
        || source.chirality != Chirality::None
        || !source.element.is_organic_subset()
        || !product.element.is_organic_subset()
    {
        return Vec::new();
    }
    inputs.input_mols[molecule_index]
        .neighbors(source_index)
        .map(|(neighbor, _)| neighbor)
        .filter(|&neighbor| {
            plain_hydrogen(inputs.input_mols[molecule_index], neighbor)
                && !inputs
                    .all_template_atoms
                    .contains(&(molecule_index, neighbor))
        })
        .collect()
}

fn should_refill_mapped_hydrogens(
    inputs: &ProductBuildInputs<'_>,
    state: &TemplateAssembly,
    index: usize,
    molecule_index: usize,
    source: &chematic_core::Atom,
    product: &chematic_core::Atom,
    explicit_h: &[AtomIdx],
) -> bool {
    explicit_h.is_empty()
        && inputs.carry_substituents
        && source.chirality == Chirality::None
        && source.hydrogen_count == Some(0)
        && source.charge == 0
        && (!source.aromatic || source.element == chematic_core::Element::C)
        && !product.aromatic
        && source.element.is_organic_subset()
        && product.element.is_organic_subset()
        && inputs.template_h(index).is_none()
        && reactant_has_explicit_h(inputs.input_mols, &state.explicit_h_cache, molecule_index)
}

#[allow(clippy::too_many_arguments)]
fn apply_folded_hydrogen_policy(
    inputs: &ProductBuildInputs<'_>,
    index: usize,
    source: &chematic_core::Atom,
    product: &mut chematic_core::Atom,
    degree_unchanged: bool,
    element_unchanged: bool,
    explicit_h: &[AtomIdx],
    refill: bool,
) {
    if explicit_h.is_empty() && !refill {
        return;
    }
    product.hydrogen_count = inputs.template_h(index);
    if product.hydrogen_count.is_none()
        && product.aromatic
        && source.aromatic
        && degree_unchanged
        && element_unchanged
        && product.charge == source.charge
    {
        product.hydrogen_count = Some(explicit_h.len() as u8);
    }
}

fn finish_product(
    inputs: &ProductBuildInputs<'_>,
    state: TemplateAssembly,
) -> (TracedProduct, Vec<Option<AtomIdx>>) {
    let TemplateAssembly {
        mut builder,
        template_idx_to_new,
        mut src_to_new,
        core_keys,
        folded_h,
        dearomatized,
        folded_cores,
        explicit_h_cache,
    } = state;
    let discoverer = carry_substituent_atoms(
        &mut builder,
        inputs.input_mols,
        inputs.all_template_atoms,
        &folded_h,
        core_keys,
        inputs.carry_substituents,
        inputs.keep_carried_maps,
        &mut src_to_new,
    );
    let DenseAtomMaps { source_of, new_of } =
        dense_atom_maps(builder.atom_count(), inputs.input_mols, &src_to_new);
    let colon_atoms = add_product_template_bonds(
        &mut builder,
        inputs.product_template,
        &template_idx_to_new,
        &source_of,
        inputs.input_mols,
    );
    carry_source_bonds(
        &mut builder,
        &source_of,
        &new_of,
        inputs.input_mols,
        inputs.all_template_atoms,
        inputs.global_map,
        inputs.reactant_bond_maps,
    );
    let mut product = correct_product_stereo(
        builder.build(),
        inputs.product_template,
        &template_idx_to_new,
        inputs.global_map,
        inputs.all_template_atoms,
        inputs.input_mols,
        &src_to_new,
    );
    let sourced_atoms = restore_carried_atom_stereo(
        &mut product,
        inputs.product_template,
        &template_idx_to_new,
        inputs.input_mols,
        &source_of,
        &discoverer,
        &src_to_new,
        inputs.keep_carried_maps,
    );
    transfer_folded_h_stereo(&mut product, &folded_h, inputs.input_mols, &new_of);
    carry_direction_stash(&mut product, inputs.input_mols, &src_to_new);
    chematic_core::ez_markers::reconcile_ez_markers(&mut product);

    let mut molecule = clear_orphaned_stereo_bonds(product);
    refill_product_hydrogens(
        &mut molecule,
        &sourced_atoms,
        inputs.input_mols,
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
        inputs.product_template,
        template_idx_to_new,
        &sourced_atoms,
    )
}
