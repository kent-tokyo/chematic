use super::*;

type Config = Option<(bool, Vec<u32>)>;

type AtomSource = Option<(usize, AtomIdx)>;
type StereoNodeConfig = Option<(bool, Vec<Node>)>;

/// One earlier traversal's bonds to a mapped atom: the traversal root's own
/// bond (if any), then carried atoms' bonds in an unknown order.
type BondGroup = (Option<Node>, Vec<Node>);

/// Bonds RDKit adds to one mapped atom while traversing the reactant from
/// mapped atoms that come earlier in the match: one group per earlier
/// atom, the atom's own bond (if any) first, then bonds from carried atoms
/// in an order that depends on the traversal.
#[derive(Default, Clone)]
struct Earlier {
    groups: Vec<BondGroup>,
    visited: Vec<bool>,
}

/// One reactant's atoms as the product sees them.
struct ReactantView {
    /// Reactant atom → product-template atom, for atoms mapped into it.
    mapped: FxHashMap<AtomIdx, usize>,
    /// Matched atoms not in this product.
    skipped: FxHashSet<AtomIdx>,
    /// Mapped atoms in match (template-atom) order.
    order: Vec<AtomIdx>,
}

type ReactantViews = (Vec<AtomSource>, Vec<ReactantView>, bool);

struct StereoApplication<'a> {
    info: &'a ReactionStereo,
    reactants: &'a [&'a Molecule],
    product_index: usize,
    template_idx_to_new: &'a [Option<AtomIdx>],
    src_to_new: &'a FxHashMap<(usize, AtomIdx), AtomIdx>,
}

impl StereoApplication<'_> {
    fn template(&self) -> &[TemplateAtom] {
        &self.info.products[self.product_index]
    }

    fn build_reactant_views(&self, per_reactant: &[FxHashMap<usize, AtomIdx>]) -> ReactantViews {
        let template = self.template();
        let mut product_by_map: FxHashMap<u16, usize> = FxHashMap::default();
        let mut duplicate = false;
        for (i, atom) in template.iter().enumerate() {
            if let Some(map) = atom.map {
                duplicate |= product_by_map.insert(map, i).is_some();
            }
        }

        let mut source = vec![None; template.len()];
        let mut views = Vec::with_capacity(self.reactants.len());
        for (reactant_index, matched) in per_reactant.iter().enumerate() {
            let mut pairs: Vec<_> = matched.iter().map(|(&t, &r)| (t, r)).collect();
            pairs.sort_unstable_by_key(|&(template_index, _)| template_index);
            let mut view = ReactantView {
                mapped: FxHashMap::default(),
                skipped: FxHashSet::default(),
                order: Vec::new(),
            };
            for (template_index, reactant_atom) in pairs {
                let target = self.info.reactants[reactant_index]
                    .get(template_index)
                    .and_then(|atom| atom.map)
                    .and_then(|map| product_by_map.get(&map).copied());
                if let Some(product_atom) = target {
                    duplicate |= source[product_atom]
                        .replace((reactant_index, reactant_atom))
                        .is_some();
                    view.mapped.insert(reactant_atom, product_atom);
                    view.order.push(reactant_atom);
                } else {
                    view.skipped.insert(reactant_atom);
                }
            }
            views.push(view);
        }
        (source, views, duplicate)
    }

    fn preserve_carried_stereo(&self, product: &mut Molecule) {
        for (&(reactant_index, reactant_atom), &product_atom) in self.src_to_new {
            if self.template_idx_to_new.contains(&Some(product_atom)) {
                continue;
            }
            let reactant = self.reactants[reactant_index];
            let atom = reactant.atom(reactant_atom);
            if !atom.chirality.is_tetrahedral() {
                continue;
            }
            let remapped: Option<Vec<u32>> = reactant
                .stereo_neighbor_order(reactant_atom)
                .and_then(|order| {
                    order
                        .iter()
                        .map(|&neighbor| {
                            if neighbor == STEREO_H_SENTINEL {
                                Some(neighbor)
                            } else {
                                self.src_to_new
                                    .get(&(reactant_index, AtomIdx(neighbor)))
                                    .map(|atom| atom.0)
                            }
                        })
                        .collect()
                });
            let degree_kept = reactant.degree(reactant_atom) == product.degree(product_atom);
            match remapped {
                Some(order) if degree_kept => {
                    product.set_chirality(product_atom, atom.chirality);
                    product.set_stereo_neighbor_order(product_atom, order);
                }
                _ => product.set_chirality(product_atom, Chirality::None),
            }
        }
    }

    fn involved_atoms(&self, source: &[Option<(usize, AtomIdx)>]) -> Vec<usize> {
        self.template()
            .iter()
            .enumerate()
            .filter_map(|(i, atom)| {
                let chiral_source = source[i]
                    .is_some_and(|(k, r)| self.reactants[k].atom(r).chirality.is_tetrahedral());
                (chiral_source || atom.tag.is_some()).then_some(i)
            })
            .collect()
    }

    fn clear_template_tetrahedral(&self, product: &mut Molecule) {
        for &atom in self.template_idx_to_new.iter().flatten() {
            if product.atom(atom).chirality.is_tetrahedral() {
                product.set_chirality(atom, Chirality::None);
            }
        }
    }

    fn collect_earlier_bonds(
        &self,
        views: &[ReactantView],
        involved: &[usize],
    ) -> FxHashMap<(usize, AtomIdx), Earlier> {
        let template = self.template();
        let template_bond = |i: usize, j: usize| template[i].order.contains(&j);
        let mut earlier = FxHashMap::default();
        for (reactant_index, view) in views.iter().enumerate() {
            if !view
                .order
                .iter()
                .any(|atom| involved.contains(&view.mapped[atom]))
            {
                continue;
            }
            let molecule = self.reactants[reactant_index];
            let mut visited = vec![false; molecule.atom_count()];
            let mut added = vec![false; molecule.atom_count()];
            let mut pending: FxHashMap<AtomIdx, Vec<BondGroup>> = FxHashMap::default();
            for &mapped_atom in &view.order {
                let template_atom = view.mapped[&mapped_atom];
                earlier.insert(
                    (reactant_index, mapped_atom),
                    Earlier {
                        groups: pending.remove(&mapped_atom).unwrap_or_default(),
                        visited: visited.clone(),
                    },
                );
                let mut round: FxHashMap<AtomIdx, BondGroup> = FxHashMap::default();
                let mut queue = std::collections::VecDeque::from([mapped_atom]);
                while let Some(local_atom) = queue.pop_front() {
                    visited[local_atom.0 as usize] = true;
                    for (neighbor, _) in molecule.neighbors(local_atom) {
                        let neighbor_index = neighbor.0 as usize;
                        if visited[neighbor_index] || view.skipped.contains(&neighbor) {
                            continue;
                        }
                        if let Some(&neighbor_template) = view.mapped.get(&neighbor) {
                            let group = round.entry(neighbor).or_default();
                            if local_atom == mapped_atom {
                                let left_map = template[template_atom].map.unwrap_or(0);
                                let right_map = template[neighbor_template].map.unwrap_or(0);
                                if !template_bond(template_atom, neighbor_template)
                                    && !self.info.reactant_bond_maps.contains(&(
                                        left_map.min(right_map),
                                        left_map.max(right_map),
                                    ))
                                {
                                    group.0 = Some(Node::Template(template_atom));
                                }
                            } else {
                                group.1.push(Node::Carried(reactant_index, local_atom));
                            }
                        } else if !added[neighbor_index] {
                            added[neighbor_index] = true;
                            queue.push_back(neighbor);
                        }
                    }
                }
                for (neighbor, group) in round {
                    if group.0.is_some() || !group.1.is_empty() {
                        pending.entry(neighbor).or_default().push(group);
                    }
                }
            }
        }
        earlier
    }
}

impl StereoApplication<'_> {
    fn atom_stereo_results(
        &self,
        template_atom: usize,
        source: &[Option<(usize, AtomIdx)>],
        views: &[ReactantView],
        earlier: &FxHashMap<(usize, AtomIdx), Earlier>,
    ) -> Result<Vec<StereoNodeConfig>, StereoUnsupported> {
        let template = self.template();
        let flag = self.info.flags[self.product_index][template_atom];
        let template_part: Vec<Node> = template[template_atom]
            .order
            .iter()
            .map(|&neighbor| Node::Template(neighbor))
            .collect();
        let Some((reactant_index, reactant_atom)) = source[template_atom] else {
            return Ok(vec![
                (flag == Some(4))
                    .then_some(template[template_atom].tag)
                    .flatten()
                    .map(|tag| (tag, template_part)),
            ]);
        };

        let molecule = self.reactants[reactant_index];
        let chiral = molecule.atom(reactant_atom).chirality.is_tetrahedral();
        let created = flag == Some(4) && template[template_atom].tag.is_some();
        if !chiral && !created {
            return Ok(vec![None]);
        }

        let view = &views[reactant_index];
        let earlier = earlier
            .get(&(reactant_index, reactant_atom))
            .cloned()
            .unwrap_or_default();
        let (prefix, rings) = reactant_bond_order(molecule, reactant_atom);
        let ring_orders =
            permutations(&rings, MAX_ORDER_OPTIONS).ok_or(StereoUnsupported::BondOrder)?;
        let mut group_options: Vec<Vec<Node>> = vec![Vec::new()];
        for (fixed, carried) in &earlier.groups {
            let permutations =
                permutations(carried, MAX_ORDER_OPTIONS).ok_or(StereoUnsupported::BondOrder)?;
            let mut next = Vec::new();
            for base in &group_options {
                for permutation in &permutations {
                    let mut option = base.clone();
                    option.extend(fixed.iter().copied());
                    option.extend(permutation.iter().copied());
                    next.push(option);
                }
            }
            if next.len() * ring_orders.len() > MAX_ORDER_OPTIONS {
                return Err(StereoUnsupported::BondOrder);
            }
            group_options = next;
        }

        let mut results = Vec::new();
        for ring_order in &ring_orders {
            let mut reactant_order = prefix.clone();
            reactant_order.extend(ring_order.iter().copied());
            let own: Vec<Node> = reactant_order
                .iter()
                .filter(|atom| {
                    !earlier
                        .visited
                        .get(atom.0 as usize)
                        .copied()
                        .unwrap_or(false)
                        && !view.skipped.contains(atom)
                })
                .filter_map(|&atom| match view.mapped.get(&atom) {
                    Some(&neighbor_template) => {
                        let left_map = template[template_atom].map.unwrap_or(0);
                        let right_map = template[neighbor_template].map.unwrap_or(0);
                        (!template[template_atom].order.contains(&neighbor_template)
                            && !self
                                .info
                                .reactant_bond_maps
                                .contains(&(left_map.min(right_map), left_map.max(right_map))))
                        .then_some(Node::Template(neighbor_template))
                    }
                    None => Some(Node::Carried(reactant_index, atom)),
                })
                .collect();
            for groups in &group_options {
                let mut list = template_part.clone();
                list.extend(groups.iter().copied());
                list.extend(own.iter().copied());
                results.push(rdkit_tag(
                    self.info,
                    self.product_index,
                    template_atom,
                    flag,
                    molecule,
                    reactant_index,
                    reactant_atom,
                    &reactant_order,
                    &list,
                    source,
                )?);
            }
        }
        Ok(results)
    }

    fn decide_atom(
        &self,
        template_atom: usize,
        source: &[Option<(usize, AtomIdx)>],
        views: &[ReactantView],
        earlier: &FxHashMap<(usize, AtomIdx), Earlier>,
        product: &mut Molecule,
    ) -> Result<Option<(AtomIdx, Vec<Config>)>, StereoUnsupported> {
        let Some(product_atom) = self
            .template_idx_to_new
            .get(template_atom)
            .copied()
            .flatten()
        else {
            return Ok(None);
        };
        if product.atom(product_atom).aromatic {
            if product.atom(product_atom).chirality != Chirality::None {
                product.set_chirality(product_atom, Chirality::None);
            }
            return Ok(None);
        }

        let mut outcomes: Vec<(Option<bool>, Config)> = Vec::new();
        for result in self.atom_stereo_results(template_atom, source, views, earlier)? {
            let converted = match result {
                None => None,
                Some((clockwise, list)) => Some(to_product_config(
                    clockwise,
                    &list,
                    product_atom,
                    product,
                    self.template_idx_to_new,
                    self.src_to_new,
                )?),
            };
            let key = converted.as_ref().map(|(clockwise, order)| {
                let mut sorted = order.clone();
                sorted.sort_unstable();
                *clockwise == permutation_parity(order, &sorted).unwrap_or(true)
            });
            if !outcomes.iter().any(|(existing, _)| *existing == key) {
                outcomes.push((key, converted));
            }
        }
        Ok(Some((
            product_atom,
            outcomes.into_iter().map(|(_, config)| config).collect(),
        )))
    }

    fn collect_decisions(
        &self,
        source: &[Option<(usize, AtomIdx)>],
        views: &[ReactantView],
        earlier: &FxHashMap<(usize, AtomIdx), Earlier>,
        product: &mut Molecule,
    ) -> Result<Vec<(AtomIdx, Vec<Config>)>, StereoUnsupported> {
        (0..self.template().len())
            .filter_map(|template_atom| {
                self.decide_atom(template_atom, source, views, earlier, product)
                    .transpose()
            })
            .collect()
    }
}

fn apply_and_validate_decisions(
    product: &mut Molecule,
    decided: &[(AtomIdx, Vec<Config>)],
) -> Result<(), StereoUnsupported> {
    let ambiguous: Vec<usize> = decided
        .iter()
        .enumerate()
        .filter_map(|(i, (_, outcomes))| (outcomes.len() > 1).then_some(i))
        .collect();
    for (atom, outcomes) in decided {
        set_config(product, *atom, &outcomes[0]);
    }
    if ambiguous.is_empty() {
        return Ok(());
    }

    let combinations = ambiguous
        .iter()
        .try_fold(1usize, |count, &i| count.checked_mul(decided[i].1.len()))
        .filter(|&count| count <= 64)
        .ok_or(StereoUnsupported::BondOrder)?;
    let reference = chematic_smiles::canonical_smiles(product);
    for combination in 1..combinations {
        let mut variant = product.clone();
        let mut rest = combination;
        for &decision in &ambiguous {
            let (atom, outcomes) = &decided[decision];
            set_config(&mut variant, *atom, &outcomes[rest % outcomes.len()]);
            rest /= outcomes.len();
        }
        if chematic_smiles::canonical_smiles(&variant) != reference {
            return Err(StereoUnsupported::BondOrder);
        }
    }
    Ok(())
}

/// Recompute the chiral tags of the product-template atoms of product `p`
/// as RDKit would, overriding the native tags on those atoms. Carried atoms
/// keep the native (geometric) stereo, which RDKit also preserves.
pub(crate) fn apply_product_stereo(
    info: &ReactionStereo,
    reactants: &[&Molecule],
    per_reactant: &[FxHashMap<usize, AtomIdx>],
    p: usize,
    template_idx_to_new: &[Option<AtomIdx>],
    src_to_new: &FxHashMap<(usize, AtomIdx), AtomIdx>,
    product: &mut Molecule,
) -> Result<(), StereoUnsupported> {
    let application = StereoApplication {
        info,
        reactants,
        product_index: p,
        template_idx_to_new,
        src_to_new,
    };
    let (source, views, duplicate) = application.build_reactant_views(per_reactant);
    application.preserve_carried_stereo(product);
    let involved = application.involved_atoms(&source);
    if involved.is_empty() {
        application.clear_template_tetrahedral(product);
        return Ok(());
    }
    if duplicate {
        return Err(StereoUnsupported::BondOrder);
    }

    let earlier = application.collect_earlier_bonds(&views, &involved);

    let decided = application.collect_decisions(&source, &views, &earlier, product)?;
    apply_and_validate_decisions(product, &decided)
}

fn set_config(product: &mut Molecule, idx: AtomIdx, config: &Config) {
    match config {
        Some((cw, order)) => {
            product.set_chirality(
                idx,
                if *cw {
                    Chirality::Clockwise
                } else {
                    Chirality::CounterClockwise
                },
            );
            product.set_stereo_neighbor_order(idx, order.clone());
        }
        None => {
            if product.atom(idx).chirality.is_tetrahedral() {
                product.set_chirality(idx, Chirality::None);
            }
        }
    }
}
