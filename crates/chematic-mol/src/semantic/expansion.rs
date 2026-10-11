use std::collections::{BTreeMap, BTreeSet};

use super::*;

struct SemanticExpansion<'a> {
    model: &'a SemanticModel,
    limits: &'a SemanticExpansionLimits,
    base: &'a Molecule,
    molecule: Molecule,
    mapping: BTreeMap<SemanticId, Vec<AtomIdx>>,
}

impl<'a> SemanticExpansion<'a> {
    fn new(
        model: &'a SemanticModel,
        base: &'a Molecule,
        limits: &'a SemanticExpansionLimits,
    ) -> Result<Self, SemanticError> {
        if limits.max_atoms == 0 || limits.max_repeat_count == 0 {
            return Err(SemanticError::InvalidExpansion {
                id: "model".into(),
                reason: "expansion limits must be greater than zero".into(),
            });
        }
        model.validate()?;
        if model.atom_ids.len() != base.atom_count() {
            return Err(SemanticError::InvalidExpansion {
                id: "model".into(),
                reason: "atom_ids must match base molecule atom count".into(),
            });
        }
        if base.atom_count() > limits.max_atoms {
            return Err(SemanticError::ExpansionLimit {
                id: "model".into(),
                resource: "atoms",
                requested: base.atom_count(),
                limit: limits.max_atoms,
            });
        }
        Ok(Self {
            model,
            limits,
            base,
            molecule: base.clone(),
            mapping: model
                .atom_ids
                .iter()
                .cloned()
                .zip((0..base.atom_count()).map(|index| vec![AtomIdx(index as u32)]))
                .collect(),
        })
    }

    fn atom_index(&self, id: &str) -> Result<AtomIdx, SemanticError> {
        self.model
            .atom_ids
            .iter()
            .position(|candidate| candidate == id)
            .map(|index| AtomIdx(index as u32))
            .ok_or_else(|| SemanticError::MissingAtom(id.into()))
    }

    fn copy_fragment(
        &mut self,
        fragment: &Molecule,
        excluded: &BTreeSet<AtomIdx>,
        id: &str,
        collected: &mut Vec<AtomIdx>,
    ) -> Result<BTreeMap<AtomIdx, AtomIdx>, SemanticError> {
        let mut remap = BTreeMap::new();
        for (index, atom) in fragment.atoms() {
            if !excluded.contains(&index) {
                let added = self.molecule.add_atom(atom.clone());
                collected.push(added);
                remap.insert(index, added);
            }
        }
        for (_, bond) in fragment.bonds() {
            let (Some(&left), Some(&right)) = (remap.get(&bond.atom1), remap.get(&bond.atom2))
            else {
                continue;
            };
            self.molecule
                .add_bond(left, right, bond.order)
                .map_err(|error| SemanticError::InvalidExpansion {
                    id: id.into(),
                    reason: error.to_string(),
                })?;
        }
        Ok(remap)
    }

    fn add_single_bond(
        &mut self,
        left: AtomIdx,
        right: AtomIdx,
        id: &str,
    ) -> Result<(), SemanticError> {
        self.molecule
            .add_bond(left, right, chematic_core::BondOrder::Single)
            .map(|_| ())
            .map_err(|error| SemanticError::InvalidExpansion {
                id: id.into(),
                reason: error.to_string(),
            })
    }

    fn expand_r_group(&mut self, group: &RGroupDefinition) -> Result<(), SemanticError> {
        let Some(selected) = group.selected_alternative else {
            return Err(SemanticError::Unsupported {
                construct: group.id.clone(),
                reason: "no explicit alternative selected".into(),
            });
        };
        let fragment = chematic_smiles::parse(&group.alternatives[selected]).map_err(|error| {
            SemanticError::InvalidExpansion {
                id: group.id.clone(),
                reason: error.to_string(),
            }
        })?;
        let wildcards = fragment
            .atoms()
            .filter_map(|(index, atom)| atom.wildcard.then_some(index))
            .collect::<Vec<_>>();
        if fragment.atom_count() < 2 || wildcards.len() != group.attachment_atoms.len() {
            return Err(SemanticError::InvalidExpansion {
                id: group.id.clone(),
                reason: format!(
                    "alternative has {} wildcard attachment markers for {} source attachments",
                    wildcards.len(),
                    group.attachment_atoms.len()
                ),
            });
        }
        ensure_atom_budget(
            &self.molecule,
            fragment.atom_count() - wildcards.len(),
            self.limits,
            &group.id,
        )?;
        let base_atoms = group
            .attachment_atoms
            .iter()
            .map(|reference| self.atom_index(&reference.atom_id))
            .collect::<Result<Vec<_>, _>>()?;
        for wildcard in &wildcards {
            if fragment.neighbors(*wildcard).count() != 1 {
                return Err(SemanticError::AmbiguousAttachment(group.id.clone()));
            }
        }
        let excluded = wildcards.iter().copied().collect();
        let mut added = Vec::new();
        let remap = self.copy_fragment(&fragment, &excluded, &group.id, &mut added)?;
        for (wildcard, base_atom) in wildcards.iter().zip(base_atoms) {
            let neighbor = fragment
                .neighbors(*wildcard)
                .next()
                .expect("validated above")
                .0;
            self.add_single_bond(base_atom, remap[&neighbor], &group.id)?;
        }
        self.mapping.insert(group.id.clone(), added);
        Ok(())
    }

    fn expand_r_groups(&mut self) -> Result<(), SemanticError> {
        let mut groups = Vec::new();
        collect_r_groups(&self.model.r_groups, &mut groups);
        for group in groups {
            self.expand_r_group(group)?;
        }
        Ok(())
    }

    fn repeat_endpoints(
        unit: &PolymerRepeatUnit,
        fragment: &Molecule,
        placeholders: bool,
    ) -> Result<(AtomIdx, AtomIdx, BTreeSet<AtomIdx>), SemanticError> {
        if placeholders {
            let left_marker = AtomIdx(0);
            let right_marker = AtomIdx(fragment.atom_count().saturating_sub(1) as u32);
            let left = fragment.neighbors(left_marker).collect::<Vec<_>>();
            let right = fragment.neighbors(right_marker).collect::<Vec<_>>();
            if fragment.atom_count() < 4 || left.len() != 1 || right.len() != 1 {
                return Err(SemanticError::Unsupported {
                    construct: unit.id.clone(),
                    reason: "repeat requires exactly one neighbor at each linkage".into(),
                });
            }
            return Ok((
                left[0].0,
                right[0].0,
                [left_marker, right_marker].into_iter().collect(),
            ));
        }
        let [left, right] =
            unit.repeat_endpoint_atoms
                .ok_or_else(|| SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: "repeat endpoints are missing".into(),
                })?;
        if left == right
            || left as usize >= fragment.atom_count()
            || right as usize >= fragment.atom_count()
        {
            return Err(SemanticError::InvalidExpansion {
                id: unit.id.clone(),
                reason: "repeat endpoint atom index is invalid".into(),
            });
        }
        Ok((AtomIdx(left), AtomIdx(right), BTreeSet::new()))
    }

    fn expand_repeat_chain(
        &mut self,
        unit: &PolymerRepeatUnit,
        pattern: &str,
        placeholders: bool,
        repeats: u32,
        endpoints: [AtomIdx; 2],
        unit_atoms: &mut Vec<AtomIdx>,
    ) -> Result<(), SemanticError> {
        let mut previous_right = None;
        for repeat_index in 0..repeats {
            let fragment = chematic_smiles::parse(pattern).map_err(|error| {
                SemanticError::InvalidExpansion {
                    id: unit.id.clone(),
                    reason: error.to_string(),
                }
            })?;
            let (fragment_left, fragment_right, excluded) =
                Self::repeat_endpoints(unit, &fragment, placeholders)?;
            ensure_atom_budget(
                &self.molecule,
                fragment.atom_count() - excluded.len(),
                self.limits,
                &unit.id,
            )?;
            let remap = self.copy_fragment(&fragment, &excluded, &unit.id, unit_atoms)?;
            self.add_single_bond(
                previous_right.unwrap_or(endpoints[0]),
                remap[&fragment_left],
                &unit.id,
            )?;
            let current_right = remap[&fragment_right];
            if repeat_index + 1 == repeats {
                self.add_single_bond(endpoints[1], current_right, &unit.id)?;
            } else {
                previous_right = Some(current_right);
            }
        }
        Ok(())
    }

    fn end_group_specs(unit: &PolymerRepeatUnit) -> Vec<(String, String)> {
        if unit.end_group_definitions.is_empty() {
            unit.end_groups
                .iter()
                .cloned()
                .enumerate()
                .map(|(side, smiles)| {
                    (
                        format!(
                            "{}.end_group_{}",
                            unit.id,
                            if side == 0 { "left" } else { "right" }
                        ),
                        smiles,
                    )
                })
                .collect()
        } else {
            unit.end_group_definitions
                .iter()
                .map(|definition| (definition.id.clone(), definition.smiles.clone()))
                .collect()
        }
    }

    fn expand_end_group(
        &mut self,
        unit: &PolymerRepeatUnit,
        side: usize,
        mapping_id: &str,
        pattern: &str,
        endpoint: AtomIdx,
        unit_atoms: &mut Vec<AtomIdx>,
    ) -> Result<(), SemanticError> {
        let fragment =
            chematic_smiles::parse(pattern).map_err(|error| SemanticError::InvalidExpansion {
                id: unit.id.clone(),
                reason: format!("invalid end-group {side}: {error}"),
            })?;
        let wildcards = fragment
            .atoms()
            .filter_map(|(index, atom)| atom.wildcard.then_some(index))
            .collect::<Vec<_>>();
        if fragment.atom_count() < 2 || wildcards.len() != 1 {
            return Err(SemanticError::InvalidExpansion {
                id: unit.id.clone(),
                reason: format!("end-group {side} requires exactly one [*] marker"),
            });
        }
        ensure_atom_budget(
            &self.molecule,
            fragment.atom_count() - 1,
            self.limits,
            &unit.id,
        )?;
        let wildcard = wildcards[0];
        let neighbors = fragment.neighbors(wildcard).collect::<Vec<_>>();
        if neighbors.len() != 1 {
            return Err(SemanticError::AmbiguousAttachment(unit.id.clone()));
        }
        let excluded = [wildcard].into_iter().collect();
        let mut added = Vec::new();
        let remap = self.copy_fragment(&fragment, &excluded, &unit.id, &mut added)?;
        self.add_single_bond(endpoint, remap[&neighbors[0].0], &unit.id)?;
        unit_atoms.extend(&added);
        self.mapping.insert(mapping_id.into(), added);
        Ok(())
    }

    fn expand_polymer_unit(&mut self, unit: &PolymerRepeatUnit) -> Result<(), SemanticError> {
        let pattern = unit
            .repeat_smiles
            .as_deref()
            .ok_or_else(|| SemanticError::Unsupported {
                construct: unit.id.clone(),
                reason: "repeat_smiles is not provided".into(),
            })?;
        let placeholders = pattern.trim().starts_with("[*]")
            && pattern.trim().ends_with("[*]")
            && pattern.matches("[*]").count() == 2;
        if !placeholders && unit.repeat_endpoint_atoms.is_none() {
            return Err(SemanticError::InvalidExpansion {
                id: unit.id.clone(),
                reason: "repeat requires two [*] markers or explicit endpoint atoms".into(),
            });
        }
        let repeats = unit
            .repeat_count
            .ok_or_else(|| SemanticError::Unsupported {
                construct: unit.id.clone(),
                reason: "repeat count must be explicit".into(),
            })?;
        if repeats > self.limits.max_repeat_count {
            return Err(SemanticError::ExpansionLimit {
                id: unit.id.clone(),
                resource: "repeats",
                requested: repeats as usize,
                limit: self.limits.max_repeat_count as usize,
            });
        }
        let endpoints = [
            self.atom_index(&unit.attachment_atoms[0].atom_id)?,
            self.atom_index(&unit.attachment_atoms[1].atom_id)?,
        ];
        let mut unit_atoms = Vec::new();
        self.expand_repeat_chain(
            unit,
            pattern,
            placeholders,
            repeats,
            endpoints,
            &mut unit_atoms,
        )?;
        for (side, (mapping_id, end_group)) in Self::end_group_specs(unit).iter().enumerate() {
            self.expand_end_group(
                unit,
                side,
                mapping_id,
                end_group,
                endpoints[usize::from(side != 0)],
                &mut unit_atoms,
            )?;
        }
        self.mapping.insert(unit.id.clone(), unit_atoms);
        Ok(())
    }

    fn expand_polymers(&mut self) -> Result<(), SemanticError> {
        for unit in &self.model.polymer_units {
            self.expand_polymer_unit(unit)?;
        }
        Ok(())
    }

    fn expand_s_groups(&mut self) -> Result<(), SemanticError> {
        for group in &self.model.s_groups {
            let members = group
                .member_atoms
                .iter()
                .map(|atom| self.atom_index(&atom.atom_id))
                .collect::<Result<Vec<_>, _>>()?;
            self.mapping.insert(group.id.clone(), members);
            if let Some(linkage) = &group.linkage {
                let endpoints = [&linkage.left_atom, &linkage.right_atom]
                    .into_iter()
                    .map(|atom| self.atom_index(&atom.atom_id))
                    .collect::<Result<Vec<_>, _>>()?;
                self.mapping.insert(linkage.id.clone(), endpoints);
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<ExpandedSemantic, SemanticError> {
        self.expand_r_groups()?;
        self.expand_polymers()?;
        self.expand_s_groups()?;
        Ok(ExpandedSemantic {
            molecule: self.molecule,
            source_to_expanded: self.mapping,
            base_molecule: self.base.clone(),
        })
    }
}

pub(super) fn expand(
    model: &SemanticModel,
    base: &Molecule,
    limits: &SemanticExpansionLimits,
) -> Result<ExpandedSemantic, SemanticError> {
    SemanticExpansion::new(model, base, limits)?.finish()
}
