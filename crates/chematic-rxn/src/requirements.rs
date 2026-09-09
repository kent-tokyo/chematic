//! Conservative prefilter requirements derived from compiled reaction queries.

use std::collections::BTreeMap;

use chematic_core::{BondOrder, Molecule};
use chematic_smarts::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

/// The bond kinds for which a query can provide a safe lower-bound count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReactionBondKind {
    Single,
    Double,
    Triple,
    Aromatic,
}

/// A conservative lower bound for bonds between two element types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactionBondLowerBound {
    pub first_atomic_number: u8,
    pub second_atomic_number: u8,
    pub bond: ReactionBondKind,
    pub count: u16,
}

/// Cheap, fail-open requirements useful for filtering reaction candidates.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReactionRequirements {
    pub min_atom_count: usize,
    pub min_bond_count: usize,
    pub element_lower_bounds: Vec<(u8, u16)>,
    pub aromatic_element_lower_bounds: Vec<(u8, u16)>,
    pub aliphatic_element_lower_bounds: Vec<(u8, u16)>,
    pub bond_lower_bounds: Vec<ReactionBondLowerBound>,
}

impl ReactionRequirements {
    pub(crate) fn from_queries(queries: &[QueryMolecule]) -> Self {
        let mut result = Self {
            min_atom_count: queries.iter().map(QueryMolecule::atom_count).sum(),
            min_bond_count: queries.iter().map(|q| q.bonds.len()).sum(),
            ..Self::default()
        };
        let mut elements = BTreeMap::new();
        let mut aromatic = BTreeMap::new();
        let mut aliphatic = BTreeMap::new();
        let mut bonds = BTreeMap::new();
        for query in queries {
            let atom_facts: Vec<_> = query
                .atoms
                .iter()
                .map(|a| {
                    (
                        positive_atomic_number(&a.query),
                        positive_aromatic(&a.query),
                    )
                })
                .collect();
            for (number, is_aromatic) in &atom_facts {
                if let Some(number) = number {
                    *elements.entry(*number).or_insert(0u16) += 1;
                    match is_aromatic {
                        Some(true) => *aromatic.entry(*number).or_insert(0u16) += 1,
                        Some(false) => *aliphatic.entry(*number).or_insert(0u16) += 1,
                        None => {}
                    }
                }
            }
            for bond in &query.bonds {
                let (Some(first), Some(second), Some(kind)) = (
                    atom_facts[bond.atom1].0,
                    atom_facts[bond.atom2].0,
                    simple_bond(&bond.query),
                ) else {
                    continue;
                };
                let (first, second) = if first <= second {
                    (first, second)
                } else {
                    (second, first)
                };
                *bonds.entry((first, second, kind)).or_insert(0u16) += 1;
            }
        }
        result.element_lower_bounds = elements.into_iter().collect();
        result.aromatic_element_lower_bounds = aromatic.into_iter().collect();
        result.aliphatic_element_lower_bounds = aliphatic.into_iter().collect();
        result.bond_lower_bounds = bonds
            .into_iter()
            .map(|((first, second, bond), count)| ReactionBondLowerBound {
                first_atomic_number: first,
                second_atomic_number: second,
                bond,
                count,
            })
            .collect();
        result
    }

    /// Return whether supplied reactants satisfy these conservative bounds.
    /// Missing requirements never cause rejection.
    pub fn could_match(&self, reactants: &[&Molecule]) -> bool {
        let atoms: usize = reactants.iter().map(|m| m.atom_count()).sum();
        let bonds_count: usize = reactants.iter().map(|m| m.bond_count()).sum();
        if atoms < self.min_atom_count || bonds_count < self.min_bond_count {
            return false;
        }
        let mut elements = BTreeMap::new();
        let mut aromatic = BTreeMap::new();
        let mut aliphatic = BTreeMap::new();
        let mut bonds = BTreeMap::new();
        for molecule in reactants {
            for (_, atom) in molecule.atoms() {
                let number = atom.element.atomic_number();
                *elements.entry(number).or_insert(0u16) += 1;
                if atom.aromatic {
                    *aromatic.entry(number).or_insert(0u16) += 1;
                } else {
                    *aliphatic.entry(number).or_insert(0u16) += 1;
                }
            }
            for (_, bond) in molecule.bonds() {
                let Some(kind) = target_bond_kind(bond.order) else {
                    continue;
                };
                let first = molecule.atom(bond.atom1).element.atomic_number();
                let second = molecule.atom(bond.atom2).element.atomic_number();
                let (first, second) = if first <= second {
                    (first, second)
                } else {
                    (second, first)
                };
                *bonds.entry((first, second, kind)).or_insert(0u16) += 1;
            }
        }
        self.element_lower_bounds
            .iter()
            .all(|(n, min)| elements.get(n).copied().unwrap_or(0) >= *min)
            && self
                .aromatic_element_lower_bounds
                .iter()
                .all(|(n, min)| aromatic.get(n).copied().unwrap_or(0) >= *min)
            && self
                .aliphatic_element_lower_bounds
                .iter()
                .all(|(n, min)| aliphatic.get(n).copied().unwrap_or(0) >= *min)
            && self.bond_lower_bounds.iter().all(|b| {
                bonds
                    .get(&(b.first_atomic_number, b.second_atomic_number, b.bond))
                    .copied()
                    .unwrap_or(0)
                    >= b.count
            })
    }
}

fn positive_atomic_number(query: &AtomQuery) -> Option<u8> {
    match query {
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
        AtomQuery::And(a, b) => positive_atomic_number(a).or_else(|| positive_atomic_number(b)),
        AtomQuery::Primitive(_) | AtomQuery::Or(_, _) | AtomQuery::Not(_) => None,
    }
}

fn positive_aromatic(query: &AtomQuery) -> Option<bool> {
    match query {
        AtomQuery::Primitive(AtomPrimitive::Aromatic(v)) => Some(*v),
        AtomQuery::And(a, b) => positive_aromatic(a).or_else(|| positive_aromatic(b)),
        AtomQuery::Primitive(_) | AtomQuery::Or(_, _) | AtomQuery::Not(_) => None,
    }
}

fn simple_bond(query: &BondQuery) -> Option<ReactionBondKind> {
    match query {
        BondQuery::Primitive(BondPrimitive::Single) => Some(ReactionBondKind::Single),
        BondQuery::Primitive(BondPrimitive::Double) => Some(ReactionBondKind::Double),
        BondQuery::Primitive(BondPrimitive::Triple) => Some(ReactionBondKind::Triple),
        BondQuery::Primitive(BondPrimitive::Aromatic) => Some(ReactionBondKind::Aromatic),
        BondQuery::Primitive(_)
        | BondQuery::And(_, _)
        | BondQuery::Or(_, _)
        | BondQuery::Not(_)
        | BondQuery::Any => None,
    }
}

fn target_bond_kind(order: BondOrder) -> Option<ReactionBondKind> {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => {
            Some(ReactionBondKind::Single)
        }
        BondOrder::Double => Some(ReactionBondKind::Double),
        BondOrder::Triple => Some(ReactionBondKind::Triple),
        BondOrder::Aromatic => Some(ReactionBondKind::Aromatic),
        BondOrder::Quadruple
        | BondOrder::Zero
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PreparedReaction;
    use chematic_smiles::parse;

    #[test]
    fn summarizes_compiled_template_and_filters_without_false_negative() {
        let prepared = PreparedReaction::new("[C:1][O:2]>>[C:1][O:2]").unwrap();
        let ethanol = parse("CCO").unwrap();
        let methane = parse("C").unwrap();
        let requirements = prepared.requirements();

        assert_eq!(requirements.min_atom_count, 2);
        assert_eq!(requirements.min_bond_count, 1);
        assert_eq!(requirements.element_lower_bounds, vec![(6, 1), (8, 1)]);
        assert_eq!(
            requirements.aliphatic_element_lower_bounds,
            vec![(6, 1), (8, 1)]
        );
        assert_eq!(
            requirements.bond_lower_bounds,
            vec![ReactionBondLowerBound {
                first_atomic_number: 6,
                second_atomic_number: 8,
                bond: ReactionBondKind::Single,
                count: 1,
            }]
        );
        assert!(!prepared.find_matches(&[&ethanol]).unwrap().is_empty());
        assert!(prepared.could_match(&[&ethanol]));
        assert!(!prepared.could_match(&[&methane]));
    }

    #[test]
    fn atomic_number_variants_fail_open_for_concrete_element_bounds() {
        let prepared = PreparedReaction::new("[#7:1][#6:2]>>[#7:1][#6:2]").unwrap();
        let requirements = prepared.requirements();
        assert_eq!(requirements.min_atom_count, 2);
        assert_eq!(requirements.min_bond_count, 1);
        assert!(requirements.element_lower_bounds.is_empty());
        assert!(requirements.bond_lower_bounds.is_empty());
        assert!(prepared.could_match(&[&parse("NC").unwrap()]));
    }
}
