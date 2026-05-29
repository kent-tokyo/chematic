//! SMILES writer: serialize a Molecule to an OpenSMILES string via DFS traversal.
//!
//! Produces valid (non-canonical) SMILES. Canonical SMILES (Morgan-rank ordering)
//! is a planned future milestone.

use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// Write a [`Molecule`] to a SMILES string.
///
/// Disconnected fragments are joined with `.`.
/// Aromatic atoms are written in lowercase.
pub fn write(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }
    SmilesWriter::new(mol).write_all()
}

struct SmilesWriter<'a> {
    mol: &'a Molecule,
    /// Bonds that are back-edges in the DFS tree (ring closures).
    ring_bonds: HashSet<BondIdx>,
    /// ring number(s) each atom must write when serialized.
    /// Both the "open" ancestor and "close" descendant of a ring store the same number.
    atom_ring_nums: HashMap<AtomIdx, Vec<(u8, BondOrder)>>,
    /// Whether each atom has been serialized in phase 2.
    written: Vec<bool>,
    next_ring: u8,
    out: String,
}

impl<'a> SmilesWriter<'a> {
    fn new(mol: &'a Molecule) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            ring_bonds: HashSet::new(),
            atom_ring_nums: HashMap::new(),
            written: vec![false; n],
            next_ring: 1,
            out: String::new(),
        }
    }

    fn write_all(mut self) -> String {
        // Phase 1: find all back-edges and assign ring-closure numbers.
        self.find_ring_closures();

        // Phase 2: DFS serialization, one fragment at a time.
        let mut first = true;
        for i in 0..self.mol.atom_count() {
            if !self.written[i] {
                if !first { self.out.push('.'); }
                first = false;
                self.write_chain(AtomIdx(i as u32), None, None);
            }
        }

        self.out
    }

    fn find_ring_closures(&mut self) {
        let n = self.mol.atom_count();
        let mut visited  = vec![false; n];
        let mut in_stack = vec![false; n];

        for start in 0..n {
            if !visited[start] {
                self.dfs_mark(AtomIdx(start as u32), None, &mut visited, &mut in_stack);
            }
        }
    }

    /// DFS that marks back-edges and assigns ring-closure numbers.
    ///
    /// `from_bond`: the bond used to arrive at `atom` (skip it to avoid re-visiting the parent).
    fn dfs_mark(
        &mut self,
        atom: AtomIdx,
        from_bond: Option<BondIdx>,
        visited: &mut Vec<bool>,
        in_stack: &mut Vec<bool>,
    ) {
        visited[atom.0 as usize] = true;
        in_stack[atom.0 as usize] = true;

        for (neighbor, bidx) in self.mol.neighbors(atom) {
            // Skip the edge we came from (undirected graph: would look like a back-edge otherwise).
            if Some(bidx) == from_bond { continue; }
            // Skip bonds already classified.
            if self.ring_bonds.contains(&bidx) { continue; }

            if !visited[neighbor.0 as usize] {
                // Tree edge: recurse.
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                // Back-edge: `atom` (descendant) closes a ring back to `neighbor` (ancestor).
                self.ring_bonds.insert(bidx);
                let rn = self.next_ring;
                self.next_ring += 1;
                let bond_order = self.mol.bond(bidx).order;
                // Both endpoints need to emit this ring number when serialized.
                self.atom_ring_nums.entry(neighbor).or_default().push((rn, bond_order)); // open
                self.atom_ring_nums.entry(atom).or_default().push((rn, bond_order));     // close
            }
        }

        in_stack[atom.0 as usize] = false;
    }

    /// Write `atom` and then recurse into its unvisited tree-edge neighbors.
    /// `incoming_bond`: the bond type on the edge leading to this atom (None for the root).
    fn write_chain(
        &mut self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        incoming_bond: Option<BondOrder>,
    ) {
        self.written[atom.0 as usize] = true;

        // Write the incoming bond (if explicit / non-default).
        if let Some(bond) = incoming_bond {
            self.out.push(bond.smiles_char());
        }

        // Write the atom symbol.
        self.emit_atom(atom);

        // Write ring-closure digits for this atom (both open and close digits).
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, bond_order) in rings {
                // Write bond type unless it is implicit.
                let atom_aromatic = self.mol.atom(atom).aromatic;
                // For ring closures we can't know the other atom's aromaticity here,
                // so we emit the bond type unless it is a plain aromatic ring bond.
                if !(bond_order == BondOrder::Aromatic && atom_aromatic) && bond_order != BondOrder::Single {
                    self.out.push(bond_order.smiles_char());
                }
                // Ring number: single digit for 1-9, `%NN` form for 10+.
                if rn >= 10 {
                    self.out.push('%');
                    self.out.push(char::from_digit((rn / 10) as u32, 10).unwrap());
                    self.out.push(char::from_digit((rn % 10) as u32, 10).unwrap());
                } else {
                    self.out.push(char::from_digit(rn as u32, 10).unwrap());
                }
            }
        }

        // Collect tree-edge children (unvisited, non-ring-closure bonds).
        let children: Vec<(AtomIdx, BondOrder)> = self.mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, bidx)| (nb, self.mol.bond(bidx).order))
            .collect();

        // Write children: all but the last one are branches (wrapped in parentheses).
        let n = children.len();
        for (i, (child, bond_order)) in children.into_iter().enumerate() {
            let is_last = i == n - 1;

            // Determine whether the bond should be written explicitly.
            let parent_arom = self.mol.atom(atom).aromatic;
            let child_arom  = self.mol.atom(child).aromatic;
            let implicit = match bond_order {
                BondOrder::Single   => !(parent_arom && child_arom), // single is implicit
                BondOrder::Aromatic => parent_arom && child_arom,    // aromatic is implicit
                _ => false,
            };
            let written_bond = if implicit { None } else { Some(bond_order) };

            if !is_last {
                self.out.push('(');
                self.write_chain(child, Some(atom), written_bond);
                self.out.push(')');
            } else {
                self.write_chain(child, Some(atom), written_bond);
            }
        }
    }

    fn emit_atom(&mut self, idx: AtomIdx) {
        let atom = self.mol.atom(idx);

        // An atom needs bracket notation when:
        //  - it has an isotope, charge, explicit H count, atom map, or
        //  - it is not in the organic subset (cannot rely on implicit-H rules).
        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogen_count.is_some()
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some();

        if needs_bracket {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            let sym = if atom.aromatic {
                atom.element.symbol().to_lowercase()
            } else {
                atom.element.symbol().to_string()
            };
            self.out.push_str(&sym);

            match atom.chirality {
                chematic_core::Chirality::CounterClockwise => self.out.push('@'),
                chematic_core::Chirality::Clockwise        => self.out.push_str("@@"),
                chematic_core::Chirality::None             => {}
            }

            if let Some(h) = atom.hydrogen_count {
                if h > 0 {
                    self.out.push('H');
                    if h > 1 { self.out.push_str(&h.to_string()); }
                }
            }

            match atom.charge {
                0 => {}
                1 => self.out.push('+'),
                -1 => self.out.push('-'),
                c if c > 0 => self.out.push_str(&format!("+{c}")),
                c => self.out.push_str(&c.to_string()),
            }

            if let Some(m) = atom.atom_map {
                self.out.push(':');
                self.out.push_str(&m.to_string());
            }

            self.out.push(']');
        } else if atom.aromatic {
            self.out.push_str(&atom.element.symbol().to_lowercase());
        } else {
            self.out.push_str(atom.element.symbol());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Parse → write → re-parse → verify atom/bond counts are preserved.
    fn roundtrip(smiles: &str) {
        let mol1 = parse(smiles).expect(smiles);
        let out  = write(&mol1);
        let mol2 = parse(&out).unwrap_or_else(|e| {
            panic!("roundtrip failed for '{}': wrote '{}', error: {e}", smiles, out)
        });
        assert_eq!(mol1.atom_count(), mol2.atom_count(),
            "atom count mismatch: input='{}' output='{}'", smiles, out);
        assert_eq!(mol1.bond_count(), mol2.bond_count(),
            "bond count mismatch: input='{}' output='{}'", smiles, out);
    }

    #[test]
    fn test_write_methane()          { assert_eq!(write(&parse("C").unwrap()), "C"); }
    #[test]
    fn test_write_ethane()           { assert_eq!(write(&parse("CC").unwrap()), "CC"); }

    #[test]
    fn test_roundtrip_propane()      { roundtrip("CCC"); }
    #[test]
    fn test_roundtrip_isobutane()    { roundtrip("CC(C)C"); }
    #[test]
    fn test_roundtrip_ethanol()      { roundtrip("CCO"); }
    #[test]
    fn test_roundtrip_acetic_acid()  { roundtrip("CC(=O)O"); }
    #[test]
    fn test_roundtrip_cyclohexane()  { roundtrip("C1CCCCC1"); }
    #[test]
    fn test_roundtrip_benzene_kekule() { roundtrip("C1=CC=CC=C1"); }
    #[test]
    fn test_roundtrip_benzene_arom() { roundtrip("c1ccccc1"); }
    #[test]
    fn test_roundtrip_pyridine()     { roundtrip("c1ccncc1"); }
    #[test]
    fn test_roundtrip_naphthalene()  { roundtrip("c1ccc2ccccc2c1"); }
    #[test]
    fn test_roundtrip_chlorobenzene() { roundtrip("c1ccccc1Cl"); }
    #[test]
    fn test_roundtrip_13c()          { roundtrip("[13C]"); }
    #[test]
    fn test_roundtrip_ammonium()     { roundtrip("[NH4+]"); }
    #[test]
    fn test_roundtrip_disconnected() { roundtrip("[Na+].[Cl-]"); }
    #[test]
    fn test_roundtrip_aspirin()      { roundtrip("CC(=O)Oc1ccccc1C(=O)O"); }
    #[test]
    fn test_roundtrip_caffeine()     { roundtrip("Cn1cnc2c1c(=O)n(c(=O)n2C)C"); }
}
