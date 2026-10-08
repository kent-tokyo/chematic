//! A read-only view of an explicit-hydrogen molecule as RDKit models it
//! after `Chem.AddHs(Chem.MolFromSmiles(s))`: RDKit atom and bond order,
//! bond begin/end atoms, per-atom bond lists, sanitization results
//! (hybridization, conjugation, aromaticity, symmetrized SSSR) and legacy
//! stereo (chiral tags relative to RDKit's bond order, double-bond stereo
//! with stereo atoms). Consumed by RDKit-exact ports that need RDKit's
//! molecule model (e.g. distance-geometry embedding).

use chematic_core::{AtomIdx, Element, Molecule};

use super::RdkitSmilesError;
use super::mol::{BondStereo, ChiralTag, Hybridization};

/// One atom of [`RdkitMolView`].
#[derive(Clone, Debug)]
pub struct RdkitViewAtom {
    pub atomic_num: u32,
    pub formal_charge: i32,
    pub is_aromatic: bool,
    /// `Atom::HybridizationType`: 0 unspecified, 1 S, 2 SP, 3 SP2, 4 SP3,
    /// 5 SP3D, 6 SP3D2 (RDKit's enum values).
    pub hybridization: u8,
    /// 0 unspecified, 1 CW, 2 CCW (RDKit's `ChiralType`).
    pub chiral_tag: u8,
    /// `getTotalNumHs(includeNeighbors=true)`.
    pub total_num_hs_with_neighbors: u32,
    /// `getTotalNumHs()`.
    pub total_num_hs: u32,
    /// `PeriodicTable::getRvdw`.
    pub rvdw: f64,
}

/// One bond of [`RdkitMolView`].
#[derive(Clone, Debug)]
pub struct RdkitViewBond {
    pub begin: usize,
    pub end: usize,
    /// RDKit `BondType` value (1 single, 2 double, 3 triple, 4 quadruple,
    /// 12 aromatic, 17 dative).
    pub bond_type: u8,
    /// `getBondTypeAsDouble`.
    pub bond_type_as_double: f64,
    pub is_aromatic: bool,
    pub is_conjugated: bool,
    /// RDKit `BondStereo` value: 0 none, 1 any, 2 Z, 3 E.
    pub stereo: u8,
    pub stereo_atoms: Vec<usize>,
}

/// See the module documentation.
#[derive(Clone, Debug)]
pub struct RdkitMolView {
    pub atoms: Vec<RdkitViewAtom>,
    pub bonds: Vec<RdkitViewBond>,
    /// Per atom: its bond indices in RDKit's adjacency order.
    pub atom_bonds: Vec<Vec<usize>>,
    /// `RingInfo::atomRings()` / `bondRings()`.
    pub atom_rings: Vec<Vec<usize>>,
    pub bond_rings: Vec<Vec<usize>>,
    /// For each RDKit bond index, the chematic bond index.
    pub chematic_bond: Vec<usize>,
}

impl RdkitMolView {
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    pub fn degree(&self, a: usize) -> usize {
        self.atom_bonds[a].len()
    }

    pub fn other_atom(&self, b: usize, a: usize) -> usize {
        let bond = &self.bonds[b];
        if bond.begin == a {
            bond.end
        } else {
            bond.begin
        }
    }

    /// `getBondBetweenAtoms`.
    pub fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.atom_bonds[a]
            .iter()
            .copied()
            .find(|&x| self.other_atom(x, a) == b)
    }

    /// Neighbours in adjacency order.
    pub fn neighbors(&self, a: usize) -> impl Iterator<Item = usize> + '_ {
        self.atom_bonds[a]
            .iter()
            .map(move |&b| self.other_atom(b, a))
    }

    pub fn num_atom_rings(&self, a: usize) -> usize {
        self.atom_rings.iter().filter(|r| r.contains(&a)).count()
    }

    pub fn num_bond_rings(&self, b: usize) -> usize {
        self.bond_rings.iter().filter(|r| r.contains(&b)).count()
    }

    pub fn is_atom_in_ring_of_size(&self, a: usize, size: usize) -> bool {
        self.atom_rings
            .iter()
            .any(|r| r.len() == size && r.contains(&a))
    }

    pub fn is_bond_in_ring_of_size(&self, b: usize, size: usize) -> bool {
        self.bond_rings
            .iter()
            .any(|r| r.len() == size && r.contains(&b))
    }

    /// `RingInfo::atomRingSizes`.
    pub fn atom_ring_sizes(&self, a: usize) -> Vec<usize> {
        self.atom_rings
            .iter()
            .filter(|r| r.contains(&a))
            .map(|r| r.len())
            .collect()
    }

    /// RDKit's connected components (`MolOps::getMolFrags(mol, mapping)`:
    /// fragments numbered in order of their lowest atom index), each as its
    /// sorted atom indices.
    pub fn fragments(&self) -> Vec<Vec<usize>> {
        let n = self.num_atoms();
        let mut comp = vec![usize::MAX; n];
        let mut frags: Vec<Vec<usize>> = Vec::new();
        for s in 0..n {
            if comp[s] != usize::MAX {
                continue;
            }
            let id = frags.len();
            comp[s] = id;
            let mut stack = vec![s];
            let mut members = vec![s];
            while let Some(a) = stack.pop() {
                for b in self.neighbors(a) {
                    if comp[b] == usize::MAX {
                        comp[b] = id;
                        stack.push(b);
                        members.push(b);
                    }
                }
            }
            members.sort_unstable();
            frags.push(members);
        }
        frags
    }

    /// The fragment `MolOps::getMolFrags(mol, sanitizeFrags=true)` returns
    /// for the sorted atom set `atoms` (one connected component): atoms and
    /// bonds keep their relative order, and ring perception is redone on the
    /// fragment as `sanitizeMol` does. `None` where RDKit's ring finder
    /// would fall back to its approximate algorithm.
    pub fn fragment(&self, atoms: &[usize]) -> Option<RdkitMolView> {
        let mut new_idx = vec![usize::MAX; self.num_atoms()];
        for (k, &a) in atoms.iter().enumerate() {
            new_idx[a] = k;
        }
        let mut new_bond = vec![usize::MAX; self.num_bonds()];
        let mut bonds = Vec::new();
        let mut chematic_bond = Vec::new();
        for (bi, b) in self.bonds.iter().enumerate() {
            if new_idx[b.begin] == usize::MAX || new_idx[b.end] == usize::MAX {
                continue;
            }
            new_bond[bi] = bonds.len();
            let mut nb = b.clone();
            nb.begin = new_idx[b.begin];
            nb.end = new_idx[b.end];
            nb.stereo_atoms = b.stereo_atoms.iter().map(|&s| new_idx[s]).collect();
            bonds.push(nb);
            chematic_bond.push(self.chematic_bond[bi]);
        }
        let atom_bonds: Vec<Vec<usize>> = atoms
            .iter()
            .map(|&a| self.atom_bonds[a].iter().map(|&b| new_bond[b]).collect())
            .collect();
        let ring_input: Vec<(usize, usize, bool)> = bonds
            .iter()
            .map(|b: &RdkitViewBond| (b.begin, b.end, b.bond_type != 17))
            .collect();
        let atom_rings = chematic_perception::rdkit_symmetrized_sssr(atoms.len(), &ring_input)?;
        let mut view = RdkitMolView {
            atoms: atoms.iter().map(|&a| self.atoms[a].clone()).collect(),
            bonds,
            atom_bonds,
            atom_rings: Vec::new(),
            bond_rings: Vec::new(),
            chematic_bond,
        };
        let bond_rings = atom_rings
            .iter()
            .map(|ring| {
                (0..ring.len())
                    .map(|k| {
                        view.bond_between(ring[k], ring[(k + 1) % ring.len()])
                            .expect("ring atoms are bonded")
                    })
                    .collect()
            })
            .collect();
        view.atom_rings = atom_rings;
        view.bond_rings = bond_rings;
        Some(view)
    }
}

/// Number of trailing hydrogen atoms appended by `add_hydrogens` (RDKit
/// `AddHs`): H atoms without a hydrogen count, bonded once.
fn added_h_start(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    let mut h0 = n;
    while h0 > 0 {
        let a = mol.atom(AtomIdx((h0 - 1) as u32));
        if !a.wildcard
            && a.element == Element::H
            && a.hydrogen_count.is_none()
            && mol.degree(AtomIdx((h0 - 1) as u32)) == 1
        {
            h0 -= 1;
        } else {
            break;
        }
    }
    h0
}

/// The RDKit view of `mol`, which must be an explicit-hydrogen molecule
/// produced by `add_hydrogens` from a molecule read from SMILES (RDKit's
/// `AddHs(MolFromSmiles(s))`).
pub fn rdkit_mol_view(mol: &Molecule) -> Result<RdkitMolView, RdkitSmilesError> {
    let h0 = added_h_start(mol);
    // RDKit numbers the parsed molecule's bonds first (chain bonds, then
    // ring closures) and appends AddHs' hydrogen bonds in creation order.
    let added = |b: chematic_core::BondIdx| {
        let bond = mol.bond(b);
        bond.atom1.0 as usize >= h0 || bond.atom2.0 as usize >= h0
    };
    let mut order: Vec<chematic_core::BondIdx> = mol
        .rdkit_bond_order()
        .into_iter()
        .filter(|&b| !added(b))
        .collect();
    let mut hs: Vec<chematic_core::BondIdx> = (0..mol.bond_count() as u32)
        .map(chematic_core::BondIdx)
        .filter(|&b| added(b))
        .collect();
    hs.sort_unstable_by_key(|b| b.0);
    order.extend(hs);

    let mut m = super::parse::from_chematic_ordered(mol, &order)?;
    super::sanitize::sanitize_keeping_hs(&mut m)?;
    super::stereo::legacy_stereo_perception(&mut m, true, true);

    let rings = m.ring_info().clone();
    let atoms = (0..m.atoms.len())
        .map(|a| {
            let at = &m.atoms[a];
            let h_nbrs = m.nbrs(a).filter(|&o| m.atoms[o].anum == 1).count() as u32;
            RdkitViewAtom {
                atomic_num: at.anum,
                formal_charge: at.charge,
                is_aromatic: at.aromatic,
                // AddHs does not set the new atoms' hybridization.
                hybridization: if a >= h0 {
                    0
                } else {
                    match at.hybrid {
                        Hybridization::Unspecified => 0,
                        Hybridization::S => 1,
                        Hybridization::Sp => 2,
                        Hybridization::Sp2 => 3,
                        Hybridization::Sp3 => 4,
                        Hybridization::Sp3d => 5,
                        Hybridization::Sp3d2 => 6,
                    }
                },
                chiral_tag: match at.chiral {
                    ChiralTag::Unspecified => 0,
                    ChiralTag::Cw => 1,
                    ChiralTag::Ccw => 2,
                },
                total_num_hs_with_neighbors: m.total_num_hs(a) + h_nbrs,
                total_num_hs: m.total_num_hs(a),
                rvdw: super::periodic::RVDW[(at.anum as usize).min(118)],
            }
        })
        .collect();
    let bonds = m
        .bonds
        .iter()
        .map(|b| RdkitViewBond {
            begin: b.begin,
            end: b.end,
            bond_type: b.bt as u8,
            bond_type_as_double: b.bt.as_double(),
            is_aromatic: b.aromatic,
            is_conjugated: b.conjugated,
            stereo: match b.stereo {
                BondStereo::None => 0,
                BondStereo::Any => 1,
                BondStereo::Z => 2,
                BondStereo::E => 3,
            },
            stereo_atoms: b.stereo_atoms.clone(),
        })
        .collect();
    Ok(RdkitMolView {
        atoms,
        bonds,
        atom_bonds: m.atom_bonds.clone(),
        atom_rings: rings.atom_rings,
        bond_rings: rings.bond_rings,
        chematic_bond: order.iter().map(|b| b.0 as usize).collect(),
    })
}
