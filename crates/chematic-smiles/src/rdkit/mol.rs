//! The parts of RDKit's molecule model (`Atom`, `Bond`, `ROMol`/`RWMol`,
//! `RingInfo`) that `MolFromSmiles` + `MolToSmiles` read and write, with
//! RDKit's indexing: atoms and bonds are numbered as RDKit numbers them and
//! every atom lists its bonds in bond-index order (RDKit's adjacency order).

use super::RdkitSmilesError;
use super::periodic;

/// `Atom::ChiralType` (the tags this pipeline carries).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ChiralTag {
    Unspecified,
    /// `CHI_TETRAHEDRAL_CW` (`@@` relative to the bond order).
    Cw,
    /// `CHI_TETRAHEDRAL_CCW` (`@` relative to the bond order).
    Ccw,
}

/// `Atom::HybridizationType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Hybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
}

/// `Bond::BondType` (numeric values as in RDKit).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BondType {
    Single = 1,
    Double = 2,
    Triple = 3,
    Quadruple = 4,
    Aromatic = 12,
    Dative = 17,
}

impl BondType {
    /// `Bond::getBondTypeAsDouble`.
    pub(crate) fn as_double(self) -> f64 {
        match self {
            BondType::Single | BondType::Dative => 1.0,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Quadruple => 4.0,
            BondType::Aromatic => 1.5,
        }
    }

    /// `getTwiceBondType`.
    pub(crate) fn twice(self) -> u32 {
        match self {
            BondType::Single | BondType::Dative => 2,
            BondType::Double => 4,
            BondType::Triple => 6,
            BondType::Quadruple => 8,
            BondType::Aromatic => 3,
        }
    }
}

/// `Bond::BondDir` (the values this pipeline carries).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BondDir {
    None,
    EndDownRight,
    EndUpRight,
}

impl BondDir {
    pub(crate) fn flipped(self) -> BondDir {
        match self {
            BondDir::EndDownRight => BondDir::EndUpRight,
            BondDir::EndUpRight => BondDir::EndDownRight,
            BondDir::None => BondDir::None,
        }
    }

    pub(crate) fn is_set(self) -> bool {
        self != BondDir::None
    }
}

/// `Bond::BondStereo` (numeric values as in RDKit).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum BondStereo {
    None = 0,
    Any = 1,
    Z = 2,
    E = 3,
}

#[derive(Clone, Debug)]
pub(crate) struct Atom {
    pub anum: u32,
    pub isotope: u32,
    pub charge: i32,
    pub num_explicit_hs: u32,
    pub no_implicit: bool,
    pub aromatic: bool,
    pub chiral: ChiralTag,
    pub radicals: u32,
    /// `molAtomMapNumber`.
    pub map: Option<u32>,
    /// `d_explicitValence` (-1: not computed).
    pub explicit_valence: i32,
    /// `d_implicitValence` (-1: not computed).
    pub implicit_valence: i32,
    pub hybrid: Hybridization,
    /// `_CIPCode` (`R` or `S`).
    pub cip_code: Option<u8>,
    /// `_ringStereoAtoms`.
    pub ring_stereo_atoms: Option<Vec<i32>>,
    /// `_ringStereochemCand`.
    pub ring_stereochem_cand: Option<bool>,
    /// `_ChiralityPossible`.
    pub chirality_possible: bool,
}

impl Atom {
    pub(crate) fn new(anum: u32) -> Atom {
        Atom {
            anum,
            isotope: 0,
            charge: 0,
            num_explicit_hs: 0,
            no_implicit: false,
            aromatic: false,
            chiral: ChiralTag::Unspecified,
            radicals: 0,
            map: None,
            explicit_valence: -1,
            implicit_valence: -1,
            hybrid: Hybridization::Unspecified,
            cip_code: None,
            ring_stereo_atoms: None,
            ring_stereochem_cand: None,
            chirality_possible: false,
        }
    }

    /// `Atom::invertChirality` for tetrahedral tags.
    pub(crate) fn invert_chirality(&mut self) {
        self.chiral = match self.chiral {
            ChiralTag::Cw => ChiralTag::Ccw,
            ChiralTag::Ccw => ChiralTag::Cw,
            ChiralTag::Unspecified => ChiralTag::Unspecified,
        };
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Bond {
    pub begin: usize,
    pub end: usize,
    pub bt: BondType,
    pub aromatic: bool,
    pub conjugated: bool,
    pub dir: BondDir,
    pub stereo: BondStereo,
    pub stereo_atoms: Vec<usize>,
}

impl Bond {
    pub(crate) fn new(begin: usize, end: usize, bt: BondType) -> Bond {
        Bond {
            begin,
            end,
            bt,
            aromatic: false,
            conjugated: false,
            dir: BondDir::None,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
        }
    }

    pub(crate) fn other(&self, atom: usize) -> usize {
        if self.begin == atom {
            self.end
        } else {
            self.begin
        }
    }

    /// `Bond::getValenceContrib`.
    pub(crate) fn valence_contrib(&self, atom: usize) -> f64 {
        if atom != self.begin && atom != self.end {
            return 0.0;
        }
        if self.bt == BondType::Dative && atom != self.end {
            return 0.0;
        }
        self.bt.as_double()
    }

    /// `canHaveDirection`.
    pub(crate) fn can_have_direction(&self) -> bool {
        matches!(self.bt, BondType::Single | BondType::Aromatic)
    }

    /// `canSetDoubleBondStereo`.
    pub(crate) fn can_set_double_bond_stereo(&self) -> bool {
        matches!(
            self.bt,
            BondType::Single | BondType::Aromatic | BondType::Dative
        )
    }
}

/// `RingInfo` with the rings of a symmetrized SSSR.
#[derive(Clone, Debug, Default)]
pub(crate) struct RingInfo {
    pub atom_rings: Vec<Vec<usize>>,
    pub bond_rings: Vec<Vec<usize>>,
    /// Ring ids each atom is in, in ring order (`atomMembers`).
    pub atom_members: Vec<Vec<usize>>,
    /// Ring ids each bond is in, in ring order.
    pub bond_members: Vec<Vec<usize>>,
}

impl RingInfo {
    pub(crate) fn num_atom_rings(&self, a: usize) -> usize {
        self.atom_members[a].len()
    }

    pub(crate) fn num_bond_rings(&self, b: usize) -> usize {
        self.bond_members[b].len()
    }

    pub(crate) fn min_bond_ring_size(&self, b: usize) -> usize {
        self.bond_members[b]
            .iter()
            .map(|&r| self.bond_rings[r].len())
            .min()
            .unwrap_or(0)
    }

    pub(crate) fn is_atom_in_ring_of_size(&self, a: usize, size: usize) -> bool {
        self.atom_members[a]
            .iter()
            .any(|&r| self.atom_rings[r].len() == size)
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Mol {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    /// Per atom: its bond indices in bond-index order.
    pub atom_bonds: Vec<Vec<usize>>,
    pub rings: Option<RingInfo>,
}

impl Mol {
    pub(crate) fn add_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atom_bonds.push(Vec::new());
        self.atoms.len() - 1
    }

    pub(crate) fn add_bond(&mut self, bond: Bond) -> usize {
        let idx = self.bonds.len();
        self.atom_bonds[bond.begin].push(idx);
        self.atom_bonds[bond.end].push(idx);
        self.bonds.push(bond);
        idx
    }

    pub(crate) fn degree(&self, a: usize) -> usize {
        self.atom_bonds[a].len()
    }

    pub(crate) fn nbrs(&self, a: usize) -> impl Iterator<Item = usize> + '_ {
        self.atom_bonds[a]
            .iter()
            .map(move |&b| self.bonds[b].other(a))
    }

    pub(crate) fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.atom_bonds[a]
            .iter()
            .copied()
            .find(|&bi| self.bonds[bi].other(a) == b)
    }

    pub(crate) fn ring_info(&self) -> &RingInfo {
        self.rings.as_ref().expect("ring info is computed")
    }

    pub(crate) fn num_atom_rings(&self, a: usize) -> usize {
        self.rings.as_ref().map_or(0, |r| r.num_atom_rings(a))
    }

    pub(crate) fn num_bond_rings(&self, b: usize) -> usize {
        self.rings.as_ref().map_or(0, |r| r.num_bond_rings(b))
    }

    /// `isAromaticAtom`.
    pub(crate) fn is_aromatic_atom(&self, a: usize) -> bool {
        if self.atoms[a].aromatic {
            return true;
        }
        self.atom_bonds[a].iter().any(|&b| {
            let bond = &self.bonds[b];
            bond.aromatic || bond.bt == BondType::Aromatic
        })
    }

    /// `Atom::getNumImplicitHs`.
    pub(crate) fn num_implicit_hs(&self, a: usize) -> u32 {
        let atom = &self.atoms[a];
        if atom.no_implicit {
            return 0;
        }
        atom.implicit_valence.max(0) as u32
    }

    /// `Atom::getTotalNumHs()`.
    pub(crate) fn total_num_hs(&self, a: usize) -> u32 {
        self.atoms[a].num_explicit_hs + self.num_implicit_hs(a)
    }

    /// `Atom::getValence(IMPLICIT)`.
    pub(crate) fn implicit_valence(&self, a: usize) -> i32 {
        let atom = &self.atoms[a];
        if atom.no_implicit {
            0
        } else {
            atom.implicit_valence
        }
    }

    /// `Atom::getTotalValence`.
    pub(crate) fn total_valence(&self, a: usize) -> i32 {
        self.atoms[a].explicit_valence + self.implicit_valence(a)
    }

    /// `Atom::getTotalDegree`.
    pub(crate) fn total_degree(&self, a: usize) -> usize {
        self.total_num_hs(a) as usize + self.degree(a)
    }

    /// `Atom::needsUpdatePropertyCache`.
    pub(crate) fn needs_update_property_cache(&self, a: usize) -> bool {
        let atom = &self.atoms[a];
        !(atom.explicit_valence >= 0 && (atom.no_implicit || atom.implicit_valence >= 0))
    }

    /// Parity of `Atom::getPerturbationOrder(probe)`: whether the bond list
    /// `probe` is an odd permutation of the atom's bond list.
    pub(crate) fn perturbation_is_odd(&self, a: usize, probe: &[usize]) -> bool {
        count_swaps(probe, &self.atom_bonds[a]) % 2 == 1
    }

    /// `getEffectiveAtomicNum(atom, false)`.
    fn effective_atomic_num(&self, a: usize) -> u32 {
        let atom = &self.atoms[a];
        (atom.anum as i32 - atom.charge).clamp(0, periodic::max_atomic_num()) as u32
    }

    /// `canBeHypervalent`.
    fn can_be_hypervalent(&self, a: usize, eff: u32) -> bool {
        let anum = self.atoms[a].anum;
        (eff > 16 && (anum == 15 || anum == 16)) || (eff > 34 && (anum == 33 || anum == 34))
    }

    /// `calculateExplicitValence(atom, strict, checkIt=false)`.
    pub(crate) fn calc_explicit_valence_value(
        &self,
        a: usize,
        strict: bool,
    ) -> Result<i32, RdkitSmilesError> {
        let atom = &self.atoms[a];
        let mut accum = 0.0f64;
        for &b in &self.atom_bonds[a] {
            accum += self.bonds[b].valence_contrib(a);
        }
        accum += f64::from(atom.num_explicit_hs);
        let ovalens = periodic::valence_list(atom.anum);
        let mut eff = atom.anum;
        if ovalens.len() > 1 || ovalens[0] != -1 {
            eff = self.effective_atomic_num(a);
        }
        let dv = periodic::default_valence(eff);
        let valens = periodic::valence_list(eff);
        // `unsigned int dv`: -1 never compares below `accum`.
        if dv >= 0 && accum > f64::from(dv) && self.is_aromatic_atom(a) {
            let mut pval = dv;
            for &val in valens {
                if val == -1 {
                    break;
                }
                if f64::from(val) > accum {
                    break;
                }
                pval = i32::from(val);
            }
            if accum - f64::from(pval) <= 1.5 {
                accum = f64::from(pval);
            }
        }
        accum += 0.1;
        let res = accum.round() as i32;
        if strict {
            let mut max_valence = i32::from(*valens.last().unwrap_or(&-1));
            let mut offset = 0;
            if self.can_be_hypervalent(a, eff) {
                max_valence = i32::from(*ovalens.last().unwrap_or(&-1));
                offset -= atom.charge;
            }
            if atom.anum == 1 && atom.charge == -1 {
                max_valence = 2;
            }
            if max_valence >= 0 && *ovalens.last().unwrap_or(&-1) >= 0 && res + offset > max_valence
            {
                return Err(RdkitSmilesError::Sanitization(format!(
                    "Explicit valence for atom # {a} {}, {res}, is greater than permitted",
                    periodic::symbol(atom.anum)
                )));
            }
        }
        Ok(res)
    }

    /// `Atom::calcExplicitValence(strict)`.
    pub(crate) fn calc_explicit_valence(
        &mut self,
        a: usize,
        strict: bool,
    ) -> Result<i32, RdkitSmilesError> {
        let v = self.calc_explicit_valence_value(a, strict)?;
        self.atoms[a].explicit_valence = v;
        Ok(v)
    }

    /// `calculateImplicitValence(atom, strict, checkIt=false)`.
    pub(crate) fn calc_implicit_valence_value(
        &self,
        a: usize,
        strict: bool,
    ) -> Result<i32, RdkitSmilesError> {
        let atom = &self.atoms[a];
        if atom.no_implicit {
            return Ok(0);
        }
        let mut explicit_valence = atom.explicit_valence;
        if explicit_valence == -1 {
            explicit_valence = self.calc_explicit_valence_value(a, strict)?;
        }
        let anum = atom.anum;
        if anum == 0 {
            return Ok(0);
        }
        let charge = atom.charge;
        let radicals = atom.radicals as i32;
        if explicit_valence == 0 && radicals == 0 && anum == 1 {
            if charge == 1 || charge == -1 {
                return Ok(0);
            } else if charge == 0 {
                return Ok(1);
            } else if strict {
                return Err(RdkitSmilesError::Sanitization(format!(
                    "Unreasonable formal charge on atom # {a}."
                )));
            } else {
                return Ok(0);
            }
        }
        // RDKit reads the stored explicit valence here.
        let mut explicit_plus_rad = atom.explicit_valence + radicals;
        let ovalens = periodic::valence_list(anum);
        let mut eff = anum;
        if ovalens.len() > 1 || ovalens[0] != -1 {
            eff = self.effective_atomic_num(a);
        }
        if eff == 0 {
            return Ok(0);
        }
        let dv = periodic::default_valence(eff);
        if dv == -1 {
            return Ok(0);
        }
        if self.can_be_hypervalent(a, eff) {
            eff = anum;
            explicit_plus_rad -= charge;
        }
        let valens = periodic::valence_list(eff);
        let mut res;
        if self.is_aromatic_atom(a) {
            if explicit_plus_rad <= dv {
                res = dv - explicit_plus_rad;
            } else {
                let mut satis = false;
                for &v in valens {
                    if v <= 0 {
                        break;
                    }
                    if explicit_plus_rad == i32::from(v) {
                        satis = true;
                        break;
                    }
                }
                if !satis && strict {
                    return Err(RdkitSmilesError::Sanitization(format!(
                        "Explicit valence for aromatic atom # {a} not equal to any accepted valence"
                    )));
                }
                res = 0;
            }
        } else {
            res = -1;
            for &v in valens {
                if v < 0 {
                    break;
                }
                let tot = i32::from(v);
                if explicit_plus_rad <= tot {
                    res = tot - explicit_plus_rad;
                    break;
                }
            }
            if res < 0 {
                if strict
                    && *valens.last().unwrap_or(&-1) != -1
                    && *ovalens.last().unwrap_or(&-1) > 0
                {
                    return Err(RdkitSmilesError::Sanitization(format!(
                        "Explicit valence for atom # {a} {} greater than permitted",
                        periodic::symbol(anum)
                    )));
                }
                res = 0;
            }
        }
        Ok(res)
    }

    /// `Atom::calcImplicitValence(strict)`.
    pub(crate) fn calc_implicit_valence(
        &mut self,
        a: usize,
        strict: bool,
    ) -> Result<i32, RdkitSmilesError> {
        if self.atoms[a].explicit_valence == -1 {
            self.calc_explicit_valence(a, strict)?;
        }
        let v = self.calc_implicit_valence_value(a, strict)?;
        self.atoms[a].implicit_valence = v;
        Ok(v)
    }

    /// `Atom::updatePropertyCache(strict)`.
    pub(crate) fn update_atom_property_cache(
        &mut self,
        a: usize,
        strict: bool,
    ) -> Result<(), RdkitSmilesError> {
        self.calc_explicit_valence(a, strict)?;
        self.calc_implicit_valence(a, strict)?;
        Ok(())
    }

    /// `ROMol::updatePropertyCache(strict)`.
    pub(crate) fn update_property_cache(&mut self, strict: bool) -> Result<(), RdkitSmilesError> {
        for a in 0..self.atoms.len() {
            self.update_atom_property_cache(a, strict)?;
        }
        Ok(())
    }

    /// `RWMol::removeAtom`: drops the atom and its bonds and renumbers the
    /// rest, keeping their order. Ring information is reset.
    pub(crate) fn remove_atom(&mut self, idx: usize) {
        let removed: Vec<usize> = self.atom_bonds[idx].clone();
        let mut keep = vec![true; self.bonds.len()];
        for &b in &removed {
            keep[b] = false;
        }
        let mut new_bond_idx = vec![usize::MAX; self.bonds.len()];
        let mut bonds = Vec::with_capacity(self.bonds.len());
        for (b, bond) in self.bonds.drain(..).enumerate() {
            if keep[b] {
                new_bond_idx[b] = bonds.len();
                bonds.push(bond);
            }
        }
        let fix_atom = |x: usize| if x > idx { x - 1 } else { x };
        for bond in &mut bonds {
            bond.begin = fix_atom(bond.begin);
            bond.end = fix_atom(bond.end);
            bond.stereo_atoms = bond.stereo_atoms.iter().map(|&x| fix_atom(x)).collect();
        }
        self.bonds = bonds;
        self.atoms.remove(idx);
        self.atom_bonds.remove(idx);
        for list in &mut self.atom_bonds {
            *list = list
                .iter()
                .filter(|&&b| keep[b])
                .map(|&b| new_bond_idx[b])
                .collect();
        }
        self.rings = None;
    }

    /// `MolOps::symmetrizeSSSR`: (re)computes the ring information.
    pub(crate) fn find_rings(&mut self) -> Result<(), RdkitSmilesError> {
        let bonds: Vec<(usize, usize, bool)> = self
            .bonds
            .iter()
            .map(|b| (b.begin, b.end, b.bt != BondType::Dative))
            .collect();
        let rings = chematic_perception::rdkit_symmetrized_sssr(self.atoms.len(), &bonds)
            .ok_or_else(|| {
                RdkitSmilesError::Unsupported(
                    "ring perception falls back to RDKit's approximate ring finder".into(),
                )
            })?;
        self.set_rings(rings);
        Ok(())
    }

    pub(crate) fn set_rings(&mut self, atom_rings: Vec<Vec<usize>>) {
        let mut info = RingInfo {
            atom_members: vec![Vec::new(); self.atoms.len()],
            bond_members: vec![Vec::new(); self.bonds.len()],
            ..RingInfo::default()
        };
        for (ri, ring) in atom_rings.iter().enumerate() {
            let mut bring = Vec::with_capacity(ring.len());
            for k in 0..ring.len() {
                let b = self
                    .bond_between(ring[k], ring[(k + 1) % ring.len()])
                    .expect("ring atoms are bonded");
                bring.push(b);
            }
            for &a in ring {
                info.atom_members[a].push(ri);
            }
            for &b in &bring {
                info.bond_members[b].push(ri);
            }
            info.bond_rings.push(bring);
        }
        info.atom_rings = atom_rings;
        self.rings = Some(info);
    }
}

/// RDKit `countSwapsToInterconvert(ref, probe)`.
pub(crate) fn count_swaps(reference: &[usize], probe: &[usize]) -> usize {
    let mut probe = probe.to_vec();
    let mut n = 0;
    for i in 0..reference.len().min(probe.len()) {
        if probe[i] != reference[i]
            && let Some(j) = (i..probe.len()).find(|&j| probe[j] == reference[i])
        {
            probe.swap(i, j);
            n += 1;
        }
    }
    n
}
