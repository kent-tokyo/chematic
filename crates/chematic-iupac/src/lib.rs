//! `chematic-iupac` — local IUPAC name generation, no network required.
//!
//! Supports:
//! - Linear alkanes and cycloalkanes
//! - Alkenes (`-ene`) and alkynes (`-yne`) with one unsaturation
//! - Simple derivatives: alcohols (`-ol`), amines (`-amine`), aldehydes (`-al`),
//!   ketones (`-one` with position locant), carboxylic acids (`-oic acid`)
//! - Esters (`alkyl alkanoate`) — linear, primary esters
//! - Primary/secondary amides (`-anamide`)
//! - Halogen substituents: fluoro-, chloro-, bromo-, iodo-
//! - Common aromatic heterocycles: benzene, pyridine, furan, thiophene,
//!   pyrrole, imidazole, pyrimidine
//!
//! Complex polycyclic systems, stereo descriptors, and structures outside
//! the above scope return [`IupacError::NotSupported`].

#![forbid(unsafe_code)]

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;
use std::collections::{HashSet, VecDeque};

// ---------------------------------------------------------------------------
// Public error type
// ---------------------------------------------------------------------------

/// Error returned by [`name`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IupacError {
    /// The molecule contains no atoms.
    Empty,
    /// The molecule is outside the supported naming scope.
    NotSupported,
}

impl core::fmt::Display for IupacError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Empty => write!(f, "empty molecule"),
            Self::NotSupported => write!(f, "IUPAC name not supported for this structure"),
        }
    }
}

impl std::error::Error for IupacError {}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Generate a local IUPAC name for `mol`.
///
/// Returns `Err(IupacError::NotSupported)` for structures outside the current
/// scope (polycyclic systems, multi-functional groups, stereocenters, etc.).
pub fn name(mol: &Molecule) -> Result<String, IupacError> {
    if mol.atom_count() == 0 {
        return Err(IupacError::Empty);
    }
    Namer::new(mol).name()
}

// ---------------------------------------------------------------------------
// Internal namer
// ---------------------------------------------------------------------------

struct Namer<'a> {
    mol: &'a Molecule,
}

impl<'a> Namer<'a> {
    fn new(mol: &'a Molecule) -> Self {
        Self { mol }
    }

    fn name(&self) -> Result<String, IupacError> {
        let mol = self.mol;

        if count_components(mol) != 1 {
            return Err(IupacError::NotSupported);
        }

        let rings = find_sssr(mol);
        let ring_atoms: HashSet<AtomIdx> = rings
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();

        let carbons:  Vec<AtomIdx> = atoms_of(mol, 6);
        let o_atoms:  Vec<AtomIdx> = atoms_of(mol, 8);
        let n_atoms:  Vec<AtomIdx> = atoms_of(mol, 7);
        let s_atoms:  Vec<AtomIdx> = atoms_of(mol, 16);
        let halogens: Vec<AtomIdx> = mol
            .atoms()
            .filter(|(_, a)| matches!(a.element.atomic_number(), 9 | 17 | 35 | 53))
            .map(|(i, _)| i)
            .collect();

        // Reject elements outside C, H, N, O, S, halogens.
        let het_elements: HashSet<u8> = mol
            .atoms()
            .filter(|(_, a)| { let an = a.element.atomic_number(); an != 6 && an != 1 })
            .map(|(_, a)| a.element.atomic_number())
            .collect();
        if het_elements.iter().any(|&an| !matches!(an, 7 | 8 | 9 | 16 | 17 | 35 | 53)) {
            return Err(IupacError::NotSupported);
        }

        let cyclic = !ring_atoms.is_empty();

        if cyclic {
            let any_aromatic = ring_atoms.iter().any(|&i| mol.atom(i).aromatic);
            if any_aromatic {
                return self.name_aromatic_ring(&ring_atoms);
            }
            // Non-aromatic ring: only unsubstituted cycloalkanes.
            if !het_elements.is_empty() {
                return Err(IupacError::NotSupported);
            }
            return self.name_cycloalkane(&ring_atoms, &carbons);
        }

        // Acyclic dispatch on heteroatom composition.
        match (o_atoms.len(), n_atoms.len(), s_atoms.len(), halogens.len()) {
            (0, 0, 0, 0) => self.name_acyclic_hydrocarbon(&carbons),
            (1, 0, 0, 0) => self.name_one_oxygen(&carbons, o_atoms[0]),
            (2, 0, 0, 0) => self.name_two_oxygens(&carbons, &o_atoms),
            (1, 1, 0, 0) => self.name_amide(&carbons, o_atoms[0], n_atoms[0]),
            (0, 1, 0, 0) => {
                // Nitrile (C≡N) takes priority over amine.
                if self.is_nitrile(n_atoms[0]) {
                    self.name_nitrile(&carbons, n_atoms[0])
                } else {
                    self.name_amine(&carbons, n_atoms[0])
                }
            }
            (0, 0, 0, _) if !halogens.is_empty() => {
                if het_elements.len() != 1 {
                    return Err(IupacError::NotSupported);
                }
                let prefix = match het_elements.iter().next().copied().unwrap() {
                    9  => "fluoro",
                    17 => "chloro",
                    35 => "bromo",
                    53 => "iodo",
                    _  => return Err(IupacError::NotSupported),
                };
                self.name_haloalkane(&carbons, &halogens, prefix)
            }
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Aromatic ring naming
    // -----------------------------------------------------------------------

    fn name_aromatic_ring(&self, ring_atoms: &HashSet<AtomIdx>) -> Result<String, IupacError> {
        let mol = self.mol;
        // All ring atoms must be aromatic.
        if !ring_atoms.iter().all(|&i| mol.atom(i).aromatic) {
            return Err(IupacError::NotSupported);
        }

        let n_n = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 7).count();
        let n_o = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 8).count();
        let n_s = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 16).count();
        let sz  = ring_atoms.len();

        // Case 1: Pure aromatic ring (no substituents).
        if ring_atoms.len() == mol.atom_count() {
            return match (sz, n_n, n_o, n_s) {
                (6, 0, 0, 0) => Ok("benzene".into()),
                (6, 1, 0, 0) => Ok("pyridine".into()),
                (6, 2, 0, 0) => Ok("pyrimidine".into()),
                (5, 0, 1, 0) => Ok("furan".into()),
                (5, 0, 0, 1) => Ok("thiophene".into()),
                (5, 1, 0, 0) => Ok("pyrrole".into()),
                (5, 2, 0, 0) => Ok("imidazole".into()),
                _            => Err(IupacError::NotSupported),
            };
        }

        // Case 2: Monosubstituted benzene (phenol, toluene, aniline, etc.)
        // Only support pure benzene ring (6 C, no N/O/S in ring).
        if sz == 6 && n_n == 0 && n_o == 0 && n_s == 0 {
            let sub_atoms: Vec<AtomIdx> = mol.atoms()
                .filter(|(i, _)| !ring_atoms.contains(i))
                .map(|(i, _)| i)
                .collect();
            return self.name_monosubstituted_benzene(ring_atoms, &sub_atoms);
        }

        Err(IupacError::NotSupported)
    }

    // -----------------------------------------------------------------------
    // Monosubstituted benzene naming
    // -----------------------------------------------------------------------

    fn name_monosubstituted_benzene(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        sub_atoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Ensure exactly one ring C is attached to the substituent.
        let attach_count = ring_atoms.iter().filter(|&&r| {
            mol.neighbors(r).any(|(nb, _)| !ring_atoms.contains(&nb))
        }).count();
        if attach_count != 1 {
            return Err(IupacError::NotSupported);
        }

        // Classify substituent by element counts + bond types.
        let mut n_c = 0usize; let mut n_n = 0usize;
        let mut n_o = 0usize; let mut n_hal = 0usize;
        let mut halogen_an = 0u8;
        for &a in sub_atoms {
            match mol.atom(a).element.atomic_number() {
                6  => n_c += 1,
                7  => n_n += 1,
                8  => n_o += 1,
                1  => {},
                an @ (9 | 17 | 35 | 53) => { n_hal += 1; halogen_an = an; }
                _  => return Err(IupacError::NotSupported),
            }
        }

        let sub_set: HashSet<AtomIdx> = sub_atoms.iter().copied().collect();
        let has_triple = mol.bonds().any(|(_, b)| {
            b.order == BondOrder::Triple
                && (sub_set.contains(&b.atom1) || sub_set.contains(&b.atom2))
        });
        let has_double = mol.bonds().any(|(_, b)| {
            b.order == BondOrder::Double
                && (sub_set.contains(&b.atom1) || sub_set.contains(&b.atom2))
        });

        match (n_c, n_n, n_o, n_hal, has_double, has_triple) {
            // Phenol: c1ccccc1O
            (0, 0, 1, 0, false, false) => Ok("phenol".into()),
            // Aniline: c1ccccc1N
            (0, 1, 0, 0, false, false) => Ok("aniline".into()),
            // Halo-benzenes
            (0, 0, 0, 1, false, false) => {
                let prefix = match halogen_an {
                    9 => "fluoro", 17 => "chloro", 35 => "bromo", 53 => "iodo", _ => return Err(IupacError::NotSupported),
                };
                Ok(format!("{prefix}benzene"))
            }
            // Toluene: c1ccccc1C (one CH3)
            (1, 0, 0, 0, false, false) => Ok("toluene".into()),
            // Benzaldehyde: c1ccccc1C=O (n_c=1, n_o=1, has_double)
            (1, 0, 1, 0, true, false) => Ok("benzaldehyde".into()),
            // Benzoic acid: c1ccccc1C(=O)O (n_c=1, n_o=2, has_double)
            (1, 0, 2, 0, true, false) => Ok("benzoic acid".into()),
            // Benzonitrile: c1ccccc1C#N (n_c=1, n_n=1, has_triple)
            (1, 1, 0, 0, false, true) => Ok("benzonitrile".into()),
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Cycloalkane naming
    // -----------------------------------------------------------------------

    fn name_cycloalkane(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
    ) -> Result<String, IupacError> {
        if ring_atoms.len() != carbons.len() {
            return Err(IupacError::NotSupported);
        }
        if carbons.iter().any(|&c| self.mol.atom(c).aromatic) {
            return Err(IupacError::NotSupported);
        }
        Ok(format!("cyclo{}", alkane_suffix(ring_atoms.len())))
    }

    // -----------------------------------------------------------------------
    // Acyclic hydrocarbon naming
    // -----------------------------------------------------------------------

    fn name_acyclic_hydrocarbon(&self, carbons: &[AtomIdx]) -> Result<String, IupacError> {
        let mol = self.mol;
        let n = carbons.len();

        let double_bonds = mol.bonds().filter(|(_, b)| b.order == BondOrder::Double).count();
        let triple_bonds = mol.bonds().filter(|(_, b)| b.order == BondOrder::Triple).count();
        if double_bonds > 1 || triple_bonds > 1 || (double_bonds > 0 && triple_bonds > 0) {
            return Err(IupacError::NotSupported);
        }

        // Check for branching.
        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        let is_branched = carbons.iter().any(|&c| {
            mol.neighbors(c).filter(|(nb, _)| c_set.contains(nb)).count() > 2
        });

        if is_branched {
            // Only saturated branched alkanes supported for now.
            if double_bonds > 0 || triple_bonds > 0 {
                return Err(IupacError::NotSupported);
            }
            return self.name_branched_alkane(carbons);
        }

        Ok(if triple_bonds == 1 {
            alkyne_suffix(n)
        } else if double_bonds == 1 {
            alkene_suffix(n)
        } else {
            alkane_suffix(n)
        })
    }

    // -----------------------------------------------------------------------
    // One-oxygen compound: alcohol / aldehyde / ketone
    // -----------------------------------------------------------------------

    fn name_one_oxygen(&self, carbons: &[AtomIdx], o_idx: AtomIdx) -> Result<String, IupacError> {
        let mol = self.mol;
        let is_double = mol.neighbors(o_idx).any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);

        if !is_double {
            // Alcohol: C–OH
            let n = carbons.len();
            return Ok(format!("{}anol", alkane_stem(n)));
        }

        // Carbonyl: find the C=O carbon.
        let carbonyl_c = mol
            .neighbors(o_idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        if implicit_hcount(mol, carbonyl_c) > 0 {
            // Aldehyde: terminal CHO.
            let n = carbons.len();
            return Ok(format!("{}anal", alkane_stem(n)));
        }

        // Ketone: internal C=O.  Determine chain length and position.
        let c_sides: Vec<AtomIdx> = mol
            .neighbors(carbonyl_c)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .collect();
        if c_sides.len() < 2 {
            return Err(IupacError::NotSupported);
        }
        let left  = count_c_chain(mol, c_sides[0], carbonyl_c);
        let right = count_c_chain(mol, c_sides[1], carbonyl_c);
        let n   = left + right + 1;
        let pos = left.min(right) + 1;
        Ok(format!("{}-{}-one", alkane_base(n), pos))
    }

    // -----------------------------------------------------------------------
    // Two-oxygen compound: carboxylic acid or ester
    // -----------------------------------------------------------------------

    fn name_two_oxygens(&self, carbons: &[AtomIdx], o_atoms: &[AtomIdx]) -> Result<String, IupacError> {
        let mol = self.mol;
        let o1 = o_atoms[0];
        let o2 = o_atoms[1];

        let o1_dbl = mol.neighbors(o1).any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);
        let o2_dbl = mol.neighbors(o2).any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);

        let (carbonyl_o, ester_o) = match (o1_dbl, o2_dbl) {
            (true, false) => (o1, o2),
            (false, true) => (o2, o1),
            _ => return Err(IupacError::NotSupported),
        };

        // Carbonyl C is bonded to the =O oxygen.
        let carbonyl_c = mol
            .neighbors(carbonyl_o)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        // Carbonyl C must also be bonded to the single-bond O.
        if !mol.neighbors(carbonyl_c).any(|(nb, _)| nb == ester_o) {
            return Err(IupacError::NotSupported);
        }

        // Is the single-bond O also bonded to another C (→ ester) or only H (→ acid)?
        let alcohol_c = mol
            .neighbors(ester_o)
            .filter(|(nb, _)| *nb != carbonyl_c && mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next();

        if let Some(alc_c) = alcohol_c {
            // Ester: alkyl alkanoate
            let acid_n    = count_c_chain(mol, carbonyl_c, ester_o);
            let alcohol_n = count_c_chain(mol, alc_c, ester_o);
            Ok(format!("{}yl {}anoate", alkane_stem(alcohol_n), alkane_stem(acid_n)))
        } else {
            // Carboxylic acid
            let n = carbons.len();
            Ok(format!("{}anoic acid", alkane_stem(n)))
        }
    }

    // -----------------------------------------------------------------------
    // Amide: C(=O)–N
    // -----------------------------------------------------------------------

    fn name_amide(
        &self,
        _carbons: &[AtomIdx],
        o_idx: AtomIdx,
        n_idx: AtomIdx,
    ) -> Result<String, IupacError> {
        let mol = self.mol;

        // O must be a carbonyl (C=O).
        if !mol.neighbors(o_idx).any(|(_, bi)| mol.bond(bi).order == BondOrder::Double) {
            return Err(IupacError::NotSupported);
        }

        let carbonyl_c = mol
            .neighbors(o_idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;

        // Carbonyl C must be bonded to N.
        if !mol.neighbors(carbonyl_c).any(|(nb, _)| nb == n_idx) {
            return Err(IupacError::NotSupported);
        }

        // Only primary/secondary amides (N has ≥ 1 H).
        if implicit_hcount(mol, n_idx) == 0 {
            return Err(IupacError::NotSupported);
        }

        // Acid chain length: all Cs reachable from carbonyl_c not through N.
        let n_carbons = count_c_chain(mol, carbonyl_c, n_idx);
        Ok(format!("{}anamide", alkane_stem(n_carbons)))
    }

    // -----------------------------------------------------------------------
    // Amine naming
    // -----------------------------------------------------------------------

    fn name_amine(&self, carbons: &[AtomIdx], n_idx: AtomIdx) -> Result<String, IupacError> {
        let mol = self.mol;
        let n_h = implicit_hcount(mol, n_idx);
        let n   = carbons.len();
        let base = alkane_stem(n);
        match n_h {
            2 => Ok(format!("{base}an-1-amine")),
            1 => Ok(format!("di{base}ylamine")),
            0 => Ok(format!("tri{base}ylamine")),
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Haloalkane naming
    // -----------------------------------------------------------------------

    fn name_haloalkane(
        &self,
        carbons: &[AtomIdx],
        halogen_atoms: &[AtomIdx],
        prefix: &str,
    ) -> Result<String, IupacError> {
        let n     = carbons.len();
        let base  = alkane_suffix(n);
        let count = halogen_atoms.len();
        let mult  = match count {
            1 => prefix.to_string(),
            2 => format!("di{prefix}"),
            3 => format!("tri{prefix}"),
            _ => return Err(IupacError::NotSupported),
        };
        Ok(format!("{mult}{base}"))
    }

    // -----------------------------------------------------------------------
    // Nitrile naming (R-C≡N → "...nitrile")
    // -----------------------------------------------------------------------

    fn is_nitrile(&self, n_idx: AtomIdx) -> bool {
        self.mol.neighbors(n_idx)
            .any(|(_, bi)| self.mol.bond(bi).order == BondOrder::Triple)
    }

    fn name_nitrile(&self, carbons: &[AtomIdx], n_idx: AtomIdx) -> Result<String, IupacError> {
        let mol = self.mol;
        // Find the C≡N carbon.
        let nitrile_c = mol.neighbors(n_idx)
            .filter(|(_, bi)| mol.bond(*bi).order == BondOrder::Triple)
            .map(|(nb, _)| nb)
            .next()
            .ok_or(IupacError::NotSupported)?;
        // Count the total C chain (nitrile C + alkyl chain).
        // count_c_chain gives all C reachable from nitrile_c without crossing N.
        let n_carbons = count_c_chain(mol, nitrile_c, n_idx);
        // n_carbons already includes the nitrile carbon itself.
        if n_carbons == 0 { return Err(IupacError::NotSupported); }
        // Verify no branching on the C chain
        let c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();
        for &c in carbons {
            if mol.neighbors(c)
                .filter(|(nb, _)| c_set.contains(nb))
                .count() > 2
            {
                return Err(IupacError::NotSupported); // branched nitrile not supported
            }
        }
        Ok(format!("{}enitrile", alkane_base(n_carbons)))
    }

    // -----------------------------------------------------------------------
    // Branched alkane naming (e.g., "2-methylpropane", "2,2-dimethylpropane")
    // -----------------------------------------------------------------------

    fn name_branched_alkane(&self, carbons: &[AtomIdx]) -> Result<String, IupacError> {
        let mol = self.mol;

        // Find the principal chain (longest C–C path).
        let chain = find_longest_c_chain(mol, carbons);
        let n = chain.len();
        if n < 2 {
            return Err(IupacError::NotSupported);
        }

        let chain_set: std::collections::HashSet<AtomIdx> = chain.iter().copied().collect();
        let all_c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();

        // Collect substituents: (chain_position_1based, alkyl_length).
        let mut subs: Vec<(usize, usize)> = Vec::new();
        for (pos0, &chain_c) in chain.iter().enumerate() {
            let position = pos0 + 1;
            for (nb, _) in mol.neighbors(chain_c) {
                if all_c_set.contains(&nb) && !chain_set.contains(&nb) {
                    // Substituent rooted at `nb`, blocked by chain_c.
                    let sub_len = count_c_chain(mol, nb, chain_c);
                    // Only support methyl (1) and ethyl (2) substituents.
                    if sub_len > 2 {
                        return Err(IupacError::NotSupported);
                    }
                    subs.push((position, sub_len));
                }
            }
        }

        if subs.is_empty() {
            return Err(IupacError::NotSupported);
        }

        // Apply IUPAC lowest-locant rule: try forward and reverse numbering.
        let subs_rev: Vec<(usize, usize)> = subs.iter()
            .map(|&(pos, len)| (n + 1 - pos, len))
            .collect();

        let first_fwd = subs.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
        let first_rev = subs_rev.iter().map(|&(p, _)| p).min().unwrap_or(usize::MAX);
        let best_subs = if first_fwd <= first_rev { subs } else { subs_rev };

        Ok(format!(
            "{}{}",
            format_substituents(&best_subs),
            alkane_suffix(n)
        ))
    }
}

// ---------------------------------------------------------------------------
// Graph helpers
// ---------------------------------------------------------------------------

fn atoms_of(mol: &Molecule, atomic_num: u8) -> Vec<AtomIdx> {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() == atomic_num)
        .map(|(i, _)| i)
        .collect()
}

/// BFS count of C atoms reachable from `start` without crossing `blocked`.
fn count_c_chain(mol: &Molecule, start: AtomIdx, blocked: AtomIdx) -> usize {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    visited.insert(start);
    queue.push_back(start);
    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if nb == blocked { continue; }
            if mol.atom(nb).element.atomic_number() == 6 && visited.insert(nb) {
                queue.push_back(nb);
            }
        }
    }
    visited.len()
}

/// Find the longest carbon chain in a C-subgraph using two-pass BFS.
///
/// Returns the sequence of AtomIdx forming the longest simple path.
/// For branched alkanes this gives the principal chain (IUPAC rule: longest chain).
fn find_longest_c_chain(mol: &Molecule, carbons: &[AtomIdx]) -> Vec<AtomIdx> {
    if carbons.is_empty() { return Vec::new(); }

    let c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();

    // BFS to find the farthest atom from a given start, returning (farthest, parents).
    let bfs_far = |start: AtomIdx| -> (AtomIdx, std::collections::HashMap<AtomIdx, AtomIdx>) {
        let mut parent: std::collections::HashMap<AtomIdx, AtomIdx> = std::collections::HashMap::new();
        let mut visited: std::collections::HashSet<AtomIdx> = std::collections::HashSet::new();
        let mut queue = VecDeque::new();
        let mut farthest = start;
        visited.insert(start);
        queue.push_back(start);
        while let Some(cur) = queue.pop_front() {
            farthest = cur;
            for (nb, _) in mol.neighbors(cur) {
                if c_set.contains(&nb) && visited.insert(nb) {
                    parent.insert(nb, cur);
                    queue.push_back(nb);
                }
            }
        }
        (farthest, parent)
    };

    let reconstruct = |end: AtomIdx, start: AtomIdx,
                        parents: &std::collections::HashMap<AtomIdx, AtomIdx>| -> Vec<AtomIdx> {
        let mut path = vec![end];
        let mut cur = end;
        while cur != start {
            cur = parents[&cur];
            path.push(cur);
        }
        path.reverse();
        path
    };

    // Pass 1: BFS from first carbon to find one endpoint of the longest chain.
    let (end1, _) = bfs_far(carbons[0]);
    // Pass 2: BFS from end1 to find the other endpoint.
    let (end2, parents) = bfs_far(end1);

    reconstruct(end2, end1, &parents)
}

/// Format substituents as an IUPAC prefix string ("2-methyl", "2,2-dimethyl", etc.).
fn format_substituents(subs: &[(usize, usize)]) -> String {
    // Group by alkyl name; sort alphabetically.
    let mut groups: std::collections::BTreeMap<&str, Vec<usize>> =
        std::collections::BTreeMap::new();
    for &(pos, len) in subs {
        let alkyl = match len {
            1 => "methyl",
            2 => "ethyl",
            3 => "propyl",
            _ => continue,
        };
        groups.entry(alkyl).or_default().push(pos);
    }

    let mut parts: Vec<String> = Vec::new();
    for (alkyl, mut positions) in groups {
        positions.sort_unstable();
        let locants = positions.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(",");
        let mult = match positions.len() {
            1 => String::new(),
            2 => "di".to_string(),
            3 => "tri".to_string(),
            _ => "?".to_string(),
        };
        parts.push(format!("{}-{}{}", locants, mult, alkyl));
    }
    parts.join("-")
}

fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    if n == 0 { return 0; }
    let mut visited = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if visited[start] { continue; }
        count += 1;
        let mut queue = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));
        visited[start] = true;
        while let Some(cur) = queue.pop_front() {
            for (nb, _) in mol.neighbors(cur) {
                if !visited[nb.0 as usize] {
                    visited[nb.0 as usize] = true;
                    queue.push_back(nb);
                }
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Naming helpers
// ---------------------------------------------------------------------------

fn alkane_stem(n: usize) -> &'static str {
    match n {
        1 => "meth", 2 => "eth",  3 => "prop", 4 => "but",
        5 => "pent", 6 => "hex",  7 => "hept", 8 => "oct",
        9 => "non",  10 => "dec", _ => "long",
    }
}

/// Stem with "an" appended — base for most suffix compounds.
fn alkane_base(n: usize) -> String {
    format!("{}an", alkane_stem(n))
}

fn alkane_suffix(n: usize) -> String {
    match n {
        1  => "methane".into(),   2  => "ethane".into(),
        3  => "propane".into(),   4  => "butane".into(),
        5  => "pentane".into(),   6  => "hexane".into(),
        7  => "heptane".into(),   8  => "octane".into(),
        9  => "nonane".into(),    10 => "decane".into(),
        11 => "undecane".into(),  12 => "dodecane".into(),
        13 => "tridecane".into(), 14 => "tetradecane".into(),
        15 => "pentadecane".into(), 16 => "hexadecane".into(),
        17 => "heptadecane".into(), 18 => "octadecane".into(),
        19 => "nonadecane".into(), 20 => "icosane".into(),
        _  => format!("{n}alkane"),
    }
}

fn alkene_suffix(n: usize) -> String { alkane_suffix(n).replace("ane", "ene") }
fn alkyne_suffix(n: usize) -> String { alkane_suffix(n).replace("ane", "yne") }

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule { parse(s).unwrap() }

    // --- Existing tests (must remain green) ---------------------------------

    #[test]
    fn test_alkanes() {
        assert_eq!(name(&mol("C")).unwrap(),      "methane");
        assert_eq!(name(&mol("CC")).unwrap(),     "ethane");
        assert_eq!(name(&mol("CCC")).unwrap(),    "propane");
        assert_eq!(name(&mol("CCCC")).unwrap(),   "butane");
        assert_eq!(name(&mol("CCCCC")).unwrap(),  "pentane");
        assert_eq!(name(&mol("CCCCCC")).unwrap(), "hexane");
    }

    #[test]
    fn test_alkenes_alkynes() {
        assert_eq!(name(&mol("C=C")).unwrap(),   "ethene");
        assert_eq!(name(&mol("CC=C")).unwrap(),  "propene");
        assert_eq!(name(&mol("C#C")).unwrap(),   "ethyne");
        assert_eq!(name(&mol("CC#C")).unwrap(),  "propyne");
    }

    #[test]
    fn test_cycloalkanes() {
        assert_eq!(name(&mol("C1CC1")).unwrap(),   "cyclopropane");
        assert_eq!(name(&mol("C1CCC1")).unwrap(),  "cyclobutane");
        assert_eq!(name(&mol("C1CCCC1")).unwrap(), "cyclopentane");
        assert_eq!(name(&mol("C1CCCCC1")).unwrap(),"cyclohexane");
    }

    #[test]
    fn test_alcohol() {
        assert_eq!(name(&mol("CO")).unwrap(),   "methanol");
        assert_eq!(name(&mol("CCO")).unwrap(),  "ethanol");
        assert_eq!(name(&mol("CCCO")).unwrap(), "propanol");
    }

    #[test]
    fn test_amine() {
        assert_eq!(name(&mol("CN")).unwrap(),  "methan-1-amine");
        assert_eq!(name(&mol("CCN")).unwrap(), "ethan-1-amine");
    }

    #[test]
    fn test_haloalkane() {
        assert_eq!(name(&mol("CCCl")).unwrap(), "chloroethane");
        assert_eq!(name(&mol("CCBr")).unwrap(), "bromoethane");
        assert_eq!(name(&mol("CF")).unwrap(),   "fluoromethane");
        assert_eq!(name(&mol("CI")).unwrap(),   "iodomethane");
    }

    #[test]
    fn test_not_supported() {
        assert!(name(&mol("CC.CC")).is_err());  // disconnected
    }

    #[test]
    fn test_empty() {
        use chematic_core::MoleculeBuilder;
        let mol = MoleculeBuilder::new().build();
        assert_eq!(name(&mol), Err(IupacError::Empty));
    }

    // --- New: benzene & aromatic heterocycles --------------------------------

    #[test]
    fn test_benzene() {
        assert_eq!(name(&mol("c1ccccc1")).unwrap(), "benzene");
    }

    #[test]
    fn test_aromatic_heterocycles() {
        assert_eq!(name(&mol("c1ccncc1")).unwrap(),   "pyridine");
        assert_eq!(name(&mol("c1ccoc1")).unwrap(),    "furan");
        assert_eq!(name(&mol("c1ccsc1")).unwrap(),    "thiophene");
        assert_eq!(name(&mol("c1cc[nH]c1")).unwrap(), "pyrrole");
        assert_eq!(name(&mol("c1cnc[nH]1")).unwrap(), "imidazole");
    }

    // --- New: ketones with position locant -----------------------------------

    #[test]
    fn test_ketones() {
        assert_eq!(name(&mol("CC(=O)C")).unwrap(),    "propan-2-one");
        assert_eq!(name(&mol("CC(=O)CC")).unwrap(),   "butan-2-one");
        assert_eq!(name(&mol("CCC(=O)CC")).unwrap(),  "pentan-3-one");
        assert_eq!(name(&mol("CCCC(=O)C")).unwrap(),  "pentan-2-one");
    }

    // --- New: carboxylic acids -----------------------------------------------

    #[test]
    fn test_carboxylic_acids() {
        assert_eq!(name(&mol("CC(=O)O")).unwrap(),  "ethanoic acid");
        assert_eq!(name(&mol("CCC(=O)O")).unwrap(), "propanoic acid");
        assert_eq!(name(&mol("C(=O)O")).unwrap(),   "methanoic acid");
    }

    // --- New: esters ---------------------------------------------------------

    #[test]
    fn test_esters() {
        assert_eq!(name(&mol("CC(=O)OC")).unwrap(),  "methyl ethanoate");
        assert_eq!(name(&mol("C(=O)OC")).unwrap(),   "methyl methanoate");
        assert_eq!(name(&mol("CC(=O)OCC")).unwrap(), "ethyl ethanoate");
    }

    // --- New: amides ---------------------------------------------------------

    #[test]
    fn test_amides() {
        assert_eq!(name(&mol("CC(=O)N")).unwrap(),   "ethanamide");
        assert_eq!(name(&mol("C(=O)N")).unwrap(),    "methanamide");
        assert_eq!(name(&mol("CCC(=O)N")).unwrap(),  "propanamide");
    }

    // ---- New: branched alkanes (v0.1.101) ------------------------------------

    #[test]
    fn test_branched_alkanes() {
        assert_eq!(name(&mol("CC(C)C")).unwrap(),    "2-methylpropane");
        assert_eq!(name(&mol("CC(C)CC")).unwrap(),   "2-methylbutane");
        assert_eq!(name(&mol("CC(C)(C)C")).unwrap(), "2,2-dimethylpropane");
        assert_eq!(name(&mol("CCCC(C)CC")).unwrap(), "3-methylhexane");
    }

    #[test]
    fn test_branched_alkane_lowest_locant() {
        // CCC(C)C = 2-methylbutane (not 3-methylbutane — lower locant wins).
        assert_eq!(name(&mol("CCC(C)C")).unwrap(), "2-methylbutane");
    }

    // ---- New: substituted benzenes (v0.1.101) --------------------------------

    #[test]
    fn test_substituted_benzenes() {
        assert_eq!(name(&mol("c1ccccc1O")).unwrap(),     "phenol");
        assert_eq!(name(&mol("c1ccccc1N")).unwrap(),     "aniline");
        assert_eq!(name(&mol("c1ccccc1Cl")).unwrap(),    "chlorobenzene");
        assert_eq!(name(&mol("c1ccccc1Br")).unwrap(),    "bromobenzene");
    }

    #[test]
    fn test_substituted_benzene_carbonyl() {
        assert_eq!(name(&mol("c1ccccc1C=O")).unwrap(),        "benzaldehyde");
        assert_eq!(name(&mol("c1ccccc1C(=O)O")).unwrap(),     "benzoic acid");
    }

    // ---- New: nitriles (v0.1.101) -------------------------------------------

    #[test]
    fn test_nitriles() {
        assert_eq!(name(&mol("CC#N")).unwrap(),  "ethanenitrile");
        assert_eq!(name(&mol("CCC#N")).unwrap(), "propanenitrile");
    }
}
