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
            (0, 1, 0, 0) => self.name_amine(&carbons, n_atoms[0]),
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
        // All atoms must be in the ring (no substituents).
        if ring_atoms.len() != mol.atom_count() {
            return Err(IupacError::NotSupported);
        }
        // All ring atoms must be aromatic.
        if !ring_atoms.iter().all(|&i| mol.atom(i).aromatic) {
            return Err(IupacError::NotSupported);
        }
        let n_n = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 7).count();
        let n_o = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 8).count();
        let n_s = ring_atoms.iter().filter(|&&i| mol.atom(i).element.atomic_number() == 16).count();
        let sz  = ring_atoms.len();
        match (sz, n_n, n_o, n_s) {
            (6, 0, 0, 0) => Ok("benzene".into()),
            (6, 1, 0, 0) => Ok("pyridine".into()),
            (6, 2, 0, 0) => Ok("pyrimidine".into()),
            (5, 0, 1, 0) => Ok("furan".into()),
            (5, 0, 0, 1) => Ok("thiophene".into()),
            (5, 1, 0, 0) => Ok("pyrrole".into()),
            (5, 2, 0, 0) => Ok("imidazole".into()),
            _            => Err(IupacError::NotSupported),
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

        // Unbranched chain: each C has ≤ 2 C neighbors.
        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        for &c in carbons {
            if mol.neighbors(c).filter(|(nb, _)| c_set.contains(nb)).count() > 2 {
                return Err(IupacError::NotSupported);
            }
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
}
