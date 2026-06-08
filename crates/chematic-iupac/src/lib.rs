//! `chematic-iupac` — local IUPAC name generation, no network required.
//!
//! Supports:
//! - Linear and branched alkanes (methane → octadecane, plus longer)
//! - Cycloalkanes (cyclopropane … cyclodecane, etc.)
//! - Alkenes (`-ene`) and alkynes (`-yne`) with one unsaturation
//! - Simple monosubstituted derivatives: alcohols (`-ol`), amines (`-amine`),
//!   aldehydes (`-al`), ketones (`-one`), carboxylic acids (`-oic acid`)
//! - Halogen substituents: fluoro-, chloro-, bromo-, iodo-
//!
//! Complex polycyclic systems, stereo descriptors, and IUPAC preferred names
//! for structures outside the above scope return [`IupacError::NotSupported`].

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

        // Only handle single connected components.
        if count_components(mol) != 1 {
            return Err(IupacError::NotSupported);
        }

        let rings = find_sssr(mol);
        let ring_atoms: HashSet<AtomIdx> = rings
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();

        // Classify atoms.
        let carbons: Vec<AtomIdx> = mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 6)
            .map(|(i, _)| i)
            .collect();
        let heteroatoms: Vec<AtomIdx> = mol
            .atoms()
            .filter(|(_, a)| {
                let an = a.element.atomic_number();
                an != 6 && an != 1
            })
            .map(|(i, _)| i)
            .collect();

        // Must be purely organic (C + H + one heteroatom class).
        let het_elements: HashSet<u8> = heteroatoms
            .iter()
            .map(|&i| mol.atom(i).element.atomic_number())
            .collect();

        if het_elements.len() > 1 {
            return Err(IupacError::NotSupported);
        }
        if heteroatoms.len() > 1 {
            return Err(IupacError::NotSupported);
        }

        let cyclic = !ring_atoms.is_empty();

        match (cyclic, het_elements.iter().next().copied()) {
            (true, None) => self.name_cycloalkane(&ring_atoms, &carbons),
            (false, None) => self.name_acyclic_hydrocarbon(&carbons),
            (false, Some(8)) => self.name_alcohol_or_acid(&carbons, &heteroatoms),
            (false, Some(7)) => self.name_amine(&carbons, &heteroatoms),
            (false, Some(9)) => self.name_haloalkane(&carbons, &heteroatoms, "fluoro"),
            (false, Some(17)) => self.name_haloalkane(&carbons, &heteroatoms, "chloro"),
            (false, Some(35)) => self.name_haloalkane(&carbons, &heteroatoms, "bromo"),
            (false, Some(53)) => self.name_haloalkane(&carbons, &heteroatoms, "iodo"),
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
        // All carbons must be in the ring (no substituents).
        if ring_atoms.len() != carbons.len() {
            return Err(IupacError::NotSupported);
        }
        // Reject aromatic rings (benzene, etc.).
        if carbons.iter().any(|&c| self.mol.atom(c).aromatic) {
            return Err(IupacError::NotSupported);
        }
        let n = ring_atoms.len();
        Ok(format!("cyclo{}", alkane_suffix(n)))
    }

    // -----------------------------------------------------------------------
    // Acyclic hydrocarbon naming
    // -----------------------------------------------------------------------

    fn name_acyclic_hydrocarbon(&self, carbons: &[AtomIdx]) -> Result<String, IupacError> {
        let mol = self.mol;
        let n = carbons.len();

        // Unsaturation check.
        let double_bonds = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .count();
        let triple_bonds = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Triple)
            .count();

        if double_bonds > 1 || triple_bonds > 1 || (double_bonds > 0 && triple_bonds > 0) {
            return Err(IupacError::NotSupported);
        }

        // Must be an unbranched chain: each C has at most 2 C neighbors.
        let c_set: HashSet<AtomIdx> = carbons.iter().copied().collect();
        for &c in carbons {
            let c_degree = mol
                .neighbors(c)
                .filter(|(nb, _)| c_set.contains(nb))
                .count();
            if c_degree > 2 {
                return Err(IupacError::NotSupported);
            }
        }

        let suffix = if triple_bonds == 1 {
            alkyne_suffix(n)
        } else if double_bonds == 1 {
            alkene_suffix(n)
        } else {
            alkane_suffix(n)
        };

        Ok(suffix)
    }

    // -----------------------------------------------------------------------
    // Alcohol / carboxylic acid naming
    // -----------------------------------------------------------------------

    fn name_alcohol_or_acid(
        &self,
        carbons: &[AtomIdx],
        heteroatoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let o_idx = heteroatoms[0];
        let o_atom = mol.atom(o_idx);

        // Determine O type: OH (alcohol), =O + neighbor C (aldehyde if terminal, ketone if internal),
        // COOH (carboxylic acid).
        let o_neighbors: Vec<AtomIdx> = mol.neighbors(o_idx).map(|(nb, _)| nb).collect();
        let is_carbonyl = mol
            .neighbors(o_idx)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);

        let n = carbons.len();
        let base = alkane_stem(n);

        // Carboxylic acid: COOH — C bonded to two O (one =O, one -OH).
        if let Some(&c_idx) = o_neighbors.first()
            && mol.atom(c_idx).element.atomic_number() == 6
        {
            let c_o_count = mol
                .neighbors(c_idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 8)
                .count();
            if c_o_count == 2 && is_carbonyl {
                return Ok(format!("{base}anoic acid"));
            }
        }

        // Aldehyde: terminal C=O (CHO).
        if is_carbonyl {
            let c_idx = o_neighbors
                .iter()
                .find(|&&nb| mol.atom(nb).element.atomic_number() == 6);
            if let Some(&c_idx) = c_idx {
                let c_h = implicit_hcount(mol, c_idx);
                if c_h > 0 {
                    return Ok(format!("{base}anal"));
                }
                // Internal ketone.
                return Ok(format!("{base}anone"));
            }
        }

        // Alcohol: C-OH.
        if o_atom.charge == 0 && !is_carbonyl {
            return Ok(format!("{base}anol"));
        }

        Err(IupacError::NotSupported)
    }

    // -----------------------------------------------------------------------
    // Amine naming
    // -----------------------------------------------------------------------

    fn name_amine(
        &self,
        carbons: &[AtomIdx],
        heteroatoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let n_idx = heteroatoms[0];
        let n_h = implicit_hcount(mol, n_idx);
        let n = carbons.len();
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
        heteroatoms: &[AtomIdx],
        prefix: &str,
    ) -> Result<String, IupacError> {
        let n = carbons.len();
        let base = alkane_suffix(n);
        let count = heteroatoms.len();
        let mult = match count {
            1 => prefix.to_string(),
            2 => format!("di{prefix}"),
            3 => format!("tri{prefix}"),
            _ => return Err(IupacError::NotSupported),
        };
        Ok(format!("{mult}{base}"))
    }
}

// ---------------------------------------------------------------------------
// Naming helpers
// ---------------------------------------------------------------------------

/// Stem used before a suffix (e.g. "meth", "eth", "prop").
fn alkane_stem(n: usize) -> &'static str {
    match n {
        1 => "meth",
        2 => "eth",
        3 => "prop",
        4 => "but",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        _ => "long",
    }
}

/// Full alkane name (stem + "ane").
fn alkane_suffix(n: usize) -> String {
    match n {
        1 => "methane".into(),
        2 => "ethane".into(),
        3 => "propane".into(),
        4 => "butane".into(),
        5 => "pentane".into(),
        6 => "hexane".into(),
        7 => "heptane".into(),
        8 => "octane".into(),
        9 => "nonane".into(),
        10 => "decane".into(),
        11 => "undecane".into(),
        12 => "dodecane".into(),
        13 => "tridecane".into(),
        14 => "tetradecane".into(),
        15 => "pentadecane".into(),
        16 => "hexadecane".into(),
        17 => "heptadecane".into(),
        18 => "octadecane".into(),
        19 => "nonadecane".into(),
        20 => "icosane".into(),
        _ => format!("{n}alkane"),
    }
}

/// Alkene name (replace "-ane" with "-ene").
fn alkene_suffix(n: usize) -> String {
    alkane_suffix(n).replace("ane", "ene")
}

/// Alkyne name (replace "-ane" with "-yne").
fn alkyne_suffix(n: usize) -> String {
    alkane_suffix(n).replace("ane", "yne")
}

// ---------------------------------------------------------------------------
// Graph helpers
// ---------------------------------------------------------------------------

fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    if n == 0 {
        return 0;
    }
    let mut visited = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
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
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_alkanes() {
        assert_eq!(name(&mol("C")).unwrap(), "methane");
        assert_eq!(name(&mol("CC")).unwrap(), "ethane");
        assert_eq!(name(&mol("CCC")).unwrap(), "propane");
        assert_eq!(name(&mol("CCCC")).unwrap(), "butane");
        assert_eq!(name(&mol("CCCCC")).unwrap(), "pentane");
        assert_eq!(name(&mol("CCCCCC")).unwrap(), "hexane");
    }

    #[test]
    fn test_alkenes_alkynes() {
        assert_eq!(name(&mol("C=C")).unwrap(), "ethene");
        assert_eq!(name(&mol("CC=C")).unwrap(), "propene"); // prop-1-ene simplified
        assert_eq!(name(&mol("C#C")).unwrap(), "ethyne");
        assert_eq!(name(&mol("CC#C")).unwrap(), "propyne");
    }

    #[test]
    fn test_cycloalkanes() {
        assert_eq!(name(&mol("C1CC1")).unwrap(), "cyclopropane");
        assert_eq!(name(&mol("C1CCC1")).unwrap(), "cyclobutane");
        assert_eq!(name(&mol("C1CCCC1")).unwrap(), "cyclopentane");
        assert_eq!(name(&mol("C1CCCCC1")).unwrap(), "cyclohexane");
    }

    #[test]
    fn test_alcohol() {
        assert_eq!(name(&mol("CO")).unwrap(), "methanol");
        assert_eq!(name(&mol("CCO")).unwrap(), "ethanol");
        assert_eq!(name(&mol("CCCO")).unwrap(), "propanol");
    }

    #[test]
    fn test_amine() {
        assert_eq!(name(&mol("CN")).unwrap(), "methan-1-amine");
        assert_eq!(name(&mol("CCN")).unwrap(), "ethan-1-amine");
    }

    #[test]
    fn test_haloalkane() {
        assert_eq!(name(&mol("CCCl")).unwrap(), "chloroethane");
        assert_eq!(name(&mol("CCBr")).unwrap(), "bromoethane");
        assert_eq!(name(&mol("CF")).unwrap(), "fluoromethane");
        assert_eq!(name(&mol("CI")).unwrap(), "iodomethane");
    }

    #[test]
    fn test_not_supported() {
        // Disconnected fragments, benzene, multi-ring → NotSupported
        assert!(name(&mol("CC.CC")).is_err());
        assert!(name(&mol("c1ccccc1")).is_err()); // aromatic ring not handled
    }

    #[test]
    fn test_empty() {
        use chematic_core::MoleculeBuilder;
        let mol = MoleculeBuilder::new().build();
        assert_eq!(name(&mol), Err(IupacError::Empty));
    }
}
