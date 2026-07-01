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

use chematic_core::{AtomIdx, Molecule};
use std::collections::HashSet;

mod acyclic;
mod helpers;
mod heteroatoms;
mod oxygen;
mod rings;
#[cfg(test)]
mod tests;

use helpers::atoms_of;

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
//
// Compound-class handlers live in sibling modules (rings/oxygen/heteroatoms/
// acyclic), each contributing an `impl<'a> Namer<'a> { ... }` block — legal
// since Rust allows multiple inherent impl blocks for the same type across
// files within one crate, and since those modules are descendants of the
// crate root they can see this private `Namer` struct.
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

        if helpers::count_components(mol) != 1 {
            return Err(IupacError::NotSupported);
        }

        let rings = chematic_perception::find_sssr(mol);
        let ring_atoms: HashSet<AtomIdx> = rings
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();

        let carbons: Vec<AtomIdx> = atoms_of(mol, 6);
        let o_atoms: Vec<AtomIdx> = atoms_of(mol, 8);
        let n_atoms: Vec<AtomIdx> = atoms_of(mol, 7);
        let s_atoms: Vec<AtomIdx> = atoms_of(mol, 16);
        let halogens: Vec<AtomIdx> = mol
            .atoms()
            .filter(|(_, a)| matches!(a.element.atomic_number(), 9 | 17 | 35 | 53))
            .map(|(i, _)| i)
            .collect();

        // Reject elements outside C, H, N, O, S, halogens.
        let het_elements: HashSet<u8> = mol
            .atoms()
            .filter(|(_, a)| {
                let an = a.element.atomic_number();
                an != 6 && an != 1
            })
            .map(|(_, a)| a.element.atomic_number())
            .collect();
        if het_elements
            .iter()
            .any(|&an| !matches!(an, 7 | 8 | 9 | 16 | 17 | 35 | 53))
        {
            return Err(IupacError::NotSupported);
        }

        let cyclic = !ring_atoms.is_empty();

        if cyclic {
            let any_aromatic = ring_atoms.iter().any(|&i| mol.atom(i).aromatic);
            if any_aromatic {
                return self.name_aromatic_ring(&ring_atoms);
            }
            let only_oxygen = het_elements.len() == 1 && het_elements.contains(&8);
            let only_nitrogen = het_elements.len() == 1 && het_elements.contains(&7);
            let n_and_o =
                het_elements.len() == 2 && het_elements.contains(&7) && het_elements.contains(&8);
            // Non-aromatic ring: supported heteroatom patterns only.
            if !het_elements.is_empty() && !only_oxygen && !only_nitrogen && !n_and_o {
                return Err(IupacError::NotSupported);
            }
            if only_oxygen {
                return self.name_cycloalkanol(&ring_atoms, &carbons, &o_atoms);
            }
            if only_nitrogen {
                return self.name_aza_ring(&ring_atoms);
            }
            if n_and_o {
                return self.name_oxaaza_ring(&ring_atoms);
            }
            // Check for spiro / bridged polycyclic systems (all-carbon only).
            if het_elements.is_empty()
                && let Ok(name) = self.name_polycyclic(&ring_atoms, &carbons)
            {
                return Ok(name);
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
            (0, 0, 1, 0) => {
                if chematic_core::implicit_hcount(self.mol, s_atoms[0]) > 0 {
                    self.name_thiol(&carbons, s_atoms[0])
                } else {
                    self.name_sulfide(&carbons, s_atoms[0])
                }
            }
            (0, 0, 0, _) if !halogens.is_empty() => {
                if het_elements.len() != 1 {
                    return Err(IupacError::NotSupported);
                }
                let prefix = match het_elements.iter().next().copied().unwrap() {
                    9 => "fluoro",
                    17 => "chloro",
                    35 => "bromo",
                    53 => "iodo",
                    _ => return Err(IupacError::NotSupported),
                };
                self.name_haloalkane(&carbons, &halogens, prefix)
            }
            _ => Err(IupacError::NotSupported),
        }
    }
}
