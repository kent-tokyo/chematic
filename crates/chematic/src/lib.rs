#![forbid(unsafe_code)]
//! # chematic — Pure-Rust Cheminformatics
//!
//! **chematic** is a comprehensive cheminformatics toolkit written entirely in
//! Rust with zero C/C++ dependencies.  It runs in browsers via WebAssembly,
//! in serverless functions, and in native Rust applications — all from the
//! same codebase.
//!
//! ## Feature highlights
//!
//! | Capability | Crate | Key API |
//! |---|---|---|
//! | SMILES parsing & writing | `smiles` | [`smiles::parse`] |
//! | SMARTS substructure search | `smarts` | [`smarts::find_matches`] |
//! | Maximum common substructure | `smarts` | [`smarts::find_mcs`] |
//! | 2D SVG depiction | `depict` | [`depict::depict_svg`] |
//! | SDF / MOL V2000 + V3000 I/O | `mol` | [`mol::parse_mol`], [`mol::parse_mol_v3000`] |
//! | CIP stereochemistry (R/S, E/Z) | `perception` | [`perception::assign_stereo_from_2d`] |
//! | Fingerprints (ECFP4/6, FCFP, MACCS, AtomPair) | `fp` | [`fp::ecfp4`], [`fp::maccs`] |
//! | Tanimoto / Dice similarity | `fp` | [`fp::BitVec2048::tanimoto`] |
//! | Molecular descriptors (LogP, TPSA, MW, …) | `chem` | [`chem::logp_crippen`], [`chem::tpsa`] |
//! | Exact / monoisotopic mass | `chem` | [`chem::exact_mass`] |
//! | Isotope distribution | `chem` | [`chem::isotope_distribution`] |
//! | Bemis–Murcko scaffold | `chem` | [`chem::murcko_scaffold`] |
//! | Structure standardization & salt stripping | `chem` | [`chem::standardize`], [`chem::largest_fragment`] |
//! | Canonical tautomers | `chem` | [`chem::canonical_tautomer`] |
//! | Reaction SMILES parsing | `rxn` | [`rxn::parse_reaction`] |
//! | SMIRKS template application | `rxn` | [`rxn::run_reactants`] |
//! | InChI / InChIKey generation | `inchi` | [`inchi::inchi`], [`inchi::inchi_key`] |
//! | 3D coordinate generation | `threed` | [`threed::generate_coords`] |
//! | XYZ / PDB I/O | `threed` | [`threed::parse_xyz`], [`threed::parse_pdb_atoms`] |
//! | UFF-derived geometry minimization | `threed` | [`threed::minimize`] |
//!
//! ## Quick start
//!
//! Add to `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! chematic = { version = "0.4", features = ["full"] }
//! ```
//!
//! Parse a molecule and compute properties:
//!
//! ```rust,ignore
//! use chematic::smiles;
//! use chematic::chem;
//! use chematic::depict;
//!
//! let mol = smiles::parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
//! println!("MW:   {:.3}", chem::molecular_weight(&mol));
//! println!("LogP: {:.3}", chem::logp_crippen(&mol));
//! println!("TPSA: {:.1}", chem::tpsa(&mol));
//! let svg = depict::depict_svg(&mol);
//! ```
//!
//! ## Feature flags
//!
//! | Flag | Includes |
//! |---|---|
//! | `smiles` | SMILES parser + core molecule graph |
//! | `perception` | Aromaticity, ring perception, stereochemistry |
//! | `mol` | SDF / MOL V2000 + V3000 read/write |
//! | `depict` | 2D layout + SVG rendering |
//! | `fp` | ECFP, FCFP, MACCS, AtomPair fingerprints |
//! | `chem` | Descriptors, LogP, TPSA, scaffolds, standardization |
//! | `smarts` | SMARTS parser, substructure search, MCS |
//! | `rxn` | Reaction SMILES, SMIRKS transforms |
//! | `inchi` | InChI and InChIKey generation |
//! | `threed` | 3D coordinates, UFF minimization, XYZ/PDB I/O |
//! | `iupac` | IUPAC nomenclature for simple molecules |
//! | `crystal` | Periodic crystal structures (lattice, PBC, neighbors, supercells, POSCAR/CONTCAR I/O) -- see `chematic-crystal`'s README |
//! | `full` | All of the above |

#[cfg(feature = "threed")]
pub use chematic_3d as threed;
#[cfg(feature = "chem")]
pub use chematic_chem as chem;
#[cfg(feature = "smiles")]
pub use chematic_core as core;
#[cfg(feature = "crystal")]
pub use chematic_crystal as crystal;
#[cfg(feature = "depict")]
pub use chematic_depict as depict;
#[cfg(feature = "fp")]
pub use chematic_fp as fp;
#[cfg(feature = "inchi")]
pub use chematic_inchi as inchi;
#[cfg(feature = "iupac")]
pub use chematic_iupac as iupac;
#[cfg(feature = "mol")]
pub use chematic_mol as mol;
#[cfg(feature = "perception")]
pub use chematic_perception as perception;
#[cfg(feature = "rxn")]
pub use chematic_rxn as rxn;
#[cfg(feature = "smarts")]
pub use chematic_smarts as smarts;
#[cfg(feature = "smiles")]
pub use chematic_smiles as smiles;

#[cfg(test)]
mod tests {
    #[test]
    fn test_crate_compiles() {
        // Umbrella crate compiles with no features enabled.
    }
}
