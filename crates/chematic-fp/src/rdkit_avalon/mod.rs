//! RDKit's Avalon fingerprint (`rdkit.Avalon.pyAvalonTools.GetAvalonFP`),
//! bit for bit.
//!
//! RDKit computes it by writing the molecule as a MOL block
//! (`MolToMolBlock(mol, includeStereo=true)`: Kekulé bonds, explicit atoms
//! only), reading that back into the Avalon toolkit's `reaccs_molecule_t`
//! and calling the toolkit's `SetFingerprintBits` twice — once with its own
//! six-ring aromaticity and once with Daylight-style aromaticity — and
//! OR-ing the two results.
//!
//! This module is a port of that toolkit code (Avalon Cheminformatics
//! Toolkit `AvalonToolkit_2.0.5-pre.3`, the version RDKit 2026.03 builds,
//! BSD-3-Clause, Copyright Novartis Pharma AG; see `THIRD_PARTY_NOTICES.md`):
//! [`AvalonMolecule`] is the `reaccs_molecule_t` subset the fingerprint
//! reads, [`read_molblock`] fills it from a MOL block the way `MolStr2Mol`
//! does, and [`avalon_fp_bytes`] runs the fingerprint. [`rdkit_avalon_fp`]
//! builds the [`AvalonMolecule`] RDKit's MOL block would give for a chematic
//! [`Molecule`].

mod convert;
mod fingerprint;
mod molfile;
mod rings;

use chematic_core::Molecule;

pub use convert::{RdkitAvalonError, avalon_molecule_from_rdkit_view};
pub use molfile::read_molblock;

pub(crate) const SINGLE: i32 = 1;
pub(crate) const DOUBLE: i32 = 2;
pub(crate) const TRIPLE: i32 = 3;
pub(crate) const AROMATIC: i32 = 4;

/// RDKit's default `bitFlags` for `GetAvalonFP` (`avalonSimilarityBits`).
pub const RDKIT_AVALON_DEFAULT_BIT_FLAGS: u32 = 15_761_407;

/// An atom of an Avalon-toolkit molecule: the fields the fingerprint reads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AvalonAtom {
    /// MOL-file atom symbol (`"C"`, `"Cl"`, `"R"`, ...).
    pub symbol: String,
    /// Formal charge.
    pub charge: i32,
    /// MDL radical code (`2` = doublet).
    pub radical: i32,
    /// Atom text of an `"R"` atom (MOL-file `A` line), else empty.
    pub atext: String,
}

/// A bond of an Avalon-toolkit molecule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AvalonBond {
    /// 1-based atom numbers.
    pub atoms: [i32; 2],
    /// MOL-file bond type (1 single, 2 double, 3 triple, 4 aromatic, ...).
    pub bond_type: i32,
}

/// The part of the Avalon toolkit's `reaccs_molecule_t` the fingerprint
/// reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AvalonMolecule {
    /// Atoms in MOL-file order.
    pub atoms: Vec<AvalonAtom>,
    /// Bonds in MOL-file order.
    pub bonds: Vec<AvalonBond>,
}

/// RDKit's `getFp` + `reaccsToFingerprint`: the fingerprint of `mol` as
/// `n_bits / 8` bytes, bit `i` at byte `i / 8`, mask `1 << (i % 8)`.
///
/// `which_bits` is RDKit's `bitFlags`
/// ([`RDKIT_AVALON_DEFAULT_BIT_FLAGS`] by default); the query variant
/// (`isQuery=True`) is not supported.
pub fn avalon_fp_bytes(mol: &AvalonMolecule, n_bits: usize, which_bits: u32) -> Vec<u8> {
    let n_bytes = n_bits / 8;
    let mut padded = n_bytes;
    while padded % 4 != 0 {
        padded += 1;
    }
    let mut out = vec![0u8; n_bytes];
    if padded == 0 {
        return out;
    }
    let ncounts = padded * 8;
    for dy in [false, true] {
        let mut counts = vec![0i32; ncounts];
        fingerprint::count_fingerprint_patterns(mol, &mut counts, which_bits, dy);
        for (i, &c) in counts.iter().enumerate() {
            if c > 0 && i / 8 < n_bytes {
                out[i / 8] |= 1 << (i % 8);
            }
        }
    }
    out
}

/// RDKit's `GetAvalonFP(mol, nBits=n_bits)` (default `bitFlags`,
/// `isQuery=False`) for a chematic molecule, as `n_bits / 8` bytes with
/// bit `i` at byte `i / 8`, mask `1 << (i % 8)`.
pub fn rdkit_avalon_fp(mol: &Molecule, n_bits: usize) -> Result<Vec<u8>, RdkitAvalonError> {
    let av = convert::avalon_molecule_from_molecule(mol)?;
    Ok(avalon_fp_bytes(&av, n_bits, RDKIT_AVALON_DEFAULT_BIT_FLAGS))
}

/// The on-bit indices of `bytes` (LSB-first within each byte).
pub fn on_bits(bytes: &[u8]) -> Vec<usize> {
    let mut v = Vec::new();
    for (i, &b) in bytes.iter().enumerate() {
        for j in 0..8 {
            if b & (1 << j) != 0 {
                v.push(i * 8 + j);
            }
        }
    }
    v
}
