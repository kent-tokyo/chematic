// Ported from the Avalon Cheminformatics Toolkit, Copyright 2001-2011 Novartis
// Pharma AG, BSD-3-Clause license; see THIRD_PARTY_NOTICES.md.
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

// The port keeps the C code's index loops and argument lists so it reads
// side by side with the toolkit source.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

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
    while !padded.is_multiple_of(4) {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// `(name, SMILES, nBits, on-bits of RDKit 2026.03.1
    /// pyAvalonTools.GetAvalonFP(Chem.MolFromSmiles(SMILES), nBits))`.
    const RDKIT_CASES: &[(&str, &str, usize, &[usize])] = &[
        (
            "benzene",
            "c1ccccc1",
            2048,
            &[343, 464, 853, 863, 1317, 1479, 1674, 1693],
        ),
        (
            "aspirin",
            "CC(=O)Oc1ccccc1C(=O)O",
            2048,
            &[
                15, 85, 157, 246, 262, 281, 322, 329, 334, 343, 365, 444, 464, 553, 579, 702, 707,
                708, 729, 734, 751, 763, 785, 799, 847, 853, 863, 980, 1040, 1139, 1240, 1250,
                1317, 1346, 1348, 1364, 1394, 1439, 1450, 1460, 1466, 1475, 1479, 1487, 1502, 1509,
                1652, 1674, 1689, 1693, 1734, 1841, 1866, 1924, 1948, 2046,
            ],
        ),
        (
            "phenanthrene",
            "c1ccc2c(c1)ccc1ccccc12",
            2048,
            &[
                23, 59, 343, 455, 464, 617, 800, 832, 863, 952, 1014, 1040, 1091, 1191, 1227, 1238,
                1317, 1339, 1395, 1402, 1577, 1674, 1693, 1841, 1852, 1953,
            ],
        ),
        (
            "pyrrole",
            "c1cc[nH]c1",
            2048,
            &[
                6, 15, 36, 60, 81, 117, 125, 181, 212, 251, 263, 325, 343, 397, 411, 414, 522, 528,
                555, 569, 657, 785, 816, 857, 882, 893, 894, 1004, 1037, 1039, 1042, 1108, 1112,
                1118, 1129, 1174, 1218, 1264, 1266, 1282, 1317, 1318, 1331, 1374, 1375, 1409, 1553,
                1564, 1674, 1693, 1699, 1708, 1712, 1732, 1756, 1761, 1794, 1813, 1817, 1881, 2024,
                2026, 2041, 2043,
            ],
        ),
        ("methyl_radical", "[CH3]", 2048, &[579]),
        (
            "sodium_acetate",
            "CC(=O)[O-].[Na+]",
            2048,
            &[60, 262, 505, 579, 734, 1004, 1355, 1652, 1960, 2016, 2046],
        ),
        (
            "deuterated_methanol",
            "[2H]C([2H])([2H])O",
            2048,
            &[334, 490, 579, 1578],
        ),
        (
            "dummy_phenyl",
            "*c1ccccc1",
            2048,
            &[343, 464, 579, 853, 863, 1317, 1479, 1674, 1693, 1841],
        ),
        (
            "cisplatin_dative",
            "[NH3]->[Pt](Cl)(Cl)<-[NH3]",
            2048,
            &[
                49, 93, 97, 125, 357, 382, 418, 505, 579, 649, 1004, 1208, 1277, 1434, 1518, 1717,
                1718, 1946, 1960, 2010,
            ],
        ),
        (
            "cubane",
            "C12C3C4C1C5C2C3C45",
            2048,
            &[
                44, 59, 70, 80, 215, 252, 325, 350, 423, 441, 453, 507, 563, 636, 800, 868, 921,
                979, 1040, 1056, 1072, 1073, 1191, 1212, 1213, 1238, 1270, 1309, 1314, 1339, 1385,
                1393, 1402, 1433, 1438, 1570, 1652, 1698, 1699, 1745, 1966, 2012, 2019, 2026,
            ],
        ),
        (
            "indole_512",
            "c1ccc2[nH]ccc2c1",
            512,
            &[
                6, 9, 10, 13, 15, 16, 18, 20, 27, 28, 29, 41, 42, 48, 55, 57, 59, 60, 63, 65, 69,
                71, 74, 77, 81, 82, 88, 91, 104, 105, 113, 117, 125, 126, 138, 144, 145, 152, 155,
                157, 159, 163, 172, 180, 181, 182, 187, 190, 194, 196, 198, 199, 203, 212, 214,
                224, 225, 240, 244, 250, 251, 257, 258, 263, 266, 270, 273, 274, 280, 281, 287,
                288, 293, 298, 304, 305, 312, 319, 320, 334, 337, 341, 343, 347, 350, 351, 352,
                353, 357, 366, 370, 371, 378, 381, 382, 385, 389, 391, 392, 396, 398, 399, 412,
                417, 418, 421, 439, 440, 442, 445, 446, 455, 464, 466, 471, 475, 476, 481, 484,
                485, 488, 490, 492, 493, 502, 505, 508,
            ],
        ),
        (
            "aspirin_1000",
            "CC(=O)Oc1ccccc1C(=O)O",
            1000,
            &[
                2, 15, 85, 115, 157, 216, 226, 246, 262, 281, 293, 322, 324, 329, 334, 340, 343,
                365, 370, 415, 426, 436, 442, 444, 451, 455, 463, 464, 478, 485, 553, 579, 628,
                650, 665, 669, 702, 707, 708, 710, 729, 734, 751, 763, 785, 799, 817, 842, 847,
                853, 863, 900, 924, 980,
            ],
        ),
    ];

    #[test]
    fn matches_rdkit_on_bits() {
        for &(name, smiles, n_bits, expected) in RDKIT_CASES {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let fp = rdkit_avalon_fp(&mol, n_bits).unwrap();
            assert_eq!(fp.len(), n_bits / 8, "{name}");
            assert_eq!(on_bits(&fp), expected, "{name}");
        }
    }

    #[test]
    fn molblock_path_matches_smiles_path() {
        // RDKit 2026.03.1 `Chem.MolToMolBlock(Chem.MolFromSmiles("c1ccccc1O"))`
        // (coordinates zeroed; the fingerprint does not read them).
        let block = "\n     RDKit          2D\n\n  7  7  0  0  0  0  0  0  0  0999 V2000\n"
            .to_string()
            + &"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n".repeat(6)
            + "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  2  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  2  0\n  6  1  1  0\n  6  7  1  0\nM  END\n";
        let av = read_molblock(&block).unwrap();
        let from_block = avalon_fp_bytes(&av, 2048, RDKIT_AVALON_DEFAULT_BIT_FLAGS);
        let from_smiles =
            rdkit_avalon_fp(&chematic_smiles::parse("c1ccccc1O").unwrap(), 2048).unwrap();
        assert_eq!(from_block, from_smiles);
    }

    #[test]
    fn charges_and_radicals_from_property_lines() {
        let block = "\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n".to_string()
            + "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            + "  1  2  1  0\nM  CHG  1   2  -1\nM  RAD  1   1   2\nM  END\n";
        let av = read_molblock(&block).unwrap();
        assert_eq!(av.atoms[1].charge, -1);
        assert_eq!(av.atoms[0].radical, 2);
    }
}
