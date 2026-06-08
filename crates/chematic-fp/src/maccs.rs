//! MACCS 166-bit structural keys fingerprint.
//!
//! Each of the 166 MACCS keys corresponds to a SMARTS query.
//! A bit is set if the query matches at least once in the molecule.
//! Key 1 (bit 0) is unused by convention; keys 2–166 cover element
//! presence (rare metals) through common functional groups.

use chematic_core::Molecule;
use chematic_smarts::{find_matches, parse_smarts};

use crate::bitvec::BitVec2048;

/// MACCS key SMARTS patterns.
/// Index i corresponds to MACCS key (i+1); bit i is set when pattern i matches.
/// Empty string = unused / always-false key.
static MACCS_SMARTS: &[&str] = &[
    "",                                               // key 1  - unused
    "[#103]",                                         // key 2  - Lr
    "[#102]",                                         // key 3  - No
    "[#101]",                                         // key 4  - Md
    "[#100]",                                         // key 5  - Fm
    "[#99]",                                          // key 6  - Es
    "[#98]",                                          // key 7  - Cf
    "[#97]",                                          // key 8  - Bk
    "[#96]",                                          // key 9  - Cm
    "[#95]",                                          // key 10 - Am
    "[#94]",                                          // key 11 - Pu
    "[#93]",                                          // key 12 - Np
    "[#92]",                                          // key 13 - U
    "[#90,#91]",                                      // key 14 - Th/Pa
    "[#89]",                                          // key 15 - Ac
    "[#88]",                                          // key 16 - Ra
    "[#87]",                                          // key 17 - Fr
    "[#85]",                                          // key 18 - At
    "[#84]",                                          // key 19 - Po
    "[#83]",                                          // key 20 - Bi
    "[#82]",                                          // key 21 - Pb
    "[#81]",                                          // key 22 - Tl
    "[#80]",                                          // key 23 - Hg
    "[#79]",                                          // key 24 - Au
    "[#77,#78]",                                      // key 25 - Ir/Pt
    "[#76]",                                          // key 26 - Os
    "[#75]",                                          // key 27 - Re
    "[#74]",                                          // key 28 - W
    "[#73]",                                          // key 29 - Ta
    "[#72]",                                          // key 30 - Hf
    "[#71]",                                          // key 31 - Lu
    "[#70]",                                          // key 32 - Yb
    "[#69]",                                          // key 33 - Tm
    "[#68]",                                          // key 34 - Er
    "[#67]",                                          // key 35 - Ho
    "[#66]",                                          // key 36 - Dy
    "[#65]",                                          // key 37 - Tb
    "[#64]",                                          // key 38 - Gd
    "[#63]",                                          // key 39 - Eu
    "[#62]",                                          // key 40 - Sm
    "[#61]",                                          // key 41 - Pm
    "[#60]",                                          // key 42 - Nd
    "[#59]",                                          // key 43 - Pr
    "[#58]",                                          // key 44 - Ce
    "[#57]",                                          // key 45 - La
    "[#55,#56]",                                      // key 46 - Cs/Ba
    "[#52,#53,#54]",                                  // key 47 - Te/I/Xe
    "[#51]",                                          // key 48 - Sb
    "[#50]",                                          // key 49 - Sn
    "[#49]",                                          // key 50 - In
    "[#47,#48]",                                      // key 51 - Ag/Cd
    "[#46]",                                          // key 52 - Pd
    "[#45]",                                          // key 53 - Rh
    "[#44]",                                          // key 54 - Ru
    "[#43]",                                          // key 55 - Tc
    "[#42]",                                          // key 56 - Mo
    "[#41]",                                          // key 57 - Nb
    "[#40]",                                          // key 58 - Zr
    "[#39]",                                          // key 59 - Y
    "[#37,#38]",                                      // key 60 - Rb/Sr
    "[#36]",                                          // key 61 - Kr
    "[#35]",                                          // key 62 - Br
    "[#34]",                                          // key 63 - Se
    "[#33]",                                          // key 64 - As
    "[#32]",                                          // key 65 - Ge
    "[#31]",                                          // key 66 - Ga
    "[#30]",                                          // key 67 - Zn
    "[#29]",                                          // key 68 - Cu
    "[#28]",                                          // key 69 - Ni
    "[#27]",                                          // key 70 - Co
    "[#26]",                                          // key 71 - Fe
    "[#25]",                                          // key 72 - Mn
    "[#24]",                                          // key 73 - Cr
    "[#23]",                                          // key 74 - V
    "[#22]",                                          // key 75 - Ti
    "[#21]",                                          // key 76 - Sc
    "[#16;R]",                                        // key 77 - S in ring
    "[#8;R]",                                         // key 78 - O in ring
    "[#7;R]",                                         // key 79 - N in ring
    "[#16]",                                          // key 80 - any S
    "[#15]",                                          // key 81 - any P
    "[#14]",                                          // key 82 - Si
    "[#6]~[#16]",                                     // key 83 - C-S bond
    "[#7]~[#6]~[#7]",                                 // key 84 - N-C-N
    "[#7]~[#7]",                                      // key 85 - N-N bond
    "[#8]~[#8]",                                      // key 86 - O-O bond
    "[#8]~[#15]",                                     // key 87 - O-P bond
    "[#16]~[#8]",                                     // key 88 - S-O bond
    "[#6]=[#16]",                                     // key 89 - C=S
    "[#16]=[#7]",                                     // key 90 - S=N
    "[#6]=[#7]",                                      // key 91 - C=N
    "[#7]~[#6]=[#8]",                                 // key 92 - N-C=O (amide-like)
    "[#8]~[#6]=[#8]",                                 // key 93 - O-C=O (ester/acid)
    "[#6]=[#6]",                                      // key 94 - C=C
    "[#6]#[#7]",                                      // key 95 - C#N (nitrile)
    "[#6]#[#6]",                                      // key 96 - C#C (alkyne)
    "[#6]~[#15]",                                     // key 97 - C-P
    "[#6]~[#8]~[#6]",                                 // key 98 - C-O-C (ether)
    "[#6]~[#7]~[#6]",                                 // key 99 - C-N-C
    "[#6]~[#16]~[#6]",                                // key 100 - C-S-C
    "[#8]~[#6]~[#8]",                                 // key 101 - O-C-O
    "[#7]~[#6]~[#8]",                                 // key 102 - N-C-O
    "[#7]~[#6]~[#16]",                                // key 103 - N-C-S
    "[#6]=[#6]~[#6]",                                 // key 104 - C=C-C
    "[#6]=[#6]~[#7]",                                 // key 105 - C=C-N
    "[#6]=[#6]~[#8]",                                 // key 106 - C=C-O
    "[#6]=[#6]~[#16]",                                // key 107 - C=C-S
    "[#6]=[#6]~[#6]=[#6]",                            // key 108 - diene
    "[#6]=[#7]~[#6]=[#8]",                            // key 109 - C=N-C=O
    "[#6]=[#7]~[#6]=[#7]",                            // key 110 - C=N-C=N
    "[#6]=[#8]~[#7]~[#6]=[#8]",                       // key 111
    "[#6]~[#6]~[#8]~[#6]=[#8]",                       // key 112 - ester chain
    "[#6]~[#6]~[#7]~[#6]=[#8]",                       // key 113 - amide chain
    "[#6]~[#8]~[#6]=[#8]",                            // key 114 - O-C=O ester
    "[#7]~[#6](=[#8])~[#7]",                          // key 115 - urea
    "[#6]=[#8]",                                      // key 116 - C=O carbonyl
    "[#6]~[#7](~[#6])~[#6]",                          // key 117 - tertiary amine
    "[#8]~[#6]~[#7]",                                 // key 118 - O-C-N
    "[!#1;!#6]~[#6]=[#8]",                            // key 119 - heteroatom adj to C=O
    "[#6]=[#8]~[#8]",                                 // key 120 - peracid
    "[#7]=[#8]",                                      // key 121 - N=O
    "[#7;R]~[#6;R]=[#7;R]",                           // key 122 - amidin in ring
    "[#6]~[#8]~[#8]~[#6]",                            // key 123 - peroxide
    "[#16]=[#8]",                                     // key 124 - S=O
    "[!#6;!#1]~[!#6;!#1]",                            // key 125 - het-het bond
    "[!#6;!#1;!#7;!#8;!#16;!#15;!#9;!#17;!#35;!#53]", // key 126 - unusual het
    "[#7]~[#6]~[#7]~[#6]~[#8]",                       // key 127
    "[#7]~[#6]~[#16]",                                // key 128 - N-C-S
    "[#7]~[#7]~[#6]",                                 // key 129 - N-N-C
    "[#7]~[#6]=[#6]~[#7]",                            // key 130 - N-C=C-N
    "[#6]=[#7]~[#7]=[#6]",                            // key 131
    "[#8]~[#16](=[#8])=[#8]",                         // key 132 - sulfate
    "[#16]~[#6]~[#16]",                               // key 133
    "[!#1;!#6]~[!#1;!#6]~[!#1;!#6]",                  // key 134 - 3 het chain
    "[#6]~[#16]~[#8]~[#6]",                           // key 135
    "[#6]~[#7]~[#8]",                                 // key 136
    "[#7]~[#7]~[#7]",                                 // key 137 - triazole/azide
    "[#6]~[#7]~[#7]~[#7]",                            // key 138
    "[#8;!R]~[#6;R]",                                 // key 139 - exocyclic O on ring C
    "[#7;R]~[#6;!R]=[#8]",                            // key 140 - exocyclic C=O on ring N
    "[#6]~[#8]~[#6]~[#8]",                            // key 141
    "[#7]~[#6](~[#8])~[#7]",                          // key 142 - urea variant
    "[!#1;!#6]~[!#1;!#6]~[!#1;!#6]~[!#1;!#6]",        // key 143
    "[#6]~[#7;R]~[#6]~[#7;R]",                        // key 144
    "[#6]~[#6]~[#8]~[#6]~[#6]",                       // key 145
    "[#7]~[#7]~[#6]=[#8]",                            // key 146 - hydrazide
    "[#6]~[#6]~[#7]~[#7]",                            // key 147
    "[#7;R]~[#6;R]~[#7;R]~[#6;R]",                    // key 148
    "[#6]~[#8]~[#6]~[#6]",                            // key 149
    "[#16;R]~[#6;R]~[#7;R]",                          // key 150
    "[#16;R]~[#6;R]~[#8;R]",                          // key 151
    "[#16;R]~[#6;R]=[#7;R]",                          // key 152
    "[#7;R]~[#6;R]=[#7;R]",                           // key 153 - imidazole N-C=N
    "[#7;R]~[#6;R]=[#8;R]",                           // key 154 - lactam N-C=O
    "[#8;R]~[#6;R]=[#8;R]",                           // key 155
    "[#8;R]~[#6;R]~[#7;R]",                           // key 156
    "[#8;R]~[#6;R]~[#8;R]",                           // key 157
    "[#8;R]~[#6;R]~[#6;R]",                           // key 158 - O-C-C in ring
    "[#7;R]~[#6;R]~[#6;R]",                           // key 159 - N-C-C in ring
    "[#6;R]~[#6;R]~[#6;R]~[#6;R]~[#6;R]~[#6;R]",      // key 160 - 6C chain in ring
    "[a]~[a]~[a]~[a]~[a]~[a]",                        // key 161 - 6-atom aromatic chain
    "[a]",                                            // key 162 - any aromatic atom
    "[!#6;a]",                                        // key 163 - heteroaromatic
    "[!#6;!#1]",                                      // key 164 - any heteroatom (non-C, non-H)
    "[#6;R]",                                         // key 165 - any ring C
    "[R]",                                            // key 166 - any ring atom
];

/// Compute the MACCS 166-bit structural keys fingerprint for `mol`.
///
/// Each of the 166 bits corresponds to a structural feature (SMARTS pattern).
/// Bit `i` (0-indexed) is set if MACCS key `i+1` matches the molecule.
pub fn maccs(mol: &Molecule) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    for (i, &pattern) in MACCS_SMARTS.iter().enumerate() {
        if pattern.is_empty() {
            continue;
        }
        if let Ok(query) = parse_smarts(pattern)
            && !find_matches(&query, mol).is_empty()
        {
            fp.set(i);
        }
        // Silently skip patterns that fail to parse or match errors
    }
    fp
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn maccs_benzene_nonzero() {
        let mol = parse("c1ccccc1").unwrap();
        let fp = maccs(&mol);
        assert!(fp.popcount() > 0, "benzene maccs should be nonzero");
    }

    #[test]
    fn maccs_ethanol_nonzero() {
        let mol = parse("CCO").unwrap();
        let fp = maccs(&mol);
        assert!(fp.popcount() > 0, "ethanol maccs should be nonzero");
    }

    #[test]
    fn maccs_benzene_has_aromatic_bit() {
        let mol = parse("c1ccccc1").unwrap();
        let fp = maccs(&mol);
        // bit 161 (0-indexed) = key 162 = any aromatic atom
        assert!(
            fp.get(161),
            "benzene should have aromatic bit (key 162, index 161) set"
        );
    }

    #[test]
    fn maccs_deterministic() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(maccs(&mol), maccs(&mol), "maccs must be deterministic");
    }

    #[test]
    fn maccs_aspirin_has_carbonyl_bit() {
        // aspirin: CC(=O)Oc1ccccc1C(=O)O — has C=O (key 116 = bit 115)
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(115),
            "aspirin should have C=O bit (key 116, index 115) set"
        );
    }

    #[test]
    fn maccs_acetonitrile_has_triple_bond_bit() {
        // CC#N — has C#N (key 95 = bit 94)
        let mol = parse("CC#N").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(94),
            "acetonitrile should have C#N bit (key 95, index 94) set"
        );
    }

    #[test]
    fn maccs_bromobenzene_has_bromine_bit() {
        // c1ccccc1Br — has Br (key 62 = bit 61)
        let mol = parse("c1ccccc1Br").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(61),
            "bromobenzene should have Br bit (key 62, index 61) set"
        );
    }
}
