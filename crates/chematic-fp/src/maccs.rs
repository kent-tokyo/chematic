//! MACCS 166-bit structural keys fingerprint.
//!
//! Bit i (0-indexed) corresponds to MACCS key (i+1).
//! SMARTS patterns match the official RDKit MACCSkeys.py definitions.
//!
//! Three keys use code-only logic (no SMARTS pattern):
//! - Key 1  (index 0):   any isotopically-labelled atom
//! - Key 125 (index 124): aromatic ring count > 1
//! - Key 166 (index 165): fragment count > 1 (disconnected molecule)
//!
//! For count-based keys the second tuple element is the minimum count:
//! the bit is set when `matches.len() > min_count`.

use chematic_core::Molecule;
use chematic_perception::find_sssr;
use chematic_smarts::{find_matches, parse_smarts};

use crate::bitvec::BitVec2048;

/// `(smarts, min_count)` pairs for all 166 MACCS keys (index = key - 1).
/// Empty string or "?" → special-case key handled in `maccs()`.
static MACCS_SMARTS: &[(&str, usize)] = &[
    ("?", 0),                                                   // key 1  - ISOTOPE (special)
    ("[#104]", 0),                                              // key 2  - Atomic num >103
    ("[#32,#33,#34,#50,#51,#52,#82,#83,#84]", 0),               // key 3  - Group IVa,Va,VIa
    ("[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]", 0),       // key 4  - Actinide
    ("[Sc,Ti,Y,Zr,Hf]", 0),                                     // key 5  - Group IIIB,IVB
    ("[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]", 0),      // key 6  - Lanthanide
    ("[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]", 0),                          // key 7  - Group VB,VIB,VIIB
    ("[!#6;!#1]1~*~*~*~1", 0),                                  // key 8  - 4M heterocycle
    ("[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]", 0),                        // key 9  - Group VIII
    ("[Be,Mg,Ca,Sr,Ba,Ra]", 0),                                 // key 10 - Group IIa
    ("*1~*~*~*~1", 0),                                          // key 11 - 4M ring
    ("[Cu,Zn,Ag,Cd,Au,Hg]", 0),                                 // key 12 - Group IB,IIB
    ("[#8]~[#7](~[#6])~[#6]", 0),                               // key 13 - ON(C)C
    ("[#16]-[#16]", 0),                                         // key 14 - S-S
    ("[#8]~[#6](~[#8])~[#8]", 0),                               // key 15 - OC(O)O
    ("[!#6;!#1]1~*~*~1", 0),                                    // key 16 - 3M heterocycle
    ("[#6]#[#6]", 0),                                           // key 17 - alkyne
    ("[#5,#13,#31,#49,#81]", 0),                                // key 18 - Group IIIA
    ("*1~*~*~*~*~*~*~1", 0),                                    // key 19 - 7M ring
    ("[#14]", 0),                                               // key 20 - Si
    ("[#6]=[#6](~[!#6;!#1])~[!#6;!#1]", 0),                     // key 21 - C=C(Q)Q
    ("*1~*~*~1", 0),                                            // key 22 - 3M ring
    ("[#7]~[#6](~[#8])~[#8]", 0),                               // key 23 - NC(O)O
    ("[#7]-[#8]", 0),                                           // key 24 - N-O
    ("[#7]~[#6](~[#7])~[#7]", 0),                               // key 25 - NC(N)N
    ("[#6;R]=[#6;R]", 0), // key 26 - ring C=C (was =;@ which fails to parse)
    ("[I]", 0),           // key 27 - I
    ("[!#6;!#1]~[CH2]~[!#6;!#1]", 0), // key 28 - QCH2Q
    ("[#15]", 0),         // key 29 - P
    ("[#6]~[!#6;!#1](~[#6])(~[#6])~*", 0), // key 30 - CQ(C)(C)A
    ("[!#6;!#1]~[F,Cl,Br,I]", 0), // key 31 - QX
    ("[#6]~[#16]~[#7]", 0), // key 32 - CSN
    ("[#7]~[#16]", 0),    // key 33 - NS
    ("[CH2]=*", 0),       // key 34 - CH2=A
    ("[Li,Na,K,Rb,Cs,Fr]", 0), // key 35 - Alkali metal
    ("[#16R]", 0),        // key 36 - S heterocycle
    ("[#7]~[#6](~[#8])~[#7]", 0), // key 37 - NC(O)N
    ("[#7]~[#6](~[#6])~[#7]", 0), // key 38 - NC(C)N
    ("[#8]~[#16](~[#8])~[#8]", 0), // key 39 - OS(O)O
    ("[#16]-[#8]", 0),    // key 40 - S-O
    ("[#6]#[#7]", 0),     // key 41 - nitrile
    ("F", 0),             // key 42 - F
    ("[!#6;!#1;!H0]~*~[!#6;!#1;!H0]", 0), // key 43 - QHAQH
    ("[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]", 0), // key 44 - OTHER
    ("[#6]=[#6]~[#7]", 0), // key 45 - C=CN
    ("Br", 0),            // key 46 - Br
    ("[#16]~*~[#7]", 0),  // key 47 - SAN
    ("[#8]~[!#6;!#1](~[#8])(~[#8])", 0), // key 48 - OQ(O)O
    ("[!+0]", 0),         // key 49 - charged atom
    ("[#6]=[#6](~[#6])~[#6]", 0), // key 50 - C=C(C)C
    ("[#6]~[#16]~[#8]", 0), // key 51 - CSO
    ("[#7]~[#7]", 0),     // key 52 - N-N
    ("[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]", 0), // key 53 - QHAAAQH
    ("[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]", 0), // key 54 - QHAAQH
    ("[#8]~[#16]~[#8]", 0), // key 55 - OSO
    ("[#8]~[#7](~[#8])~[#6]", 0), // key 56 - ON(O)C
    ("[#8R]", 0),         // key 57 - O heterocycle
    ("[!#6;!#1]~[#16]~[!#6;!#1]", 0), // key 58 - QSQ
    ("[#16]!:*:*", 0),    // key 59 - S-nonarom-arom chain
    ("[#16]=[#8]", 0),    // key 60 - S=O
    ("*~[#16](~*)~*", 0), // key 61 - AS(A)A
    ("*@*!@*@*", 0),      // key 62 - ring-open-ring
    ("[#7]=[#8]", 0),     // key 63 - N=O
    ("*@*!@[#16]", 0),    // key 64 - ring-open-S
    ("c:n", 0),           // key 65 - aromatic C-N
    ("[#6]~[#6](~[#6])(~[#6])~*", 0), // key 66 - quaternary C
    ("[!#6;!#1]~[#16]", 0), // key 67 - QS
    ("[!#6;!#1;!H0]~[!#6;!#1;!H0]", 0), // key 68 - QHQH
    ("[!#6;!#1]~[!#6;!#1;!H0]", 0), // key 69 - QQH
    ("[!#6;!#1]~[#7]~[!#6;!#1]", 0), // key 70 - QNQ
    ("[#7]~[#8]", 0),     // key 71 - N-O (any bond)
    ("[#8]~*~*~[#8]", 0), // key 72 - O..O (2 bonds)
    ("[#16]=*", 0),       // key 73 - S=A
    ("[CH3]~*~[CH3]", 0), // key 74 - CH3-A-CH3
    ("*!@[#7]@*", 0),     // key 75 - open-N-ring
    ("[#6]=[#6](~*)~*", 0), // key 76 - C=C(A)A
    ("[#7]~*~[#7]", 0),   // key 77 - N-A-N
    ("[#6]=[#7]", 0),     // key 78 - C=N
    ("[#7]~*~*~[#7]", 0), // key 79 - N-A-A-N
    ("[#7]~*~*~*~[#7]", 0), // key 80 - N-A-A-A-N
    ("[#16]~*(~*)~*", 0), // key 81 - S(A)A
    ("*~[CH2]~[!#6;!#1;!H0]", 0), // key 82 - A-CH2-QH
    ("[!#6;!#1]1~*~*~*~*~1", 0), // key 83 - 5M heterocycle
    ("[NH2]", 0),         // key 84 - NH2
    ("[#6]~[#7](~[#6])~[#6]", 0), // key 85 - tertiary N
    ("[C;H2,H3][!#6;!#1][C;H2,H3]", 0), // key 86 - CH2/CH3 Q CH2/CH3
    ("[F,Cl,Br,I]!@*@*", 0), // key 87 - X-open-ring
    ("[#16]", 0),         // key 88 - any S
    ("[#8]~*~*~*~[#8]", 0), // key 89 - O..O (3 bonds)
    (
        "[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]",
        0,
    ), // key 90
    (
        "[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]",
        0,
    ), // key 91
    ("[#8]~[#6](~[#7])~[#6]", 0), // key 92 - OC(N)C
    ("[!#6;!#1]~[CH3]", 0), // key 93 - Q-CH3
    ("[!#6;!#1]~[#7]", 0), // key 94 - Q-N
    ("[#7]~*~*~[#8]", 0), // key 95 - N-A-A-O
    ("*1~*~*~*~*~1", 0),  // key 96 - 5M ring
    ("[#7]~*~*~*~[#8]", 0), // key 97 - N-A-A-A-O
    ("[!#6;!#1]1~*~*~*~*~*~1", 0), // key 98 - 6M heterocycle
    ("[#6]=[#6]", 0),     // key 99 - C=C
    ("*~[CH2]~[#7]", 0),  // key 100 - A-CH2-N
    (
        "[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]",
        0,
    ), // key 101 - 8-14M ring
    ("[!#6;!#1]~[#8]", 0), // key 102 - Q-O
    ("Cl", 0),            // key 103 - Cl
    ("[!#6;!#1;!H0]~*~[CH2]~*", 0), // key 104 - QH-A-CH2-A
    ("*@*(@*)@*", 0),     // key 105 - ring branch
    ("[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]", 0), // key 106 - QA(Q)Q
    ("[F,Cl,Br,I]~*(~*)~*", 0), // key 107 - X-A(A)A
    ("[CH3]~*~*~*~[CH2]~*", 0), // key 108 - CH3-A-A-A-CH2-A
    ("*~[CH2]~[#8]", 0),  // key 109 - A-CH2-O
    ("[#7]~[#6]~[#8]", 0), // key 110 - NCO
    ("[#7]~*~[CH2]~*", 0), // key 111 - N-A-CH2-A
    ("*~*(~*)(~*)~*", 0), // key 112 - branched atom
    ("[#8]!:*:*", 0),     // key 113 - O-nonarom-arom
    ("[CH3]~[CH2]~*", 0), // key 114 - CH3-CH2-A
    ("[CH3]~*~[CH2]~*", 0), // key 115 - CH3-A-CH2-A
    ("[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]", 0), // key 116 - CH3-A-A-CH2
    ("[#7]~*~[#8]", 0),   // key 117 - N-A-O
    ("[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", 1), // key 118 - ACH2CH2A > 1
    ("[#7]=*", 0),        // key 119 - N=A
    ("[!#6;R]", 1),       // key 120 - het in ring > 1
    ("[#7;R]", 0),        // key 121 - N heterocycle
    ("*~[#7](~*)~*", 0),  // key 122 - AN(A)A
    ("[#8]~[#6]~[#8]", 0), // key 123 - OCO
    ("[!#6;!#1]~[!#6;!#1]", 0), // key 124 - QQ
    ("?", 0),             // key 125 - aromatic ring > 1 (special)
    ("*!@[#8]!@*", 0),    // key 126 - A-open-O-open-A
    ("*@*!@[#8]", 1),     // key 127 - ring-open-O > 1
    (
        "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]",
        0,
    ), // key 128
    (
        "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1),$(*~[CH2]~*1~*~[CH2]1)]",
        0,
    ), // key 129
    ("[!#6;!#1]~[!#6;!#1]", 1), // key 130 - QQ > 1
    ("[!#6;!#1;!H0]", 1), // key 131 - QH > 1
    ("[#8]~*~[CH2]~*", 0), // key 132 - O-A-CH2-A
    ("*@*!@[#7]", 0),     // key 133 - ring-open-N
    ("[F,Cl,Br,I]", 0),   // key 134 - halogen
    ("[#7]!:*:*", 0),     // key 135 - N-nonarom-arom
    ("[#8]=*", 1),        // key 136 - O=A > 1
    ("[!C;!c;R]", 0),     // key 137 - non-C in ring
    ("[!#6;!#1]~[CH2]~*", 1), // key 138 - QCH2A > 1
    ("[O;!H0]", 0),       // key 139 - OH
    ("[#8]", 3),          // key 140 - O > 3
    ("[CH3]", 2),         // key 141 - CH3 > 2
    ("[#7]", 1),          // key 142 - N > 1
    ("*@*!@[#8]", 0),     // key 143 - ring-open-O
    ("*!:*:*!:*", 0),     // key 144 - nonarom-arom-nonarom
    ("*1~*~*~*~*~*~1", 1), // key 145 - 6M ring > 1
    ("[#8]", 2),          // key 146 - O > 2
    ("[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]", 0), // key 147 - ACH2CH2A
    ("*~[!#6;!#1](~*)~*", 0), // key 148 - AQ(A)A
    ("[C;H3,H4]", 1),     // key 149 - CH3 > 1
    ("*!@*@*!@*", 0),     // key 150 - open-ring-open
    ("[#7;!H0]", 0),      // key 151 - NH
    ("[#8]~[#6](~[#6])~[#6]", 0), // key 152 - OC(C)C
    ("[!#6;!#1]~[CH2]~*", 0), // key 153 - QCH2A
    ("[#6]=[#8]", 0),     // key 154 - C=O
    ("*!@[CH2]!@*", 0),   // key 155 - open-CH2-open
    ("[#7]~*(~*)~*", 0),  // key 156 - NA(A)A
    ("[#6]-[#8]", 0),     // key 157 - C-O
    ("[#6]-[#7]", 0),     // key 158 - C-N
    ("[#8]", 1),          // key 159 - O > 1
    ("[C;H3,H4]", 0),     // key 160 - CH3
    ("[#7]", 0),          // key 161 - N
    ("[a]", 0),           // key 162 - aromatic atom
    ("*1~*~*~*~*~*~1", 0), // key 163 - 6M ring
    ("[#8]", 0),          // key 164 - O
    ("[R]", 0),           // key 165 - ring atom
    ("?", 0),             // key 166 - fragments > 1 (special)
];

fn count_aromatic_rings(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atom(idx).aromatic))
        .count()
}

/// Compute the MACCS 166-bit structural keys fingerprint for `mol`.
///
/// Bit `i` (0-indexed) corresponds to MACCS key `i+1`. Keys 1, 125, and 166
/// use code-only logic; all other keys use SMARTS matching.
/// For count-based keys (14 total), the bit is set only when the match count
/// exceeds the threshold stored in `MACCS_SMARTS`.
pub fn maccs(mol: &Molecule) -> BitVec2048 {
    let mut fp = BitVec2048::new();

    for (i, &(pattern, min_count)) in MACCS_SMARTS.iter().enumerate() {
        if pattern == "?" || pattern.is_empty() {
            continue;
        }
        if let Ok(query) = parse_smarts(pattern)
            && find_matches(&query, mol).len() > min_count
        {
            fp.set(i);
        }
        // Patterns that fail to parse silently leave the bit unset.
    }

    // Key 1 (index 0): isotope — any atom carrying an explicit mass label
    if mol.atoms().any(|(_, a)| a.isotope.is_some()) {
        fp.set(0);
    }

    // Key 125 (index 124): aromatic ring count > 1
    if count_aromatic_rings(mol) > 1 {
        fp.set(124);
    }

    // Key 166 (index 165): fragment count > 1 (disconnected molecule)
    if !mol.is_connected() {
        fp.set(165);
    }

    fp
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn maccs_benzene_has_aromatic_bits() {
        let mol = parse("c1ccccc1").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(161),
            "key 162 (any aromatic atom) should be set for benzene"
        );
        assert!(fp.get(162), "key 163 (6M ring) should be set for benzene");
        assert!(fp.get(164), "key 165 (ring atom) should be set for benzene");
    }

    #[test]
    fn maccs_benzene_aromatic_ring_gt1_not_set() {
        // Benzene has exactly 1 aromatic ring → key 125 should NOT be set
        let mol = parse("c1ccccc1").unwrap();
        let fp = maccs(&mol);
        assert!(
            !fp.get(124),
            "benzene has only 1 aromatic ring; key 125 must be 0"
        );
    }

    #[test]
    fn maccs_naphthalene_aromatic_ring_gt1_set() {
        // Naphthalene has 2 aromatic rings → key 125 should be set
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(124),
            "naphthalene has 2 aromatic rings; key 125 must be 1"
        );
    }

    #[test]
    fn maccs_ethanol_has_oxygen_bits() {
        let mol = parse("CCO").unwrap();
        let fp = maccs(&mol);
        assert!(fp.get(163), "key 164 (O) should be set for ethanol");
        assert!(fp.get(156), "key 157 (C-O) should be set for ethanol");
        assert!(fp.get(138), "key 139 (OH) should be set for ethanol");
    }

    #[test]
    fn maccs_isotope_key() {
        // [13C] has an isotope label → key 1 should be set
        let labeled = parse("[13CH4]").unwrap();
        let unlabeled = parse("C").unwrap();
        assert!(maccs(&labeled).get(0), "key 1 should be set for [13C]");
        assert!(!maccs(&unlabeled).get(0), "key 1 should not be set for CH4");
    }

    #[test]
    fn maccs_fragment_key() {
        // NaCl ([Na+].[Cl-]) is disconnected → key 166 should be set
        let salt = parse("[Na+].[Cl-]").unwrap();
        let ethanol = parse("CCO").unwrap();
        assert!(
            maccs(&salt).get(165),
            "key 166 should be set for [Na+].[Cl-]"
        );
        assert!(
            !maccs(&ethanol).get(165),
            "key 166 should not be set for CCO"
        );
    }

    #[test]
    fn maccs_oxygen_count_keys() {
        // Citric acid has 6 O → keys 159/146/140 (O>1, O>2, O>3) should be set
        let citric = parse("OC(CC(O)(C(=O)O)C(=O)O)(C(=O)O)").unwrap();
        let fp = maccs(&citric);
        assert!(fp.get(158), "key 159 (O>1) should be set for citric acid");
        assert!(fp.get(145), "key 146 (O>2) should be set for citric acid");
        assert!(fp.get(139), "key 140 (O>3) should be set for citric acid");
    }

    #[test]
    fn maccs_nitrile_key() {
        let mol = parse("CC#N").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(40),
            "key 41 (C#N nitrile) should be set for acetonitrile"
        );
    }

    #[test]
    fn maccs_bromine_key() {
        let mol = parse("c1ccccc1Br").unwrap();
        let fp = maccs(&mol);
        assert!(fp.get(45), "key 46 (Br) should be set for bromobenzene");
    }

    #[test]
    fn maccs_carbonyl_key() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let fp = maccs(&mol);
        assert!(fp.get(153), "key 154 (C=O) should be set for aspirin");
    }

    #[test]
    fn maccs_deterministic() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(maccs(&mol), maccs(&mol));
    }

    #[test]
    fn maccs_n_gt1_key() {
        // Ethylenediamine has 2 N → key 142 (N>1) should be set
        let mol = parse("NCCN").unwrap();
        let fp = maccs(&mol);
        assert!(
            fp.get(141),
            "key 142 (N>1) should be set for ethylenediamine"
        );
    }

    #[test]
    fn maccs_key26_ring_double_bond() {
        // Cyclohexene has a ring C=C → key 26 (index 25) should be set
        let mol = parse("C1=CCCCC1").unwrap();
        assert!(
            maccs(&mol).get(25),
            "cyclohexene should set key 26 (ring C=C)"
        );
        // Ethylene has C=C but not in ring → key 26 should NOT be set
        let mol2 = parse("C=C").unwrap();
        assert!(
            !maccs(&mol2).get(25),
            "ethylene: key 26 (ring C=C) should be 0"
        );
    }

    #[test]
    fn all_maccs_smarts_parse() {
        let mut failed = Vec::new();
        for (i, &(pattern, _)) in MACCS_SMARTS.iter().enumerate() {
            if pattern == "?" || pattern.is_empty() {
                continue;
            }
            if parse_smarts(pattern).is_err() {
                failed.push((i + 1, pattern));
            }
        }
        assert!(
            failed.is_empty(),
            "MACCS SMARTS parse failures: {:?}",
            failed
        );
    }
}
