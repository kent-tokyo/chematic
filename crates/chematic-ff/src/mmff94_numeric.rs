//! MMFF94 numeric atom types (1–99) and faithful partial-charge calculation.
//!
//! ## Charge formula (Halgren 1996, MMFF.V, equation 15)
//!
//! For atom i with neighbors j:
//!   q_i = (1 - M·v)·q0_i + v·Σq_j + Σ bci(j→i)
//!
//! where:
//! - v = fcadj(type_i) from PBCI table (0.0 for most organic atoms)
//! - M = coordination number
//! - bci(j→i) comes from CHG table or PBCI fallback
//!
//! ## CHG sign convention
//!
//! Entry (bond_type, a, b, bci): `bci` is the charge added to the SECOND atom (b).
//! For atom i bonded to j:
//!   - Found as (bt, j, i): atom i is second → sign=+1, contrib = +bci
//!   - Found as (bt, i, j): atom i is first  → sign=−1, contrib = −bci
//!   - Not found: contrib = PBCI(type_i) − PBCI(type_j)

#![allow(clippy::approx_constant)]

use chematic_core::{AtomIdx, BondOrder, Element, Molecule, implicit_hcount};

// ── Error type ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NumericTypeError(pub String);

impl std::fmt::Display for NumericTypeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "MMFF94 numeric type error: {}", self.0)
    }
}

// ── PBCI table: (atom_type, pbci, fcadj) ────────────────────────────────────
// Source: RDKit Code/ForceField/MMFF/Params.cpp — defaultMMFFPBCI
// 99 entries, one per MMFF94 atom type.
//
// Issue #227 Phase 2 Step 6: a table-level cross-check against a fresh
// `Params.cpp` fetch (pinned commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`)
// confirmed every (pbci, fcadj) VALUE below is byte-identical to RDKit's real
// `defaultMMFFPBCI` -- the table itself was never wrong. Five of the
// trailing `//` comment symbols/descriptions were, however (types 32, 34,
// 45, 47, 53, corrected below with a note on each row) -- these are
// human-written labels only, not read by any lookup, so the mislabeling
// never affected `pbci_for`'s actual behavior. It did, before this
// investigation, seed a false lead in this file's own audit trail (an
// earlier hypothesis assumed the wrong-looking comments meant wrong
// values).
static MMFF94_PBCI: &[(u8, f64, f64)] = &[
    (1, 0.000, 0.000),   // CR       alkyl carbon
    (2, -0.135, 0.000),  // C=C      vinylic carbon
    (3, -0.095, 0.000),  // C=O      carbonyl carbon
    (4, -0.200, 0.000),  // CSP      acetylenic carbon
    (5, -0.023, 0.000),  // HC       H on carbon
    (6, -0.243, 0.000),  // OR       single-bond oxygen (ether/alcohol)
    (7, -0.687, 0.000),  // O=C      carbonyl oxygen
    (8, -0.253, 0.000),  // NR       sp3 amine nitrogen
    (9, -0.306, 0.000),  // N=C      imine nitrogen
    (10, -0.244, 0.000), // NC=O     amide nitrogen
    (11, -0.317, 0.000), // F
    (12, -0.304, 0.000), // CL
    (13, -0.238, 0.000), // BR
    (14, -0.208, 0.000), // I
    (15, -0.236, 0.000), // S        thiol/sulfide
    (16, -0.475, 0.000), // S=C      S doubly bonded to C
    (17, -0.191, 0.000), // SO       sulfoxide S
    (18, -0.118, 0.000), // SO2      sulfone S
    (19, 0.094, 0.000),  // SI       silicon
    (20, -0.019, 0.000), // P        phosphorus
    (21, 0.157, 0.000),  // =P       phosphorus =P
    (22, -0.095, 0.000), // P=O      phosphoryl P
    (23, 0.193, 0.000),  // HNR      H on amine N
    (24, 0.257, 0.000),  // HOCO     H on O in acid/alcohol
    (25, 0.012, 0.000),  // PO4      phosphate P
    (26, -0.142, 0.000), // P=S      P=S
    (27, 0.094, 0.000),  // HN=C     H on imine N
    (28, 0.058, 0.000),  // HNCO     H on amide N
    (29, 0.207, 0.000),  // HOCO     H on O in ester
    (30, -0.166, 0.000), // N2OX     N-oxide N
    (31, 0.161, 0.000),  // HOH      H in water
    (32, -0.732, 0.500), // O2CM     carboxylate/nitro-nitrate/sulfate/phosphate O (fcadj=0.5) -- was mislabeled "NR+", value already correct
    (33, 0.257, 0.000),  // HOX      H on O in N-oxide
    (34, -0.491, 0.000), // NR+      quaternary/protonated nitrogen -- was mislabeled "O-", value already correct
    (35, -0.456, 0.500), // OM       oxide oxygen (fcadj=0.5)
    (36, -0.031, 0.000), // HNR+     H on protonated N
    (37, -0.127, 0.000), // C5A      aromatic C in 5-ring alpha to N
    (38, -0.437, 0.000), // C5B      aromatic C in 5-ring
    (39, -0.104, 0.000), // C5       generic 5-ring aromatic C
    (40, -0.264, 0.000), // N5A      aromatic N in 5-ring (NH type)
    (41, 0.052, 0.000),  // N5B      aromatic N in 5-ring (sp2 type)
    (42, -0.757, 0.000), // N5+      protonated aromatic N
    (43, -0.326, 0.000), // O5       aromatic O (furan)
    (44, -0.237, 0.000), // S5       aromatic S (thiophene)
    (45, -0.260, 0.000), // NO2/NO3  nitro/nitrate nitrogen -- was mislabeled "N5", value already correct
    (46, -0.429, 0.000), // NO2      nitro N
    (47, -0.418, 0.000), // NAZT     terminal N in azido/diazo group -- was mislabeled "NO3", value already correct
    (48, -0.525, 0.000), // O2NO     nitro O
    (49, -0.283, 0.000), // O3NO     nitrate O
    (50, 0.284, 0.000),  // OP       phosphate O
    (51, -1.046, 0.000), // O2P      phosphonate =O
    (52, -0.546, 0.000), // O3P      bridging phosphate O
    (53, -0.048, 0.000), // =N=      central N in C=N=N or N=N=N (azide) -- was mislabeled "O4P", value already correct
    (54, -0.424, 0.000), // O4CL     perchlorate O
    (55, -0.476, 0.000), // CLO4     perchlorate Cl
    (56, -0.438, 0.000), // C=ON     C in amide (alternative)
    (57, -0.105, 0.000), // CORR     corrected aromatic C
    (58, -0.488, 0.000), // N5A+     aromatic N in 5-ring (N=C side)
    (59, -0.337, 0.000), // N5B+     aromatic N in 5-ring
    (60, -0.635, 0.000), // NC=O2    urea N
    (61, -0.265, 0.000), // NC=O3    carbamate N
    (62, -0.125, 0.250), // NM       anionic N (fcadj=0.25)
    (63, -0.180, 0.000), // CB       aromatic C in 6-ring (benzene)
    (64, -0.181, 0.000), // C_6ring  aromatic C in 6-ring (variant)
    (65, -0.475, 0.000), // N5       aromatic 5-ring N (pyrrole NH)
    (66, -0.467, 0.000), // N5+      aromatic 5-ring N (imidazole N=C)
    (67, -0.099, 0.000), // N6A      aromatic 6-ring N (pyridine)
    (68, -0.135, 0.000), // N6B      aromatic 6-ring N (pyrimidine)
    (69, -0.099, 0.000), // NPYD     pyridinium N
    (70, -0.269, 0.000), // NPYR     pyridyl N
    (71, -0.071, 0.000), // HS       H on sulfur
    (72, -0.580, 0.500), // SO4      sulfonate/sulfate O (fcadj=0.5)
    (73, -0.200, 0.000), // S2CM     thiocarboxylate
    (74, -0.301, 0.000), // SCN      thiocyanate S
    (75, -0.255, 0.000), // NSO2     sulfonamide N
    (76, -0.568, 0.250), // =N-      anionic N=C (fcadj=0.25)
    (77, -0.282, 0.000), // NC=S     thioamide N
    (78, -0.168, 0.000), // NRNH     N=N-H nitrogen
    (79, -0.471, 0.000), // OXN      N-oxide O (alternative)
    (80, -0.144, 0.000), // C6+      cationic aromatic C
    (81, -0.514, 0.000), // N6+      pyridinium N+ (protonated pyridine)
    (82, -0.099, 0.000), // CB6      benzimidazole bridging C
    (83, 0.000, 0.000),  // OXO      placeholder
    (84, 0.000, 0.000),  // placeholder
    (85, 0.000, 0.000),  // placeholder
    (86, 0.000, 0.000),  // placeholder
    (87, 2.000, 0.000),  // NA+2     doubly protonated N (q=+2)
    (88, 3.000, 0.000),  // FE+3     tripositive iron
    (89, -1.000, 0.000), // F-       fluoride anion
    (90, -1.000, 0.000), // CL-      chloride anion
    (91, -1.000, 0.000), // BR-      bromide anion
    (92, 1.000, 0.000),  // LI+      lithium cation
    (93, 1.000, 0.000),  // NA+      sodium cation
    (94, 1.000, 0.000),  // K+       potassium cation
    (95, 2.000, 0.000),  // CA2+     calcium dication
    (96, 2.000, 0.000),  // MG2+     magnesium dication
    (97, 1.000, 0.000),  // ZN2+     zinc dication (formal +1 per ligand)
    (98, 2.000, 0.000),  // ZN+2     zinc dication
    (99, 2.000, 0.000),  // CU+2     copper dication
];

// ── CHG table: (bond_type, a, b, bci) ────────────────────────────────────────
// Convention: bci = charge added to the SECOND atom (b) from bond with first (a).
// Source: RDKit Code/ForceField/MMFF/Params.cpp — defaultMMFFChg
static MMFF94_CHG: &[(u8, u8, u8, f64)] = &[
    (0, 1, 1, 0.0),
    (0, 1, 2, -0.1382),
    (0, 1, 3, -0.061),
    (0, 1, 4, -0.2),
    (0, 1, 5, 0.0),
    (0, 1, 6, -0.28),
    (0, 1, 8, -0.27),
    (0, 1, 9, -0.246),
    (0, 1, 10, -0.3001),
    (0, 1, 11, -0.34),
    (0, 1, 12, -0.29),
    (0, 1, 13, -0.23),
    (0, 1, 14, -0.19),
    (0, 1, 15, -0.23),
    (0, 1, 17, -0.1935),
    (0, 1, 18, -0.1052),
    (0, 1, 19, 0.0805),
    (0, 1, 20, 0.0),
    (0, 1, 22, -0.095),
    (0, 1, 25, 0.0),
    (0, 1, 26, -0.1669),
    (0, 1, 34, -0.503),
    (0, 1, 35, -0.4274),
    (0, 1, 37, -0.1435),
    (0, 1, 39, -0.2556),
    (0, 1, 40, -0.3691),
    (0, 1, 41, 0.106),
    (0, 1, 43, -0.3557),
    (0, 1, 45, -0.2402),
    (0, 1, 46, -0.3332),
    (0, 1, 54, -0.3461),
    (0, 1, 55, -0.4895),
    (0, 1, 56, -0.3276),
    (0, 1, 57, -0.105),
    (0, 1, 58, -0.488),
    (0, 1, 61, -0.2657),
    (0, 1, 62, -0.2),
    (0, 1, 63, -0.18),
    (0, 1, 64, -0.181),
    (0, 1, 67, -0.099),
    (0, 1, 68, -0.256),
    (0, 1, 72, -0.55),
    (0, 1, 73, -0.0877),
    (0, 1, 75, -0.255),
    (0, 1, 78, -0.168),
    (0, 1, 80, -0.144),
    (0, 1, 81, -0.514),
    (0, 2, 2, 0.0),
    (1, 2, 2, 0.0),
    (1, 2, 3, -0.0144),
    (0, 2, 4, -0.065),
    (1, 2, 4, -0.065),
    (0, 2, 5, 0.15),
    (0, 2, 6, -0.0767),
    (1, 2, 9, -0.171),
    (0, 2, 10, -0.109),
    (0, 2, 11, -0.1495),
    (0, 2, 12, -0.14),
    (0, 2, 13, -0.11),
    (0, 2, 14, -0.09),
    (0, 2, 15, -0.101),
    (0, 2, 17, -0.056),
    (0, 2, 18, 0.017),
    (0, 2, 19, 0.229),
    (0, 2, 20, 0.116),
    (0, 2, 22, 0.04),
    (0, 2, 25, 0.147),
    (0, 2, 30, -0.031),
    (0, 2, 34, -0.356),
    (0, 2, 35, -0.35),
    (1, 2, 37, 0.0284),
    (1, 2, 39, 0.031),
    (0, 2, 40, -0.1),
    (0, 2, 41, 0.25),
    (0, 2, 43, -0.191),
    (0, 2, 45, -0.2044),
    (0, 2, 46, -0.294),
    (0, 2, 55, -0.341),
    (0, 2, 56, -0.303),
    (0, 2, 62, -0.05),
    (1, 2, 63, -0.045),
    (1, 2, 64, -0.046),
    (1, 2, 67, 0.036),
    (0, 2, 72, -0.45),
    (1, 2, 81, -0.379),
    (1, 3, 3, 0.0),
    (1, 3, 4, -0.105),
    (0, 3, 5, 0.06),
    (0, 3, 6, -0.15),
    (0, 3, 7, -0.57),
    (0, 3, 9, -0.45),
    (1, 3, 9, -0.211),
    (0, 3, 10, -0.06),
    (0, 3, 11, -0.222),
    (0, 3, 12, -0.209),
    (0, 3, 15, -0.141),
    (0, 3, 16, -0.38),
    (0, 3, 17, -0.096),
    (0, 3, 18, -0.023),
    (0, 3, 20, 0.053),
    (0, 3, 22, 0.0),
    (0, 3, 25, 0.107),
    (1, 3, 30, -0.071),
    (0, 3, 35, -0.361),
    (1, 3, 37, 0.0862),
    (1, 3, 39, -0.009),
    (0, 3, 40, -0.05),
    (0, 3, 41, 0.147),
    (0, 3, 43, -0.2363),
    (0, 3, 45, -0.165),
    (0, 3, 48, -0.43),
    (0, 3, 51, -0.95),
    (0, 3, 53, -0.0134),
    (0, 3, 54, -0.4),
    (1, 3, 54, -0.329),
    (0, 3, 55, -0.381),
    (0, 3, 56, -0.343),
    (1, 3, 57, -0.01),
    (1, 3, 58, -0.393),
    (0, 3, 62, -0.03),
    (1, 3, 63, -0.085),
    (1, 3, 64, -0.086),
    (0, 3, 67, -0.004),
    (0, 3, 74, -0.319),
    (0, 3, 75, -0.2474),
    (1, 3, 78, -0.073),
    (1, 3, 80, -0.049),
    (0, 4, 5, 0.177),
    (0, 4, 6, -0.043),
    (0, 4, 7, -0.487),
    (0, 4, 9, -0.3),
    (1, 4, 9, -0.106),
    (0, 4, 10, -0.044),
    (0, 4, 15, -0.036),
    (0, 4, 20, 0.181),
    (0, 4, 22, 0.105),
    (0, 4, 30, 0.034),
    (1, 4, 37, 0.073),
    (0, 4, 40, -0.064),
    (0, 4, 42, -0.5571),
    (0, 4, 43, -0.126),
    (1, 4, 63, 0.02),
    (1, 4, 64, 0.019),
    (0, 5, 19, 0.2),
    (0, 5, 20, 0.0),
    (0, 5, 22, -0.1),
    (0, 5, 30, -0.15),
    (0, 5, 37, -0.15),
    (0, 5, 41, 0.2203),
    (0, 5, 57, -0.15),
    (0, 5, 63, -0.15),
    (0, 5, 64, -0.15),
    (0, 5, 78, -0.15),
    (0, 5, 80, -0.15),
    (0, 6, 6, 0.0),
    (0, 6, 8, -0.1),
    (0, 6, 9, -0.063),
    (0, 6, 10, 0.0355),
    (0, 6, 15, 0.007),
    (0, 6, 17, 0.052),
    (0, 6, 18, 0.1837),
    (0, 6, 19, 0.2974),
    (0, 6, 20, 0.2579),
    (0, 6, 21, 0.4),
    (0, 6, 22, 0.148),
    (0, 6, 24, 0.5),
    (0, 6, 25, 0.2712),
    (0, 6, 26, 0.101),
    (0, 6, 29, 0.45),
    (0, 6, 30, 0.077),
    (0, 6, 33, 0.5),
    (0, 6, 37, 0.0825),
    (0, 6, 39, 0.139),
    (0, 6, 40, -0.021),
    (0, 6, 41, 0.295),
    (0, 6, 43, -0.083),
    (0, 6, 45, -0.009),
    (0, 6, 54, -0.181),
    (0, 6, 55, -0.233),
    (0, 6, 57, 0.138),
    (0, 6, 58, -0.245),
    (0, 6, 63, 0.063),
    (0, 6, 64, 0.062),
    (0, 7, 17, 0.5),
    (0, 7, 46, 0.1618),
    (0, 7, 74, 0.5),
    (0, 8, 8, 0.0),
    (0, 8, 9, -0.053),
    (0, 8, 10, 0.009),
    (0, 8, 12, -0.051),
    (0, 8, 15, 0.017),
    (0, 8, 17, 0.062),
    (0, 8, 19, 0.347),
    (0, 8, 20, 0.2096),
    (0, 8, 22, 0.158),
    (0, 8, 23, 0.36),
    (0, 8, 25, 0.2679),
    (0, 8, 26, 0.111),
    (0, 8, 34, -0.238),
    (0, 8, 39, 0.149),
    (0, 8, 40, -0.011),
    (0, 8, 43, -0.073),
    (0, 8, 45, -0.007),
    (0, 8, 46, -0.176),
    (0, 8, 55, -0.223),
    (0, 8, 56, -0.185),
    (0, 9, 9, 0.0),
    (0, 9, 10, 0.062),
    (0, 9, 12, 0.002),
    (0, 9, 15, 0.07),
    (0, 9, 18, 0.188),
    (0, 9, 19, 0.4),
    (0, 9, 20, 0.287),
    (0, 9, 25, 0.318),
    (0, 9, 27, 0.4),
    (0, 9, 34, -0.185),
    (0, 9, 35, -0.15),
    (1, 9, 37, 0.179),
    (1, 9, 39, 0.202),
    (0, 9, 40, 0.042),
    (0, 9, 41, 0.358),
    (0, 9, 45, 0.046),
    (0, 9, 53, 0.3179),
    (0, 9, 54, -0.118),
    (0, 9, 55, -0.17),
    (0, 9, 56, -0.132),
    (1, 9, 57, 0.201),
    (0, 9, 62, 0.181),
    (1, 9, 63, 0.126),
    (1, 9, 64, 0.125),
    (0, 9, 67, 0.207),
    (1, 9, 78, 0.138),
    (1, 9, 81, -0.208),
    (0, 10, 10, 0.0),
    (0, 10, 13, 0.006),
    (0, 10, 14, 0.036),
    (0, 10, 15, 0.008),
    (0, 10, 17, 0.053),
    (0, 10, 20, 0.225),
    (0, 10, 22, 0.149),
    (0, 10, 25, 0.256),
    (0, 10, 26, 0.102),
    (0, 10, 28, 0.37),
    (0, 10, 34, -0.247),
    (0, 10, 35, -0.212),
    (0, 10, 37, 0.117),
    (0, 10, 39, 0.14),
    (0, 10, 40, -0.02),
    (0, 10, 41, 0.296),
    (0, 10, 45, -0.016),
    (0, 10, 63, 0.064),
    (0, 10, 64, 0.063),
    (0, 11, 20, 0.298),
    (0, 11, 22, 0.2317),
    (0, 11, 25, 0.329),
    (0, 11, 26, 0.175),
    (0, 11, 37, 0.19),
    (0, 11, 40, 0.053),
    (0, 12, 15, 0.068),
    (0, 12, 18, 0.186),
    (0, 12, 19, 0.3701),
    (0, 12, 20, 0.29),
    (0, 12, 22, 0.2273),
    (0, 12, 25, 0.316),
    (0, 12, 26, 0.2112),
    (0, 12, 37, 0.177),
    (0, 12, 40, 0.04),
    (0, 12, 57, 0.199),
    (0, 12, 63, 0.124),
    (0, 12, 64, 0.123),
    (0, 13, 20, 0.219),
    (0, 13, 22, 0.143),
    (0, 13, 37, 0.111),
    (0, 13, 64, 0.057),
    (0, 14, 20, 0.189),
    (0, 14, 37, 0.081),
    (0, 15, 15, 0.0),
    (0, 15, 18, 0.118),
    (0, 15, 19, 0.33),
    (0, 15, 20, 0.217),
    (0, 15, 22, 0.141),
    (0, 15, 25, 0.248),
    (0, 15, 26, 0.094),
    (0, 15, 30, 0.07),
    (0, 15, 37, 0.1015),
    (0, 15, 40, -0.028),
    (0, 15, 43, -0.09),
    (0, 15, 57, 0.131),
    (0, 15, 63, 0.056),
    (0, 15, 64, 0.055),
    (0, 15, 71, 0.18),
    (0, 16, 16, 0.0),
    (0, 17, 17, 0.0),
    (0, 17, 20, 0.172),
    (0, 17, 22, 0.096),
    (0, 17, 37, 0.064),
    (0, 17, 43, -0.135),
    (0, 18, 18, 0.0),
    (0, 18, 20, 0.099),
    (0, 18, 22, 0.023),
    (0, 18, 32, -0.65),
    (0, 18, 37, -0.009),
    (0, 18, 39, 0.014),
    (0, 18, 43, -0.138),
    (0, 18, 48, -0.5895),
    (0, 18, 55, -0.358),
    (0, 18, 58, -0.37),
    (0, 18, 62, 0.2099),
    (0, 18, 63, -0.062),
    (0, 18, 64, -0.063),
    (0, 18, 80, -0.026),
    (0, 19, 19, 0.0),
    (0, 19, 20, -0.113),
    (0, 19, 37, -0.221),
    (0, 19, 40, -0.358),
    (0, 19, 63, -0.274),
    (0, 19, 75, -0.349),
    (0, 20, 20, 0.0),
    (0, 20, 22, -0.076),
    (0, 20, 25, 0.031),
    (0, 20, 26, -0.123),
    (0, 20, 30, -0.138),
    (0, 20, 34, -0.472),
    (0, 20, 37, -0.108),
    (0, 20, 40, -0.245),
    (0, 20, 41, 0.071),
    (0, 20, 43, -0.307),
    (0, 20, 45, -0.241),
    (0, 22, 22, 0.0),
    (0, 22, 30, -0.071),
    (0, 22, 34, -0.396),
    (0, 22, 37, -0.032),
    (0, 22, 40, -0.169),
    (0, 22, 41, 0.147),
    (0, 22, 43, -0.231),
    (0, 22, 45, -0.165),
    (0, 23, 39, -0.27),
    (0, 23, 62, -0.4),
    (0, 23, 67, -0.292),
    (0, 23, 68, -0.36),
    (0, 25, 25, 0.0),
    (0, 25, 32, -0.7),
    (0, 25, 37, -0.139),
    (0, 25, 39, -0.116),
    (0, 25, 40, -0.276),
    (0, 25, 43, -0.338),
    (0, 25, 57, -0.117),
    (0, 25, 63, -0.192),
    (0, 25, 71, -0.0362),
    (0, 25, 72, -0.6773),
    (0, 26, 26, 0.0),
    (0, 26, 34, -0.349),
    (0, 26, 37, 0.015),
    (0, 26, 40, -0.122),
    (0, 26, 71, 0.096),
    (0, 28, 40, -0.4),
    (0, 28, 43, -0.42),
    (0, 28, 48, -0.4),
    (0, 30, 30, 0.0),
    (0, 30, 40, -0.098),
    (1, 30, 67, 0.067),
    (0, 31, 70, -0.43),
    (0, 32, 41, 0.65),
    (0, 32, 45, 0.52),
    (0, 32, 67, 0.633),
    (0, 32, 68, 0.75),
    (0, 32, 69, 0.75),
    (0, 32, 73, 0.35),
    (0, 32, 77, 0.45),
    (0, 32, 82, 0.633),
    (0, 34, 36, 0.45),
    (0, 34, 37, 0.364),
    (0, 34, 43, 0.165),
    (0, 35, 37, 0.329),
    (0, 35, 63, 0.276),
    (0, 36, 54, -0.4),
    (0, 36, 55, -0.45),
    (0, 36, 56, -0.45),
    (0, 36, 58, -0.457),
    (4, 36, 58, -0.45),
    (0, 36, 81, -0.45),
    (0, 37, 37, 0.0),
    (1, 37, 37, 0.0),
    (0, 37, 38, -0.31),
    (0, 37, 39, 0.023),
    (1, 37, 39, 0.023),
    (0, 37, 40, -0.1),
    (0, 37, 41, 0.179),
    (0, 37, 43, -0.199),
    (0, 37, 45, -0.133),
    (0, 37, 46, -0.302),
    (0, 37, 55, -0.349),
    (0, 37, 56, -0.311),
    (1, 37, 57, 0.022),
    (0, 37, 58, -0.361),
    (1, 37, 58, -0.361),
    (4, 37, 58, -0.35),
    (0, 37, 61, -0.138),
    (0, 37, 62, 0.002),
    (0, 37, 63, 0.0),
    (1, 37, 63, -0.053),
    (0, 37, 64, 0.0),
    (1, 37, 64, -0.054),
    (1, 37, 67, 0.028),
    (0, 37, 69, -0.0895),
    (0, 37, 78, -0.041),
    (0, 37, 81, -0.387),
    (1, 37, 81, -0.387),
    (0, 38, 38, 0.0),
    (0, 38, 63, 0.257),
    (0, 38, 64, 0.256),
    (0, 38, 69, 0.338),
    (0, 38, 78, 0.269),
    (1, 39, 39, 0.0),
    (0, 39, 40, -0.16),
    (0, 39, 45, -0.156),
    (0, 39, 63, -0.1516),
    (1, 39, 63, -0.076),
    (0, 39, 64, -0.077),
    (1, 39, 64, -0.077),
    (0, 39, 65, -0.418),
    (0, 39, 78, -0.064),
    (0, 40, 40, 0.0),
    (0, 40, 45, 0.004),
    (0, 40, 46, -0.165),
    (0, 40, 54, -0.16),
    (0, 40, 63, 0.084),
    (0, 40, 64, 0.083),
    (0, 40, 78, 0.096),
    (0, 41, 41, 0.0),
    (0, 41, 55, -0.528),
    (0, 41, 62, -0.177),
    (0, 41, 72, -0.5),
    (0, 41, 80, -0.196),
    (0, 42, 61, 0.492),
    (0, 43, 43, 0.0),
    (0, 43, 45, 0.066),
    (0, 43, 64, 0.145),
    (0, 44, 63, 0.04),
    (0, 44, 65, -0.2207),
    (0, 44, 78, 0.069),
    (0, 44, 80, 0.093),
    (0, 45, 63, 0.08),
    (0, 45, 64, 0.079),
    (0, 45, 78, 0.092),
    (0, 47, 53, 0.37),
    (0, 49, 50, 0.5673),
    (0, 51, 52, 0.5),
    (0, 55, 57, 0.3544),
    (0, 55, 62, 0.351),
    (0, 55, 64, 0.295),
    (0, 55, 80, 0.332),
    (0, 56, 57, 0.4),
    (0, 56, 63, 0.258),
    (0, 56, 80, 0.27),
    (4, 57, 58, -0.4),
    (1, 57, 63, -0.075),
    (1, 57, 64, -0.076),
    (0, 58, 63, 0.308),
    (0, 58, 64, 0.307),
    (0, 59, 63, 0.14),
    (0, 59, 65, -0.1209),
    (0, 59, 78, 0.169),
    (0, 59, 80, 0.193),
    (0, 59, 82, 0.238),
    (0, 60, 61, 0.37),
    (0, 62, 63, -0.055),
    (0, 62, 64, -0.056),
    (0, 63, 63, 0.0),
    (1, 63, 63, 0.0),
    (0, 63, 64, 0.0),
    (0, 63, 66, -0.3381),
    (0, 63, 72, -0.4),
    (0, 63, 78, 0.012),
    (0, 63, 81, -0.334),
    (0, 64, 64, 0.0),
    (0, 64, 65, -0.2888),
    (0, 64, 66, -0.2272),
    (0, 64, 78, 0.013),
    (0, 64, 81, -0.333),
    (0, 64, 82, 0.082),
    (0, 65, 66, 0.0),
    (0, 65, 78, 0.307),
    (0, 65, 81, -0.039),
    (0, 65, 82, 0.376),
    (0, 66, 66, 0.0),
    (0, 66, 78, 0.299),
    (0, 66, 81, -0.047),
    (0, 71, 75, -0.0958),
    (0, 72, 73, 0.45),
    (0, 76, 76, 0.0),
    (0, 76, 78, 0.4),
    (0, 78, 78, 0.0),
    (1, 78, 78, 0.0),
    (0, 78, 79, -0.303),
    (0, 78, 81, -0.35),
    (0, 79, 81, -0.043),
    (0, 80, 81, -0.4),
];

// ── Lookup helpers ───────────────────────────────────────────────────────────

/// Returns (pbci, fcadj) for the given MMFF94 numeric atom type.
pub fn pbci_for(atom_type: u8) -> (f64, f64) {
    for &(t, pbci, fcadj) in MMFF94_PBCI {
        if t == atom_type {
            return (pbci, fcadj);
        }
    }
    (0.0, 0.0)
}

/// Look up the BCI contribution to `type_i` from a bond with `type_j`.
///
/// Returns `Some(contribution)` if found in CHG table, `None` for PBCI fallback.
/// Convention: entry (bt, a, b, bci) means b gains `bci` from bond with a.
/// - If type_i is b (second): contribution = +bci
/// - If type_i is a (first): contribution = −bci
fn lookup_chg_contribution(bond_type: u8, type_i: u8, type_j: u8) -> Option<f64> {
    // Try type_i as b (second) — entry is (bt, type_j, type_i, bci)
    for &(bt, a, b, bci) in MMFF94_CHG {
        if bt == bond_type && a == type_j && b == type_i {
            return Some(bci); // type_i is the recipient
        }
    }
    // Try type_i as a (first) — entry is (bt, type_i, type_j, bci), negate
    for &(bt, a, b, bci) in MMFF94_CHG {
        if bt == bond_type && a == type_i && b == type_j {
            return Some(-bci); // type_i is the donor, negate
        }
    }
    None
}

// A local `bond_type_for(order: BondOrder) -> u8` used to live here,
// mapping Single/Double/Triple/Aromatic to 0/1/2/4. Removed (issue #227
// Phase 2, BCI investigation): it does not implement RDKit's real MMFF94
// "bond type" concept at all. RDKit's `getMMFFBondType(bond)` (the same
// function `getMMFFBondStretchParams`/`computeMMFFCharges` both call, per a
// direct read of `AtomTyper.cpp:2457-2475,3462-3474` at the pinned commit)
// returns 0 unless the bond is formally SINGLE *and* both atom types are
// flagged `sbmb`/`arom` (RDKit's "single bond between two conjugation-
// capable atoms" special case) -- never a function of bond multiplicity.
// `MMFF94_CHG`'s own bond-type column is `{0, 1, 4}` (verified by a full
// scan, matching a fresh parse of RDKit's pinned `defaultMMFFChg`
// byte-for-byte, 498 rows both sides) -- the 3 `bond_type=4` rows are for
// a single atom-type pair (58/36, 58/37, 58/57) that `getMMFFBondType`
// itself can never produce (confirmed: `getMMFFChgParams` is called from
// exactly one call site in the whole pinned RDKit source, always with
// `getMMFFBondType`'s 0/1 result) -- vestigial, unreachable data RDKit
// itself never queries either. The correct bond-type formula already
// exists, oracle-validated, as [`crate::mmff94_minimizer::bond_type_for`]
// (`(ti, tj, order) -> u8`) -- reused directly below instead of
// reimplementing it a second time.

// ── Atom type assignment ─────────────────────────────────────────────────────

/// Assign MMFF94 numeric atom types (1–99) to all atoms in the molecule.
///
/// This implements the core atom type perception rules for organic chemistry.
/// For atoms not handled, returns `Err`.
///
/// Thin wrapper over [`assign_mmff94_numeric_types_with_view`] that discards
/// the re-perceived molecule -- correct for callers that only need numeric
/// types (charges, vdW, atom-level reporting). Any caller that ALSO does
/// bond-order-dependent classification (`bond_type_for`/`angle_type_for`/
/// `torsion_type_for`/`stretch_bend_type_for`, all of which read
/// [`chematic_core::BondOrder`] directly) must call
/// [`assign_mmff94_numeric_types_with_view`] instead and use its returned
/// molecule for that purpose, not its own input `mol` -- see that function's
/// doc for why (issue #227 Phase 1, torsion parameter gap root cause).
pub fn assign_mmff94_numeric_types(mol: &Molecule) -> Result<Vec<u8>, NumericTypeError> {
    assign_mmff94_numeric_types_with_view(mol).map(|(types, _)| types)
}

/// Like [`assign_mmff94_numeric_types`], but also returns the MMFF-specific
/// re-perceived molecule ([`compute_mmff94_aromatic_view`]'s output, same
/// atom count/bond topology as `mol`, only [`chematic_core::BondOrder`]
/// values on ring bonds can differ) used to derive those types.
///
/// Root cause this exists to fix (issue #227 Phase 1, torsion parameter gap,
/// 2026-08-15): `bond_type_for`'s "`BondOrder::Aromatic` forces
/// `bond_type=0`" rule is itself correct (confirmed against a live RDKit
/// oracle for benzene's own ring bond), but every one of chematic-ff's
/// bond-order-dependent classification call sites
/// (`bond_type_for`/`angle_type_for`/`torsion_type_for`/
/// `stretch_bend_type_for`, both in production energy/gradient code and in
/// coverage-gate diagnostics) was being fed the CALLER'S original `mol`, not
/// this re-perceived view -- even though the numeric TYPES those same call
/// sites use were already correctly derived from it. For a ring system where
/// chematic's general/SMILES aromaticity perception and MMFF94's own
/// stricter, Kekule-based perception (`setMMFFAromaticity`) disagree (e.g.
/// caffeine's pyrimidinedione ring: chematic's general model treats the
/// whole fused bicyclic system as one delocalized aromatic ring, RDKit's
/// real sanitizer Kekulizes that ring to alternating single/double bonds
/// while leaving only the fused imidazole ring aromatic -- oracle-confirmed
/// via `MolFromSmiles(...).GetBondBetweenAtoms(5,6).GetIsAromatic() ==
/// False`, and independently confirmed against chematic's own
/// pre-existing, already-oracle-validated
/// `validation/results/mmff94_aromaticity_bond_parity_227_oracle.json` dump,
/// which already recorded `bond_aromatic["5-6"]: false` for caffeine before
/// this fix), this meant the classification code was computed from the
/// WRONG bond order, landing on a torsion/bond/angle-type code with no table
/// row even though chematic's own, unmodified parameter table already
/// carries the correct row at the code RDKit's real (Kekulized) bond type
/// resolves to. Measured on the 265-molecule Wave 1 corpus: 254/254 of the
/// `torsions_missing` instances with `present_at_different_classification =
/// Some` (`crates/chematic-3d/examples/mmff94_term_coverage_audit.rs`) have
/// this exact shape -- oracle-validated (all 254, not a sample) via a live
/// `GetMMFFTorsionParams` cross-check: the oracle's returned value matches
/// EXACTLY one of chematic's own pre-existing rows at a different
/// classification code in 254/254 cases, and RDKit's real bond object is
/// non-aromatic (`GetIsAromatic() == False`) in 254/254 cases. No Halgren
/// empirical-rule implementation was needed for any of them (see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry for the two
/// falsified alternative hypotheses this ruled out first).
pub fn assign_mmff94_numeric_types_with_view(
    mol: &Molecule,
) -> Result<(Vec<u8>, Molecule), NumericTypeError> {
    let n = mol.atom_count();
    let mut types = vec![0u8; n];
    // MMFF94's aromaticity loop follows RDKit's symmetrized SSSR semantics;
    // the extra same-size representatives matter for degenerate fused ring
    // systems. The general perception APIs continue to use the Horton SSSR.
    let rings = chematic_perception::find_symmetrized_sssr(mol)
        .rings()
        .to_vec();
    // MMFF94 has its own, stricter, Kekule-based aromaticity perception
    // (RDKit's `setMMFFAromaticity`), distinct from chematic's own general
    // Huckel model -- most visibly for "mancude" ring systems where a ring
    // atom carries a genuine exocyclic multiple bond (e.g. caffeine's
    // carbonyl carbons, which chematic's own model still treats as
    // aromatic). All typing below reads aromaticity/bond-order state from
    // this re-perceived view, never from `mol` directly, so every
    // `assign_*_type`/`aromatic_*_type` function's existing `atom.aromatic`
    // and `BondOrder::Aromatic` checks stay correct for MMFF94 purposes
    // without needing their own signatures changed.
    let mmff_mol = compute_mmff94_aromatic_view(mol, &rings)?;

    for (i, ty) in types.iter_mut().enumerate().take(n) {
        let idx = AtomIdx(i as u32);
        let atom = mmff_mol.atom(idx);
        let t = match atom.element {
            Element::C => assign_c_type(&mmff_mol, &rings, idx)?,
            Element::N => assign_n_type(&mmff_mol, &rings, idx)?,
            Element::O => assign_o_type(&mmff_mol, &rings, idx)?,
            Element::S => assign_s_type(&mmff_mol, idx)?,
            Element::P => assign_p_type(&mmff_mol, idx)?,
            Element::SI => 19,
            Element::F => 11,
            Element::CL => 12,
            Element::BR => 13,
            Element::I => 14,
            Element::H => assign_h_type(&mmff_mol, idx)?,
            _ => {
                return Err(NumericTypeError(format!(
                    "unsupported element {:?} at atom {i}",
                    atom.element
                )));
            }
        };
        *ty = t;
    }

    // Construction-time semantic-compatibility invariant (issue #227,
    // Phase 1B-0): a numeric type assigned to an atom must represent the
    // same element the registry says it does. This can only fire if a
    // `assign_*_type` function above has a bug (assigns a type belonging to
    // the wrong element's registry entry) -- it is a defense-in-depth
    // backstop, not the primary correctness mechanism, but it converts any
    // such bug from a silent wrong-parameter collision (issue #227's
    // `furan` finding: a real-but-semantically-wrong table row silently
    // accepted because the lookup only ever checked `Some`/`None`) into a
    // typed, fail-closed `Err` here instead.
    for (i, &assigned) in types.iter().enumerate().take(n) {
        let idx = AtomIdx(i as u32);
        let actual_element = mol.atom(idx).element;
        match chematic_ff_numeric_type_registry_lookup(assigned) {
            Some(info) if info.element == actual_element => {}
            Some(info) => {
                return Err(NumericTypeError(format!(
                    "internal error: atom {i} ({actual_element:?}) was assigned MMFF94 \
                     numeric type {assigned} ({}), whose registry element is {:?} -- \
                     semantic-compatibility invariant violated",
                    info.symbol, info.element
                )));
            }
            None => {
                return Err(NumericTypeError(format!(
                    "internal error: atom {i} ({actual_element:?}) was assigned MMFF94 \
                     numeric type {assigned}, which does not exist in the numeric type \
                     registry"
                )));
            }
        }
    }

    Ok((types, mmff_mol))
}

fn chematic_ff_numeric_type_registry_lookup(
    id: u8,
) -> Option<&'static crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo> {
    crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(id)
}

// ── Helper: bond iteration ───────────────────────────────────────────────────

struct BondInfo {
    neighbor: AtomIdx,
    order: BondOrder,
}

fn bonds_of(mol: &Molecule, idx: AtomIdx) -> Vec<BondInfo> {
    mol.bonds()
        .filter_map(|(_, b)| {
            if b.atom1 == idx {
                Some(BondInfo {
                    neighbor: b.atom2,
                    order: b.order,
                })
            } else if b.atom2 == idx {
                Some(BondInfo {
                    neighbor: b.atom1,
                    order: b.order,
                })
            } else {
                None
            }
        })
        .collect()
}

fn count_bond_order(mol: &Molecule, idx: AtomIdx, order: BondOrder) -> usize {
    bonds_of(mol, idx)
        .iter()
        .filter(|b| b.order == order)
        .count()
}

fn is_bonded_to(mol: &Molecule, idx: AtomIdx, elem: Element, order: BondOrder) -> bool {
    bonds_of(mol, idx)
        .iter()
        .any(|b| mol.atom(b.neighbor).element == elem && b.order == order)
}

// ── Aromatic ring helpers ────────────────────────────────────────────────────
//
// Faithful port of RDKit's `RingMembershipSize`/`isAtomNOxide`
// (`Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`, tag
// `Release_2026_03_3`, commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` --
// see `scripts/mmff94_provenance/PROVENANCE.md`), not a re-derivation. An
// "aromatic ring of size N" requires every consecutive bond in the ring to
// literally be `BondOrder::Aromatic` (chematic's SMILES parser already
// represents aromatic-SMILES ring bonds this way, matching RDKit's
// `Bond::AROMATIC` requirement in `isRingAromatic`), not just that the
// endpoint atoms carry `aromatic = true`.

fn ring_is_fully_aromatic(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        matches!(
            mol.bond_between(a, b).map(|(_, e)| e.order),
            Some(BondOrder::Aromatic)
        )
    })
}

fn atom_in_aromatic_ring_of_size(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    atom: AtomIdx,
    size: usize,
) -> bool {
    rings
        .iter()
        .any(|r| r.len() == size && r.contains(&atom) && ring_is_fully_aromatic(mol, r))
}

/// RDKit's `RingInfo::isAtomInRingOfSize` -- plain ring membership, no
/// aromaticity requirement (used for CR3R/CR4R strained-aliphatic-ring
/// typing, unlike `atom_in_aromatic_ring_of_size` above).
fn atom_in_ring_of_size(rings: &[Vec<AtomIdx>], atom: AtomIdx, size: usize) -> bool {
    rings.iter().any(|r| r.len() == size && r.contains(&atom))
}

fn atoms_share_aromatic_ring_of_size(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    a: AtomIdx,
    b: AtomIdx,
    size: usize,
) -> bool {
    rings.iter().any(|r| {
        r.len() == size && r.contains(&a) && r.contains(&b) && ring_is_fully_aromatic(mol, r)
    })
}

/// RDKit's `getTotalDegree()`: heavy-atom bond count plus implicit/explicit
/// hydrogens.
fn total_degree(mol: &Molecule, idx: AtomIdx) -> usize {
    bonds_of(mol, idx).len() + implicit_hcount(mol, idx) as usize
}

/// Sum of Kekule bond orders plus implicit H count -- RDKit's
/// `getValence(Atom::ValenceType::EXPLICIT) + getNumImplicitHs()`. Requires
/// real (non-`Aromatic`) bond orders on `mol`.
fn total_valence(mol: &Molecule, idx: AtomIdx) -> u32 {
    let explicit: u32 = bonds_of(mol, idx)
        .iter()
        .map(|b| b.order.order_int() as u32)
        .sum();
    explicit + implicit_hcount(mol, idx) as u32
}

/// Partial, behaviorally-calibrated port of RDKit's `MolOps::
/// setMMFFAromaticity` (`Code/GraphMol/Aromaticity.cpp`, pinned commit --
/// see `scripts/mmff94_provenance/PROVENANCE.md`). "Partial" specifically
/// because one gate (the hybridization check just below) is approximated
/// rather than ported; every other rule (ring-by-ring pi-electron counting,
/// the exocyclic-double-bond and NOS lone-pair-bonus rules, the multi-pass
/// resolution loop) is a direct, line-cited port, not a re-derivation.
///
/// This is a **separate, stricter** aromaticity model than chematic's own
/// general Huckel engine (`chematic_perception::apply_aromaticity`) --
/// confirmed empirically, not assumed (issue #227 Priority 1A): re-kekulizing
/// caffeine and re-running chematic's own Huckel model on the result still
/// marks its two carbonyl-bearing ring carbons aromatic (matching general
/// chemistry convention, and RDKit's own *default* aromaticity model), but
/// RDKit's MMFF94-specific perception does not. The reason is the
/// `is_nos_in_ring && !exo_double_bond` rule below: a ring's heteroatom
/// lone-pair pi bonus is withheld from the *entire ring* whenever any ring
/// atom carries a genuine exocyclic multiple bond, which chematic's general
/// model instead treats as a legitimate 0-pi-contributing ring member (still
/// part of the aromatic system, just contributing no electrons) -- a
/// different, also-legitimate aromaticity convention, but not MMFF94's own.
/// This can't be solved by reusing the existing general aromaticity engine;
/// it needs its own pass.
///
/// Operates on a freshly-Kekulized copy of `mol` (needs real Single/Double
/// bond orders, not `BondOrder::Aromatic`). `rings` is the caller's SSSR --
/// atom indices are unaffected by Kekulization, only bond orders change, so
/// the same ring atom-index lists remain valid. Returns a per-atom-index
/// bool, true where MMFF94 considers the atom aromatic.
///
/// Not ported: `getHybridization() != Atom::SP2` (`Aromaticity.cpp` line
/// 1023 -- RDKit's real, general hybridization perception, computed during
/// standard sanitization, not itself MMFF-specific). Approximated as
/// `total_degree(atom) > 3` for ring C/N atoms -- a saturated (4-connected)
/// ring atom can't be SP2, which is the real-world failure mode this RDKit
/// gate exists to catch; a full hybridization inference engine doesn't
/// otherwise exist in chematic to port this faithfully.
///
/// Measured gap (issue #227 Phase 0.3, not assumed): on the 265-molecule
/// Wave 1 corpus, 4,172 ring C/N heavy atoms are subject to this gate.
/// `scripts/mmff94_hybridization_gate_gap_227_report.py` (ground truth:
/// RDKit's real `atom.GetHybridization()` per atom, joined against
/// chematic's own `total_degree` dump) classifies all 4,172 into four
/// exclusive buckets, 0 unclassified:
///   - same_decision: 4,128 (98.9%).
///   - rdkit_rejects_chematic_accepts: 44. All degree-3 ring nitrogens RDKit
///     resolves as pyramidal SP3 (e.g. an N-substituted amine nitrogen
///     fused into a ring) that this approximation's threshold does not
///     catch, since a 3-connected atom can legitimately be either SP2 or
///     SP3 and only real hybridization perception can tell them apart --
///     the approximation under-triggers in exactly this one direction,
///     never the other.
///   - rdkit_accepts_chematic_rejects: 0. The approximation never wrongly
///     rejects a genuinely SP2 ring atom on this corpus.
///   - oracle_unavailable: 0.
///
/// See `validation/results/mmff94_hybridization_gate_gap_227_report.txt`
/// for the full data and example atoms.
///
/// Returns a **re-perceived molecule**, not just a bool vector: Kekulized
/// bond orders everywhere, except that bonds belonging to an MMFF-aromatic
/// ring are promoted to `BondOrder::Aromatic` and their atoms' `.aromatic`
/// flag is set true (every other atom's `.aromatic` flag is explicitly
/// false, overriding whatever the input molecule's own aromaticity model
/// said). This lets every existing `atom.aromatic`/`atom_in_aromatic_ring_
/// of_size`/`bonds_of`-double-bond-counting helper below keep working
/// completely unchanged -- they just operate on this view instead of the
/// caller's original molecule for the MMFF-aromaticity-sensitive decisions.
/// Known limitation shared with the rest of this module's ring handling:
/// inherits whatever `rings` (the caller's SSSR) contains, including its
/// existing fused-ring envelope-vs-component-ring behavior (see
/// `chematic_perception::augmented_ring_set`'s doc comment) -- not a new
/// limitation introduced here.
///
/// Fail-closed on Kekulization failure (issue #227 Priority 1A-1, Phase
/// 0.1): this re-perception fundamentally requires real Single/Double bond
/// orders (see above), so a molecule chematic cannot Kekulize has no valid
/// MMFF view to compute -- returning `Err` here, never silently reusing
/// `mol` itself as a stand-in "MMFF view" (that molecule's atoms/bonds still
/// carry whatever *general* aromaticity perception produced them, which is
/// exactly the model this function exists to NOT use for MMFF94 typing, per
/// the doc above). The caller (`assign_mmff94_numeric_types`) propagates
/// this as a typed `NumericTypeError`, which `chematic-3d`'s
/// `From<NumericTypeError> for ForceFieldBridgeError` already converts
/// structurally (every `NumericTypeError` variant, regardless of message
/// content) to `UnsupportedAtomType` -- never by matching on the error
/// string.
///
/// `pub` (not part of the crate's primary typing API, which only ever needs
/// numeric types) specifically so issue #227's corpus-wide MMFF-aromaticity
/// parity audits (`crates/chematic-3d/examples/
/// mmff94_aromaticity_corpus_parity_dump_227.rs`) can dump this
/// intermediate view's atom/bond flags directly against a live RDKit
/// oracle, rather than only being able to check the final numeric types
/// this view feeds into.
pub fn compute_mmff94_aromatic_view(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
) -> Result<Molecule, NumericTypeError> {
    let n = mol.atom_count();
    if rings.is_empty() {
        return Ok(mol.clone());
    }
    let kmol = match chematic_core::kekulize(mol) {
        Ok(kek) if kek.is_empty() => mol.clone(),
        Ok(kek) => chematic_core::apply_kekule(mol, &kek),
        Err(e) => {
            return Err(NumericTypeError(format!(
                "MMFF94 aromaticity re-perception failed at the Kekulization stage: {} \
                 -- refusing to type-assign this molecule rather than silently reusing its \
                 pre-Kekulization (possibly-aromatic-perceived) form as the MMFF view",
                e.detail
            )));
        }
    };

    let atom_in_any_ring = |a: AtomIdx| -> bool { rings.iter().any(|r| r.contains(&a)) };

    let mut resolved = vec![false; n]; // aromBitVect
    let mut is_arom = vec![false; n];
    // Directly records which *rings* (by index into `rings`) themselves
    // passed the Huckel check below -- issue #227 Phase 0.2. Bond promotion
    // must read this, never be reverse-derived from `is_arom` (see the
    // fused-ring correctness note at its use site further down): in a fused
    // system a ring R can fail its own Huckel check yet have every one of
    // its atoms independently marked aromatic via *other*, unrelated
    // accepted rings that happen to share those atoms -- `is_arom` alone
    // cannot distinguish "this ring passed" from "this ring's atoms all
    // happen to belong to other rings that passed."
    let mut ring_accepted = vec![false; rings.len()];

    let mut old_n_resolved: i64 = -1;
    let mut n_resolved: i64 = 0;
    let max_passes = rings.len() + 2; // defensive bound; progress is monotonic

    for _pass in 0..max_passes {
        if n_resolved <= old_n_resolved {
            break;
        }
        old_n_resolved = n_resolved;

        for (ring_idx, ring) in rings.iter().enumerate() {
            let len = ring.len();
            let mut pi_e: i32 = 0;
            let mut move_to_next_ring = false;
            let mut is_nos_in_ring = false;
            let mut exo_double_bond = false;

            for (j, &atom_idx) in ring.iter().enumerate() {
                if move_to_next_ring {
                    break;
                }
                let atom = kmol.atom(atom_idx);
                let is_divalent_s =
                    atom.element == Element::S && total_degree(&kmol, atom_idx) == 2;
                if atom.element == Element::N || atom.element == Element::O || is_divalent_s {
                    is_nos_in_ring = true;
                }

                let next_idx = ring[(j + 1) % len];
                let ring_bond_order = kmol
                    .bond_between(atom_idx, next_idx)
                    .map(|(_, b)| b.order)
                    .unwrap_or(BondOrder::Single);

                if ring_bond_order == BondOrder::Double {
                    pi_e += 2;
                    continue;
                }

                let is_candidate = atom.element == Element::C
                    || (atom.element == Element::N && total_valence(&kmol, atom_idx) == 4);
                if !is_candidate {
                    continue;
                }

                for nb in bonds_of(&kmol, atom_idx) {
                    if ring.contains(&nb.neighbor) {
                        continue; // looking for exocyclic neighbors only
                    }
                    if nb.order == BondOrder::Single {
                        continue;
                    }
                    if atom_in_any_ring(nb.neighbor) && !resolved[nb.neighbor.0 as usize] {
                        move_to_next_ring = true;
                        break;
                    }
                    if nb.order == BondOrder::Double {
                        if is_arom[nb.neighbor.0 as usize] {
                            pi_e += 1;
                        } else {
                            exo_double_bond = true;
                        }
                    }
                }
            }

            if move_to_next_ring {
                continue;
            }

            let mut can_be_aromatic = true;
            for &atom_idx in ring {
                resolved[atom_idx.0 as usize] = true;
                let atom = kmol.atom(atom_idx);
                if matches!(atom.element, Element::C | Element::N)
                    && total_degree(&kmol, atom_idx) > 3
                {
                    can_be_aromatic = false;
                }
            }
            if !can_be_aromatic {
                continue;
            }

            if is_nos_in_ring && !exo_double_bond && (len % 2 == 1) {
                pi_e += 2;
            }

            if pi_e > 2 && (pi_e - 2) % 4 == 0 {
                ring_accepted[ring_idx] = true;
                for &atom_idx in ring {
                    is_arom[atom_idx.0 as usize] = true;
                }
            }
        }

        n_resolved = rings
            .iter()
            .map(|r| r.iter().filter(|&&a| resolved[a.0 as usize]).count() as i64)
            .sum();

        // RDKit's setMMFFAromaticity also stops once every atom belonging to
        // an SSSR ring has been resolved. This is distinct from the progress
        // plateau check above: a pass may resolve all ring atoms while the
        // aggregate counter would otherwise permit another pass. Keep this
        // condition explicit so termination does not depend on the chosen
        // fixed-point counter representation.
        let arom_rings_all_set = rings
            .iter()
            .all(|ring| ring.iter().all(|&a| resolved[a.0 as usize]));
        if arom_rings_all_set {
            break;
        }
    }

    // A ring's bonds are MMFF-aromatic iff *that ring itself* passed the
    // Huckel check above (`ring_accepted`) -- issue #227 Phase 0.2. This
    // must NOT be reconstructed by checking whether every atom in the ring
    // independently carries `is_arom == true`: in a fused/polycyclic
    // system, a non-aromatic ring's atoms can each separately belong to a
    // different accepted aromatic ring, which would make the ring look
    // "all atoms aromatic" without that ring ever having passed its own
    // Huckel check -- wrongly promoting that non-aromatic ring's own bonds
    // to `BondOrder::Aromatic`. Reading `ring_accepted` (set at the exact
    // moment each ring's own check passed, above) makes this structurally
    // impossible instead of relying on a derived invariant that fused rings
    // can break.
    let mut aromatic_bonds: std::collections::HashSet<(AtomIdx, AtomIdx)> =
        std::collections::HashSet::new();
    for (ring_idx, ring) in rings.iter().enumerate() {
        if ring_accepted[ring_idx] {
            let len = ring.len();
            for (j, &a) in ring.iter().enumerate() {
                let b = ring[(j + 1) % len];
                aromatic_bonds.insert((a.min(b), a.max(b)));
            }
        }
    }

    let mut builder = chematic_core::MoleculeBuilder::new();
    for (idx, atom) in kmol.atoms() {
        let mut atom = atom.clone();
        atom.aromatic = is_arom[idx.0 as usize];
        builder.add_atom(atom);
    }
    for (_, bond) in kmol.bonds() {
        let key = (bond.atom1.min(bond.atom2), bond.atom1.max(bond.atom2));
        let order = if aromatic_bonds.contains(&key) {
            BondOrder::Aromatic
        } else {
            bond.order
        };
        builder
            .add_bond(bond.atom1, bond.atom2, order)
            .expect("duplicate bond during MMFF94 aromaticity re-perception");
    }
    Ok(builder.build())
}

/// RDKit's `isAtomNOxide`: a >=3-connected nitrogen with a terminal
/// (degree-1) oxygen neighbor.
fn is_atom_n_oxide(mol: &Molecule, idx: AtomIdx) -> bool {
    if mol.atom(idx).element != Element::N || total_degree(mol, idx) < 3 {
        return false;
    }
    bonds_of(mol, idx)
        .iter()
        .any(|b| mol.atom(b.neighbor).element == Element::O && total_degree(mol, b.neighbor) == 1)
}

/// The alpha/beta heteroatom classification RDKit's `setMMFFHeavyAtomType`
/// computes once per aromatic 5-ring C/N atom and reuses across several
/// branches. `alpha`/`beta` neighbors are ring-adjacent O/S/(non-N-oxide,
/// 3-connected) N in the *same* 5-membered aromatic ring as `atom`.
struct AlphaBetaHeteroatoms {
    alpha: Vec<AtomIdx>,
    beta: Vec<AtomIdx>,
    is_alpha_os: bool,
    is_beta_os: bool,
    alpha_or_beta_in_same_ring: bool,
}

fn is_alpha_beta_heteroatom_candidate(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    matches!(atom.element, Element::O | Element::S)
        || (atom.element == Element::N && total_degree(mol, idx) == 3 && !is_atom_n_oxide(mol, idx))
}

fn find_alpha_beta_heteroatoms(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    atom: AtomIdx,
) -> AlphaBetaHeteroatoms {
    let mut alpha = Vec::new();
    let mut beta = Vec::new();

    if matches!(mol.atom(atom).element, Element::C | Element::N) {
        for nb in bonds_of(mol, atom) {
            let nbr = nb.neighbor;
            if !atom_in_aromatic_ring_of_size(mol, rings, nbr, 5) {
                continue;
            }
            if atoms_share_aromatic_ring_of_size(mol, rings, atom, nbr, 5)
                && is_alpha_beta_heteroatom_candidate(mol, nbr)
            {
                alpha.push(nbr);
            }
            for nb2 in bonds_of(mol, nbr) {
                let nbr2 = nb2.neighbor;
                if nbr2 == atom {
                    continue;
                }
                if !atom_in_aromatic_ring_of_size(mol, rings, nbr2, 5) {
                    continue;
                }
                if atoms_share_aromatic_ring_of_size(mol, rings, atom, nbr2, 5)
                    && is_alpha_beta_heteroatom_candidate(mol, nbr2)
                {
                    beta.push(nbr2);
                }
            }
        }
    }

    let is_alpha_os = alpha
        .iter()
        .any(|&a| matches!(mol.atom(a).element, Element::O | Element::S));
    let is_beta_os = beta
        .iter()
        .any(|&b| matches!(mol.atom(b).element, Element::O | Element::S));
    let alpha_or_beta_in_same_ring = !alpha.is_empty()
        && !beta.is_empty()
        && alpha.iter().any(|&a| {
            beta.iter()
                .any(|&b| atoms_share_aromatic_ring_of_size(mol, rings, a, b, 5))
        });

    AlphaBetaHeteroatoms {
        alpha,
        beta,
        is_alpha_os,
        is_beta_os,
        alpha_or_beta_in_same_ring,
    }
}

// ── C type assignment ────────────────────────────────────────────────────────

fn assign_c_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);

    // Aromatic: detect 5-ring vs 6-ring (RDKit AtomTyper.cpp, aromatic block)
    if atom.aromatic
        && let Some(t) = aromatic_c_type(mol, rings, idx)
    {
        return Ok(t);
    }
    // Aromatic-flagged but not in a fully-aromatic-bonded 5- or 6-ring (e.g.
    // a Kekule-only representation slipping through, or an aromatic ring
    // size MMFF94 doesn't special-case) -- fall through to the aliphatic
    // rules below rather than silently mis-assign a wrong-element type,
    // matching RDKit's own fallthrough when its aromatic switch doesn't set
    // atomType.

    let double_bonds = count_bond_order(mol, idx, BondOrder::Double);

    // sp carbon: RDKit's real CSP (type 4) rule for a non-aromatic carbon
    // is simply `getTotalDegree() == 2` (`AtomTyper.cpp`, pinned commit
    // `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, lines ~954-960 -- the
    // "2 neighbors" branch reached once the earlier degree-4 and degree-3
    // branches don't match), with no check on which elements the two
    // bonds go to. It covers both true acetylenic carbons (one triple
    // bond leaves exactly one remaining substituent, so degree is always
    // 2) *and* cumulated-double-bond ("allenic") carbons such as the
    // central C of an aryl isothiocyanate's N=C=S (two double bonds, zero
    // remaining substituents, also degree 2). Chematic previously only
    // special-cased `triple_bonds > 0`, so a degree-2 carbon reached via
    // two double bonds instead fell into the `double_bonds > 0`
    // "double-bonded to N/O/P/S" branch below and was mistyped 3 (generic
    // carbonyl-family) instead of 4 -- confirmed live against RDKit on
    // `chembl_tier_b_0071`/`_0082`'s isothiocyanate carbon (issue #337).
    // `total_degree(mol, idx) == 2` is a strict superset of the old
    // `triple_bonds > 0` gate (a carbon triple bond always consumes 3 of
    // its 4 valence units, leaving exactly one more substituent), not a
    // narrower replacement, so this cannot regress any previously-correct
    // triple-bond CSP assignment.
    if total_degree(mol, idx) == 2 {
        return Ok(4); // CSP: acetylenic or cumulated-double-bond ("allenic") carbon
    }

    // sp2 carbon
    if double_bonds > 0 {
        // Issue #227 Priority 1A-2: RDKit's real type-3 "C=O" row is an
        // umbrella covering a 3-connected carbon double-bonded to N, O, P,
        // *or* S (`AtomTyper.cpp` lines 907-943 at the pinned commit,
        // `doubleBondedElement ∈ {7,8,15,16}`), not literal C=O alone --
        // it also names C=N (imine), C=P, and thio- variants explicitly in
        // its own symbol list. Previously only O/S were checked here, so a
        // carbon double-bonded to nitrogen (e.g. the ring carbon of a
        // cyclic hydrazide/dione tautomer, or an amidine carbon) fell
        // through to the generic vinylic type below -- verified this was
        // the single root cause of a 39-atom residual (chematic: C=C(2),
        // RDKit: C=O(3)) via a live RDKit oracle re-measurement, not
        // assumed.
        if is_bonded_to(mol, idx, Element::O, BondOrder::Double)
            || is_bonded_to(mol, idx, Element::S, BondOrder::Double)
            || is_bonded_to(mol, idx, Element::N, BondOrder::Double)
            || is_bonded_to(mol, idx, Element::P, BondOrder::Double)
        {
            return Ok(3); // C=O / C=N / C=P / C=S family (generic carbonyl-like)
        }
        // Otherwise generic sp2/vinylic carbon (C=C).
        return Ok(2); // C=C vinylic
    }

    // sp3: small-ring strain context (RDKit AtomTyper.cpp aliphatic-carbon
    // block, `getTotalDegree() == 4` gate -- 3-membered ring checked before
    // 4-membered, matching RDKit's own if/if (not if/else if) order, though
    // a carbon can't be in both a 3- and 4-ring simultaneously in practice).
    if total_degree(mol, idx) == 4 {
        if atom_in_ring_of_size(rings, idx, 3) {
            return Ok(22); // CR3R
        }
        if atom_in_ring_of_size(rings, idx, 4) {
            return Ok(20); // CR4R
        }
    }
    Ok(1) // CR alkyl carbon
}

/// Faithful port of RDKit's aromatic-carbon cases (`AtomTyper.cpp`
/// `setMMFFHeavyAtomType`, `case 6:` under both the 5-ring and 6-ring
/// aromatic switches). Returns `None` if the atom is flagged aromatic but
/// isn't actually in a fully-aromatic-bonded 5- or 6-membered ring.
///
/// Not ported: CIM+ (type 80, aromatic C between two imidazolium N's) --
/// rare in the corpus this targets and not yet needed to close the
/// dominant gap; falls through to C5/C5A/C5B below instead of being
/// misclassified as a different element (still element-correct, just not
/// maximally specific).
fn aromatic_c_type(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> Option<u8> {
    if atom_in_aromatic_ring_of_size(mol, rings, idx, 5) {
        let het = find_alpha_beta_heteroatoms(mol, rings, idx);

        // General C5: no alpha/beta heteroatoms but ring not all benzene-like,
        // or alpha+beta present but in different rings / neither is O/S.
        if het.alpha.len() == het.beta.len() {
            let surrounded_by_benzene_c = bonds_of(mol, idx).iter().all(|nb| {
                mol.atom(nb.neighbor).element == Element::C
                    && atom_in_aromatic_ring_of_size(mol, rings, nb.neighbor, 6)
            });
            let surrounded_by_arom = bonds_of(mol, idx).iter().all(|nb| {
                !atoms_share_aromatic_ring_of_size(mol, rings, idx, nb.neighbor, 5)
                    || mol.atom(nb.neighbor).aromatic
            });
            if (het.alpha.is_empty()
                && het.beta.is_empty()
                && !surrounded_by_benzene_c
                && surrounded_by_arom)
                || (!het.alpha.is_empty()
                    && !het.beta.is_empty()
                    && (!het.alpha_or_beta_in_same_ring || (!het.is_alpha_os && !het.is_beta_os)))
            {
                return Some(78); // C5: general 5-ring aromatic C
            }
        }
        if !het.alpha.is_empty() && (het.beta.is_empty() || het.is_alpha_os) {
            return Some(63); // C5A: alpha to N/O/S
        }
        if !het.beta.is_empty() && (het.alpha.is_empty() || het.is_beta_os) {
            return Some(64); // C5B: beta to N/O/S
        }
    }

    if atom_in_aromatic_ring_of_size(mol, rings, idx, 6) {
        return Some(37); // CB: benzene/pyridine-ring carbon
    }

    None
}

// ── N type assignment ────────────────────────────────────────────────────────

/// Per-nitrogen-atom aggregate of the structural facts RDKit's `case 7:`
/// 3-connected branch (`AtomTyper.cpp` lines 1093-1325 at the pinned
/// commit) computes across all of a nitrogen's carbon neighbors, used to
/// pick NC=C (40) vs NC=O (10) vs plain NR (8).
struct N3CarbonContext {
    has_carbon_neighbor: bool,
    /// `isNCOorNCS`: at least one carbon neighbor carries its own real
    /// (non-aromatic) double bond to O or S.
    is_carbonyl_like: bool,
    /// The `elementTripleBondedToC == 7` contribution to
    /// `isNSO2orNSO3orNCN`: at least one carbon neighbor is triple-bonded
    /// to a nitrogen (cyano). The sulfonamide (P/S-with->=2-terminal-O)
    /// contribution to the same RDKit flag is a separate, ipso-N-level
    /// check this port does not implement (type 43, a distinct tiny
    /// residual bucket -- see `mmff94_hybridization_gate_gap_227_report.py`-
    /// style audit tooling for issue #227 Priority 1A-2).
    is_cyano_like: bool,
    /// True if *any* carbon neighbor independently qualifies for the
    /// NC=C-family trigger (see [`nc_eq_c_carbon_neighbor_qualifies`]).
    any_carbon_qualifies_nc_eq_c: bool,
}

/// Source-grounded deterministic port of RDKit's `case 7:` 3-connected-
/// nitrogen carbon-neighbor scan (`AtomTyper.cpp` lines 1093-1325 at the
/// pinned commit -- see `scripts/mmff94_provenance/PROVENANCE.md`), used by
/// `assign_n_type` to
/// pick NC=C (type 40: enamine/aniline/amidine/N-C%C nitrogen with a
/// delocalized lone pair) over generic NR (8) or NC=O (10).
///
/// For each carbon neighbor of the ipso nitrogen, RDKit checks (per
/// [`nc_eq_c_carbon_neighbor_qualifies`]) whether that carbon: carries its
/// own real C=O/C=S double bond (excludes NC=C, routes toward NC=O
/// instead); is triple-bonded to another nitrogen (cyano, excludes both);
/// or otherwise qualifies via being an aromatic 6-ring carbon with no
/// attached aromatic O/S, or via its own genuine double bond to
/// carbon/nitrogen/phosphorus (an *aromatic* bond to a plain carbon, or to
/// a nitrogen in exactly one ring, also counts per RDKit's own rule at
/// lines 1124-1130), or via its own triple bond to carbon.
///
/// Not ported: the charged amidinium (`NCN+`, type 55) / guanidinium
/// (`NGD+`, type 56) sub-cases (lines 1197-1206, 1273-1284) -- moot for
/// `assign_n_type`'s actual callers, since every positively-charged
/// nitrogen is already routed to type 34 before reaching this function at
/// all (issue #227's construction-time semantic-compatibility work);
/// verified 0/129 corpus atoms in the NR-vs-NC=C residual this port closes
/// carry a charge, so this omission is not a silent gap for this fix's
/// actual target population.
///
/// Structural divergence point from RDKit's literal C++ (not a
/// simplification of the underlying chemistry, a difference in how a
/// nitrogen with *several* differently-behaved carbon neighbors is
/// resolved): RDKit's loop uses shared mutable variables for
/// `elementDoubleBondedToC`/`isNbrBenzeneC`/`nObondedToC`/`nSbondedToC`
/// that are declared once before the carbon-neighbor loop and are not all
/// reset the same way per iteration -- `nObondedToC`/`nSbondedToC` reset to
/// 0 at the start of each carbon neighbor's processing (so after the loop
/// they only reflect the *last* carbon neighbor visited), while
/// `isNbrBenzeneC` never resets (true if *any* carbon neighbor was a
/// benzene carbon) and `elementDoubleBondedToC` is only overwritten by a
/// later neighbor that itself qualifies. In principle this makes the C++
/// code's outcome depend on adjacency-list iteration order for a nitrogen
/// with several differently-behaved carbon neighbors -- an implementation
/// incidental, not an intentional MMFF rule, and not something chematic's
/// own neighbor iteration order is guaranteed to reproduce anyway. This
/// port instead evaluates each carbon neighbor independently and triggers
/// NC=C if *any one* qualifies (after the shared carbonyl/cyano exclusion
/// gates, which genuinely are OR-accumulated across all neighbors in
/// RDKit's own code too, and are ported as such below) -- well-defined and
/// order-independent by construction.
///
/// Empirically, this theoretical divergence does not appear to be
/// observable for legal-valence neutral organic molecules: 8 constructed
/// multi-carbon-context molecules (nitrogen bonded to two structurally
/// distinct carbons, each combination of qualifying/non-qualifying/
/// carbonyl-blocked) were probed against a live RDKit oracle across 32
/// `Chem.RenumberAtoms` atom orderings apiece (256 trials total) and
/// RDKit's own output never varied by order in any of them. The
/// mechanism appears to be that any aromatic-6-ring carbon neighbor's own
/// ring bonds set `elementDoubleBondedToC` into `{6,7}` regardless of which
/// neighbor is processed last (a neutral aromatic 6-ring's atoms are always
/// C or N), which happens to make the order-sensitive branch practically
/// unreachable in that population. The corresponding test,
/// `nc_eq_c_multi_carbon_context_is_order_independent_and_matches_rdkit`,
/// therefore pins these as exact RDKit-parity regressions (both chematic's
/// own determinism *and* RDKit agreement), not merely a documented
/// divergence -- and this port is identical to RDKit whenever a nitrogen
/// has only one structurally-relevant carbon neighbor, which is the
/// overwhelming majority of real molecules regardless.
fn classify_n_c3_carbon_context(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> N3CarbonContext {
    let mut ctx = N3CarbonContext {
        has_carbon_neighbor: false,
        is_carbonyl_like: false,
        is_cyano_like: false,
        any_carbon_qualifies_nc_eq_c: false,
    };
    for nb in bonds_of(mol, idx) {
        if mol.atom(nb.neighbor).element != Element::C {
            continue;
        }
        ctx.has_carbon_neighbor = true;
        let (carbonyl, cyano, qualifies) =
            nc_eq_c_carbon_neighbor_qualifies(mol, rings, nb.neighbor);
        ctx.is_carbonyl_like |= carbonyl;
        ctx.is_cyano_like |= cyano;
        ctx.any_carbon_qualifies_nc_eq_c |= qualifies;
    }
    ctx
}

/// Evaluates a single carbon neighbor of a 3-connected nitrogen against
/// RDKit's per-carbon structural tests (`AtomTyper.cpp` lines 1093-1188,
/// 1266-1298 at the pinned commit). Returns
/// `(has_own_carbonyl_or_thiocarbonyl, is_cyano_carbon,
/// qualifies_for_nc_eq_c)`. See [`classify_n_c3_carbon_context`]'s doc for
/// how the three are combined and this function's role in that port.
fn nc_eq_c_carbon_neighbor_qualifies(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    c_idx: AtomIdx,
) -> (bool, bool, bool) {
    let is_benzene_c = mol.atom(c_idx).aromatic && atom_in_ring_of_size(rings, c_idx, 6);
    let mut has_own_carbonyl_or_thiocarbonyl = false;
    let mut is_cyano_carbon = false;
    let mut has_aromatic_o_or_s_neighbor = false;
    let mut double_bonded_c_n_or_p = false;
    let mut triple_bonded_c = false;

    for nb in bonds_of(mol, c_idx) {
        let nbr_elem = mol.atom(nb.neighbor).element;
        match nb.order {
            BondOrder::Double => match nbr_elem {
                Element::O | Element::S => has_own_carbonyl_or_thiocarbonyl = true,
                Element::C | Element::N | Element::P => double_bonded_c_n_or_p = true,
                _ => {}
            },
            BondOrder::Triple => {
                if nbr_elem == Element::N {
                    is_cyano_carbon = true;
                }
                if nbr_elem == Element::C {
                    triple_bonded_c = true;
                }
            }
            BondOrder::Aromatic => {
                let counts_as_double = nbr_elem == Element::C
                    || (nbr_elem == Element::N
                        && rings.iter().filter(|r| r.contains(&nb.neighbor)).count() == 1);
                if counts_as_double {
                    double_bonded_c_n_or_p = true;
                }
            }
            _ => {}
        }
        if mol.atom(nb.neighbor).aromatic && matches!(nbr_elem, Element::O | Element::S) {
            has_aromatic_o_or_s_neighbor = true;
        }
    }

    let qualifies = (is_benzene_c && !has_aromatic_o_or_s_neighbor)
        || double_bonded_c_n_or_p
        || triple_bonded_c;

    (has_own_carbonyl_or_thiocarbonyl, is_cyano_carbon, qualifies)
}

fn assign_n_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);

    // Aromatic nitrogen
    if atom.aromatic
        && let Some(t) = aromatic_n_type(mol, rings, idx)
    {
        return Ok(t);
    }

    let double_bonds = count_bond_order(mol, idx, BondOrder::Double);
    let triple_bonds = count_bond_order(mol, idx, BondOrder::Triple);
    let nbrs = bonds_of(mol, idx);
    let degree = total_degree(mol, idx);

    // Terminal (degree-1) nitrogen: nitrile/isocyanide (NSP, type 42) or the
    // terminal nitrogen of an azide/diazo group (NAZT, type 47).
    // Source-grounded port of RDKit's degree-1 nitrogen branch
    // (`AtomTyper.cpp` lines ~1454-1481 at the pinned commit, issue #227) --
    // must run before the generic `triple_bonds > 0 -> 9` fallback below,
    // which previously caught every real nitrile nitrogen (always degree-1
    // in a legal structure) before this more specific check could fire.
    if degree == 1
        && let Some(nb) = nbrs.first()
    {
        if nb.order == BondOrder::Triple {
            return Ok(42); // NSP: nitrile/isocyanide nitrogen
        }
        if mol.atom(nb.neighbor).element == Element::N && total_degree(mol, nb.neighbor) == 2 {
            // ipso is bonded to a 2-connected nitrogen (the azide/diazo
            // center) -- NAZT iff that center's OTHER neighbor is itself a
            // 2-connected nitrogen or a 3-connected carbon.
            let is_azt = bonds_of(mol, nb.neighbor).iter().any(|b2| {
                b2.neighbor != idx
                    && ((mol.atom(b2.neighbor).element == Element::N
                        && total_degree(mol, b2.neighbor) == 2)
                        || (mol.atom(b2.neighbor).element == Element::C
                            && total_degree(mol, b2.neighbor) == 3))
            });
            if is_azt {
                return Ok(47); // NAZT: terminal azide/diazo nitrogen
            }
        }
    }

    // Central, charged, 2-connected cumulated nitrogen (the "=N=" center of
    // an azide/diazo group, type 53). Must run before the generic
    // `charge > 0 -> 34` fallback below, which would otherwise mask it
    // (issue #227's "azide/diazo typing" gap).
    if atom.charge > 0
        && degree == 2
        && nbrs
            .iter()
            .all(|b| mol.atom(b.neighbor).element == Element::N)
    {
        return Ok(53); // =N=: central cumulated nitrogen (azide/diazo)
    }

    // Iminium nitrogen (N+=C, type 54). RDKit's MMFF atom typer checks this
    // before the generic positive-N fallback: a three-connected, positively
    // charged N with total bond order at least four and a real N=C/C=N
    // double bond is the iminium class, unless it is a terminal-oxygen
    // environment handled by the dedicated oxygen/nitro cases below.
    let double_bonded_to_c_or_n = nbrs.iter().any(|b| {
        b.order == BondOrder::Double
            && matches!(mol.atom(b.neighbor).element, Element::C | Element::N)
    });
    let total_bond_order: u32 = nbrs.iter().map(|b| b.order.order_int() as u32).sum();
    let iminium_terminal_o_count = nbrs
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::O && bonds_of(mol, b.neighbor).len() == 1
        })
        .count();
    if atom.charge > 0
        && degree == 3
        && total_bond_order >= 4
        && double_bonded_to_c_or_n
        && iminium_terminal_o_count == 0
    {
        return Ok(54); // N+=C: iminium nitrogen
    }

    // Nitro nitrogen (NO2/NO3, type 45). Also must run before the generic
    // `charge > 0 -> 34` fallback: a nitro N is only ever written
    // charge-separated ([N+](=O)[O-]) by a sanitizable structure, so the
    // generic charge check would otherwise mask it every time (issue #227's
    // "charge-shortcut masking nitro-N" gap).
    let terminal_o_count = nbrs
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::O && bonds_of(mol, b.neighbor).len() == 1
        })
        .count();
    if atom.charge > 0 && terminal_o_count >= 2 {
        return Ok(45); // NO2 / NO3
    }

    // Formal charge: quaternary ammonium / protonated N.
    // Registry-verified: type 34 is NR+ (N+, QUATERNARY N); type 32 is
    // O2CM (O, CARBOXYLATE ANION), an oxygen-only type -- the previous
    // `32` here was exactly the silent element-collision the numeric
    // type registry's construction-time invariant now catches instead
    // of allowing through as a false "success".
    if atom.charge > 0 {
        return Ok(34); // NR+
    }

    // Sulfonamide/sulfonate/phosphonamide nitrogen (NSO2/NSO3, type 43):
    // ipso attached to a P or S bonded to >=2 terminal oxygens. Source-
    // grounded port of the S/P-neighbor half of RDKit's `isNSO2orNSO3orNCN`
    // (`AtomTyper.cpp` lines ~985-1000 at the pinned commit, issue #227) --
    // the cyanamide (N-C%N) half of the same RDKit flag is handled by
    // `ctx.is_cyano_like` in the 3-connected branch below, which already
    // existed but wasn't wired to return 43 until now.
    if nbrs.iter().any(|b| {
        let e = mol.atom(b.neighbor).element;
        (e == Element::P || e == Element::S) && count_terminal_o_neighbors(mol, b.neighbor) >= 2
    }) {
        return Ok(43); // NSO2 / NSO3
    }

    // Nitrile / isocyanide (N≡C). Unreachable for any real (degree-1)
    // nitrile now that the branch above handles it -- kept as a
    // conservative fallback for a hypothetical non-degree-1 triple-bonded N
    // this port hasn't observed in practice.
    if triple_bonds > 0 {
        return Ok(9); // N=C (close approximation for nitrile)
    }

    // N=C or N=N (imine, hydrazone, etc.)
    if double_bonds > 0 {
        return Ok(9); // N=C imine
    }

    // sp3 N, 3-connected: enamine/aniline (NC=C) vs amide (NC=O) vs
    // cyanamide (NC%N) vs plain (NR). Source-grounded deterministic port of
    // RDKit's `case 7:` 3-connected branch (`AtomTyper.cpp` lines 1093-1325
    // at the pinned commit) -- see `classify_n_c3_carbon_context`'s doc for
    // the exact condition and its one documented, empirically-unobserved
    // structural divergence from RDKit's literal C++.
    if total_degree(mol, idx) == 3 {
        let ctx = classify_n_c3_carbon_context(mol, rings, idx);
        if ctx.has_carbon_neighbor {
            if ctx.is_cyano_like {
                return Ok(43); // NC%N: nitrogen attached to a cyano carbon
            }
            if !ctx.is_carbonyl_like && ctx.any_carbon_qualifies_nc_eq_c {
                return Ok(40); // NC=C / NC=N / NC=P / NC%C: deloc. lone pair
            }
            if ctx.is_carbonyl_like {
                return Ok(10); // NC=O / NC=S amide/thioamide nitrogen
            }
        }
    }

    // sp3 N — check if amide (bonded to carbonyl C). Fallback for N atoms
    // not covered by the 3-connected-with-carbon-neighbor branch above
    // (e.g. degree != 3, or a degree-3 N with no carbon neighbor at all).
    let is_amide = nbrs.iter().any(|b| {
        let nbr = mol.atom(b.neighbor);
        nbr.element == Element::C && {
            // Check if that C has a C=O double bond
            bonds_of(mol, b.neighbor).iter().any(|bb| {
                bb.order == BondOrder::Double && mol.atom(bb.neighbor).element == Element::O
            })
        }
    });

    if is_amide {
        return Ok(10); // NC=O amide nitrogen
    }

    Ok(8) // NR plain amine
}

/// Faithful port of RDKit's aromatic-nitrogen cases (`AtomTyper.cpp`
/// `setMMFFHeavyAtomType`, `case 7:` under both the 5-ring and 6-ring
/// aromatic switches). Returns `None` if the atom is flagged aromatic but
/// isn't actually in a fully-aromatic-bonded 5- or 6-membered ring.
///
/// Not ported: the 5-ring N-oxide alpha/beta/other sub-distinction (N5AX
/// vs N5BX vs N5OX) -- RDKit itself collapses all three to numeric type 82,
/// so this does too, faithfully, not as a simplification of RDKit's own
/// behavior.
fn aromatic_n_type(mol: &Molecule, rings: &[Vec<AtomIdx>], idx: AtomIdx) -> Option<u8> {
    if atom_in_aromatic_ring_of_size(mol, rings, idx, 5) {
        if is_atom_n_oxide(mol, idx) {
            return Some(82); // N5AX/N5BX/N5OX, collapsed like RDKit's own code
        }
        let het = find_alpha_beta_heteroatoms(mol, rings, idx);
        if het.alpha.is_empty() && het.beta.is_empty() {
            if total_degree(mol, idx) == 3 {
                return Some(39); // NPYL: pyrrole-type N with pi lone pair
            }
            return Some(76); // N5M: anionic 5-ring aromatic N
        }
        if total_degree(mol, idx) == 3 && het.alpha.len() != het.beta.len() {
            return Some(81); // NIM+/N5A+/N5B+/N5+: positively charged 5-ring N
        }
        if !het.alpha.is_empty() && (het.beta.is_empty() || het.is_alpha_os) {
            return Some(65); // N5A: alpha to N/O/S
        }
        if !het.beta.is_empty() && (het.alpha.is_empty() || het.is_beta_os) {
            return Some(66); // N5B: beta to N/O/S
        }
        if !het.alpha.is_empty() && !het.beta.is_empty() {
            return Some(79); // N5: general 5-ring aromatic N
        }
    }

    if atom_in_aromatic_ring_of_size(mol, rings, idx, 6) {
        if is_atom_n_oxide(mol, idx) {
            return Some(69); // NPOX: pyridinium N-oxide
        }
        if total_degree(mol, idx) == 3 {
            return Some(58); // NPD+: protonated pyridinium N
        }
        return Some(38); // NPYD: neutral pyridine-type N
    }

    None
}

// ── O type assignment ────────────────────────────────────────────────────────

/// A genuinely terminal oxygen in RDKit's sense (`AtomTyper.cpp`'s
/// `atom->getDegree() <= 1` gate at the pinned commit's line 1554, reached
/// only after the 3- and (total-degree-)2-neighbor branches have already
/// been ruled out): at most one heavy neighbor *and* no implicit H. This
/// excludes -OH (1 heavy neighbor + 1 implicit H, `getTotalDegree()==2`)
/// and bridging ether/ester -O- (2 heavy neighbors) alike, leaving only
/// =O / -O⁻-style oxygens whose valence is already fully satisfied by
/// their one explicit bond.
fn is_terminal_o(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.atom(idx).element == Element::O
        && bonds_of(mol, idx).len() <= 1
        && implicit_hcount(mol, idx) == 0
}

/// RDKit's `nObondedToCorNorS`, O contribution: count of terminal oxygens
/// bonded to `central` (`AtomTyper.cpp` lines 1597-1600). Note this counts
/// `central`'s terminal-O neighbors generally, so when called with
/// `central` = the ipso oxygen's own neighbor, it naturally includes the
/// ipso oxygen itself alongside any sibling terminal oxygens -- exactly as
/// RDKit's loop does (it never excludes the atom it was reached from).
fn count_terminal_o_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| is_terminal_o(mol, b.neighbor))
        .count()
}

/// RDKit's `nNbondedToCorNorS`: count of degree-2 (`getTotalDegree()==2`)
/// nitrogens bonded to `central` (`AtomTyper.cpp` lines 1593-1596).
fn count_deg2_n_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::N && total_degree(mol, b.neighbor) == 2
        })
        .count()
}

/// RDKit's `nSbondedToCorNorS`: count of terminal (`getTotalDegree()==1`)
/// sulfurs bonded to `central` (`AtomTyper.cpp` lines 1601-1604).
fn count_terminal_s_neighbors(mol: &Molecule, central: AtomIdx) -> usize {
    bonds_of(mol, central)
        .iter()
        .filter(|b| {
            mol.atom(b.neighbor).element == Element::S && total_degree(mol, b.neighbor) == 1
        })
        .count()
}

/// Source-grounded deterministic port of RDKit's `case 8:` degree-≤1
/// ("terminal") oxygen branch (`AtomTyper.cpp` lines 1554-1737 at the
/// pinned commit -- issue #227 Priority 1A-3). Distinguishes OM (35, oxide
/// oxygen on an unremarkable sp3 C/N or bonded to literal H), O2CM (32, the
/// carboxylate / nitro-nitrate / genuine N-oxide / thiosulfinate /
/// sulfate-sulfonate-sulfonamide-sulfone / phosphate-phosphonate-
/// phosphine-oxide / perchlorate terminal-oxygen union), and O=C (7,
/// generic carbonyl/nitroso/sulfoxide oxygen) for a genuinely terminal
/// oxygen. Returns `None` if `idx` is not terminal or if none of RDKit's
/// conditions fire (falls back to `assign_o_type`'s pre-existing generic
/// handling below).
///
/// This is not a "central element implies type" shortcut: routing depends
/// on this bond's order and on the formal degree/valence of the *central*
/// atom and how many *other* terminal O/S atoms (or degree-2 N atoms) share
/// it, exactly as RDKit computes it -- e.g. a carboxylic acid's neutral
/// C=O is *not* O2CM (its carbon has only one terminal oxygen), only a
/// carboxylate *anion*'s two terminal oxygens are (both share a carbon
/// with two terminal oxygens), matching real MMFF94/RDKit output.
///
/// One RDKit quirk ported faithfully, not smoothed over: for a phosphorus
/// or chlorine central atom, `isPhosphateOrPerchlorateO` really is
/// unconditional in RDKit's own source (`atomicNum==15 || atomicNum==17`,
/// no bond-order or valence check at all) -- this is RDKit's own
/// documented behavior (its type-32 symbol list explicitly names "OP,
/// Oxygen in phosphine oxide"), not a chematic simplification, so a
/// neutral phosphine oxide's O also becomes O2CM here, same as RDKit.
fn classify_terminal_o(mol: &Molecule, idx: AtomIdx) -> Option<u8> {
    let central_bond = bonds_of(mol, idx).into_iter().next()?;
    let central = central_bond.neighbor;
    let bond_order = central_bond.order;
    let central_elem = mol.atom(central).element;

    let n_o = count_terminal_o_neighbors(mol, central);
    let n_n = count_deg2_n_neighbors(mol, central);
    let n_s = count_terminal_s_neighbors(mol, central);

    let mut is_oxide_on_c_or_n = false;
    let mut is_o2cm = false;
    let mut is_carbonyl_like = false;

    match central_elem {
        Element::C => {
            if n_o == 2 {
                is_o2cm = true; // isCarboxylateO
            }
            if bond_order == BondOrder::Double {
                is_carbonyl_like = true; // isCarbonylO
            } else if bond_order == BondOrder::Single && n_o == 1 {
                is_oxide_on_c_or_n = true; // isOxideOBondedToC
            }
        }
        Element::N => {
            if n_o >= 2 {
                is_o2cm = true; // isNitroO (nitro/nitrate)
            }
            if bond_order == BondOrder::Double {
                is_carbonyl_like = true; // isNitrosoO
            } else if bond_order == BondOrder::Single && n_o == 1 {
                // RDKit distinguishes a genuine N-oxide (isNOxideO, ipso
                // -> O2CM) from a plain oxide oxygen on an otherwise-normal
                // trivalent nitrogen (isOxideOBondedToN, ipso -> OM) by
                // whether the *central* N's real bond-order-sum valence is
                // 4 vs 2-or-3. Summing bond orders directly is unsafe here:
                // `central` may be a ring atom whose bonds are still
                // `BondOrder::Aromatic` in this post-`compute_mmff94_aromatic_view`
                // molecule (order_int()==1, undercounting a ring bond that's
                // really a Kekule double bond -- e.g. pyridine N-oxide's
                // ring N would wrongly read total_valence==3, not 4). Formal
                // charge is a bond-order-representation-independent proxy
                // that holds for every standard, RDKit-sanitizable N-oxide
                // depiction (charge-separated N+/O-, the only form RDKit
                // itself produces/accepts): a genuinely oxidized nitrogen
                // carries a +1 formal charge; a merely-oxide-substituted
                // normal trivalent nitrogen does not. A non-charge-separated
                // "neutral pentavalent N=O" spelling, if it parses at all,
                // is out of scope for this port (RDKit's own sanitizer does
                // not produce or expect that form).
                if mol.atom(central).charge > 0 {
                    is_o2cm = true; // isNOxideO: genuine N-oxide
                } else {
                    is_oxide_on_c_or_n = true; // isOxideOBondedToN
                }
            }
        }
        Element::S => {
            if n_s == 1 {
                is_o2cm = true; // isThioSulfinateO
            }
            if bond_order == BondOrder::Single
                || (bond_order == BondOrder::Double && (n_o + n_n) > 1)
            {
                is_o2cm = true; // isSulfateO: sulfate/sulfonate/sulfonamide/sulfone
            } else if bond_order == BondOrder::Double && (n_o + n_n) == 1 {
                is_carbonyl_like = true; // isSulfoxideO
            }
        }
        Element::P | Element::CL => {
            is_o2cm = true; // isPhosphateOrPerchlorateO -- unconditional, see doc above
        }
        Element::H => {
            is_oxide_on_c_or_n = true; // isOxideOBondedToH
        }
        _ => return None,
    }

    if is_oxide_on_c_or_n {
        return Some(35); // OM
    }
    if is_o2cm {
        return Some(32); // O2CM
    }
    if is_carbonyl_like {
        return Some(7); // O=C / O=N / O=S generic
    }
    None
}

fn assign_o_type(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    idx: AtomIdx,
) -> Result<u8, NumericTypeError> {
    // Aromatic 5-ring oxygen (furan-type) -- RDKit's aromatic switch has no
    // 6-ring case for oxygen (no common neutral 6-ring aromatic O in the
    // chemistry this typer covers), so only the 5-ring case is ported.
    if mol.atom(idx).aromatic && atom_in_aromatic_ring_of_size(mol, rings, idx, 5) {
        return Ok(59); // OFUR
    }

    // Terminal oxygen (=O, -O⁻, ...): issue #227 Priority 1A-3's OM(35) vs
    // O2CM(32) vs O=C(7) disambiguation. Falls through to the generic
    // handling below for non-terminal O (ether/ester/alcohol) or the rare
    // terminal O whose central atom isn't C/N/S/P/Cl/H.
    if is_terminal_o(mol, idx)
        && let Some(t) = classify_terminal_o(mol, idx)
    {
        return Ok(t);
    }

    // Double bond to C or N → carbonyl/similar oxygen (type 7)
    if count_bond_order(mol, idx, BondOrder::Double) > 0 {
        return Ok(7); // O=C
    }

    // Anionic O (formal charge -1) not resolved above: phenoxide and
    // similar oxide oxygens on an unremarkable carbon/nitrogen.
    // Registry-verified: type 35 is OM (OXIDE OXYGEN ON SP3 C); type 34
    // is NR+, a nitrogen-only type -- the previous `34` here was the
    // same class of silent element collision fixed in `assign_n_type`
    // above (that function's protonated-N branch had the mirror-image
    // bug: it returned 32, O2CM's id, instead of 34).
    if mol.atom(idx).charge < 0 {
        return Ok(35); // OM
    }

    // Single-bond O (ether, alcohol, ester, amide O)
    Ok(6) // OR
}

// ── S type assignment ────────────────────────────────────────────────────────

/// Source-grounded deterministic port of RDKit's `case 16:` (sulfur)
/// degree-based branch (`AtomTyper.cpp` lines ~1815-1917 at the pinned
/// commit, issue #227). The previous version only counted *explicit double
/// bonds* to oxygen, which missed every charge-separated sulfoxide/sulfone
/// spelling (e.g. `[S+]([O-])`, MMFF94's only valid form for a charged
/// sulfinyl/sulfonyl sulfur) -- this port instead follows RDKit's own
/// degree/terminal-neighbor-counting rule, which is bond-order-separation
/// agnostic by construction (a terminal O counts whether it's reached via a
/// double bond or a charge-separated single bond).
fn assign_s_type(mol: &Molecule, idx: AtomIdx) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);
    if atom.aromatic {
        return Ok(44); // S5 aromatic sulfur (thiophene)
    }

    let degree = total_degree(mol, idx);
    let nbrs = bonds_of(mol, idx);

    // 3- or 4-connected sulfur (sulfoxide/sulfone/thiosulfinate family).
    if degree == 3 || degree == 4 {
        let n_o_or_n_bonded = nbrs
            .iter()
            .filter(|b| {
                (mol.atom(b.neighbor).element == Element::O && bonds_of(mol, b.neighbor).len() == 1)
                    || (mol.atom(b.neighbor).element == Element::N
                        && total_degree(mol, b.neighbor) == 2)
            })
            .count();
        let n_s_bonded = nbrs
            .iter()
            .filter(|b| {
                mol.atom(b.neighbor).element == Element::S && bonds_of(mol, b.neighbor).len() == 1
            })
            .count();
        let c_double_bonded = nbrs
            .iter()
            .any(|b| b.order == BondOrder::Double && mol.atom(b.neighbor).element == Element::C);

        if (degree == 3 && n_o_or_n_bonded == 2 && c_double_bonded) || degree == 4 {
            return Ok(18); // SO2: sulfone sulfur
        }
        if (n_o_or_n_bonded > 0 && n_s_bonded > 0) || (n_o_or_n_bonded == 2 && !c_double_bonded) {
            return Ok(73); // SSOM: anionic thiosulfinate sulfur
        }
        return Ok(17); // S=O / >S=N: sulfoxide sulfur
    }

    // 2-connected sulfur.
    if degree == 2 {
        let o_double_bonded = nbrs
            .iter()
            .any(|b| b.order == BondOrder::Double && mol.atom(b.neighbor).element == Element::O);
        if o_double_bonded {
            return Ok(74); // =S=O: sulfinyl sulfur (e.g. C=S=O)
        }
        return Ok(15); // S: thiol, sulfide, or disulfide
    }

    // 1-connected (terminal) sulfur.
    if degree == 1
        && let Some(nb) = nbrs.first()
    {
        let n_term_s_on_nbr = bonds_of(mol, nb.neighbor)
            .iter()
            .filter(|bb| {
                mol.atom(bb.neighbor).element == Element::S && bonds_of(mol, bb.neighbor).len() == 1
            })
            .count();
        let c_double_bonded =
            nb.order == BondOrder::Double && mol.atom(nb.neighbor).element == Element::C;
        if c_double_bonded && n_term_s_on_nbr != 2 {
            return Ok(16); // S=C: sulfur doubly bonded to carbon
        }
        return Ok(72); // S-P / SM / SSMO: other terminal sulfur
    }

    Ok(15) // conservative fallback (degree 0, e.g. an isolated ion)
}

// ── P type assignment ────────────────────────────────────────────────────────

fn assign_p_type(mol: &Molecule, idx: AtomIdx) -> Result<u8, NumericTypeError> {
    // P with =O → phosphoryl (type 25)
    if is_bonded_to(mol, idx, Element::O, BondOrder::Double) {
        return Ok(25); // PO4
    }
    Ok(20) // P generic sp3
}

// ── H type assignment ────────────────────────────────────────────────────────

fn assign_h_type(mol: &Molecule, idx: AtomIdx) -> Result<u8, NumericTypeError> {
    let nbrs = bonds_of(mol, idx);
    if nbrs.is_empty() {
        return Ok(5); // H_C fallback
    }
    let nbr_atom = mol.atom(nbrs[0].neighbor);

    Ok(match nbr_atom.element {
        Element::C => 5,  // HC  H on carbon
        Element::O => 24, // HOCO H on O in acid/alcohol
        Element::S => 71, // HS  H on sulfur
        Element::N => {
            // Distinguish: amide NH (type 28) vs amine NH (type 23) vs imine=NH (type 27)
            let n_idx = nbrs[0].neighbor;
            let n_atom = mol.atom(n_idx);
            if n_atom.aromatic {
                return Ok(23); // treat as HNR for aromatic NH
            }
            let n_is_amide = bonds_of(mol, n_idx).iter().any(|b| {
                b.order == BondOrder::Single
                    && mol.atom(b.neighbor).element == Element::C
                    && bonds_of(mol, b.neighbor).iter().any(|bb| {
                        bb.order == BondOrder::Double && mol.atom(bb.neighbor).element == Element::O
                    })
            });
            if n_is_amide {
                28 // HNCO H on amide N
            } else if count_bond_order(mol, n_idx, BondOrder::Double) > 0 {
                27 // HN=C H on imine N
            } else {
                23 // HNR  H on amine N
            }
        }
        _ => 5,
    })
}

// ── Partial charge calculation ───────────────────────────────────────────────

/// RDKit's real `computeMMFFCharges` does NOT feed the molecule's raw/literal
/// formal charge (`atom->getFormalCharge()`) into equation 15 -- it first
/// computes a separate, MMFF-atom-TYPE-derived formal charge ("MMFFFormalCharge"
/// in RDKit's own naming) via a dedicated per-type switch statement that runs
/// BEFORE the main charge loop (`AtomTyper.cpp` lines ~3095-3350, "We need to
/// set formal charges upfront"), and uses THAT value everywhere a formal
/// charge is needed: as this atom's own q0 in `(1-M*v)*q0`, and as the
/// NEIGHBOR formal charge source for both the `v*sumFormalCharge` term and
/// the anionic-neighbor-leak adjustment. For most types (every type not
/// named in the switch) the derived charge defaults to 0.0 -- RDKit's own
/// pre-switch initial value -- even when the atom's raw SMILES formal charge
/// is nonzero (e.g. nitro nitrogen type 45, azide types 47/53, sulfoxide
/// type 17: none of these are switch cases, so RDKit treats them as
/// formal-charge-neutral for equation 15's purposes despite how the input
/// structure happened to place the literal charge).
///
/// Issue #227 Phase 2 Step 6: this was the root cause of the 67/6693-atom
/// residual left after the BCI bond-type fix (PR #331) -- see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI section for the
/// full molecule-by-molecule evidence.
///
/// **Faithful but intentionally partial port.** Implemented and
/// independently verified against a live RDKit 2026.03.4 oracle query (not
/// merely re-derived from this port's own output): the unconditional
/// +1/+2/+3/-1 "simple" type groups; the O2CM/SM (types 32/72)
/// neighbor-counting redistribution for its carbon-neighbor
/// (carboxylate/thiocarboxylate), nitro/nitrate-nitrogen-neighbor (type 45,
/// 3-terminal-O case only), and sulfone/sulfonate/sulfonamide-sulfur-neighbor
/// (type 18) branches (synthetic fixtures: acetate, nitrate ion,
/// methanesulfonate, dimethyl sulfone; see the `o2cm_sm_*`/`nitrate_ion_*`/
/// `sulfone_and_sulfonate_*` tests below). Implemented but **not**
/// independently oracle-verified (zero corpus exposure either way, so
/// nothing to falsify against -- see below): type 62's (NM) extra
/// "subtract half of each positively-charged neighbor's derived charge"
/// adjustment (`AtomTyper.cpp` lines ~3378-3383, applied in
/// `mmff94_charges_numeric`'s Step 1, not here, since -- unlike the
/// `isDoubleZero(v)` leak -- this adjustment happens BEFORE the `(1-M*v)`
/// multiplication and therefore cannot be expressed as a separate additive
/// term when `v != 0`, which it does for type 62, `fcadj=0.25`).
///
/// **NOT ported** (falls through to the 0.0 default, same as RDKit's own
/// fallback when no switch condition matches):
/// - O2CM/SM's phosphorus-neighbor (type 25, phosphate/phosphonate O),
///   thiosulfinate-sulfur-neighbor (type 73), and perchlorate-chlorine-
///   neighbor (type 77) branches.
/// - Type 76 (N5M): formal charge shared across N5M nitrogens co-membered in
///   a 5-ring (needs ring perception).
/// - Types 55/56/81 (NIM+/N5A+/N5B+): formal charge averaged across a
///   conjugated cationic-nitrogen network reached via alternating type-57/80
///   carbons (needs BFS over a subset of the molecule graph).
/// - Type 61's diazonium special case (+1 bump when bonded to a type-42
///   diazonium nitrogen).
///
/// None of these gaps are guesses papered over with an untested formula --
/// zero atoms of types 62/76/55/56/61/81, and zero O2CM/SM atoms with a
/// phosphorus/type-73/type-77 neighbor, appear anywhere in the 264-molecule
/// Wave 1 corpus this port was measured against (confirmed by a dedicated,
/// committed, independently re-runnable survey,
/// `crates/chematic-3d/examples/mmff94_fchg_type_exposure_survey_227.rs` --
/// issue #227 Phase 2 Step 6), so none of this can be silently masking a
/// corpus-visible bug. Flagged as a follow-up for whoever next touches
/// MMFF94 formal-charge handling on a corpus that does exercise these
/// types.
///
/// **Blast radius, stated explicitly in both directions**: within the
/// 264-molecule corpus, this fix changes the computed charge for exactly
/// 5/6,693 atoms (see the per-atom-join measurement in
/// `scripts/mmff94_provenance/PROVENANCE.md`) -- small, and the reason no
/// downstream 3D-pipeline re-measurement was run for this step (contrast
/// the prior BCI bond-type fix, which moved 1,620 atoms and produced one
/// genuine new stereo violation). Outside this corpus, the behavioral
/// change is broader than that number suggests: it applies to *any*
/// molecule containing a carboxylate, sulfonate/sulfamate, nitrate,
/// nitro, azide, sulfoxide, or quaternary-ammonium group -- the corpus
/// simply happens to contain none of the anionic O2CM/SM forms (every
/// type-32 atom in it has a sulfone/nitro/sulfoxide neighbor, never a
/// carbon or phosphorus neighbor) and only 3 molecules combining nitro/
/// azide/sulfoxide with the absent-from-switch types this fix's Step
/// 1/Step 3 change directly affects.
fn mmff_derived_formal_charge(mol: &Molecule, types: &[u8], idx: AtomIdx) -> f64 {
    match types[idx.0 as usize] {
        // "Non-complicated" +1/+2/+3/-1 atom types (`AtomTyper.cpp`'s
        // `computeMMFFCharges` switch, cases with a single hardcoded `fChg`
        // assignment) -- independent of the atom's own raw formal charge.
        34 | 49 | 51 | 54 | 58 | 92 | 93 | 94 | 97 => 1.0,
        87 | 95 | 96 | 98 | 99 => 2.0,
        88 => 3.0,
        35 | 62 | 89 | 90 | 91 => -1.0,
        // O2CM (32) / SM (72): formal charge shared/localized across
        // terminal O/S atoms bonded to a common neighbor.
        32 | 72 => o2cm_sm_formal_charge(mol, types, idx),
        _ => 0.0,
    }
}

/// O2CM/SM formal-charge redistribution (`AtomTyper.cpp` lines ~3095-3168,
/// `case 32: case 72:`). Reuses this module's existing terminal-O/S/deg-2-N
/// neighbor counters (`count_terminal_o_neighbors`/`count_terminal_s_neighbors`/
/// `count_deg2_n_neighbors`), the same helpers `classify_terminal_o` already
/// uses to *assign* type 32 in the first place, rather than re-deriving the
/// same counts a second way -- with one known, pre-existing divergence from
/// RDKit's real `nSecNbondedToNbr` inherited by this reuse:
/// `count_deg2_n_neighbors` omits RDKit's `!nbr2Atom->getIsAromatic()`
/// condition (`AtomTyper.cpp` line ~3116), so a degree-2 AROMATIC nitrogen
/// would be counted here where RDKit's real algorithm would not (this would
/// flip the sulfonamide fixup and the type-18 branch's `total` by 1 for such
/// a case). Not changed here (the shared helper's existing, already-shipped
/// type-ASSIGNMENT behavior in `classify_terminal_o` is out of scope for a
/// charge-calculation fix); the full-corpus, zero-regression per-atom join
/// (`scripts/mmff94_provenance/PROVENANCE.md`) is the corpus-level evidence
/// this reuse is safe for every type-18-neighbor atom actually measured.
fn o2cm_sm_formal_charge(mol: &Molecule, types: &[u8], idx: AtomIdx) -> f64 {
    for nbr_bond in bonds_of(mol, idx) {
        let nbr = nbr_bond.neighbor;
        let nbr_elem = mol.atom(nbr).element;
        let nbr_type = types[nbr.0 as usize];
        let n_term_os = count_terminal_o_neighbors(mol, nbr) + count_terminal_s_neighbors(mol, nbr);
        let mut n_sec_n = count_deg2_n_neighbors(mol, nbr);
        // Deprotonated-sulfonamide fixup: a sulfur with 2 terminal O/S and 1
        // secondary N is not treated as having a "replaceable" secondary N.
        if nbr_elem == Element::S && n_term_os == 2 && n_sec_n == 1 {
            n_sec_n = 0;
        }
        if nbr_elem == Element::C && n_term_os > 0 {
            return if n_term_os == 1 {
                -1.0
            } else {
                -((n_term_os - 1) as f64) / (n_term_os as f64)
            };
        }
        if nbr_type == 45 && n_term_os == 3 {
            return -1.0 / 3.0;
        }
        if nbr_type == 18 && n_term_os > 0 {
            let total = n_sec_n + n_term_os;
            return if total == 2 {
                0.0
            } else {
                -((total as f64) - 2.0) / (n_term_os as f64)
            };
        }
        // NOT ported: type-25 (phosphate/phosphonate/phosphine-oxide P) and
        // type-77 (perchlorate Cl) neighbor branches, and type-73
        // (thiosulfinate S) -- see `mmff_derived_formal_charge`'s doc.
    }
    0.0
}

/// Compute MMFF94 partial charges using the full PBCI+CHG tables (Halgren 1996).
///
/// Implements equation 15 from MMFF.V paper. For most neutral organic atoms
/// (fcadj=0, no formal charge), this reduces to:
///   q_i = Σ_{j bonded} bci(j→i)
///
/// Returns per-atom partial charges in units of elementary charge.
///
/// Issue #227 Phase 2 (BCI investigation): the BCI (bond-charge-increment)
/// step below reads bond order from [`assign_mmff94_numeric_types_with_view`]'s
/// re-perceived molecule, the same MMFF-specific Kekulized view Phase 1
/// already threads through `bond_type_for`/`angle_type_for`/
/// `torsion_type_for`/`stretch_bend_type_for` -- not this function's own
/// `mol` argument. RDKit's real `computeMMFFCharges` calls the identical
/// `getMMFFBondType(bond)` its bond-stretch code calls, on the identical
/// (sanitized/Kekulized) `mol` object (`AtomTyper.cpp:3071-3488`, pinned
/// commit) -- there is no separate "charge bond order" in RDKit's own
/// algorithm, so there must not be one here either.
pub fn mmff94_charges_numeric(mol: &Molecule) -> Result<Vec<f64>, NumericTypeError> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let n = mol.atom_count();
    let mut charges = vec![0.0f64; n];

    // Derived MMFF formal charges ("MMFFFormalCharge" in RDKit's own naming)
    // -- NOT the molecule's raw/literal `atom.charge`. See
    // `mmff_derived_formal_charge`'s doc for why this distinction is load-
    // bearing (issue #227 Phase 2 Step 6 BCI residual fix).
    let fchg: Vec<f64> = (0..n)
        .map(|i| mmff_derived_formal_charge(mol, &types, AtomIdx(i as u32)))
        .collect();

    // Step 1: formal charge contribution (scaled by fcadj)
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (_, fcadj) = pbci_for(types[i]);
        let mut q0 = fchg[i];
        // Type 62 (NM, anionic divalent N) special case (`AtomTyper.cpp`
        // lines ~3378-3383): subtract half of each positively-charged
        // neighbor's derived formal charge from q0, purely locally -- this
        // adjustment is never written back to `fchg`, so it is invisible to
        // any OTHER atom reading this atom's charge as a neighbor value
        // (matching RDKit's own `getMMFFFormalCharge`, which always returns
        // the switch-only stored value). Must happen here, before the
        // `(1-M*v)` multiplication below, not as a separate additive term
        // in Step 3 -- unlike the `isDoubleZero(v)` leak, this fires
        // whenever `v != 0` (fcadj(62) = 0.25), so it cannot be expressed
        // as a no-op-when-v-is-zero additive term the way that leak can.
        if types[i] == 62 {
            for b in bonds_of(mol, idx) {
                let nbr_fchg = fchg[b.neighbor.0 as usize];
                if nbr_fchg > 0.0 {
                    q0 -= nbr_fchg / 2.0;
                }
            }
        }
        // (1 - M*v)*q0 simplified for fcadj=0 (most atoms): charge[i] = 0
        // For charged atoms with fcadj > 0:
        let m = bonds_of(mol, idx).len() as f64;
        charges[i] = (1.0 - m * fcadj) * q0;
    }

    // Step 2: BCI contributions from each bond. `mmff_mol` (not `mol`) is
    // the bond-order source -- same reperceived view as Step 1's atom
    // types, and the same view Phase 1 already requires for every other
    // bond-order-dependent MMFF94 classification. `mmff_mol` has the same
    // atom count/bond topology as `mol` (only ring-bond `BondOrder` values
    // can differ), so atom-index-based `types[i]`/`types[j]` stay valid.
    for (_, bond) in mmff_mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let ti = types[i];
        let tj = types[j];
        let bt = crate::mmff94_minimizer::bond_type_for(ti, tj, bond.order);

        // Contribution to atom i
        let ci =
            lookup_chg_contribution(bt, ti, tj).unwrap_or_else(|| pbci_for(ti).0 - pbci_for(tj).0);

        // Contribution to atom j
        let cj =
            lookup_chg_contribution(bt, tj, ti).unwrap_or_else(|| pbci_for(tj).0 - pbci_for(ti).0);

        charges[i] += ci;
        charges[j] += cj;
    }

    // Step 3: formal charge redistribution -- equation 15's `v*sumFormalCharge`
    // term and RDKit's `isDoubleZero(v)` anionic-neighbor-leak adjustment
    // (`AtomTyper.cpp` lines ~3384-3399). These are RDKit's own two
    // MUTUALLY EXCLUSIVE branches on whether the ipso atom's own fcadj is
    // (numerically) zero, not two independent unconditional additions: the
    // leak only fires when fcadj_i == 0, and `v*sumFormalCharge` is
    // multiplied by v so it is a no-op exactly when fcadj_i == 0 anyway --
    // the two forms are equivalent, but computing them as an if/else makes
    // the "leak only applies when v==0" condition explicit rather than
    // relying on a multiply-by-zero coincidence. Both branches read the
    // derived `fchg` (not raw `atom.charge`) for the same reason Step 1
    // does. `ISDOUBLEZERO_EPS` mirrors RDKit's own `isDoubleZero` helper
    // (`Params.h`: `x < 1e-10 && x > -1e-10`), not an arbitrary tolerance.
    const ISDOUBLEZERO_EPS: f64 = 1.0e-10;
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (_, fcadj_i) = pbci_for(types[i]);
        if fcadj_i.abs() < ISDOUBLEZERO_EPS {
            for b in bonds_of(mol, idx) {
                let nbr_fchg = fchg[b.neighbor.0 as usize];
                if nbr_fchg < 0.0 {
                    let deg = bonds_of(mol, b.neighbor).len() as f64;
                    charges[i] += nbr_fchg / (2.0 * deg);
                }
            }
        } else {
            let sum_fc: f64 = bonds_of(mol, idx)
                .iter()
                .map(|b| fchg[b.neighbor.0 as usize])
                .sum();
            charges[i] += fcadj_i * sum_fc;
        }
    }

    Ok(charges)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    #![allow(clippy::needless_range_loop)]

    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    // ── Type assignment tests ────────────────────────────────────────────────

    #[test]
    fn glycine_types_match_mmff94_reference() {
        // AGLYSL01 reference: C1=type1, C2=type3, N1=type8, O5=type6, O6=type7
        // H on C = type5, H on N = type23, H on O = type24
        // SMILES: NCC(=O)O  (heavy atoms: N, C, C, O, O; explicit H via parse)
        let m = mol("NCC(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();

        // Collect heavy atom types by element
        let mut n_types: Vec<u8> = Vec::new();
        let mut c_types: Vec<u8> = Vec::new();
        let mut o_types: Vec<u8> = Vec::new();

        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            match a.element {
                Element::N => n_types.push(types[i]),
                Element::C => c_types.push(types[i]),
                Element::O => o_types.push(types[i]),
                _ => {}
            }
        }

        // Amine N → type 8
        assert!(
            n_types.iter().all(|&t| t == 8),
            "N should be type 8 (NR), got {:?}",
            n_types
        );
        // Should have sp3 C (type 1) and carbonyl C (type 3)
        assert!(
            c_types.contains(&1),
            "should have sp3 C (type 1), got {:?}",
            c_types
        );
        assert!(
            c_types.contains(&3),
            "should have carbonyl C (type 3), got {:?}",
            c_types
        );
        // O=C (type 7) and O-H (type 6)
        assert!(
            o_types.contains(&6),
            "should have OR oxygen (type 6), got {:?}",
            o_types
        );
        assert!(
            o_types.contains(&7),
            "should have O=C oxygen (type 7), got {:?}",
            o_types
        );
    }

    #[test]
    fn benzene_aromatic_c_is_type_37() {
        // Issue #227 Phase 1B-0: was asserting the old wrong value (63,
        // which the real Halgren/RDKit numbering assigns to a 5-ring alpha
        // carbon, not benzene). Verified against a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties` on benzene, this PR's audit).
        let m = mol("c1ccccc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::C {
                assert_eq!(types[i], 37, "benzene C should be type 37 (CB)");
            }
        }
    }

    #[test]
    fn pyridine_n_is_type_38() {
        // Issue #227 Phase 1B-0: was asserting the old wrong value (67,
        // which is not a valid MMFF94 numeric type in the real Halgren/
        // RDKit numbering at all -- chematic invented it). Verified against
        // a live RDKit oracle.
        let m = mol("c1ccncc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 38, "pyridine N should be type 38 (NPYD)");
            }
        }
    }

    #[test]
    fn pyridinium_n_is_type_58() {
        let m = mol("c1cc[nH+]cc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(
                    types[i], 58,
                    "protonated pyridine N should be type 58 (NPD+)"
                );
            }
        }
    }

    #[test]
    fn iminium_n_is_type_54_when_not_aromatic() {
        // Kekulized iminium spelling: the positive, three-connected N has
        // total bond order 4 and a real N=C bond, so RDKit uses N+=C (54).
        let m = mol("C[N+]1=CC=CC=C1");
        let (n_idx, _) = m.atoms().find(|(_, a)| a.element == Element::N).unwrap();
        assert_eq!(assign_n_type(&m, &[], n_idx).unwrap(), 54);
    }

    #[test]
    fn pyrrole_n_is_type_39() {
        let m = mol("c1cc[nH]c1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 39, "pyrrole N should be type 39 (NPYL)");
            }
        }
    }

    #[test]
    fn protonated_amine_n_is_type_34_not_the_o2cm_oxygen_row() {
        // Regression for the mirror-image bug to the furan C-C/N-row
        // collision: `assign_n_type`'s formal-charge branch used to
        // return `32` (O2CM, an OXYGEN type) for a positively-charged
        // nitrogen. The construction-time semantic-compatibility
        // invariant now makes that class of bug impossible to ship
        // silently -- it must resolve to type 34 (NR+), the registry's
        // only nitrogen entry among {32, 34, 35}.
        let m = mol("C[NH3+]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 34, "protonated amine N should be type 34 (NR+)");
            }
        }
    }

    #[test]
    fn carboxylate_anionic_o_is_type_32_not_the_nr_plus_nitrogen_row() {
        // Issue #227 Priority 1A-3 correction: this test used to pin
        // acetate's anionic O to 35 (OM), which was itself an instance of
        // the very O2CM residual this Priority closes -- a carboxylate
        // anion's carbon has *two* terminal oxygen neighbors, so
        // `classify_terminal_o`'s `isCarboxylateO` condition fires and both
        // oxygens resolve to O2CM (32), confirmed against a live RDKit
        // oracle. The original mirror-image bug this test was written to
        // guard against still applies and is re-asserted below:
        // `assign_o_type`'s anionic-oxygen branch used to return `34`
        // (NR+, a NITROGEN type) for a negatively-charged oxygen.
        let m = mol("CC(=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut saw_anionic_o = false;
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::O && a.charge < 0 {
                assert_eq!(
                    types[i], 32,
                    "carboxylate anion's O should be type 32 (O2CM)"
                );
                assert_ne!(
                    types[i], 34,
                    "must never collide with NR+, a nitrogen-only type"
                );
                saw_anionic_o = true;
            }
        }
        assert!(saw_anionic_o, "test fixture must contain an anionic oxygen");

        // Negative control: a simple oxide oxygen on an unremarkable
        // aromatic carbon (that carbon has only ONE terminal oxygen, so
        // the carboxylate condition does not fire) must still resolve to
        // plain OM (35), not O2CM -- confirms the fix didn't just make
        // every anionic oxygen 32.
        let phenoxide = mol("c1ccccc1[O-]");
        let types = assign_mmff94_numeric_types(&phenoxide).unwrap();
        let mut saw_phenoxide_o = false;
        for i in 0..phenoxide.atom_count() {
            let a = phenoxide.atom(AtomIdx(i as u32));
            if a.element == Element::O {
                assert_eq!(
                    types[i], 35,
                    "phenoxide O (no sibling terminal O) should stay type 35 (OM)"
                );
                saw_phenoxide_o = true;
            }
        }
        assert!(saw_phenoxide_o, "test fixture must contain the phenoxide O");
    }

    #[test]
    fn imine_carbon_is_type_3_not_generic_vinylic() {
        // Issue #227 Priority 1A-2: closes a 39-atom residual (chematic:
        // C=C(2), RDKit: C=O(3)) -- RDKit's type-3 row is a "generic
        // carbonyl-like" umbrella for a 3-connected carbon double-bonded
        // to N, O, P, *or* S, not literal C=O alone. Expected values
        // copied verbatim from a live RDKit oracle.
        let ketimine = mol("CC(C)=NC"); // N-isopropylidene methylamine
        let types = assign_mmff94_numeric_types(&ketimine).unwrap();
        assert_eq!(
            types[1], 3,
            "imine carbon (C=N) should be type 3 (C=O family)"
        );

        // Negative control: plain alkene carbon (C=C) must stay type 2.
        let propene = mol("CC=C");
        let types = assign_mmff94_numeric_types(&propene).unwrap();
        assert_eq!(types[1], 2, "alkene carbon should stay type 2 (C=C)");
        assert_eq!(types[2], 2, "alkene carbon should stay type 2 (C=C)");
    }

    #[test]
    fn enamine_and_aniline_nitrogens_are_type_40_not_generic_nr() {
        // Issue #227 Priority 1A-2: closes the 129-atom NR(8)-vs-NC=C(40)
        // residual from Priority 1A. All expected values below are copied
        // verbatim from a live RDKit oracle query
        // (`AllChem.MMFFGetMoleculeProperties`), not derived from this
        // fix's own output.

        // Aniline: N attached to a benzene-ring carbon with no aromatic
        // O/S neighbor of its own -- the "isNbrBenzeneC" case.
        let aniline = mol("Nc1ccccc1");
        let types = assign_mmff94_numeric_types(&aniline).unwrap();
        assert_eq!(types[0], 40, "aniline N should be type 40 (NC=C)");

        // N-methyl vinylamine (enamine): N attached to a carbon genuinely
        // double-bonded to another carbon -- the "elementDoubleBondedToC
        // == 6" case.
        let enamine = mol("C=CNC");
        let types = assign_mmff94_numeric_types(&enamine).unwrap();
        assert_eq!(types[2], 40, "enamine N should be type 40 (NC=C)");

        // N-methylformamidine's exocyclic NH2: attached to a carbon
        // double-bonded to *another nitrogen* -- the "elementDoubleBondedToC
        // == 7" (amidine-like, non-charged) case.
        let amidine_like = mol("C(=NC)N");
        let types = assign_mmff94_numeric_types(&amidine_like).unwrap();
        assert_eq!(
            types[3], 40,
            "amidine-like exocyclic NH2 should be type 40 (NC=N)"
        );

        // Negative control: plain aliphatic amine, no qualifying carbon
        // neighbor at all -- must stay generic NR (8), not regress to 40.
        let propylamine = mol("CCCN");
        let types = assign_mmff94_numeric_types(&propylamine).unwrap();
        assert_eq!(
            types[3], 8,
            "plain aliphatic amine N should stay type 8 (NR)"
        );

        // Negative control: amide N (carbon neighbor has its own real
        // C=O) must still resolve to NC=O (10), never NC=C (40) -- the
        // carbonyl exclusion gate must fire correctly.
        let amide = mol("CC(=O)NC");
        let types = assign_mmff94_numeric_types(&amide).unwrap();
        assert_eq!(types[3], 10, "amide N should stay type 10 (NC=O)");
    }

    // ── NC=C multi-carbon-context order-independence (issue #227,
    //    Priority 1A-2 reinforcement) ────────────────────────────────────

    /// Deterministic xorshift64-based permutation generator. Not a source of
    /// true randomness -- just a fixed, reproducible way to generate many
    /// distinct atom/bond orderings without hand-authoring each one.
    fn deterministic_permutation(n: usize, seed: u64) -> Vec<usize> {
        let mut perm: Vec<usize> = (0..n).collect();
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        for i in (1..n).rev() {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let j = (state as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    /// Rebuilds `mol` with atoms reordered per `atom_perm` (`new_idx ->
    /// old_idx`) and bonds re-added in the order given by `bond_perm`
    /// (indices into the *original* bond list). The two orderings vary
    /// independently, so this probes both atom-relabeling and
    /// bond-insertion-order sensitivity.
    fn rebuild_with_order(mol: &Molecule, atom_perm: &[usize], bond_perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; atom_perm.len()];
        for (new_idx, &old_idx) in atom_perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in atom_perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        let bonds: Vec<_> = mol.bonds().collect();
        for &bidx in bond_perm {
            let (_, bond) = bonds[bidx];
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            builder
                .add_bond(a, b, bond.order)
                .expect("relabeling a valid molecule's own bonds cannot fail");
        }
        builder.build()
    }

    /// Returns the MMFF94 numeric type of the sole nitrogen atom in `mol`.
    /// All fixtures below are constructed to contain exactly one N, so
    /// identifying it by element (a content key) rather than by index
    /// survives atom relabeling -- an index-based lookup would silently
    /// compare the wrong atoms across variants.
    fn sole_nitrogen_type(mol: &Molecule) -> u8 {
        let types = assign_mmff94_numeric_types(mol).unwrap();
        let mut found = None;
        for i in 0..mol.atom_count() {
            if mol.atom(AtomIdx(i as u32)).element == Element::N {
                assert!(found.is_none(), "fixture must have exactly one N atom");
                found = Some(types[i]);
            }
        }
        found.expect("fixture must contain a nitrogen atom")
    }

    /// Verifies `expected_type` for the sole N atom is reproduced
    /// identically across: the original parse, full atom-order reversal, 16
    /// deterministic atom relabelings, a bond-insertion-order shuffle, and
    /// an independently hand-written alternate SMILES spelling of the same
    /// molecule. This is the order-independence evidence
    /// `classify_n_c3_carbon_context`'s doc comment claims for the NC=C
    /// port (issue #227 Priority 1A-2 reinforcement).
    fn assert_n_type_is_order_independent(
        smiles: &str,
        alt_smiles: &str,
        expected_type: u8,
        label: &str,
    ) {
        let base = mol(smiles);
        assert_eq!(
            sole_nitrogen_type(&base),
            expected_type,
            "{label}: original SMILES parse"
        );

        let n = base.atom_count();
        let bond_count = base.bonds().count();
        let identity_bonds: Vec<usize> = (0..bond_count).collect();

        let reversed: Vec<usize> = (0..n).rev().collect();
        let reversed_mol = rebuild_with_order(&base, &reversed, &identity_bonds);
        assert_eq!(
            sole_nitrogen_type(&reversed_mol),
            expected_type,
            "{label}: reversed atom order"
        );

        for seed in 0..16u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            assert_eq!(
                sole_nitrogen_type(&variant),
                expected_type,
                "{label}: atom relabeling seed {seed}"
            );
        }

        let identity_atoms: Vec<usize> = (0..n).collect();
        let bond_perm = deterministic_permutation(bond_count, 0xB0AD);
        let bond_order_variant = rebuild_with_order(&base, &identity_atoms, &bond_perm);
        assert_eq!(
            sole_nitrogen_type(&bond_order_variant),
            expected_type,
            "{label}: bond insertion order shuffle"
        );

        let alt = mol(alt_smiles);
        assert_eq!(
            sole_nitrogen_type(&alt),
            expected_type,
            "{label}: alternate equivalent SMILES"
        );
    }

    #[test]
    fn nc_eq_c_multi_carbon_context_is_order_independent_and_matches_rdkit() {
        // Issue #227 Priority 1A-2 reinforcement: `classify_n_c3_carbon_context`
        // evaluates each carbon neighbor of a 3-connected N independently,
        // a deliberate divergence from RDKit's literal C++ (which threads
        // shared mutable state across its carbon-neighbor loop -- see that
        // function's doc comment for the full source-level analysis). The
        // only case where the two designs *could* observably differ is a
        // nitrogen with >=2 structurally distinct carbon neighbors, so
        // that's exactly what these 4 fixtures exercise. All expected
        // types are copied verbatim from a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties`), cross-checked across 32
        // `Chem.RenumberAtoms` atom orderings per fixture on the RDKit side
        // too (0/32 showed order sensitivity in any fixture) -- so these
        // pin exact parity, not a documented divergence.
        assert_n_type_is_order_independent(
            "CN(C)c1ccccc1",
            "c1ccc(cc1)N(C)C",
            40,
            "N bonded to two non-qualifying methyls + one qualifying aromatic C (aniline)",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C(C)(C)C",
            "c1ccc(cc1)N(C)C(C)(C)C",
            40,
            "N bonded to one qualifying aromatic C + one non-qualifying sp3 tert-butyl C",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C=C",
            "C=CN(C)c1ccccc1",
            40,
            "N bonded to two independently-qualifying carbons (aromatic ring + C=C)",
        );
        assert_n_type_is_order_independent(
            "CN(c1ccccc1)C(=O)C",
            "CC(=O)N(C)c1ccccc1",
            10,
            "N bonded to a qualifying aromatic C + a carbonyl C -- the global \
             isNCOorNCS gate must block NC=C regardless of neighbor order",
        );
    }

    #[test]
    fn caffeine_matches_rdkit_mmff_aromaticity_per_atom() {
        // Regression for issue #227 Priority 1A: caffeine's SMILES writes
        // its whole fused purine-dione system in lowercase aromatic
        // notation, including the two carbonyl-bearing ring carbons. This
        // is legitimate under chematic's own general Huckel model (which
        // still treats a ring carbon with an exocyclic C=O as a valid,
        // 0-pi-contributing aromatic-ring member), but RDKit's separate,
        // stricter MMFF94-specific aromaticity perception is not: the
        // exocyclic double bond blocks that ring's heteroatom lone-pair pi
        // bonus, so the whole 6-membered pyrimidinedione ring fails 4n+2 --
        // while the fused 5-membered imidazole ring (whose own atoms carry
        // no exocyclic multiple bonds) independently still passes. Expected
        // values below are copied verbatim from a live RDKit oracle query
        // (`validation/results/mmff94_rdkit_type_oracle.jsonl`, molecule
        // "caffeine"), not derived from this fix's own output.
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected = [
            (0, 1),   // C, methyl
            (1, 39),  // N, NPYL -- imidazole ring, aromatic
            (2, 63),  // C, C5A  -- imidazole ring, aromatic
            (3, 66),  // N, N5B  -- imidazole ring, aromatic
            (4, 64),  // C, C5B  -- fused atom, aromatic via imidazole ring only
            (5, 63),  // C, C5A  -- fused atom, aromatic via imidazole ring only
            (6, 3),   // C, C=O  -- pyrimidinedione ring, NOT aromatic
            (7, 7),   // O, O=C
            (8, 10),  // N, NC=O -- pyrimidinedione ring, NOT aromatic
            (9, 1),   // C, methyl
            (10, 3),  // C, C=O  -- pyrimidinedione ring, NOT aromatic
            (11, 7),  // O, O=C
            (12, 10), // N, NC=O -- pyrimidinedione ring, NOT aromatic
            (13, 1),  // C, methyl
        ];
        for (i, expected_type) in expected {
            assert_eq!(
                types[i], expected_type,
                "caffeine atom {i}: expected MMFF94 type {expected_type} (RDKit oracle), got {}",
                types[i]
            );
        }
    }

    // ── assign_mmff94_numeric_types_with_view (issue #227 Phase 1) ──────────

    #[test]
    fn caffeine_reperceived_view_kekulizes_the_dione_ring_bond_to_single() {
        // Root-cause regression: the SAME bond the aromaticity test above
        // already showed gets non-aromatic (type 3/type 63, not two
        // aromatic-designated types) numeric types must ALSO carry a
        // non-Aromatic BondOrder in the returned view -- otherwise
        // bond_type_for/torsion_type_for (which read BondOrder, not the
        // numeric type's own registry `arom` flag) still see the wrong
        // thing even though typing itself is correct. Oracle-confirmed via
        // `MolFromSmiles(...).GetBondBetweenAtoms(5,6).GetIsAromatic() ==
        // False` (Release 2026.03.4).
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types[5], 63);
        assert_eq!(types[6], 3);
        let order = view
            .bond_between(AtomIdx(5), AtomIdx(6))
            .expect("atoms 5-6 must be bonded")
            .1
            .order;
        assert_ne!(
            order,
            BondOrder::Aromatic,
            "caffeine's ring-6 C5A-C=O bond must be Kekulized to a real \
             Single/Double order in the MMFF view, not left Aromatic"
        );
        // Original `m` is untouched (`&Molecule`, never mutated) -- the two
        // molecules may legitimately disagree.
        assert_eq!(
            m.bond_between(AtomIdx(5), AtomIdx(6)).unwrap().1.order,
            BondOrder::Aromatic,
            "the caller's original molecule must be left exactly as parsed"
        );
    }

    #[test]
    fn assign_mmff94_numeric_types_is_a_thin_wrapper_over_the_view_variant() {
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let types_only = assign_mmff94_numeric_types(&m).unwrap();
        let (types_with_view, _) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types_only, types_with_view);
    }

    /// Finds the bond order (in the reperceived MMFF view) of the unique
    /// bond between a type-`ta` atom and a type-`tb` atom -- a CONTENT-based
    /// identification, not index-based, so it survives atom relabeling
    /// (types themselves are already independently proven order-independent
    /// by `nc_eq_c_multi_carbon_context_is_order_independent_and_matches_rdkit`
    /// and friends above). Panics if the pair is not found or is ambiguous
    /// (more than one such bond) -- both would make this an unreliable
    /// identity key for the fixture it's called on.
    fn bond_order_between_unique_type_pair(
        mol: &Molecule,
        types: &[u8],
        ta: u8,
        tb: u8,
    ) -> BondOrder {
        let mut found: Option<BondOrder> = None;
        for (_, bond) in mol.bonds() {
            let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            let matches = (types[i] == ta && types[j] == tb) || (types[i] == tb && types[j] == ta);
            if matches {
                assert!(
                    found.is_none(),
                    "type pair ({ta},{tb}) must identify a UNIQUE bond in this fixture"
                );
                found = Some(bond.order);
            }
        }
        found.unwrap_or_else(|| panic!("no bond found between type {ta} and type {tb}"))
    }

    #[test]
    fn caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering() {
        // Issue #227 Phase 1 reviewer follow-up: the deterministic-repeat
        // test above only proves no hidden randomness on a FIXED atom
        // order -- it says nothing about whether `chematic_core::kekulize`
        // (a blossom-matching solver, which CAN have genuine ties between
        // equally-valid Kekule structures for a symmetric ring) might pick
        // a different-but-equally-valid alternating bond pattern depending
        // on atom traversal order, and if so, whether that changes which
        // classification code the C5A(63)-C=O(3) ring-fusion bond this
        // whole fix depends on resolves to. Renumber caffeine's atoms 32
        // ways (deterministic_permutation, already used elsewhere in this
        // file for the same purpose on atom TYPES) and confirm the
        // reperceived BOND ORDER for that specific bond -- identified by
        // its unique (type 63, type 3) content signature, not by index, so
        // relabeling can't fool the check -- is identical every time.
        let base = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let (base_types, base_view) = assign_mmff94_numeric_types_with_view(&base).unwrap();
        let reference = bond_order_between_unique_type_pair(&base_view, &base_types, 63, 3);
        assert_ne!(
            reference,
            BondOrder::Aromatic,
            "sanity: the reference value itself must be the fix's expected outcome"
        );

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let (types, view) = assign_mmff94_numeric_types_with_view(&variant).unwrap();
            let order = bond_order_between_unique_type_pair(&view, &types, 63, 3);
            assert_eq!(
                order, reference,
                "seed {seed}: reperceived bond order for the C5A-C=O ring-fusion \
                 bond must not depend on atom renumbering"
            );
        }
    }

    #[test]
    fn benzene_reperceived_ring_bond_orders_are_invariant_under_atom_renumbering() {
        // Companion to the caffeine test above, for the textbook case of a
        // GENUINE Kekule tie (two equally-valid alternating single/double
        // patterns related by a bond-order swap) rather than caffeine's
        // substituent-constrained ring. Confirms that whichever choice
        // `chematic_core::kekulize` makes, the OUTPUT is uniformly
        // `Aromatic` on all 6 ring bonds (RDKit's real MMFF-aromaticity
        // promotion for an accepted ring, not a residual single/double
        // pattern) and that this holds identically across 32 renumberings --
        // i.e. the Kekule tie, even if the solver's internal choice varies
        // with traversal order, never leaks into the final classification
        // input.
        let base = mol("c1ccccc1");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let (_, view) = assign_mmff94_numeric_types_with_view(&variant).unwrap();
            for (_, bond) in view.bonds() {
                assert_eq!(
                    bond.order,
                    BondOrder::Aromatic,
                    "seed {seed}: every benzene ring bond must resolve to Aromatic \
                     regardless of atom renumbering or the Kekulizer's internal tie-break"
                );
            }
        }
    }

    #[test]
    fn assign_mmff94_numeric_types_with_view_is_deterministic() {
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let (types1, view1) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (types2, view2) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        assert_eq!(types1, types2);
        for (_, b) in view1.bonds() {
            let other = view2.bond_between(b.atom1, b.atom2).unwrap().1;
            assert_eq!(b.order, other.order);
        }
    }

    #[test]
    fn small_ring_sp3_carbons_are_cr3r_and_cr4r() {
        // Regression for issue #227 Priority 1A: RDKit AtomTyper.cpp types
        // any 4-total-degree ring carbon in a 3- or 4-membered ring as
        // CR3R/CR4R, not the generic aliphatic CR -- previously chematic
        // always fell through to type 1 regardless of ring size.
        let cyclopropane = mol("C1CC1");
        let types = assign_mmff94_numeric_types(&cyclopropane).unwrap();
        for i in 0..cyclopropane.atom_count() {
            assert_eq!(types[i], 22, "cyclopropane carbon {i} should be CR3R (22)");
        }

        let cyclobutane = mol("C1CCC1");
        let types = assign_mmff94_numeric_types(&cyclobutane).unwrap();
        for i in 0..cyclobutane.atom_count() {
            assert_eq!(types[i], 20, "cyclobutane carbon {i} should be CR4R (20)");
        }

        let cyclohexane = mol("C1CCCCC1");
        let types = assign_mmff94_numeric_types(&cyclohexane).unwrap();
        for i in 0..cyclohexane.atom_count() {
            assert_eq!(
                types[i], 1,
                "cyclohexane carbon {i} (6-ring, not 3- or 4-) should stay generic CR (1)"
            );
        }
    }

    #[test]
    fn furan_ring_carbons_are_type_63_and_64() {
        let m = mol("c1ccoc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert_eq!(
            c_types,
            vec![63, 64],
            "furan ring carbons should be C5A (alpha, 63) and C5B (beta, 64)"
        );
    }

    #[test]
    fn thiophene_ring_carbons_are_type_63_and_64() {
        let m = mol("c1ccsc1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert_eq!(c_types, vec![63, 64]);
    }

    #[test]
    fn indole_has_both_5ring_and_6ring_aromatic_carbon_types() {
        // Fused bicyclic (5-ring pyrrole-like + 6-ring benzo) -- exercises
        // atom_in_aromatic_ring_of_size across a shared ring-fusion bond,
        // not just a single isolated ring.
        let m = mol("c1ccc2[nH]ccc2c1");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut c_types: Vec<u8> = m
            .bonds()
            .flat_map(|(_, b)| [b.atom1, b.atom2])
            .filter(|&a| m.atom(a).element == Element::C)
            .map(|a| types[a.0 as usize])
            .collect();
        c_types.sort_unstable();
        c_types.dedup();
        assert!(
            c_types.contains(&37),
            "indole must have some 6-ring (CB) carbons: {c_types:?}"
        );
        assert!(
            c_types.iter().any(|t| [63, 64, 78].contains(t)),
            "indole must have some 5-ring aromatic carbons: {c_types:?}"
        );
    }

    #[test]
    fn halogens_map_correctly() {
        let m = mol("CF");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            match a.element {
                Element::F => assert_eq!(types[i], 11),
                Element::C => assert_eq!(types[i], 1),
                _ => {}
            }
        }
        let m2 = mol("CCl");
        let types2 = assign_mmff94_numeric_types(&m2).unwrap();
        for i in 0..m2.atom_count() {
            if m2.atom(AtomIdx(i as u32)).element == Element::CL {
                assert_eq!(types2[i], 12);
            }
        }
    }

    #[test]
    fn amide_n_is_type_10() {
        // Acetamide: NC(=O)C
        let m = mol("NC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::N {
                assert_eq!(types[i], 10, "amide N should be type 10 (NC=O)");
            }
        }
    }

    #[test]
    fn sulfoxide_is_type_17_sulfone_is_type_18() {
        let m_so = mol("CS(=O)C"); // DMSO
        let types_so = assign_mmff94_numeric_types(&m_so).unwrap();
        let m_s2 = mol("CS(=O)(=O)C"); // DMSO2
        let types_s2 = assign_mmff94_numeric_types(&m_s2).unwrap();

        for i in 0..m_so.atom_count() {
            if m_so.atom(AtomIdx(i as u32)).element == Element::S {
                assert_eq!(types_so[i], 17, "DMSO S should be type 17 (SO)");
            }
        }
        for i in 0..m_s2.atom_count() {
            if m_s2.atom(AtomIdx(i as u32)).element == Element::S {
                assert_eq!(types_s2[i], 18, "DMSO2 S should be type 18 (SO2)");
            }
        }
    }

    // ── Issue #227: nitrile/sulfonamide/nitro/azide/charged-sulfoxide
    // typing gaps (5 small pre-existing gaps named in PR #239's own
    // summary, root-caused via `scripts/mmff94_angle_bond_gap_classify.py`'s
    // live-RDKit-oracle atom-type cross-check on the failing mmff94_strict
    // corpus). All 5 expected types verified against a live RDKit 2026.03.3
    // oracle (`AllChem.MMFFGetMoleculeProperties`), not hand-derived.

    #[test]
    fn nitrile_nitrogen_is_type_42_not_generic_imine() {
        // Oracle: CC#N -> C(1), C(4), N(42).
        let m = mol("CC#N");
        assert_eq!(
            sole_nitrogen_type(&m),
            42,
            "nitrile N should be type 42 (NSP)"
        );
    }

    #[test]
    fn sulfonamide_nitrogen_is_type_43() {
        // Oracle: CS(=O)(=O)N -> C(1), S(18), O(32), O(32), N(43).
        let m = mol("CS(=O)(=O)N");
        assert_eq!(
            sole_nitrogen_type(&m),
            43,
            "sulfonamide N should be type 43 (NSO2)"
        );
    }

    #[test]
    fn nitro_nitrogen_is_type_45_not_generic_charged() {
        // Oracle: C[N+](=O)[O-] -> C(1), N(45), O(32), O(32).
        let m = mol("C[N+](=O)[O-]");
        assert_eq!(
            sole_nitrogen_type(&m),
            45,
            "nitro N should be type 45 (NO2)"
        );
    }

    #[test]
    fn azide_terminal_and_central_nitrogens_are_type_47_and_53() {
        // Oracle: CN=[N+]=[N-] -> C(1), N(9, attachment N, unaffected),
        // N(53, central cumulated), N(47, terminal).
        let m = mol("CN=[N+]=[N-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut n_types: Vec<u8> = m
            .atoms()
            .filter(|(_, a)| a.element == Element::N)
            .map(|(i, _)| types[i.0 as usize])
            .collect();
        n_types.sort_unstable();
        assert_eq!(
            n_types,
            vec![9, 47, 53],
            "azide N's should be [9 (attachment), 47 (terminal, NAZT), 53 (central, =N=)]"
        );
    }

    #[test]
    fn charged_sulfoxide_sulfur_is_type_17_not_generic_thioether() {
        // Oracle: C[S+](C)[O-] -> C(1), S(17), C(1), O(32). The charge-
        // separated single-bond S-O spelling is MMFF94's only valid form
        // for a charged sulfoxide/sulfonium-oxide; the previous
        // `assign_s_type` only counted explicit double bonds to O and so
        // fell through to the generic S (15) for this spelling.
        let m = mol("C[S+](C)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let (s_idx, _) = m.atoms().find(|(_, a)| a.element == Element::S).unwrap();
        assert_eq!(
            types[s_idx.0 as usize], 17,
            "charge-separated sulfoxide S should be type 17 (S=O)"
        );
    }

    // ── Charge calculation tests ─────────────────────────────────────────────

    #[test]
    fn charge_sum_equals_formal_charge_methane() {
        let m = mol("C");
        let q = mmff94_charges_numeric(&m).unwrap();
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 0.1, "methane net charge = {total:.4}");
    }

    #[test]
    fn charge_sum_equals_formal_charge_glycine() {
        let m = mol("NCC(=O)O");
        let q = mmff94_charges_numeric(&m).unwrap();
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 0.15, "glycine net charge = {total:.4}");
    }

    #[test]
    fn carbonyl_oxygen_is_most_negative() {
        // Acetone: CC(=O)C — carbonyl O should be the most negative atom
        let m = mol("CC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let q = mmff94_charges_numeric(&m).unwrap();
        let (o_idx, _) = m.atoms().find(|(_, a)| a.element == Element::O).unwrap();
        let o_charge = q[o_idx.0 as usize];
        assert!(
            o_charge < -0.3,
            "ketone O charge = {o_charge:.3}, expected < -0.3"
        );
        // Also verify O is type 7
        assert_eq!(types[o_idx.0 as usize], 7, "ketone O should be type 7");
    }

    #[test]
    fn amine_n_is_negative() {
        let m = mol("CCN");
        let q = mmff94_charges_numeric(&m).unwrap();
        let n_charge = m
            .atoms()
            .find(|(_, a)| a.element == Element::N)
            .map(|(i, _)| q[i.0 as usize])
            .unwrap();
        assert!(
            n_charge < -0.1,
            "amine N charge = {n_charge:.3}, expected negative"
        );
    }

    #[test]
    fn h_on_nitrogen_is_positive() {
        // Explicit H in SMILES so they appear as atoms
        let m = mol("C[NH2]");
        let q = mmff94_charges_numeric(&m).unwrap();
        let types = assign_mmff94_numeric_types(&m).unwrap();
        // Find H atoms bonded to N (type 23)
        let h_charges: Vec<f64> = m
            .atoms()
            .filter(|(i, a)| a.element == Element::H && types[i.0 as usize] == 23)
            .map(|(i, _)| q[i.0 as usize])
            .collect();
        if h_charges.is_empty() {
            // If no explicit H-N atoms, just verify N is negative
            let n_charge = m
                .atoms()
                .find(|(_, a)| a.element == Element::N)
                .map(|(i, _)| q[i.0 as usize])
                .unwrap();
            assert!(
                n_charge < 0.0,
                "amine N charge = {n_charge:.3}, expected negative"
            );
        } else {
            for &hq in &h_charges {
                assert!(hq > 0.05, "H-N charge = {hq:.3}, expected positive");
            }
        }
    }

    #[test]
    fn pbci_table_has_99_entries() {
        assert_eq!(MMFF94_PBCI.len(), 99);
    }

    #[test]
    fn chg_table_has_498_entries() {
        assert_eq!(MMFF94_CHG.len(), 498);
    }

    // ── BCI bond-type-source fix (issue #227 Phase 2) ────────────────────────

    #[test]
    fn acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix() {
        // Regression pin for the Phase 2 BCI fix: acetone's C=O bond (type
        // 3 - type 7, `BondOrder::Double`) has a `MMFF94_CHG` row ONLY at
        // `bond_type=0` (`(0, 3, 7, -0.57)`) -- there is no `(1, 3, 7, ...)`
        // row. The old, removed local `bond_type_for(order)` mapped
        // `Double -> 1`, so this bond used to MISS the table entirely and
        // fall back to the generic `pbci_for(3) - pbci_for(7)` difference,
        // silently wrong (RDKit's own algorithm never sets `bondType=1` for
        // a bond that isn't formally SINGLE, per `getMMFFBondType`,
        // `AtomTyper.cpp:2457-2475`). The fixed `bond_type_for(ti, tj,
        // order)` (`crate::mmff94_minimizer`) maps every Double/Triple/
        // Aromatic bond to 0 unconditionally, landing on the real `(0, 3,
        // 7, -0.57)` row directly. Expected values below are copied
        // verbatim from a live RDKit oracle query
        // (`AllChem.MMFFGetMoleculeProperties(Chem.MolFromSmiles("CC(=O)C")).GetMMFFPartialCharge(i)`,
        // `rdkit==2026.03.4`), not derived from this fix's own output.
        let m = mol("CC(=O)C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 3, 7, 1], "sanity: acetone atom types");
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.061, 0.448, -0.57, 0.061];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "acetone atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn mmff94_charges_numeric_uses_reperceived_view_not_callers_bond_order() {
        // Directly proves the fix's mechanism (not just its net effect
        // above): caffeine's pyrimidinedione-ring C5A-C=O bond is
        // `BondOrder::Aromatic` in the caller's original molecule but
        // Kekulized to a real Single/Double order in
        // `assign_mmff94_numeric_types_with_view`'s returned `mmff_mol`
        // (already pinned by
        // `caffeine_reperceived_view_kekulizes_the_dione_ring_bond_to_single`
        // above). If `mmff94_charges_numeric` were still reading `mol`'s
        // own (Aromatic) bond order, every BCI contribution across that
        // bond would use `bond_type_for(ti, tj, Aromatic)` = 0 -- which
        // happens to coincide with the fixed formula's answer for THIS
        // particular pair, so this test instead checks the mechanism
        // directly: charges must be finite and identical across two calls
        // (determinism), and must differ from a hand-computed
        // "old-formula" charge vector for at least one BCI-sensitive
        // molecule with a real Double/Triple bond -- acetone
        // (`acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix`
        // above) already demonstrates that divergence with an oracle pin;
        // this test only re-confirms determinism survives the view lookup.
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let q1 = mmff94_charges_numeric(&m).unwrap();
        let q2 = mmff94_charges_numeric(&m).unwrap();
        assert_eq!(q1, q2, "mmff94_charges_numeric must be deterministic");
        assert!(
            q1.iter().all(|c| c.is_finite()),
            "all charges must be finite"
        );
    }

    #[test]
    fn mmff94_charges_numeric_is_invariant_under_atom_renumbering() {
        // Issue #227 Phase 2 (mirrors Phase 1's reviewer-requested
        // renumbering-invariance test for the same reperceived-view
        // mechanism, `caffeine_reperceived_bond_order_is_invariant_under_atom_renumbering`
        // above): the fix's BCI lookup now depends on `mmff_mol`'s
        // Kekulized bond order, computed by `chematic_core::kekulize` (a
        // blossom-matching solver that CAN have genuine ties for a
        // symmetric ring). A per-atom charge is not itself a stable
        // renumbering-invariant key (atom indices change), but the
        // molecule's own SORTED multiset of per-atom charges must be --
        // renumbering relabels which slot each charge sits in, never the
        // set of values a real, order-independent computation produces.
        let base = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C"); // caffeine
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let mut reference = mmff94_charges_numeric(&base).unwrap();
        reference.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let mut charges = mmff94_charges_numeric(&variant).unwrap();
            charges.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(charges.len(), reference.len());
            for (c, r) in charges.iter().zip(reference.iter()) {
                assert!(
                    (c - r).abs() < 1e-9,
                    "seed {seed}: sorted charge multiset must match the \
                     original ordering's (up to float-summation-order \
                     noise), got {c} vs {r}"
                );
            }
        }
    }

    // ── Derived-formal-charge fix (issue #227 Phase 2 Step 6) ────────────────
    // Root cause: `mmff94_charges_numeric` was feeding the molecule's raw,
    // literal SMILES formal charge into equation 15 (both as this atom's own
    // q0 and as the neighbor-formal-charge source for the redistribution
    // terms). RDKit's real algorithm uses a separate, MMFF-atom-TYPE-derived
    // formal charge instead (see `mmff_derived_formal_charge`'s doc for the
    // full citation). These tests pin the 3/11 residual molecules this gap
    // actually explains (the other 8/11 are unrelated atom-TYPE-assignment
    // bugs, out of scope here -- see `scripts/mmff94_provenance/PROVENANCE.md`),
    // plus the new O2CM/SM redistribution branches on synthetic fixtures the
    // 264-molecule corpus itself does not exercise.

    #[test]
    fn chembl_tier_b_0080_azide_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Azide N types 47 (NAZT, terminal) and 53 (=N=, central) are absent
        // from RDKit's derived-formal-charge switch entirely, so RDKit
        // treats both as formal-charge-neutral (0.0) for equation 15 despite
        // their raw SMILES charges being -1/+1 -- chematic previously used
        // those raw charges directly. Expected values copied verbatim from
        // the already-committed live RDKit oracle dump
        // (`validation/results/mmff94_bci_charges_227_rdkit_oracle.jsonl`,
        // `chembl_tier_b_0080`, `rdkit==2026.03.4`), not derived from this
        // fix's own output.
        let m = mol("COc1cc2nc(N3CCN(/C(S)=N/c4ccc(N=[N+]=[N-])cc4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 40, 3, 15, 9, 37, 37, 37, 37, 9, 53, 47, 37, 37, 1,
            1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0080 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3691, -0.7882, 0.641,
            -0.141, -0.629, 0.179, 0.0, 0.0, 0.179, -0.4969, 0.6879, -0.37, 0.0, 0.0, 0.3691,
            0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0080 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0159_nitro_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Nitro N (type 45) is also absent from RDKit's derived-formal-charge
        // switch (only the 3-terminal-oxygen NITRATE case is a switch
        // condition, not nitro's 2-oxygen case), and its two O2CM (type 32)
        // oxygens hit no branch of the O2CM/SM redistribution either (their
        // shared neighbor is a type-45 N with only 2 terminal oxygens, not
        // 3) -- so RDKit's derived formal charge is 0.0 for the N and both
        // O's, not the raw +1/0/-1 chematic previously used. Expected values
        // copied verbatim from the live RDKit oracle dump
        // (`chembl_tier_b_0159`), not derived from this fix's own output.
        let m = mol("N[C@@H](CCC(=O)Nc1ccc([N+](=O)[O-])cc1)C(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            8, 1, 1, 1, 3, 7, 10, 37, 37, 37, 37, 45, 32, 32, 37, 37, 3, 7, 6,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0159 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            -0.27, 0.331, 0.0, 0.061, 0.569, -0.57, -0.177, 0.117, 0.0, 0.0, 0.133, 0.907, -0.52,
            -0.52, 0.0, 0.0, 0.659, -0.57, -0.15,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0159 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0161_sulfoxide_charges_match_rdkit_oracle_after_derived_formal_charge_fix() {
        // Sulfoxide S (type 17) is absent from RDKit's derived-formal-charge
        // switch, so its derived charge is 0.0, not the raw +1 chematic
        // previously used (its O2CM neighbor's O also gets 0.0: the
        // sulfoxide S is type 17, not one of the O2CM branch's recognized
        // neighbor types). Expected values copied verbatim from the live
        // RDKit oracle dump (`chembl_tier_b_0161`), not derived from this
        // fix's own output.
        let m = mol("CON(C)C(=O)/C=C/CC1C(=O)N2[C@@H]1[S+]([O-])C(C)(C)[C@@H]2C(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 10, 1, 3, 7, 2, 2, 1, 20, 3, 7, 10, 20, 17, 32, 1, 1, 1, 1, 3, 7, 6,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0161 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3155, -0.3246, 0.3001, 0.6156, -0.57, 0.0144, -0.1382, 0.1382, 0.053, 0.577,
            -0.57, -0.5851, 0.397, 0.1755, -0.541, 0.1935, 0.0, 0.0, 0.3611, 0.659, -0.57, -0.15,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0161 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn nitrobenzene_nitro_group_charges_match_rdkit_oracle() {
        // Minimal, isolated reproduction of the `chembl_tier_b_0159` nitro
        // mechanism above (type-45 N, type-32 O's, no O2CM branch fires).
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`, `MMFFGetMoleculeProperties`), independent of
        // the 264-molecule corpus dump.
        let m = mol("c1ccccc1[N+](=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![37, 37, 37, 37, 37, 37, 45, 32, 32]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.133, 0.907, -0.52, -0.52];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "nitrobenzene atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn nitrate_ion_o2cm_three_oxygen_branch_matches_rdkit_oracle() {
        // New O2CM/SM branch (nitro/nitrate-nitrogen-neighbor, `nbr_type ==
        // 45 && n_term_os == 3`): the 264-molecule corpus contains only
        // 2-terminal-oxygen (nitro) type-45 neighbors, never the
        // 3-terminal-oxygen (nitrate) case, so this fixture independently
        // exercises the branch the corpus can't. All 3 oxygens are
        // symmetry-equivalent and must get an identical shared charge.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("[O-][N+](=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![32, 45, 32, 32]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            -0.6866666666666666,
            1.06,
            -0.6866666666666666,
            -0.6866666666666666,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "nitrate ion atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn sulfone_and_sulfonate_o2cm_type18_branch_matches_rdkit_oracle() {
        // O2CM/SM's sulfone/sulfonate/sulfonamide-sulfur-neighbor branch
        // (`nbr_type == 18`) has two arms: dimethyl sulfone's 2 terminal
        // oxygens (already implicitly corpus-validated -- 34 such atoms
        // across 17 corpus molecules already matched the oracle both before
        // and after this fix, since raw charge 0 and derived charge 0
        // coincide for a plain neutral sulfone) hit the `total == 2 -> 0.0`
        // arm; methanesulfonate's 3 terminal oxygens (its raw charge is
        // split -1/0/0 across otherwise-equivalent atoms, unlike the
        // symmetric derived charge) hit the `total != 2 -> fractional` arm,
        // which the corpus does not exercise. Expected values from a fresh
        // live RDKit oracle query (`rdkit==2026.03.4`).
        let sulfone = mol("CS(=O)(=O)C");
        let sulfone_types = assign_mmff94_numeric_types(&sulfone).unwrap();
        assert_eq!(sulfone_types, vec![1, 18, 32, 32, 1]);
        let sulfone_q = mmff94_charges_numeric(&sulfone).unwrap();
        let sulfone_expected = [0.1052, 1.0896, -0.65, -0.65, 0.1052];
        for (i, exp) in sulfone_expected.iter().enumerate() {
            assert!(
                (sulfone_q[i] - exp).abs() < 1e-6,
                "dimethyl sulfone atom {i}: expected charge {exp}, got {}",
                sulfone_q[i]
            );
        }

        let sulfonate = mol("CS(=O)(=O)[O-]");
        let sulfonate_types = assign_mmff94_numeric_types(&sulfonate).unwrap();
        assert_eq!(sulfonate_types, vec![1, 18, 32, 32, 32]);
        let sulfonate_q = mmff94_charges_numeric(&sulfonate).unwrap();
        let sulfonate_expected = [
            0.1052,
            1.3448,
            -0.8166666666666667,
            -0.8166666666666667,
            -0.8166666666666667,
        ];
        for (i, exp) in sulfonate_expected.iter().enumerate() {
            assert!(
                (sulfonate_q[i] - exp).abs() < 1e-6,
                "methanesulfonate atom {i}: expected charge {exp}, got {}",
                sulfonate_q[i]
            );
        }
    }

    #[test]
    fn o2cm_carboxylate_carbon_neighbor_branch_shares_formal_charge_evenly() {
        // O2CM/SM's carbon-neighbor branch (carboxylate/thiocarboxylate) is
        // not exercised anywhere in the 264-molecule corpus. Unlike the
        // type-45/type-18 branches above, a full end-to-end oracle
        // comparison on a real carboxylate (e.g. acetate, `CC(=O)[O-]`) is
        // confounded by a SEPARATE, pre-existing, out-of-scope gap in
        // `assign_mmff94_numeric_types`: RDKit assigns the carboxylate
        // carbon a dedicated type (41, CO2M/CS2M, `AtomTyper.cpp` lines
        // ~885-895) chematic does not yet implement (chematic assigns the
        // generic type 3 instead), which shifts the BCI bond contribution to
        // the oxygens by a constant, uniform offset -- NOT a bug in this
        // fix's formal-charge redistribution (both oxygens are still
        // affected identically, which is exactly the property being tested
        // here). Documented as a follow-up, not fixed in this PR (per the
        // stop condition: fixing it means touching atom-type assignment).
        //
        // This test instead directly exercises `mmff_derived_formal_charge`
        // (this module's own new function, in scope for a direct test) and
        // checks its output against Halgren's cited formula
        // (`-(n_term_os - 1) / n_term_os` for n_term_os == 2) by hand
        // arithmetic, not RDKit's own output -- for 2 terminal oxygens
        // sharing one carbon, each must get -(2-1)/2 = -0.5.
        let m = mol("CC(=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        // idx 2 = carbonyl O, idx 3 = the explicit [O-] -- both terminal,
        // both bonded to the same carboxylate carbon.
        assert_eq!(m.atom(AtomIdx(2)).element, Element::O);
        assert_eq!(m.atom(AtomIdx(3)).element, Element::O);
        let fchg_2 = mmff_derived_formal_charge(&m, &types, AtomIdx(2));
        let fchg_3 = mmff_derived_formal_charge(&m, &types, AtomIdx(3));
        assert!(
            (fchg_2 - (-0.5)).abs() < 1e-9,
            "carboxylate O (idx 2): expected shared formal charge -0.5, got {fchg_2}"
        );
        assert!(
            (fchg_3 - (-0.5)).abs() < 1e-9,
            "carboxylate O (idx 3): expected shared formal charge -0.5, got {fchg_3}"
        );
        assert_eq!(
            fchg_2, fchg_3,
            "both terminal oxygens on the same carboxylate carbon must share \
             the formal charge identically, regardless of which one the \
             input SMILES happened to write the literal '-' charge on"
        );
    }

    #[test]
    fn mmff94_charges_numeric_derived_formal_charge_is_invariant_under_atom_renumbering() {
        // Mirrors `mmff94_charges_numeric_is_invariant_under_atom_renumbering`
        // above, but with a fixture that actually exercises the new
        // derived-formal-charge code path (caffeine, used there, has no
        // type-32/45/47/53/17 atoms and never touches
        // `mmff_derived_formal_charge`'s non-default arms). Nitrobenzene's
        // nitro group does: type 45 N (fcadj-relevant leak source) and type
        // 32 O's (O2CM redistribution, here landing on the "no branch
        // matches" default).
        let base = mol("c1ccccc1[N+](=O)[O-]");
        let n = base.atom_count();
        let identity_bonds: Vec<usize> = (0..base.bonds().count()).collect();

        let mut reference = mmff94_charges_numeric(&base).unwrap();
        reference.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for seed in 0..32u64 {
            let perm = deterministic_permutation(n, seed);
            let variant = rebuild_with_order(&base, &perm, &identity_bonds);
            let mut charges = mmff94_charges_numeric(&variant).unwrap();
            charges.sort_by(|a, b| a.partial_cmp(b).unwrap());
            assert_eq!(charges.len(), reference.len());
            for (c, r) in charges.iter().zip(reference.iter()) {
                assert!(
                    (c - r).abs() < 1e-9,
                    "seed {seed}: sorted charge multiset must match the \
                     original ordering's (up to float-summation-order \
                     noise), got {c} vs {r}"
                );
            }
        }
    }

    #[test]
    fn glycine_h_types_correct() {
        // H on C → 5, H on N → 23, H on O → 24
        // Use explicit SMILES with H: [NH2]CC(=O)O
        // But standard parse leaves H implicit. Let's check that H type assignment works:
        let m = mol("[NH2]CC(=O)O");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::H {
                let t = types[i];
                assert!(
                    matches!(t, 5 | 23 | 24),
                    "H type should be 5/23/24, got {t}"
                );
            }
        }
    }

    #[test]
    fn furan_o_is_type_59() {
        // Issue #227 Phase 1B-0: this test used to hedge between 43 and the
        // old (wrong) fallback value 6, because aromatic-O detection wasn't
        // implemented at all. Now implemented (RDKit's real type for
        // aromatic 5-ring O is 59, OFUR -- 43 was never correct; that's
        // NSO2, sulfonamide nitrogen, an unrelated type). Verified against
        // a live RDKit oracle.
        let m = mol("o1cccc1"); // furan
        let types = assign_mmff94_numeric_types(&m).expect("furan must type successfully");
        for i in 0..m.atom_count() {
            if m.atom(AtomIdx(i as u32)).element == Element::O {
                assert_eq!(types[i], 59, "furan O should be type 59 (OFUR)");
            }
        }
    }

    #[test]
    fn o2cm_terminal_oxygen_fixture_matrix_matches_rdkit() {
        // Issue #227 Priority 1A-3: closes the 37-atom
        // `terminal_oxygen_o2cm_umbrella_gap` residual left after PR #239
        // (`classify_terminal_o`, ported from `AtomTyper.cpp` `case 8:`,
        // lines 1554-1737 at the pinned commit). All expected O types below
        // are copied verbatim from a live RDKit oracle
        // (`AllChem.MMFFGetMoleculeProperties`), covering every disjunct of
        // the real O2CM union plus the OM(35)/O=C(7) alternatives it must
        // not swallow.
        let fixtures: &[(&str, &str, &[u8])] = &[
            // Anionic carboxylate: both O's -> O2CM (its C has 2 terminal O).
            ("acetate", "CC(=O)[O-]", &[32, 32]),
            // Neutral carboxylic acid: =O generic carbonyl, -OH stays OR (6)
            // -- the carboxylate condition must NOT fire (only 1 terminal O
            // on that carbon; -OH has an implicit H so isn't terminal at
            // all).
            ("acetic_acid", "CC(=O)O", &[7, 6]),
            // Ester: carbonyl =O stays 7, bridging -O- (non-terminal) stays 6.
            ("ester", "CC(=O)OC", &[7, 6]),
            // Nitro: central N has 2 terminal O -> both O2CM, regardless of
            // which O is drawn with the formal negative charge.
            ("nitrobenzene", "c1ccc(cc1)[N+](=O)[O-]", &[32, 32]),
            // Genuine (charge-separated) N-oxide -> O2CM, not OM -- this is
            // the fixture that exercises the aromatic-ring-N charge-based
            // discriminator (see `classify_terminal_o`'s N branch doc).
            ("pyridine_n_oxide", "c1ccc[n+]([O-])c1", &[32]),
            // Sulfoxide: exactly 1 terminal O on S -> generic O=S (7), not
            // O2CM -- must not be swept in by the sulfone/sulfonate rule.
            ("dmso", "CS(=O)C", &[7]),
            // Sulfone: 2 terminal O's on S -> both O2CM.
            ("dimethyl_sulfone", "CS(=O)(=O)C", &[32, 32]),
            ("sulfonamide", "CS(=O)(=O)N", &[32, 32]),
            ("sulfonate", "CS(=O)(=O)[O-]", &[32, 32, 32]),
            // Sulfate ester: the 2 terminal =O -> O2CM; the 2 bridging
            // -O-C ester oxygens (non-terminal) stay OR (6).
            ("sulfate", "COS(=O)(=O)OC", &[6, 32, 32, 6]),
            // Phosphate ester: terminal =O -> O2CM unconditionally (RDKit's
            // own quirk, ported faithfully); bridging -O-C esters stay 6.
            ("phosphate", "COP(=O)(OC)OC", &[6, 32, 6, 6]),
            ("phosphonate", "CP(=O)(O)O", &[32, 6, 6]),
            // Perchlorate anion: all 4 O's on Cl -> O2CM unconditionally.
            ("perchlorate", "[O-]Cl(=O)(=O)=O", &[32, 32, 32, 32]),
            // Non-terminal / unrelated O's: unaffected negative controls.
            ("ether", "COC", &[6]),
            ("alcohol", "CO", &[6]),
            ("ketone", "CC(=O)C", &[7]),
            ("amide_carbonyl", "CC(=O)NC", &[7]),
        ];

        for (name, smiles, expected) in fixtures {
            let m = mol(smiles);
            let types = assign_mmff94_numeric_types(&m)
                .unwrap_or_else(|e| panic!("{name} ({smiles}) failed to type: {e}"));
            let actual: Vec<u8> = (0..m.atom_count())
                .filter(|&i| m.atom(AtomIdx(i as u32)).element == Element::O)
                .map(|i| types[i])
                .collect();
            assert_eq!(
                &actual, expected,
                "{name} ({smiles}): expected O types {expected:?}, got {actual:?}"
            );
        }
    }

    // ── Phase 0.1 (issue #227, Priority 1A pre-merge fix): Kekulization
    // failure must fail closed, never silently reuse `mol` as a stand-in
    // MMFF view ──────────────────────────────────────────────────────────

    /// Odd-membered, all-carbon, neutral aromatic ring: every atom is a
    /// `atom_must_be_matched` carbon (no lone pair, no exocyclic multiple
    /// bond), so a perfect matching over 5 atoms is impossible by simple
    /// parity -- guaranteed `Err` from `chematic_core::kekulize` regardless
    /// of which of its four matching passes runs, not a flaky/order-
    /// dependent failure.
    fn unkekulizable_five_ring() -> Molecule {
        let mut b = chematic_core::MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| b.add_atom(chematic_core::Atom::aromatic(Element::C)))
            .collect();
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// A different failure mechanism from the parity case above, and one
    /// that specifically exercises a real SSSR-detected ring (unlike a bare
    /// 2-atom edge, which `find_sssr` correctly reports as ring-free and
    /// which would make `compute_mmff94_aromatic_view` take its
    /// no-rings-to-reperceive early-return instead of ever reaching
    /// Kekulization -- not the code path this test targets): a 3-membered
    /// all-aromatic-notation ring, C-O-O, where the lone must-match atom
    /// (C: no lone pair, no exocyclic multiple bond) has only lone-pair-
    /// donating heteroatom neighbors (both O), so it has zero eligible
    /// double-bond partners in the matching graph -- unmatchable regardless
    /// of which of the four matching passes runs, not a parity artifact.
    /// Models a malformed valence/aromatic-notation combination (an
    /// aromatic ring notation whose bonding pattern cannot support any
    /// consistent Kekule structure).
    fn unkekulizable_isolated_must_match_carbon() -> Molecule {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(chematic_core::Atom::aromatic(Element::C));
        let o1 = b.add_atom(chematic_core::Atom::aromatic(Element::O));
        let o2 = b.add_atom(chematic_core::Atom::aromatic(Element::O));
        b.add_bond(c, o1, BondOrder::Aromatic).unwrap();
        b.add_bond(o1, o2, BondOrder::Aromatic).unwrap();
        b.add_bond(o2, c, BondOrder::Aromatic).unwrap();
        b.build()
    }

    #[test]
    fn kekulization_impossible_molecule_fails_closed_not_silently_reused() {
        for m in [
            unkekulizable_five_ring(),
            unkekulizable_isolated_must_match_carbon(),
        ] {
            // Sanity: confirm the fixture actually exercises the failure
            // path this test targets, not some unrelated typing error.
            assert!(
                chematic_core::kekulize(&m).is_err(),
                "test fixture must be genuinely unkekulizable"
            );

            let result = assign_mmff94_numeric_types(&m);
            let err = result.expect_err(
                "a molecule chematic cannot Kekulize must fail MMFF94 typing, \
                 never silently proceed using the un-re-perceived molecule",
            );

            // The failure must name the Kekulization stage and carry the
            // underlying cause, not be an opaque/generic message -- so a
            // caller (and `chematic-3d`'s typed `ForceFieldBridgeError`
            // conversion) can tell this apart from a different typing
            // failure by more than string content alone if it ever needs
            // to.
            assert!(
                err.0.contains("Kekulization"),
                "error must name the Kekulization stage, got: {}",
                err.0
            );
            assert!(
                err.0.contains("cannot be assigned a double bond"),
                "error must carry the underlying KekuleError cause, got: {}",
                err.0
            );

            // Determinism: repeated calls on the same input produce the
            // identical typed error, not a flaky/order-dependent result.
            let err2 = assign_mmff94_numeric_types(&m).expect_err("must fail deterministically");
            assert_eq!(err, err2, "Kekulization failure must be deterministic");

            // Never a wrong type: there is no `Ok` branch here at all --
            // the only way this test could regress into "returns a wrong
            // type" is by silently succeeding, which the above already
            // rules out.
        }
    }

    // ── Phase 0.2 (issue #227, Priority 1A pre-merge fix): a ring's bonds
    // must only be promoted to `BondOrder::Aromatic` when that ring itself
    // passed the Huckel check -- never reconstructed from whether every
    // atom in the ring happens to be aromatic via *other* accepted rings
    // ──────────────────────────────────────────────────────────────────

    #[test]
    fn pyracylene_peri_fused_artifact_ring_bonds_are_not_wrongly_promoted() {
        // Empirically confirmed hazard case (not theoretical): pyracylene's
        // SSSR contains a genuine 4-membered "envelope" ring, atoms
        // {6,7,8,13}, that is neither of the molecule's two real
        // (RDKit-agreeing) accepted aromatic rings -- verified identical on
        // both engines' own SSSR, not assumed. Every one of its 4 atoms is
        // independently aromatic via the two real accepted rings it
        // straddles, so the pre-fix reconstruction
        // (`ring.iter().all(|a| is_arom[a])`) would have wrongly promoted
        // this artifact ring's own bonds to `BondOrder::Aromatic` too.
        // RDKit ground truth (`Chem.SetAromaticity(mol,
        // Chem.AROMATICITY_MMFF94)` on `c1cc2cccc3c2c2c1cccc23`, this PR's
        // audit): all 14 atoms aromatic, but this specific 4-ring's own
        // bonds are NOT (its own Huckel check never independently passes:
        // an all-carbon 4-membered ring has no heteroatom lone-pair bonus
        // available and can supply at most 2 ring-double-bond pi electrons
        // either way, never satisfying `pi_e > 2 && (pi_e-2)%4==0`).
        let m = mol("c1cc2cccc3c2c2c1cccc23");
        let rings = chematic_perception::find_sssr(&m).rings().to_vec();
        let artifact_ring = rings
            .iter()
            .find(|r| {
                let mut set: Vec<u32> = r.iter().map(|a| a.0).collect();
                set.sort_unstable();
                set == [6, 7, 8, 13]
            })
            .expect("pyracylene's SSSR must contain the {6,7,8,13} artifact ring");
        assert_eq!(
            artifact_ring.len(),
            4,
            "sanity: the artifact ring must be the 4-membered one"
        );

        let view = compute_mmff94_aromatic_view(&m, &rings)
            .expect("pyracylene kekulizes and re-perceives successfully");

        let real_rings: Vec<&Vec<AtomIdx>> = rings
            .iter()
            .filter(|r| {
                let mut set: Vec<u32> = r.iter().map(|a| a.0).collect();
                set.sort_unstable();
                set != [6, 7, 8, 13]
            })
            .collect();
        assert_eq!(real_rings.len(), 3, "pyracylene has 4 SSSR rings total");

        // Two of the artifact ring's four edges (7-8, 8-13) are *also*
        // edges of one of the two real accepted 6-rings, so they are
        // legitimately aromatic for that unrelated reason -- not evidence
        // either way about the artifact ring itself. Only the two edges
        // exclusive to the artifact ring (6-7, 6-13) isolate the bug: they
        // must NOT be aromatic, since no real accepted ring ever claims
        // them.
        let is_edge_of_a_real_ring = |a: AtomIdx, b: AtomIdx| -> bool {
            real_rings.iter().any(|r| {
                let len = r.len();
                (0..len).any(|j| {
                    let (x, y) = (r[j], r[(j + 1) % len]);
                    (x == a && y == b) || (x == b && y == a)
                })
            })
        };

        let len = artifact_ring.len();
        let mut checked_exclusive_edge = false;
        for (j, &a) in artifact_ring.iter().enumerate() {
            let b = artifact_ring[(j + 1) % len];
            if is_edge_of_a_real_ring(a, b) {
                continue; // shared edge -- aromatic for an unrelated, legitimate reason
            }
            checked_exclusive_edge = true;
            let (_, bond) = view
                .bond_between(a, b)
                .unwrap_or_else(|| panic!("bond {}-{} must exist", a.0, b.0));
            assert_ne!(
                bond.order,
                BondOrder::Aromatic,
                "artifact-ring-exclusive bond {}-{} must NOT be promoted to aromatic -- this \
                 ring never independently passed the Huckel check, even though all 4 of its \
                 atoms are aromatic via the two other, real accepted rings",
                a.0,
                b.0
            );
        }
        assert!(
            checked_exclusive_edge,
            "sanity: the artifact ring must have at least one edge not shared with a real ring"
        );

        // Positive control in the same molecule: the real rings must still
        // be correctly promoted, proving this isn't just "nothing gets
        // promoted" degenerate behavior.
        let mut any_promoted = false;
        for ring in &real_rings {
            let len = ring.len();
            for (j, &a) in ring.iter().enumerate() {
                let b = ring[(j + 1) % len];
                if let Some((_, bond)) = view.bond_between(a, b)
                    && bond.order == BondOrder::Aromatic
                {
                    any_promoted = true;
                }
            }
        }
        assert!(
            any_promoted,
            "at least one of pyracylene's real rings must still be promoted aromatic"
        );
    }

    #[test]
    fn mmff94_aromaticity_fixture_matrix_matches_rdkit_atom_and_bond_flags() {
        // Issue #227 Phase 0.2: 12-fixture matrix, atom-level AND bond-level
        // MMFF aromaticity parity against a live RDKit oracle
        // (`Chem.SetAromaticity(mol, Chem.AROMATICITY_MMFF94)` on a
        // freshly-Kekulized copy -- RDKit's real `setMMFFAromaticity`,
        // independent of `MMFFGetMoleculeProperties`'s own internal usage of
        // it; verified on caffeine against the already-established pattern
        // (see `caffeine_matches_rdkit_mmff_aromaticity_per_atom` above) before
        // trusting it for the other 11. Values copied verbatim from
        // `validation/results/mmff94_aromaticity_bond_parity_227_oracle.json`
        // (generator: `scripts/mmff94_aromaticity_bond_parity_227.py`), not
        // derived from this fix's own output. Fixtures chosen to stress the
        // exact fused-ring bond-promotion hazard Phase 0.2 fixes: tetralin
        // (fused aromatic+non-aromatic), a spiro system, and a bridged
        // system are included specifically because a naive
        // reconstruct-from-atom-flags approach could wrongly promote their
        // non-aromatic ring's bonds. Bond-type comparison is not separately
        // re-derived here: `bond_type_for` is a pure function of exactly
        // (atom type, bond order), both already verified per-fixture below.
        struct Fixture {
            name: &'static str,
            smiles: &'static str,
            atom_aromatic: &'static [bool],
            bond_aromatic: &'static [(u32, u32, bool)],
            atom_types: &'static [u8],
        }

        // --- benzene ---
        const BENZENE_ATOM_AROM: [bool; 6] = [true, true, true, true, true, true];
        const BENZENE_BOND_AROM: [(u32, u32, bool); 6] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (4, 5, true),
        ];
        const BENZENE_ATOM_TYPES: [u8; 6] = [37, 37, 37, 37, 37, 37];

        // --- naphthalene ---
        const NAPHTHALENE_ATOM_AROM: [bool; 10] =
            [true, true, true, true, true, true, true, true, true, true];
        const NAPHTHALENE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 9, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 8, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
        ];
        const NAPHTHALENE_ATOM_TYPES: [u8; 10] = [37, 37, 37, 37, 37, 37, 37, 37, 37, 37];

        // --- anthracene ---
        const ANTHRACENE_ATOM_AROM: [bool; 14] = [
            true, true, true, true, true, true, true, true, true, true, true, true, true, true,
        ];
        const ANTHRACENE_BOND_AROM: [(u32, u32, bool); 16] = [
            (0, 1, true),
            (0, 13, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 12, true),
            (4, 5, true),
            (5, 6, true),
            (5, 10, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
            (9, 10, true),
            (10, 11, true),
            (11, 12, true),
            (12, 13, true),
        ];
        const ANTHRACENE_ATOM_TYPES: [u8; 14] =
            [37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37];

        // --- indole ---
        const INDOLE_ATOM_AROM: [bool; 9] = [true, true, true, true, true, true, true, true, true];
        const INDOLE_BOND_AROM: [(u32, u32, bool); 10] = [
            (0, 1, true),
            (0, 8, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 7, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
        ];
        const INDOLE_ATOM_TYPES: [u8; 9] = [37, 37, 37, 63, 39, 63, 64, 64, 37];

        // --- quinoline ---
        const QUINOLINE_ATOM_AROM: [bool; 10] =
            [true, true, true, true, true, true, true, true, true, true];
        const QUINOLINE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 9, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 8, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
        ];
        const QUINOLINE_ATOM_TYPES: [u8; 10] = [37, 37, 37, 37, 38, 37, 37, 37, 37, 37];

        // --- purine ---
        const PURINE_ATOM_AROM: [bool; 9] = [true, true, true, true, true, true, true, true, true];
        const PURINE_BOND_AROM: [(u32, u32, bool); 10] = [
            (0, 1, true),
            (0, 8, true),
            (1, 2, true),
            (2, 3, true),
            (2, 6, true),
            (3, 4, true),
            (4, 5, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
        ];
        const PURINE_ATOM_TYPES: [u8; 9] = [37, 38, 63, 39, 63, 66, 64, 37, 38];

        // --- caffeine ---
        const CAFFEINE_ATOM_AROM: [bool; 14] = [
            false, true, true, true, true, true, false, false, false, false, false, false, false,
            false,
        ];
        const CAFFEINE_BOND_AROM: [(u32, u32, bool); 15] = [
            (0, 1, false),
            (1, 2, true),
            (1, 5, true),
            (2, 3, true),
            (3, 4, true),
            (4, 5, true),
            (4, 12, false),
            (5, 6, false),
            (6, 7, false),
            (6, 8, false),
            (8, 9, false),
            (8, 10, false),
            (10, 11, false),
            (10, 12, false),
            (12, 13, false),
        ];
        const CAFFEINE_ATOM_TYPES: [u8; 14] = [1, 39, 63, 66, 64, 63, 3, 7, 10, 1, 3, 7, 10, 1];

        // --- azulene ---
        const AZULENE_ATOM_AROM: [bool; 10] = [
            false, false, false, false, false, false, false, false, false, false,
        ];
        const AZULENE_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, false),
            (0, 9, false),
            (1, 2, false),
            (2, 3, false),
            (3, 4, false),
            (3, 9, false),
            (4, 5, false),
            (5, 6, false),
            (6, 7, false),
            (7, 8, false),
            (8, 9, false),
        ];
        const AZULENE_ATOM_TYPES: [u8; 10] = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

        // --- tetralin_fused_nonaromatic ---
        const TETRALIN_FUSED_NONAROMATIC_ATOM_AROM: [bool; 10] = [
            true, true, true, true, true, true, false, false, false, false,
        ];
        const TETRALIN_FUSED_NONAROMATIC_BOND_AROM: [(u32, u32, bool); 11] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 9, false),
            (4, 5, true),
            (4, 6, false),
            (6, 7, false),
            (7, 8, false),
            (8, 9, false),
        ];
        const TETRALIN_FUSED_NONAROMATIC_ATOM_TYPES: [u8; 10] =
            [37, 37, 37, 37, 37, 37, 1, 1, 1, 1];

        // --- carbostyril_exocyclic_carbonyl ---
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_AROM: [bool; 11] = [
            false, false, false, false, true, true, true, true, true, true, false,
        ];
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_BOND_AROM: [(u32, u32, bool); 12] = [
            (0, 1, false),
            (1, 2, false),
            (1, 10, false),
            (2, 3, false),
            (3, 4, false),
            (4, 5, true),
            (4, 9, true),
            (5, 6, true),
            (6, 7, true),
            (7, 8, true),
            (8, 9, true),
            (9, 10, false),
        ];
        const CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_TYPES: [u8; 11] =
            [7, 3, 2, 2, 37, 37, 37, 37, 37, 37, 10];

        // --- spiro_indane_cyclohexane ---
        const SPIRO_INDANE_CYCLOHEXANE_ATOM_AROM: [bool; 14] = [
            true, true, true, true, true, true, false, false, false, false, false, false, false,
            false,
        ];
        const SPIRO_INDANE_CYCLOHEXANE_BOND_AROM: [(u32, u32, bool); 16] = [
            (0, 1, true),
            (0, 5, true),
            (1, 2, true),
            (2, 3, true),
            (3, 4, true),
            (3, 13, false),
            (4, 5, true),
            (4, 6, false),
            (6, 7, false),
            (7, 8, false),
            (7, 12, false),
            (7, 13, false),
            (8, 9, false),
            (9, 10, false),
            (10, 11, false),
            (11, 12, false),
        ];
        const SPIRO_INDANE_CYCLOHEXANE_ATOM_TYPES: [u8; 14] =
            [37, 37, 37, 37, 37, 37, 1, 1, 1, 1, 1, 1, 1, 1];

        // --- norbornane_bridged ---
        const NORBORNANE_BRIDGED_ATOM_AROM: [bool; 7] =
            [false, false, false, false, false, false, false];
        const NORBORNANE_BRIDGED_BOND_AROM: [(u32, u32, bool); 8] = [
            (0, 1, false),
            (0, 5, false),
            (1, 2, false),
            (2, 3, false),
            (2, 6, false),
            (3, 4, false),
            (4, 5, false),
            (5, 6, false),
        ];
        const NORBORNANE_BRIDGED_ATOM_TYPES: [u8; 7] = [1, 1, 1, 1, 1, 1, 1];

        let fixtures: Vec<Fixture> = vec![
            Fixture {
                name: "benzene",
                smiles: "c1ccccc1",
                atom_aromatic: &BENZENE_ATOM_AROM,
                bond_aromatic: &BENZENE_BOND_AROM,
                atom_types: &BENZENE_ATOM_TYPES,
            },
            Fixture {
                name: "naphthalene",
                smiles: "c1ccc2ccccc2c1",
                atom_aromatic: &NAPHTHALENE_ATOM_AROM,
                bond_aromatic: &NAPHTHALENE_BOND_AROM,
                atom_types: &NAPHTHALENE_ATOM_TYPES,
            },
            Fixture {
                name: "anthracene",
                smiles: "c1ccc2cc3ccccc3cc2c1",
                atom_aromatic: &ANTHRACENE_ATOM_AROM,
                bond_aromatic: &ANTHRACENE_BOND_AROM,
                atom_types: &ANTHRACENE_ATOM_TYPES,
            },
            Fixture {
                name: "indole",
                smiles: "c1ccc2[nH]ccc2c1",
                atom_aromatic: &INDOLE_ATOM_AROM,
                bond_aromatic: &INDOLE_BOND_AROM,
                atom_types: &INDOLE_ATOM_TYPES,
            },
            Fixture {
                name: "quinoline",
                smiles: "c1ccc2ncccc2c1",
                atom_aromatic: &QUINOLINE_ATOM_AROM,
                bond_aromatic: &QUINOLINE_BOND_AROM,
                atom_types: &QUINOLINE_ATOM_TYPES,
            },
            Fixture {
                name: "purine",
                smiles: "c1nc2[nH]cnc2cn1",
                atom_aromatic: &PURINE_ATOM_AROM,
                bond_aromatic: &PURINE_BOND_AROM,
                atom_types: &PURINE_ATOM_TYPES,
            },
            Fixture {
                name: "caffeine",
                smiles: "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                atom_aromatic: &CAFFEINE_ATOM_AROM,
                bond_aromatic: &CAFFEINE_BOND_AROM,
                atom_types: &CAFFEINE_ATOM_TYPES,
            },
            Fixture {
                name: "azulene",
                smiles: "c1ccc2cccccc12",
                atom_aromatic: &AZULENE_ATOM_AROM,
                bond_aromatic: &AZULENE_BOND_AROM,
                atom_types: &AZULENE_ATOM_TYPES,
            },
            Fixture {
                name: "tetralin_fused_nonaromatic",
                smiles: "c1ccc2c(c1)CCCC2",
                atom_aromatic: &TETRALIN_FUSED_NONAROMATIC_ATOM_AROM,
                bond_aromatic: &TETRALIN_FUSED_NONAROMATIC_BOND_AROM,
                atom_types: &TETRALIN_FUSED_NONAROMATIC_ATOM_TYPES,
            },
            Fixture {
                name: "carbostyril_exocyclic_carbonyl",
                smiles: "O=c1ccc2ccccc2[nH]1",
                atom_aromatic: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_AROM,
                bond_aromatic: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_BOND_AROM,
                atom_types: &CARBOSTYRIL_EXOCYCLIC_CARBONYL_ATOM_TYPES,
            },
            Fixture {
                name: "spiro_indane_cyclohexane",
                smiles: "c1ccc2c(c1)CC3(CCCCC3)C2",
                atom_aromatic: &SPIRO_INDANE_CYCLOHEXANE_ATOM_AROM,
                bond_aromatic: &SPIRO_INDANE_CYCLOHEXANE_BOND_AROM,
                atom_types: &SPIRO_INDANE_CYCLOHEXANE_ATOM_TYPES,
            },
            Fixture {
                name: "norbornane_bridged",
                smiles: "C1CC2CCC1C2",
                atom_aromatic: &NORBORNANE_BRIDGED_ATOM_AROM,
                bond_aromatic: &NORBORNANE_BRIDGED_BOND_AROM,
                atom_types: &NORBORNANE_BRIDGED_ATOM_TYPES,
            },
        ];

        for fx in &fixtures {
            let m = mol(fx.smiles);
            assert_eq!(
                m.atom_count(),
                fx.atom_aromatic.len(),
                "{}: atom count / element-sequence alignment sanity",
                fx.name
            );
            let rings = chematic_perception::find_sssr(&m).rings().to_vec();
            let view = compute_mmff94_aromatic_view(&m, &rings)
                .unwrap_or_else(|e| panic!("{}: re-perception failed: {e}", fx.name));

            for i in 0..m.atom_count() {
                let idx = AtomIdx(i as u32);
                assert_eq!(
                    view.atom(idx).aromatic,
                    fx.atom_aromatic[i],
                    "{}: atom {i} aromatic-flag mismatch vs RDKit oracle",
                    fx.name
                );
            }

            for &(a, b, expected_arom) in fx.bond_aromatic {
                let (_, bond) = view
                    .bond_between(AtomIdx(a), AtomIdx(b))
                    .unwrap_or_else(|| panic!("{}: bond {a}-{b} must exist", fx.name));
                let is_arom = bond.order == BondOrder::Aromatic;
                assert_eq!(
                    is_arom, expected_arom,
                    "{}: bond {a}-{b} aromatic-flag mismatch vs RDKit oracle (got {is_arom}, want {expected_arom})",
                    fx.name
                );
            }

            if !fx.atom_types.is_empty() {
                let types = assign_mmff94_numeric_types(&m)
                    .unwrap_or_else(|e| panic!("{}: numeric typing failed: {e}", fx.name));
                assert_eq!(
                    types, fx.atom_types,
                    "{}: numeric atom types mismatch vs RDKit oracle",
                    fx.name
                );
            }
        }
    }

    #[test]
    fn kekulization_failure_never_reaches_the_uff_or_mmff94_minimizer_bridge() {
        // Phase 0.1's "must not proceed to parameter lookup" requirement,
        // verified at this crate's own public boundary: every
        // `assign_mmff94_numeric_types` caller in this crate (e.g.
        // `crates/chematic-ff/src/mmff94_minimizer.rs`) uses `?` on this
        // call before ever touching a bond/angle/torsion/oop parameter
        // table, so a `NumericTypeError` here structurally cannot be
        // followed by a parameter-lookup attempt on the same call. This
        // test pins that `?`-propagation contract at the numeric-typing
        // boundary itself, independent of any specific downstream caller.
        let m = unkekulizable_five_ring();
        let types_result = assign_mmff94_numeric_types(&m);
        assert!(types_result.is_err());
        // If this were `Ok`, a caller doing `let types = ...?;` would
        // proceed to index `types` for parameter lookup; confirming `Err`
        // here is exactly what makes that `?` short-circuit for every
        // caller, without needing to duplicate every caller's own test.
    }

    // ── Cumulated-double-bond CSP fix (issue #337 sub-bug 2) ─────────────────
    // Root cause: `assign_c_type`'s sp-carbon check only fired on
    // `triple_bonds > 0`, so a carbon reached via *two* double bonds (a
    // cumulated diene / "allenic" carbon, e.g. the central C of an aryl
    // isothiocyanate's N=C=S) fell through to the `double_bonds > 0`
    // "double-bonded to N/O/P/S" branch and got the generic carbonyl-family
    // type 3 instead of RDKit's real CSP type 4. RDKit's actual rule
    // (`AtomTyper.cpp`, pinned commit
    // `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, lines ~954-960) is simply
    // `getTotalDegree() == 2`, unconditional on which elements the two bonds
    // go to -- see `assign_c_type`'s doc comment for the full citation.
    // These tests pin the 2/8 issue #337 molecules this fix resolves (the
    // other 6/8, the pyridinium-conjugated-exocyclic-amine sub-bug, are an
    // RDKit Kekulization/aromaticity-perception artifact, not an
    // atom-typing rule -- left as an honestly-disclosed residual, see
    // `scripts/mmff94_provenance/PROVENANCE.md`), plus isolated synthetic
    // fixtures pinning the exact discriminating condition independent of
    // the corpus molecules' other complexity.

    #[test]
    fn propyne_alkyne_carbons_still_type_csp_after_degree_based_fix() {
        // No-regression pin: a real triple bond always yields
        // `total_degree == 2` for carbon (the triple bond alone consumes 3
        // of its 4 valence units), so this must keep matching RDKit exactly
        // as it did before this fix broadened the condition. Expected
        // values from a fresh live RDKit oracle query (`rdkit==2026.03.4`,
        // `MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")` on the
        // implicit-H `Chem.MolFromSmiles` result -- MMFF atom
        // typing/charges are purely topological, no embedding needed, same
        // precedent as `scripts/mmff94_bci_charges_oracle_227.py`).
        let m = mol("CC#C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 4, 4]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.2, -0.2, 0.0];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "propyne atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn allene_central_carbon_types_csp_not_generic_vinylic() {
        // Broader-than-corpus pin: plain carbon allene (no heteroatoms at
        // all) is not in the 264-molecule corpus, but is the clearest
        // demonstration that the real discriminating condition is
        // `total_degree == 2`, not "cumulated bond to a heteroatom" --
        // before this fix, the central carbon was mistyped 2 (generic
        // vinylic C=C), same failure mode as the isothiocyanate carbons.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("C=C=C");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![2, 4, 2]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.065, -0.13, 0.065];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "allene atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn methyl_isothiocyanate_minimal_ncs_fixture_matches_rdkit_oracle() {
        // Minimal isolated reproduction of the corpus mechanism below,
        // independent of the aromatic ring and amide/piperazine complexity
        // both `chembl_tier_b_0071`/`_0082` also carry. Expected values
        // from a fresh live RDKit oracle query (`rdkit==2026.03.4`).
        let m = mol("CN=C=S");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![1, 9, 4, 16]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.246, -0.546, 0.575, -0.275];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "methyl isothiocyanate atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn phenyl_isothiocyanate_aryl_ncs_fixture_matches_rdkit_oracle() {
        // Same motif, aryl instead of methyl (closer to the corpus
        // molecules' own aryl isothiocyanate substructure). Also confirms
        // the fix doesn't disturb the adjacent aromatic ring's own typing.
        // Expected values from a fresh live RDKit oracle query
        // (`rdkit==2026.03.4`).
        let m = mol("c1ccccc1N=C=S");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        assert_eq!(types, vec![37, 37, 37, 37, 37, 37, 9, 4, 16]);
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "phenyl isothiocyanate atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0071_aryl_isothiocyanate_matches_rdkit_oracle_after_csp_fix() {
        // Corpus regression pin. Expected values cross-checked against the
        // already-committed live RDKit oracle dumps
        // (`validation/results/mmff94_rdkit_type_oracle.jsonl` and
        // `mmff94_bci_charges_227_rdkit_oracle.jsonl`, `rdkit==2026.03.4`)
        // via the corpus-wide per-atom join reported in
        // `scripts/mmff94_provenance/PROVENANCE.md`'s issue #337
        // follow-up, not derived from this fix's own output.
        let m = mol("COc1cc2nc(N3CCN(C(=O)c4ccc(N=C=S)c(I)c4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 10, 3, 7, 37, 37, 37, 37, 9, 4, 16, 37, 14, 37, 1,
            1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0071 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3001, -0.6602,
            0.5438, -0.57, 0.0862, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275, 0.081, -0.081, 0.0,
            0.3001, 0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0071 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }

    #[test]
    fn chembl_tier_b_0082_aryl_isothiocyanate_matches_rdkit_oracle_after_csp_fix() {
        // Corpus regression pin, same mechanism as `_0071` above via an
        // enone rather than a plain benzamide linker. Expected values
        // cross-checked against the already-committed live RDKit oracle
        // dumps via the same corpus-wide per-atom join, not derived from
        // this fix's own output.
        let m = mol("COc1cc2nc(N3CCN(C(=O)/C=C/c4ccc(N=C=S)cc4)CC3)nc(N)c2cc1OC");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let expected_types = [
            1, 6, 37, 37, 37, 38, 37, 40, 1, 1, 10, 3, 7, 2, 2, 37, 37, 37, 37, 9, 4, 16, 37, 37,
            1, 1, 38, 37, 40, 37, 37, 37, 6, 1,
        ];
        assert_eq!(
            types, expected_types,
            "sanity: chembl_tier_b_0082 atom types"
        );
        let q = mmff94_charges_numeric(&m).unwrap();
        let expected = [
            0.28, -0.3625, 0.0825, 0.0, 0.31, -0.62, 0.72, -0.8382, 0.3691, 0.3001, -0.6602,
            0.6156, -0.57, 0.0144, -0.0284, 0.0284, 0.0, 0.0, 0.179, -0.479, 0.575, -0.275, 0.0,
            0.0, 0.3001, 0.3691, -0.62, 0.41, -0.1, 0.0, 0.0, 0.0825, -0.3625, 0.28,
        ];
        for (i, exp) in expected.iter().enumerate() {
            assert!(
                (q[i] - exp).abs() < 1e-6,
                "chembl_tier_b_0082 atom {i}: expected charge {exp}, got {}",
                q[i]
            );
        }
    }
}
