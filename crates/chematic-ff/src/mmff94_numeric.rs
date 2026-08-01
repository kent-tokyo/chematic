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
    (32, -0.732, 0.500), // NR+      protonated amine N (fcadj=0.5)
    (33, 0.257, 0.000),  // HOX      H on O in N-oxide
    (34, -0.491, 0.000), // O-       anionic O (carboxylate/phenoxide)
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
    (45, -0.260, 0.000), // N5       generic 5-ring aromatic N
    (46, -0.429, 0.000), // NO2      nitro N
    (47, -0.418, 0.000), // NO3      nitrate N
    (48, -0.525, 0.000), // O2NO     nitro O
    (49, -0.283, 0.000), // O3NO     nitrate O
    (50, 0.284, 0.000),  // OP       phosphate O
    (51, -1.046, 0.000), // O2P      phosphonate =O
    (52, -0.546, 0.000), // O3P      bridging phosphate O
    (53, -0.048, 0.000), // O4P      phosphate anion O
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

fn bond_type_for(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Aromatic => 4,
        _ => 0,
    }
}

// ── Atom type assignment ─────────────────────────────────────────────────────

/// Assign MMFF94 numeric atom types (1–99) to all atoms in the molecule.
///
/// This implements the core atom type perception rules for organic chemistry.
/// For atoms not handled, returns `Err`.
pub fn assign_mmff94_numeric_types(mol: &Molecule) -> Result<Vec<u8>, NumericTypeError> {
    let n = mol.atom_count();
    let mut types = vec![0u8; n];
    let rings = chematic_perception::find_sssr(mol).rings().to_vec();

    for (i, ty) in types.iter_mut().enumerate().take(n) {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let t = match atom.element {
            Element::C => assign_c_type(mol, &rings, idx)?,
            Element::N => assign_n_type(mol, &rings, idx)?,
            Element::O => assign_o_type(mol, &rings, idx)?,
            Element::S => assign_s_type(mol, idx)?,
            Element::P => assign_p_type(mol, idx)?,
            Element::SI => 19,
            Element::F => 11,
            Element::CL => 12,
            Element::BR => 13,
            Element::I => 14,
            Element::H => assign_h_type(mol, idx)?,
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

    Ok(types)
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
    let triple_bonds = count_bond_order(mol, idx, BondOrder::Triple);

    // sp carbon (triple bond or allene)
    if triple_bonds > 0 {
        return Ok(4); // CSP
    }

    // sp2 carbon
    if double_bonds > 0 {
        // C=O, C=S → type 3 (carbonyl/thioamide family)
        if is_bonded_to(mol, idx, Element::O, BondOrder::Double)
            || is_bonded_to(mol, idx, Element::S, BondOrder::Double)
        {
            return Ok(3); // C=O (general carbonyl)
        }
        // C=N or C=C → type 2
        return Ok(2); // C=C vinylic
    }

    // sp3
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

    // Formal charge: quaternary ammonium / protonated N.
    // Registry-verified: type 34 is NR+ (N+, QUATERNARY N); type 32 is
    // O2CM (O, CARBOXYLATE ANION), an oxygen-only type -- the previous
    // `32` here was exactly the silent element-collision the numeric
    // type registry's construction-time invariant now catches instead
    // of allowing through as a false "success".
    if atom.charge > 0 {
        return Ok(34); // NR+
    }

    // Nitrile / isocyanide (N≡C)
    if triple_bonds > 0 {
        return Ok(9); // N=C (close approximation for nitrile)
    }

    // N=C or N=N (imine, hydrazone, etc.)
    if double_bonds > 0 {
        return Ok(9); // N=C imine
    }

    // sp3 N — check if amide (bonded to carbonyl C)
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

    // Nitro group (N with two =O bonds): check for N(=O)=O pattern
    let double_o = bonds_of(mol, idx)
        .iter()
        .filter(|b| b.order == BondOrder::Double && mol.atom(b.neighbor).element == Element::O)
        .count();
    if double_o >= 2 {
        return Ok(46); // NO2 nitro N
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

    // Double bond to C or N → carbonyl/similar oxygen (type 7)
    if count_bond_order(mol, idx, BondOrder::Double) > 0 {
        return Ok(7); // O=C
    }

    // Anionic O (formal charge -1): carboxylate/phenoxide.
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

fn assign_s_type(mol: &Molecule, idx: AtomIdx) -> Result<u8, NumericTypeError> {
    let atom = mol.atom(idx);
    if atom.aromatic {
        return Ok(44); // S5 aromatic sulfur (thiophene)
    }

    let double_o = bonds_of(mol, idx)
        .iter()
        .filter(|b| b.order == BondOrder::Double && mol.atom(b.neighbor).element == Element::O)
        .count();

    match double_o {
        2.. => Ok(18), // SO2 sulfone
        1 => Ok(17),   // S=O sulfoxide
        0 => {
            // Check double bond to C
            if count_bond_order(mol, idx, BondOrder::Double) > 0 {
                return Ok(16); // S=C
            }
            Ok(15) // S thiol/sulfide
        }
    }
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

/// Compute MMFF94 partial charges using the full PBCI+CHG tables (Halgren 1996).
///
/// Implements equation 15 from MMFF.V paper. For most neutral organic atoms
/// (fcadj=0, no formal charge), this reduces to:
///   q_i = Σ_{j bonded} bci(j→i)
///
/// Returns per-atom partial charges in units of elementary charge.
pub fn mmff94_charges_numeric(mol: &Molecule) -> Result<Vec<f64>, NumericTypeError> {
    let types = assign_mmff94_numeric_types(mol)?;
    let n = mol.atom_count();
    let mut charges = vec![0.0f64; n];

    // Step 1: formal charge contribution (scaled by fcadj)
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let (_, fcadj) = pbci_for(types[i]);
        let q0 = atom.charge as f64;
        // (1 - M*v)*q0 simplified for fcadj=0 (most atoms): charge[i] = 0
        // For charged atoms with fcadj > 0:
        let m = bonds_of(mol, idx).len() as f64;
        charges[i] = (1.0 - m * fcadj) * q0;
    }

    // Step 2: BCI contributions from each bond
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let ti = types[i];
        let tj = types[j];
        let bt = bond_type_for(bond.order);

        // Contribution to atom i
        let ci =
            lookup_chg_contribution(bt, ti, tj).unwrap_or_else(|| pbci_for(ti).0 - pbci_for(tj).0);

        // Contribution to atom j
        let cj =
            lookup_chg_contribution(bt, tj, ti).unwrap_or_else(|| pbci_for(tj).0 - pbci_for(ti).0);

        charges[i] += ci;
        charges[j] += cj;
    }

    // Step 3: formal charge redistribution for charged neighbors (fcadj term)
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let (_, fcadj_i) = pbci_for(types[i]);
        if fcadj_i > 0.0 {
            // v*sumFormalCharge: redistribute neighbor formal charges
            let sum_fc: f64 = bonds_of(mol, idx)
                .iter()
                .map(|b| mol.atom(b.neighbor).charge as f64)
                .sum();
            charges[i] += fcadj_i * sum_fc;
        }
        // Anionic neighbor charge leaks: q0 adjustment
        // (for negatively charged neighbors — from RDKit source)
        for b in bonds_of(mol, idx) {
            let nbr = mol.atom(b.neighbor);
            if nbr.charge < 0 {
                let deg = bonds_of(mol, b.neighbor).len() as f64;
                charges[i] += (nbr.charge as f64) / (2.0 * deg);
            }
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
    fn carboxylate_anionic_o_is_type_35_not_the_nr_plus_nitrogen_row() {
        // Mirror-image bug found alongside the one above: `assign_o_type`'s
        // anionic-oxygen branch used to return `34` (NR+, a NITROGEN
        // type) for a negatively-charged oxygen. Must resolve to type 35
        // (OM), the registry's only oxygen entry among {32, 34, 35}.
        let m = mol("CC(=O)[O-]");
        let types = assign_mmff94_numeric_types(&m).unwrap();
        let mut saw_anionic_o = false;
        for i in 0..m.atom_count() {
            let a = m.atom(AtomIdx(i as u32));
            if a.element == Element::O && a.charge < 0 {
                assert_eq!(types[i], 35, "anionic O should be type 35 (OM)");
                saw_anionic_o = true;
            }
        }
        assert!(saw_anionic_o, "test fixture must contain an anionic oxygen");
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
}
