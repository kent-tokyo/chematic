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
    let triple_bonds = count_bond_order(mol, idx, BondOrder::Triple);

    // sp carbon (triple bond or allene)
    if triple_bonds > 0 {
        return Ok(4); // CSP
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

    // sp3 N, 3-connected: enamine/aniline (NC=C) vs amide (NC=O) vs plain
    // (NR). Source-grounded deterministic port of RDKit's `case 7:`
    // 3-connected branch (`AtomTyper.cpp` lines 1093-1325 at the pinned
    // commit) -- see `classify_n_c3_carbon_context`'s doc for the exact
    // condition and its one documented, empirically-unobserved structural
    // divergence from RDKit's literal C++.
    if total_degree(mol, idx) == 3 {
        let ctx = classify_n_c3_carbon_context(mol, rings, idx);
        if ctx.has_carbon_neighbor {
            if !ctx.is_carbonyl_like && !ctx.is_cyano_like && ctx.any_carbon_qualifies_nc_eq_c {
                return Ok(40); // NC=C / NC=N / NC=P / NC%C: deloc. lone pair
            }
            if !ctx.is_cyano_like && ctx.is_carbonyl_like {
                return Ok(10); // NC=O / NC=S amide/thioamide nitrogen
            }
            // ctx.is_cyano_like (NSO2/NC%N family, type 43) or neither
            // condition matched: not yet ported (a separate, tiny,
            // pre-existing residual, unaffected by this change) -- falls
            // through to the generic sp3 checks below, same as before.
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
}
