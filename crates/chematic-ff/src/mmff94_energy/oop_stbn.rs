//! MMFF94 Out-of-Plane and Stretch-Bend parameters (Halgren MMFF.V/VI).

/// MMFF94 Out-of-Plane parameters (117 entries, Halgren MMFF.VI)
/// Format: (type_i, type_j, type_k, type_l, koop)
/// type_i = central sp2 atom; j/k/l = three attached atoms (0=wildcard)
/// E_oop = (0.043844 × koop / 2) × χ² (χ in degrees, Halgren eq. 6)
pub static MMFF94_OOP: &[(u8, u8, u8, u8, f64)] = &[
    (0, 2, 0, 0, 0.0200),
    (0, 3, 0, 0, 0.1300),
    (0, 8, 0, 0, 0.0000),
    (0, 10, 0, 0, -0.0200),
    (0, 17, 0, 0, 0.0000),
    (0, 26, 0, 0, 0.0000),
    (0, 30, 0, 0, 0.0100),
    (0, 37, 0, 0, 0.0350),
    (0, 39, 0, 0, 0.0200),
    (0, 40, 0, 0, -0.0050),
    (0, 41, 0, 0, 0.1800),
    (0, 43, 0, 0, 0.0000),
    (0, 45, 0, 0, 0.1500),
    (0, 49, 0, 0, 0.0000),
    (0, 54, 0, 0, 0.0200),
    (0, 55, 0, 0, 0.0200),
    (0, 56, 0, 0, 0.0200),
    (0, 57, 0, 0, 0.0800),
    (0, 58, 0, 0, 0.0250),
    (0, 63, 0, 0, 0.0500),
    (0, 64, 0, 0, 0.0400),
    (0, 67, 0, 0, 0.0700),
    (0, 69, 0, 0, 0.0700),
    (0, 73, 0, 0, 0.0000),
    (0, 78, 0, 0, 0.0450),
    (0, 80, 0, 0, 0.0800),
    (0, 81, 0, 0, 0.0250),
    (0, 82, 0, 0, 0.0000),
    (1, 2, 1, 2, 0.0300),
    (1, 2, 2, 2, 0.0270),
    (1, 2, 2, 3, 0.0260),
    (1, 2, 2, 5, 0.0130),
    (1, 2, 2, 37, 0.0320),
    (1, 3, 1, 7, 0.1460),
    (1, 3, 2, 7, 0.1380),
    (1, 3, 3, 7, 0.1340),
    (1, 3, 5, 7, 0.1220),
    (1, 3, 6, 7, 0.1410),
    (1, 3, 7, 10, 0.1290),
    (1, 3, 7, 37, 0.1380),
    (1, 10, 1, 3, -0.0200),
    (1, 10, 3, 6, -0.0330),
    (1, 10, 3, 28, -0.0200),
    (1, 37, 37, 37, 0.0400),
    (1, 39, 63, 63, 0.0120),
    (1, 40, 28, 37, -0.0060),
    (1, 41, 32, 32, 0.1780),
    (1, 54, 3, 36, 0.0160),
    (1, 55, 36, 57, 0.0200),
    (1, 56, 36, 57, 0.0200),
    (2, 2, 2, 5, 0.0130),
    (2, 2, 3, 5, 0.0120),
    (2, 2, 5, 5, 0.0060),
    (2, 2, 5, 6, 0.0270),
    (2, 2, 5, 37, 0.0170),
    (2, 2, 5, 40, 0.0120),
    (2, 2, 5, 41, 0.0080),
    (2, 3, 5, 7, 0.1130),
    (2, 3, 5, 9, 0.0810),
    (2, 3, 6, 7, 0.1270),
    (2, 3, 7, 10, 0.1160),
    (2, 37, 37, 37, 0.0310),
    (2, 40, 28, 28, -0.0070),
    (2, 41, 32, 32, 0.1610),
    (3, 3, 5, 7, 0.1130),
    (3, 3, 6, 7, 0.1270),
    (3, 10, 3, 28, -0.0300),
    (3, 10, 28, 28, -0.0190),
    (3, 37, 37, 37, 0.0270),
    (3, 40, 28, 28, -0.0070),
    (3, 54, 36, 36, 0.0180),
    (5, 3, 5, 7, 0.1030),
    (5, 3, 5, 9, 0.0740),
    (5, 3, 5, 54, 0.0780),
    (5, 3, 6, 7, 0.1190),
    (5, 3, 7, 10, 0.1020),
    (5, 3, 9, 40, 0.0670),
    (5, 30, 20, 30, 0.0080),
    (5, 37, 37, 37, 0.0150),
    (5, 37, 37, 38, 0.0460),
    (5, 37, 37, 63, 0.0080),
    (5, 37, 37, 64, 0.0120),
    (5, 37, 37, 69, 0.0160),
    (5, 37, 38, 38, 0.0840),
    (5, 41, 32, 32, 0.1580),
    (5, 57, 55, 55, 0.0380),
    (5, 63, 39, 64, 0.0190),
    (5, 63, 39, 66, 0.0680),
    (5, 63, 44, 64, 0.0140),
    (5, 63, 44, 66, 0.0550),
    (5, 63, 59, 64, 0.0330),
    (5, 63, 59, 66, 0.0850),
    (5, 64, 63, 64, 0.0060),
    (5, 64, 63, 66, 0.0430),
    (5, 64, 64, 65, 0.0520),
    (5, 64, 65, 66, 0.0940),
    (5, 78, 78, 81, 0.0460),
    (5, 80, 81, 81, 0.0570),
    (6, 3, 7, 37, 0.1270),
    (6, 37, 37, 37, 0.0480),
    (7, 3, 10, 10, 0.1130),
    (7, 3, 20, 20, 0.1510),
    (9, 3, 40, 40, 0.0570),
    (15, 37, 37, 37, 0.0250),
    (23, 39, 63, 63, -0.0140),
    (23, 39, 63, 65, 0.0210),
    (23, 39, 65, 65, 0.0620),
    (28, 40, 28, 37, 0.0040),
    (32, 69, 37, 37, 0.0670),
    (36, 55, 36, 57, 0.0200),
    (36, 56, 36, 57, 0.0200),
    (36, 81, 78, 80, 0.0160),
    (37, 37, 37, 40, 0.0460),
    (37, 63, 39, 64, 0.0100),
    (37, 64, 63, 64, -0.0110),
    (50, 49, 50, 50, 0.0000),
    (56, 57, 56, 56, 0.1580),
];

/// MMFF94 Stretch-Bend parameters (282 entries, Halgren MMFF.V)
/// Format: (stretch_bend_type, type_i, type_j, type_k, kba_ijk, kba_kji)
/// E_sb = 2.51210 × (kba_ijk × Δr_ij + kba_kji × Δr_kj) × Δθ (degrees)
///
/// The leading key column is RDKit's **stretch-bend type** (0-11, from
/// `getMMFFStretchBendType` — see [`crate::mmff94_minimizer::stretch_bend_type_for`]),
/// NOT the angle type (0-8, from [`crate::mmff94_minimizer::angle_type_for`]).
/// A prior version of this doc comment (and of `mmff94_stbn`/
/// `mmff94_stbn_type_only`'s own parameter) mislabeled this column as
/// `angle_type`, which fed the wrong key straight into every lookup below —
/// fixed here (issue #227 Priority 2C). Self-consistency proof (no oracle
/// needed): the one key-5 row below, `(5, 22, 22, 22, ...)`, is an all-CR3R
/// (cyclopropane ring carbon) triple; under `angle_type_for`'s own table a
/// 3-ring with `bt_sum=0` is angle type 3, never 5 — a row keyed 5 would be
/// unreachable garbage if this column really were `angle_type`. Under
/// `getMMFFStretchBendType`, angle type 3 (3-ring, bt_sum=0) maps to
/// stretch-bend type 5 — exactly this row.
pub static MMFF94_STBN: &[(u8, u8, u8, u8, f64, f64)] = &[
    (0, 1, 1, 1, 0.2060, 0.2060),
    (0, 1, 1, 2, 0.1360, 0.1970),
    (0, 1, 1, 3, 0.2110, 0.0920),
    (0, 1, 1, 5, 0.2270, 0.0700),
    (0, 1, 1, 6, 0.1730, 0.4170),
    (0, 1, 1, 8, 0.1360, 0.2820),
    (0, 1, 1, 10, 0.1870, 0.3380),
    (0, 1, 1, 11, 0.2090, 0.6330),
    (0, 1, 1, 12, 0.1760, 0.3860),
    (0, 1, 1, 15, 0.1390, 0.2170),
    (0, 1, 1, 34, 0.2360, 0.4360),
    (0, 1, 1, 37, 0.1520, 0.2600),
    (0, 1, 1, 39, 0.1440, 0.5950),
    (0, 1, 1, 41, 0.1220, 0.0510),
    (0, 1, 1, 56, 0.2620, 0.4510),
    (0, 1, 1, 68, 0.1860, 0.1250),
    (0, 1, 2, 1, 0.2500, 0.2500),
    (0, 1, 2, 2, 0.2030, 0.2070),
    (0, 1, 2, 5, 0.2150, 0.1280),
    (0, 1, 3, 1, 0.3580, 0.3580),
    (0, 1, 3, 5, 0.3210, 0.1830),
    (0, 1, 3, 6, 0.3380, 0.7320),
    (0, 1, 3, 7, 0.1540, 0.8560),
    (0, 1, 3, 10, 0.2230, 0.7320),
    (0, 1, 6, 1, 0.3090, 0.3090),
    (0, 1, 6, 2, 0.1570, 0.3750),
    (0, 1, 6, 3, -0.1530, 0.2520),
    (0, 1, 6, 21, 0.2560, 0.1430),
    (0, 1, 6, 37, 0.1630, 0.3750),
    (0, 1, 8, 1, 0.3120, 0.3120),
    (0, 1, 8, 6, 0.2120, 0.3540),
    (0, 1, 8, 23, 0.3090, 0.1350),
    (0, 1, 9, 3, 0.3260, 0.5800),
    (0, 1, 10, 1, 0.0630, 0.0630),
    (0, 1, 10, 3, -0.0210, 0.3400),
    (0, 1, 10, 6, -0.0240, 0.3740),
    (0, 1, 10, 28, 0.1550, -0.0510),
    (0, 1, 15, 1, 0.1250, 0.1250),
    (0, 1, 15, 15, 0.0120, 0.2380),
    (0, 1, 15, 37, 0.0480, 0.2290),
    (0, 1, 15, 71, 0.0800, -0.0120),
    (0, 1, 18, 1, 0.0230, 0.0230),
    (0, 1, 18, 6, 0.0030, 0.2130),
    (0, 1, 18, 32, -0.0910, 0.3900),
    (0, 1, 18, 43, -0.0080, 0.6070),
    (0, 1, 20, 5, 0.2900, 0.0980),
    (0, 1, 20, 20, 0.1790, 0.0040),
    (0, 1, 22, 5, 0.0670, 0.1740),
    (0, 1, 22, 22, 0.1990, 0.0390),
    (0, 1, 34, 1, 0.2020, 0.2020),
    (0, 1, 34, 36, 0.1600, -0.0090),
    (0, 1, 37, 37, 0.4850, 0.3110),
    (0, 1, 39, 63, 0.3130, 0.5000),
    (0, 1, 40, 28, 0.2380, 0.0910),
    (0, 1, 40, 37, 0.1530, 0.5900),
    (0, 1, 41, 32, 0.5030, 0.9430),
    (0, 1, 54, 3, 0.1920, -0.0510),
    (0, 1, 54, 36, 0.2400, 0.0790),
    (0, 1, 55, 36, 0.1890, 0.0330),
    (0, 1, 55, 57, 0.1660, 0.2110),
    (0, 1, 56, 36, 0.2110, -0.0400),
    (0, 1, 56, 57, 0.0260, 0.3860),
    (0, 1, 68, 1, 0.2170, 0.2170),
    (0, 1, 68, 23, 0.2850, 0.0500),
    (0, 1, 68, 32, -0.0470, 0.5030),
    (0, 2, 1, 2, 0.2820, 0.2820),
    (0, 2, 1, 3, 0.2060, 0.0220),
    (0, 2, 1, 5, 0.2340, 0.0880),
    (0, 2, 1, 6, 0.1830, 0.3870),
    (0, 2, 1, 8, 0.2140, 0.3630),
    (0, 2, 2, 5, 0.2070, 0.1570),
    (0, 2, 2, 6, 0.1180, 0.5760),
    (0, 2, 2, 40, 0.2890, 0.3900),
    (0, 2, 2, 41, 0.1910, -0.0470),
    (0, 2, 6, 3, -0.2280, 0.0520),
    (0, 2, 6, 29, 0.2590, 0.1630),
    (0, 2, 40, 28, 0.3420, 0.1560),
    (0, 2, 41, 32, 0.5940, 0.9690),
    (0, 3, 1, 5, 0.1570, 0.1150),
    (0, 3, 1, 6, -0.0360, 0.4560),
    (0, 3, 1, 10, 0.0380, 0.1950),
    (0, 3, 6, 24, 0.2150, 0.0640),
    (0, 3, 6, 37, -0.2250, -0.3200),
    (0, 3, 9, 27, 0.4640, 0.2220),
    (0, 3, 10, 3, -0.2190, -0.2190),
    (0, 3, 10, 6, 0.4970, 0.5130),
    (0, 3, 10, 28, 0.1370, 0.0660),
    (0, 3, 20, 5, -0.0490, 0.1710),
    (0, 3, 40, 28, 0.2280, 0.1040),
    (0, 3, 54, 36, 0.0050, 0.1270),
    (0, 5, 1, 5, 0.1150, 0.1150),
    (0, 5, 1, 6, 0.0130, 0.4360),
    (0, 5, 1, 8, 0.0270, 0.3580),
    (0, 5, 1, 9, 0.0400, 0.4180),
    (0, 5, 1, 10, 0.0430, 0.2610),
    (0, 5, 1, 11, 0.0030, 0.4520),
    (0, 5, 1, 12, -0.0180, 0.3800),
    (0, 5, 1, 15, 0.0180, 0.2550),
    (0, 5, 1, 18, 0.1210, 0.2180),
    (0, 5, 1, 20, 0.0690, 0.3270),
    (0, 5, 1, 22, 0.0550, 0.2670),
    (0, 5, 1, 34, -0.0030, 0.3420),
    (0, 5, 1, 37, 0.0740, 0.2870),
    (0, 5, 1, 39, 0.0920, 0.6070),
    (0, 5, 1, 40, 0.0230, 0.3350),
    (0, 5, 1, 41, 0.0930, 0.1180),
    (0, 5, 1, 54, 0.0160, 0.3430),
    (0, 5, 1, 55, 0.0300, 0.3970),
    (0, 5, 1, 56, 0.0310, 0.3840),
    (0, 5, 1, 68, 0.0410, 0.2160),
    (0, 5, 2, 5, 0.1400, 0.1400),
    (0, 5, 2, 6, 0.2130, 0.5020),
    (0, 5, 2, 40, 0.0700, 0.4630),
    (0, 5, 2, 41, 0.1910, 0.0050),
    (0, 5, 3, 5, 0.1260, 0.1260),
    (0, 5, 3, 6, 0.1740, 0.7340),
    (0, 5, 3, 7, 0.0320, 0.8050),
    (0, 5, 3, 9, 0.0370, 0.6690),
    (0, 5, 3, 10, 0.1690, 0.6190),
    (0, 5, 3, 40, 0.0870, 0.6850),
    (0, 5, 3, 54, 0.0980, 0.2100),
    (0, 5, 20, 5, 0.1820, 0.1820),
    (0, 5, 20, 6, 0.0510, 0.3120),
    (0, 5, 20, 8, 0.0720, 0.2260),
    (0, 5, 20, 12, 0.0140, 0.5970),
    (0, 5, 20, 20, 0.1010, 0.0790),
    (0, 5, 20, 30, 0.1080, 0.1230),
    (0, 5, 22, 5, 0.2540, 0.2540),
    (0, 5, 22, 22, 0.1810, 0.1080),
    (0, 5, 26, 5, -0.1210, -0.1210),
    (0, 5, 30, 20, 0.2510, 0.0070),
    (0, 5, 30, 30, 0.2670, 0.0540),
    (0, 5, 37, 37, 0.2790, 0.2500),
    (0, 5, 37, 38, 0.2670, 0.3890),
    (0, 5, 37, 63, 0.2160, 0.4340),
    (0, 5, 37, 64, 0.1670, 0.3640),
    (0, 5, 37, 69, 0.2730, 0.3910),
    (0, 5, 41, 32, 0.2760, 0.8520),
    (0, 5, 57, 55, 0.0430, 0.4200),
    (0, 5, 63, 39, 0.0090, 0.6540),
    (0, 5, 63, 44, -0.0150, 0.4460),
    (0, 5, 63, 59, 0.0670, 0.5880),
    (0, 5, 63, 64, 0.0550, 0.3700),
    (0, 5, 63, 66, 0.1100, 0.4640),
    (0, 5, 64, 63, 0.0860, 0.3450),
    (0, 5, 64, 64, 0.0850, 0.3690),
    (0, 5, 64, 65, 0.0510, 0.4360),
    (0, 5, 64, 66, 0.1130, 0.4520),
    (0, 5, 78, 78, 0.2790, 0.2500),
    (0, 5, 78, 81, 0.0830, 0.2500),
    (0, 5, 80, 81, -0.1010, 0.6910),
    (0, 6, 1, 6, 0.3200, 0.3200),
    (0, 6, 1, 37, 0.3100, 0.1600),
    (0, 6, 3, 7, 0.4940, 0.5780),
    (0, 6, 8, 23, 0.4180, 0.0200),
    (0, 6, 18, 6, 0.0880, 0.0880),
    (0, 6, 18, 32, 0.1230, 0.3690),
    (0, 6, 37, 37, 0.8300, 0.3390),
    (0, 7, 3, 10, 0.7710, 0.3530),
    (0, 7, 3, 20, 0.8650, -0.1810),
    (0, 8, 6, 21, 0.3040, 0.0550),
    (0, 9, 3, 40, 0.6800, 0.2600),
    (0, 10, 3, 10, 1.0500, 1.0500),
    (0, 10, 6, 21, 0.4190, 0.1580),
    (0, 11, 1, 11, 0.5860, 0.5860),
    (0, 12, 1, 12, 0.5080, 0.5080),
    (0, 12, 20, 20, 0.3100, 0.0000),
    (0, 15, 15, 71, 0.1720, -0.0680),
    (0, 15, 37, 37, 0.6500, 0.2590),
    (0, 18, 6, 33, 0.3090, 0.1200),
    (0, 18, 43, 23, 0.3770, 0.0570),
    (0, 20, 8, 23, 0.1280, 0.1220),
    (0, 23, 8, 23, 0.1900, 0.1900),
    (0, 23, 39, 63, -0.1310, 0.4220),
    (0, 23, 39, 65, -0.1220, 0.2810),
    (0, 23, 43, 23, 0.0820, 0.0820),
    (0, 23, 68, 23, 0.1450, 0.1450),
    (0, 23, 68, 32, -0.1820, 0.5040),
    (0, 28, 10, 28, 0.0810, 0.0810),
    (0, 28, 40, 28, 0.0940, 0.0940),
    (0, 28, 40, 37, 0.1860, 0.4230),
    (0, 29, 6, 37, 0.1300, 0.2410),
    (0, 31, 6, 31, 0.2270, 0.2270),
    (0, 31, 70, 31, 0.2100, 0.2100),
    (0, 32, 18, 32, 0.4040, 0.4040),
    (0, 32, 18, 43, 0.3840, 0.2810),
    (0, 32, 41, 32, 0.6520, 0.6520),
    (0, 32, 69, 37, 1.0180, 0.4180),
    (0, 36, 34, 36, 0.0870, 0.0870),
    (0, 36, 54, 36, 0.1480, 0.1480),
    (0, 36, 55, 36, 0.1060, 0.1060),
    (0, 36, 55, 57, 0.0930, 0.0800),
    (0, 36, 56, 36, 0.1010, 0.1010),
    (0, 36, 56, 57, 0.1080, 0.0680),
    (0, 36, 81, 78, 0.0210, 0.3680),
    (0, 36, 81, 80, 0.0180, 0.4220),
    (0, 37, 15, 71, 0.1870, -0.0270),
    (0, 37, 37, 37, -0.4110, -0.4110),
    (0, 37, 37, 38, -0.4240, -0.4660),
    (0, 37, 37, 40, 0.4290, 0.9010),
    (0, 37, 37, 63, -0.1730, -0.2150),
    (0, 37, 37, 64, -0.2290, -0.2290),
    (0, 37, 37, 69, -0.2440, -0.5550),
    (0, 37, 38, 37, -0.3420, -0.3420),
    (0, 37, 38, 38, -0.1640, -1.1300),
    (0, 37, 63, 39, 0.1780, 0.5230),
    (0, 37, 63, 64, -0.0450, 0.4970),
    (0, 37, 64, 63, 0.0590, 0.2990),
    (0, 37, 64, 64, 0.2770, 0.3770),
    (0, 37, 69, 37, -0.1690, -0.1690),
    (0, 38, 37, 38, -0.5160, -0.5160),
    (0, 39, 63, 64, 0.4220, 0.4090),
    (0, 39, 63, 66, 0.4360, 0.5250),
    (0, 39, 65, 64, 0.5280, 0.6440),
    (0, 39, 65, 66, 0.3970, 0.2580),
    (0, 40, 3, 40, 0.4820, 0.4820),
    (0, 43, 18, 43, 0.4280, 0.4280),
    (0, 44, 63, 64, 0.5810, 0.4260),
    (0, 44, 63, 66, 0.5420, 0.3650),
    (0, 44, 65, 64, 0.8160, 0.5430),
    (0, 50, 49, 50, 0.0720, 0.0720),
    (0, 55, 57, 55, 0.1250, 0.1250),
    (0, 56, 57, 56, 0.4310, 0.4310),
    (0, 58, 57, 58, 0.7320, 0.7320),
    (0, 59, 63, 64, 0.8520, 0.3320),
    (0, 59, 63, 66, 0.7750, 0.3000),
    (0, 59, 65, 64, 1.1770, 0.5940),
    (0, 63, 39, 63, 0.4690, 0.4690),
    (0, 63, 39, 65, 0.7410, 0.5060),
    (0, 63, 44, 63, 0.5910, 0.5910),
    (0, 63, 44, 65, 0.8570, 0.9780),
    (0, 63, 59, 63, 0.4970, 0.4970),
    (0, 63, 59, 65, 0.7230, 0.8740),
    (0, 63, 64, 64, 0.2060, 0.0300),
    (0, 63, 64, 66, 0.1710, 0.0780),
    (0, 63, 66, 64, 0.2130, -0.1730),
    (0, 63, 66, 66, 0.2340, 0.0770),
    (0, 64, 64, 65, 0.0790, 0.4030),
    (0, 64, 66, 65, -0.1490, 0.3830),
    (0, 65, 39, 65, 0.7060, 0.7060),
    (0, 65, 64, 66, 0.4060, 0.0660),
    (0, 65, 66, 66, 0.1990, 0.1010),
    (0, 71, 15, 71, 0.0450, 0.0450),
    (0, 78, 78, 81, -0.3980, 0.3140),
    (0, 78, 81, 80, 0.3660, 0.4190),
    (0, 81, 80, 81, 0.7320, 0.7320),
    (1, 2, 2, 2, 0.2500, 0.2190),
    (1, 2, 2, 5, 0.2670, 0.1590),
    (1, 2, 3, 5, 0.4070, 0.1590),
    (1, 2, 3, 6, 0.4290, 0.4730),
    (1, 2, 3, 7, 0.2140, 0.7940),
    (1, 2, 3, 9, 0.2270, 0.6100),
    (1, 2, 3, 10, 0.2980, 0.6000),
    (1, 2, 37, 37, 0.3210, 0.2350),
    (1, 3, 2, 5, 0.2640, 0.1560),
    (1, 3, 3, 5, 0.2510, 0.1330),
    (1, 3, 3, 6, 0.0660, 0.6680),
    (1, 3, 3, 7, -0.0930, 0.8660),
    (1, 3, 37, 37, 0.1790, 0.2170),
    (2, 1, 2, 2, 0.2220, 0.2690),
    (2, 1, 2, 3, 0.2440, 0.2920),
    (2, 1, 2, 37, 0.2460, 0.2600),
    (2, 1, 3, 2, 0.2460, 0.4090),
    (2, 1, 3, 3, 0.3030, 0.1450),
    (2, 1, 3, 37, 0.2170, 0.2070),
    (2, 2, 2, 3, 0.1550, 0.1120),
    (2, 2, 2, 37, 0.1430, 0.1720),
    (2, 5, 2, 37, 0.1530, 0.2880),
    (2, 6, 3, 37, 0.3500, 0.1750),
    (2, 7, 3, 37, 0.7070, 0.0070),
    (4, 3, 6, 20, 0.4560, 0.3790),
    (4, 3, 20, 20, 0.6070, 0.4370),
    (4, 6, 3, 20, 1.1790, 0.7520),
    (4, 6, 20, 20, 0.8230, 0.3960),
    (4, 8, 20, 20, 0.7010, 0.3690),
    (4, 20, 3, 20, 0.5360, 0.5360),
    (4, 20, 6, 20, 0.7390, 0.7390),
    (4, 20, 8, 20, 0.6530, 0.6530),
    (4, 20, 20, 20, 0.2830, 0.2830),
    (4, 20, 20, 30, 0.3400, 0.5290),
    (4, 20, 30, 30, 0.4130, 0.7050),
    (5, 22, 22, 22, 0.0000, 0.0000),
];

/// Look up OOP bending parameter for central sp2 atom j with neighbors i, k, l.
/// Wildcard matching (0) tried as fallback.
pub fn mmff94_oop(type_j: u8, type_i: u8, type_k: u8, type_l: u8) -> Option<f64> {
    // Normalize: sort (i, k, l) except j stays central; try all orderings via wildcard
    for &(ti, tk, tl) in &[
        (type_i, type_k, type_l),
        (type_k, type_i, type_l),
        (type_l, type_k, type_i),
        (type_i, type_l, type_k),
        (type_k, type_l, type_i),
        (type_l, type_i, type_k),
    ] {
        if let Some(koop) = search_oop(type_j, ti, tk, tl) {
            return Some(koop);
        }
    }
    // Wildcard fallback
    for &(ti, tk, tl) in &[(0, 0, 0), (type_i, 0, 0), (0, type_k, 0), (0, 0, type_l)] {
        if let Some(koop) = search_oop(type_j, ti, tk, tl) {
            return Some(koop);
        }
    }
    None
}

fn search_oop(type_j: u8, type_i: u8, type_k: u8, type_l: u8) -> Option<f64> {
    MMFF94_OOP
        .binary_search_by_key(&(type_i, type_j, type_k, type_l), |&(i, j, k, l, _)| {
            (i, j, k, l)
        })
        .ok()
        .map(|idx| MMFF94_OOP[idx].4)
}

/// Look up Stretch-Bend parameters for angle i-j-k by MMFF *type* alone —
/// no element/periodic-row fallback (see [`mmff94_stbn`] for that).
///
/// Returns (kba_ijk, kba_kji). Both orderings (i,j,k) and (k,j,i) tried at
/// the requested `stretch_bend_type`, then — if `stretch_bend_type` isn't 0
/// and no row exists there — the *specific* (ti,tj,tk) triple is retried at
/// type 0 before finally falling back to the fully generic `(0, 0, type_j, 0)`
/// wildcard. `MMFF94_STBN` is overwhelmingly type-0 (246/282 rows), so
/// without this intermediate step, correctly classifying a term as a
/// non-zero type that this table doesn't happen to cover would silently
/// drop straight to the least specific fallback instead of the
/// specific-triple type-0 row a hardcoded `stretch_bend_type=0` caller would
/// have found.
///
/// `stretch_bend_type` is RDKit's `getMMFFStretchBendType` output (0-11),
/// computed via [`crate::mmff94_minimizer::stretch_bend_type_for`] — **not**
/// the angle type (0-8, [`crate::mmff94_minimizer::angle_type_for`]). Issue
/// #227 Priority 2C: this parameter used to be the angle type directly
/// (mislabeled as `angle_type` here), which used the wrong lookup key —
/// `MMFF94_STBN`'s own leading column is keyed by stretch-bend type, not
/// angle type (see the table's own doc comment for the self-consistency
/// proof). **Breaking change**: callers passing a raw angle type now get
/// silently wrong results; see `CHANGELOG.md` migration notes.
///
/// `pub` (not just used internally) because diagnostic tooling
/// (`mmff94_term_coverage_audit.rs`) specifically wants "does a *different
/// classification code* have a row for this exact type triple" —
/// independent of any element-based fallback, which [`mmff94_stbn`]'s
/// Dfsb tier is not (it doesn't vary with `stretch_bend_type` at all).
pub fn mmff94_stbn_type_only(
    stretch_bend_type: u8,
    type_i: u8,
    type_j: u8,
    type_k: u8,
) -> Option<(f64, f64)> {
    let search = |sbt: u8, ti: u8, tj: u8, tk: u8| {
        MMFF94_STBN
            .binary_search_by_key(&(sbt, ti, tj, tk), |&(a, i, j, k, _, _)| (a, i, j, k))
            .ok()
            .map(|idx| (MMFF94_STBN[idx].4, MMFF94_STBN[idx].5))
    };
    search(stretch_bend_type, type_i, type_j, type_k)
        .or_else(|| search(stretch_bend_type, type_k, type_j, type_i).map(|(a, b)| (b, a)))
        .or_else(|| {
            if stretch_bend_type != 0 {
                search(0, type_i, type_j, type_k)
                    .or_else(|| search(0, type_k, type_j, type_i).map(|(a, b)| (b, a)))
            } else {
                None
            }
        })
        .or_else(|| search(0, 0, type_j, 0))
}

/// RDKit's periodic-table-row default stretch-bend constants
/// (`defaultMMFFDfsb`, `Code/ForceField/MMFF/Params.cpp`), RDKit's own
/// residual fallback for stretch-bend once the specific/generic MMFF-type
/// table (`MMFF94_STBN`) has no row at all. 29 rows, verbatim from
/// `scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt` (programmatically
/// extracted from the pinned RDKit commit, not hand-transcribed — see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s "Stretch-bend" row).
/// `(periodic_row_i, periodic_row_j, periodic_row_k, f_ijk, f_kji)` with
/// `periodic_row_i <= periodic_row_k`, matching RDKit's own
/// `MMFFDfsbCollection::getMMFFDfsbParams` canonicalization.
static MMFF94_DFSB: &[(u8, u8, u8, f64, f64)] = &[
    (0, 1, 0, 0.15, 0.15),
    (0, 1, 1, 0.10, 0.30),
    (0, 1, 2, 0.05, 0.35),
    (0, 1, 3, 0.05, 0.35),
    (0, 1, 4, 0.05, 0.35),
    (0, 2, 0, 0.00, 0.00),
    (0, 2, 1, 0.00, 0.15),
    (0, 2, 2, 0.00, 0.15),
    (0, 2, 3, 0.00, 0.15),
    (0, 2, 4, 0.00, 0.15),
    (1, 1, 1, 0.30, 0.30),
    (1, 1, 2, 0.30, 0.50),
    (1, 1, 3, 0.30, 0.50),
    (1, 1, 4, 0.30, 0.50),
    (2, 1, 2, 0.50, 0.50),
    (2, 1, 3, 0.50, 0.50),
    (2, 1, 4, 0.50, 0.50),
    (3, 1, 3, 0.50, 0.50),
    (3, 1, 4, 0.50, 0.50),
    (4, 1, 4, 0.50, 0.50),
    (1, 2, 1, 0.30, 0.30),
    (1, 2, 2, 0.25, 0.25),
    (1, 2, 3, 0.25, 0.25),
    (1, 2, 4, 0.25, 0.25),
    (2, 2, 2, 0.25, 0.25),
    (2, 2, 3, 0.25, 0.25),
    (2, 2, 4, 0.25, 0.25),
    (3, 2, 3, 0.25, 0.25),
    (3, 2, 4, 0.25, 0.25),
];

/// RDKit's periodic-table-row bucketing (`getPeriodicTableRow`,
/// `Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp`, pinned commit —
/// see `scripts/mmff94_provenance/PROVENANCE.md`): atomic number 1-2 (H,
/// He) -> row 0, 3-10 -> row 1, 11-18 -> row 2, 19-36 -> row 3, 37-54 ->
/// row 4, anything heavier -> row 0 (RDKit's own default, not a chematic
/// omission — no row in `MMFF94_DFSB` is centered on an atom heavier than
/// Xe anyway, since stretch-bend terms only arise for organic-chemistry
/// central atoms in practice).
fn mmff94_periodic_table_row(atomic_number: u8) -> u8 {
    match atomic_number {
        3..=10 => 1,
        11..=18 => 2,
        19..=36 => 3,
        37..=54 => 4,
        _ => 0,
    }
}

/// RDKit's `MMFFDfsbCollection::getMMFFDfsbParams` — the periodic-row
/// default stretch-bend fallback. Only meaningful once
/// `mmff94_stbn_type_only` has already missed (matches RDKit's own order:
/// `MMFFMolProperties::getMMFFStretchBendParams` tries the specific/generic
/// MMFF-type lookup first, Dfsb only on failure). RDKit's own
/// `isDoubleZero(kbaIJK) && isDoubleZero(kbaKJI)` exclusion (a resolved-but
/// -both-zero row still counts as "not resolved") is replicated: the
/// table's one all-zero row, `(0, 2, 0)`, returns `None` here, not
/// `Some((0.0, 0.0))`.
fn mmff94_dfsb_stbn(atomic_num_i: u8, atomic_num_j: u8, atomic_num_k: u8) -> Option<(f64, f64)> {
    let (mut row_i, row_j, mut row_k) = (
        mmff94_periodic_table_row(atomic_num_i),
        mmff94_periodic_table_row(atomic_num_j),
        mmff94_periodic_table_row(atomic_num_k),
    );
    let swapped = row_i > row_k;
    if swapped {
        std::mem::swap(&mut row_i, &mut row_k);
    }
    let (f_ijk, f_kji) = MMFF94_DFSB
        .iter()
        .find(|&&(r1, r2, r3, ..)| r1 == row_i && r2 == row_j && r3 == row_k)
        .map(|&(.., f_ijk, f_kji)| (f_ijk, f_kji))?;
    if f_ijk == 0.0 && f_kji == 0.0 {
        return None;
    }
    Some(if swapped {
        (f_kji, f_ijk)
    } else {
        (f_ijk, f_kji)
    })
}

/// Look up Stretch-Bend parameters for angle i-j-k.
///
/// Returns (kba_ijk, kba_kji). Tries the type-based table first
/// ([`mmff94_stbn_type_only`]); if that misses entirely, falls back to
/// RDKit's own periodic-table-row default (Priority 2B, issue #227 — see
/// `mmff94_dfsb_stbn`'s doc). Verified against the pinned RDKit commit's
/// real `MMFFMolProperties::getMMFFStretchBendParams`: this is RDKit's
/// *complete* residual fallback chain for stretch-bend — no
/// equivalence-class step exists in RDKit's real stretch-bend path at all
/// (that mechanism is angle/torsion/OOP-only in RDKit), so this function
/// now matches RDKit's own coverage, not a partial subset of it.
///
/// `stretch_bend_type` (0-11) is RDKit's `getMMFFStretchBendType` output —
/// compute it via [`crate::mmff94_minimizer::stretch_bend_type_for`], **not**
/// the angle type ([`crate::mmff94_minimizer::angle_type_for`], 0-8). Issue
/// #227 Priority 2C: this parameter used to be the angle type directly
/// (mislabeled `angle_type`), which fed the wrong key into `MMFF94_STBN` —
/// see that table's own doc comment for the self-consistency proof this was
/// a real bug. **Breaking change**: existing callers passing a raw angle
/// type now get silently wrong parameters; see `CHANGELOG.md` migration
/// notes for the required call-site update.
pub fn mmff94_stbn(
    stretch_bend_type: u8,
    type_i: u8,
    type_j: u8,
    type_k: u8,
    atomic_num_i: u8,
    atomic_num_j: u8,
    atomic_num_k: u8,
) -> Option<(f64, f64)> {
    mmff94_stbn_type_only(stretch_bend_type, type_i, type_j, type_k)
        .or_else(|| mmff94_dfsb_stbn(atomic_num_i, atomic_num_j, atomic_num_k))
}

#[cfg(test)]
mod dfsb_tests {
    use super::*;

    #[test]
    fn dfsb_table_has_29_rows() {
        // Regression guard against accidental corruption of the ported
        // table -- 29 is the exact row count of
        // scripts/mmff94_provenance/rdkit_defaultMMFFDfsb.txt (verified at
        // port time via a Python regex-extraction script, not hand-counted).
        assert_eq!(MMFF94_DFSB.len(), 29);
    }

    #[test]
    fn dfsb_table_rows_are_canonical_row_i_le_row_k() {
        for &(row_i, _row_j, row_k, ..) in MMFF94_DFSB {
            assert!(
                row_i <= row_k,
                "row ({row_i}, _, {row_k}) violates row_i <= row_k canonicalization"
            );
        }
    }
}
