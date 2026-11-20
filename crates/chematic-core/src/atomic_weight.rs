//! Atomic weights and nuclide masses shared by descriptors and CIP.

use crate::Element;
use crate::rdkit_isotope_mass_table::RDKIT_ISOTOPE_MASS_TABLE;

/// Average atomic weights in daltons, indexed by atomic number minus one.
/// Historical chematic weights are preserved. Other entries come from the
/// RDKit 2026.03.3 periodic table.
static AVG_MASS_TABLE: [f64; 118] = [
    1.0080,   // 1 H (pre-existing)
    4.0030,   // 2 He (pre-existing)
    6.9410,   // 3 Li (pre-existing)
    9.0120,   // 4 Be (pre-existing)
    10.8110,  // 5 B (pre-existing)
    12.0110,  // 6 C (pre-existing)
    14.0070,  // 7 N (pre-existing)
    15.9990,  // 8 O (pre-existing)
    18.9980,  // 9 F (pre-existing)
    20.1800,  // 10 Ne (pre-existing)
    22.9900,  // 11 Na (pre-existing)
    24.3050,  // 12 Mg (pre-existing)
    26.9820,  // 13 Al (pre-existing)
    28.0860,  // 14 Si (pre-existing)
    30.9740,  // 15 P (pre-existing)
    32.0650,  // 16 S (pre-existing)
    35.4530,  // 17 Cl (pre-existing)
    39.9480,  // 18 Ar (pre-existing)
    39.0980,  // 19 K (pre-existing)
    40.0780,  // 20 Ca (pre-existing)
    44.9560,  // 21 Sc (RDKit)
    47.8670,  // 22 Ti (RDKit)
    50.9440,  // 23 V (RDKit)
    51.9960,  // 24 Cr (RDKit)
    54.9380,  // 25 Mn (RDKit)
    55.8450,  // 26 Fe (RDKit)
    58.9330,  // 27 Co (RDKit)
    58.6930,  // 28 Ni (RDKit)
    63.5460,  // 29 Cu (RDKit)
    65.3900,  // 30 Zn (RDKit)
    69.7230,  // 31 Ga (RDKit)
    72.6100,  // 32 Ge (RDKit)
    74.9220,  // 33 As (pre-existing)
    78.9710,  // 34 Se (pre-existing)
    79.9040,  // 35 Br (pre-existing)
    83.8000,  // 36 Kr (RDKit)
    85.4680,  // 37 Rb (RDKit)
    87.6200,  // 38 Sr (RDKit)
    88.9060,  // 39 Y (RDKit)
    91.2240,  // 40 Zr (RDKit)
    92.9060,  // 41 Nb (RDKit)
    95.9400,  // 42 Mo (RDKit)
    98.0000,  // 43 Tc (RDKit)
    101.0700, // 44 Ru (RDKit)
    102.9060, // 45 Rh (RDKit)
    106.4200, // 46 Pd (RDKit)
    107.8680, // 47 Ag (RDKit)
    112.4120, // 48 Cd (RDKit)
    114.8180, // 49 In (RDKit)
    118.7110, // 50 Sn (RDKit)
    121.7600, // 51 Sb (RDKit)
    127.6000, // 52 Te (RDKit)
    126.9040, // 53 I (pre-existing)
    131.2900, // 54 Xe (RDKit)
    132.9050, // 55 Cs (RDKit)
    137.3280, // 56 Ba (RDKit)
    138.9060, // 57 La (RDKit)
    140.1160, // 58 Ce (RDKit)
    140.9080, // 59 Pr (RDKit)
    144.2400, // 60 Nd (RDKit)
    145.0000, // 61 Pm (RDKit)
    150.3600, // 62 Sm (RDKit)
    151.9640, // 63 Eu (RDKit)
    157.2500, // 64 Gd (RDKit)
    158.9250, // 65 Tb (RDKit)
    162.5000, // 66 Dy (RDKit)
    164.9300, // 67 Ho (RDKit)
    167.2600, // 68 Er (RDKit)
    168.9340, // 69 Tm (RDKit)
    173.0400, // 70 Yb (RDKit)
    174.9670, // 71 Lu (RDKit)
    178.4900, // 72 Hf (RDKit)
    180.9480, // 73 Ta (RDKit)
    183.8400, // 74 W (RDKit)
    186.2070, // 75 Re (RDKit)
    190.2300, // 76 Os (RDKit)
    192.2170, // 77 Ir (RDKit)
    195.0780, // 78 Pt (RDKit)
    196.9670, // 79 Au (RDKit)
    200.5900, // 80 Hg (RDKit)
    204.3830, // 81 Tl (RDKit)
    207.2000, // 82 Pb (RDKit)
    208.9800, // 83 Bi (RDKit)
    209.0000, // 84 Po (RDKit)
    210.0000, // 85 At (RDKit)
    222.0000, // 86 Rn (RDKit)
    223.0000, // 87 Fr (RDKit)
    226.0000, // 88 Ra (RDKit)
    227.0000, // 89 Ac (RDKit)
    232.0380, // 90 Th (RDKit)
    231.0360, // 91 Pa (RDKit)
    238.0290, // 92 U (RDKit)
    237.0000, // 93 Np (RDKit)
    244.0000, // 94 Pu (RDKit)
    243.0000, // 95 Am (RDKit)
    247.0000, // 96 Cm (RDKit)
    247.0000, // 97 Bk (RDKit)
    251.0000, // 98 Cf (RDKit)
    252.0000, // 99 Es (RDKit)
    257.0000, // 100 Fm (RDKit)
    258.0000, // 101 Md (RDKit)
    259.0000, // 102 No (RDKit)
    262.0000, // 103 Lr (RDKit)
    267.0000, // 104 Rf (RDKit)
    268.0000, // 105 Db (RDKit)
    269.0000, // 106 Sg (RDKit)
    270.0000, // 107 Bh (RDKit)
    269.0000, // 108 Hs (RDKit)
    278.0000, // 109 Mt (RDKit)
    281.0000, // 110 Ds (RDKit)
    281.0000, // 111 Rg (RDKit)
    285.0000, // 112 Cn (RDKit)
    284.0000, // 113 Nh (RDKit)
    289.0000, // 114 Fl (RDKit)
    288.0000, // 115 Mc (RDKit)
    293.0000, // 116 Lv (RDKit)
    292.0000, // 117 Ts (RDKit)
    294.0000, // 118 Og (RDKit)
];

/// The average atomic weight used by native molecular-weight descriptors.
pub fn average_atomic_weight(element: Element) -> f64 {
    AVG_MASS_TABLE[usize::from(element.atomic_number()) - 1]
}

/// The average atomic weight used by the RDKit compatibility profile.
pub fn rdkit_average_atomic_weight(element: Element) -> f64 {
    match element.atomic_number() {
        5 => 10.812,
        16 => 32.067,
        34 => 78.960,
        _ => average_atomic_weight(element),
    }
}

/// RDKit periodic-table nuclide mass in daltons. Unknown isotope labels
/// use their integer mass number, matching the descriptor and CIP callers.
pub fn rdkit_isotope_mass(element: Element, isotope: u16) -> f64 {
    let key = (element.atomic_number(), isotope);
    RDKIT_ISOTOPE_MASS_TABLE
        .binary_search_by_key(&key, |&(z, a, _)| (z, a))
        .map_or(f64::from(isotope), |index| {
            RDKIT_ISOTOPE_MASS_TABLE[index].2
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generated_nuclide_table_is_sorted_and_complete() {
        assert_eq!(RDKIT_ISOTOPE_MASS_TABLE.len(), 3111);
        for pair in RDKIT_ISOTOPE_MASS_TABLE.windows(2) {
            let (left_z, left_a, _) = pair[0];
            let (right_z, right_a, _) = pair[1];
            assert!((left_z, left_a) < (right_z, right_a));
        }
        for &(z, a, mass) in &RDKIT_ISOTOPE_MASS_TABLE {
            let element = Element::from_atomic_number(z).expect("valid atomic number");
            assert_eq!(rdkit_isotope_mass(element, a), mass);
        }
    }
}
