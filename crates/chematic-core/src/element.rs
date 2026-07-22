//! Periodic-table elements indexed by atomic number (1–118).

/// An element represented as its atomic number, wrapped in a newtype.
///
/// Zero-copy: implements `Copy`; comparison and hashing are O(1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Element(u8);

impl Element {
    // -- commonly used element constants ----------------------------------
    pub const H: Element = Element(1);
    pub const HE: Element = Element(2);
    pub const LI: Element = Element(3);
    pub const BE: Element = Element(4);
    pub const B: Element = Element(5);
    pub const C: Element = Element(6);
    pub const N: Element = Element(7);
    pub const O: Element = Element(8);
    pub const F: Element = Element(9);
    pub const NE: Element = Element(10);
    pub const NA: Element = Element(11);
    pub const MG: Element = Element(12);
    pub const AL: Element = Element(13);
    pub const SI: Element = Element(14);
    pub const P: Element = Element(15);
    pub const S: Element = Element(16);
    pub const CL: Element = Element(17);
    pub const AR: Element = Element(18);
    pub const K: Element = Element(19);
    pub const CA: Element = Element(20);
    pub const SC: Element = Element(21);
    pub const TI: Element = Element(22);
    pub const V: Element = Element(23);
    pub const CR: Element = Element(24);
    pub const MN: Element = Element(25);
    pub const FE: Element = Element(26);
    pub const CO: Element = Element(27);
    pub const NI: Element = Element(28);
    pub const CU: Element = Element(29);
    pub const ZN: Element = Element(30);
    pub const GA: Element = Element(31);
    pub const GE: Element = Element(32);
    pub const AS: Element = Element(33);
    pub const SE: Element = Element(34);
    pub const BR: Element = Element(35);
    pub const KR: Element = Element(36);
    pub const RB: Element = Element(37);
    pub const SR: Element = Element(38);
    pub const Y: Element = Element(39);
    pub const ZR: Element = Element(40);
    pub const NB: Element = Element(41);
    pub const MO: Element = Element(42);
    pub const TC: Element = Element(43);
    pub const RU: Element = Element(44);
    pub const RH: Element = Element(45);
    pub const PD: Element = Element(46);
    pub const AG: Element = Element(47);
    pub const CD: Element = Element(48);
    pub const IN: Element = Element(49);
    pub const SN: Element = Element(50);
    pub const SB: Element = Element(51);
    pub const TE: Element = Element(52);
    pub const I: Element = Element(53);
    pub const XE: Element = Element(54);
    pub const CS: Element = Element(55);
    pub const BA: Element = Element(56);
    pub const LA: Element = Element(57);
    pub const CE: Element = Element(58);
    pub const PR: Element = Element(59);
    pub const ND: Element = Element(60);
    pub const PM: Element = Element(61);
    pub const SM: Element = Element(62);
    pub const EU: Element = Element(63);
    pub const GD: Element = Element(64);
    pub const TB: Element = Element(65);
    pub const DY: Element = Element(66);
    pub const HO: Element = Element(67);
    pub const ER: Element = Element(68);
    pub const TM: Element = Element(69);
    pub const YB: Element = Element(70);
    pub const LU: Element = Element(71);
    pub const HF: Element = Element(72);
    pub const TA: Element = Element(73);
    pub const W: Element = Element(74);
    pub const RE: Element = Element(75);
    pub const OS: Element = Element(76);
    pub const IR: Element = Element(77);
    pub const PT: Element = Element(78);
    pub const AU: Element = Element(79);
    pub const HG: Element = Element(80);
    pub const TL: Element = Element(81);
    pub const PB: Element = Element(82);
    pub const BI: Element = Element(83);
    pub const PO: Element = Element(84);
    pub const AT: Element = Element(85);
    pub const RN: Element = Element(86);
    pub const FR: Element = Element(87);
    pub const RA: Element = Element(88);
    pub const AC: Element = Element(89);
    pub const TH: Element = Element(90);
    pub const PA: Element = Element(91);
    pub const U: Element = Element(92);
    pub const NP: Element = Element(93);
    pub const PU: Element = Element(94);
    pub const AM: Element = Element(95);
    pub const CM: Element = Element(96);
    pub const BK: Element = Element(97);
    pub const CF: Element = Element(98);
    pub const ES: Element = Element(99);
    pub const FM: Element = Element(100);
    pub const MD: Element = Element(101);
    pub const NO: Element = Element(102);
    pub const LR: Element = Element(103);
    pub const RF: Element = Element(104);
    pub const DB: Element = Element(105);
    pub const SG: Element = Element(106);
    pub const BH: Element = Element(107);
    pub const HS: Element = Element(108);
    pub const MT: Element = Element(109);
    pub const DS: Element = Element(110);
    pub const RG: Element = Element(111);
    pub const CN: Element = Element(112);
    pub const NH: Element = Element(113);
    pub const FL: Element = Element(114);
    pub const MC: Element = Element(115);
    pub const LV: Element = Element(116);
    pub const TS: Element = Element(117);
    pub const OG: Element = Element(118);

    // -- constructors -------------------------------------------------------

    /// Create an Element from an atomic number. Returns `None` if `n` is outside 1–118.
    #[inline]
    pub const fn from_atomic_number(n: u8) -> Option<Self> {
        if n >= 1 && n <= 118 {
            Some(Self(n))
        } else {
            None
        }
    }

    /// Create an Element from its symbol (title-case, e.g. "C", "Cl").
    /// Case-sensitive: "C" is carbon, "c" is not accepted.
    pub fn from_symbol(s: &str) -> Option<Self> {
        match s {
            "H" => Some(Self::H),
            "He" => Some(Self::HE),
            "Li" => Some(Self::LI),
            "Be" => Some(Self::BE),
            "B" => Some(Self::B),
            "C" => Some(Self::C),
            "N" => Some(Self::N),
            "O" => Some(Self::O),
            "F" => Some(Self::F),
            "Ne" => Some(Self::NE),
            "Na" => Some(Self::NA),
            "Mg" => Some(Self::MG),
            "Al" => Some(Self::AL),
            "Si" => Some(Self::SI),
            "P" => Some(Self::P),
            "S" => Some(Self::S),
            "Cl" => Some(Self::CL),
            "Ar" => Some(Self::AR),
            "K" => Some(Self::K),
            "Ca" => Some(Self::CA),
            "Sc" => Some(Self::SC),
            "Ti" => Some(Self::TI),
            "V" => Some(Self::V),
            "Cr" => Some(Self::CR),
            "Mn" => Some(Self::MN),
            "Fe" => Some(Self::FE),
            "Co" => Some(Self::CO),
            "Ni" => Some(Self::NI),
            "Cu" => Some(Self::CU),
            "Zn" => Some(Self::ZN),
            "Ga" => Some(Self::GA),
            "Ge" => Some(Self::GE),
            "As" => Some(Self::AS),
            "Se" => Some(Self::SE),
            "Br" => Some(Self::BR),
            "Kr" => Some(Self::KR),
            "Rb" => Some(Self::RB),
            "Sr" => Some(Self::SR),
            "Y" => Some(Self::Y),
            "Zr" => Some(Self::ZR),
            "Nb" => Some(Self::NB),
            "Mo" => Some(Self::MO),
            "Tc" => Some(Self::TC),
            "Ru" => Some(Self::RU),
            "Rh" => Some(Self::RH),
            "Pd" => Some(Self::PD),
            "Ag" => Some(Self::AG),
            "Cd" => Some(Self::CD),
            "In" => Some(Self::IN),
            "Sn" => Some(Self::SN),
            "Sb" => Some(Self::SB),
            "Te" => Some(Self::TE),
            "I" => Some(Self::I),
            "Xe" => Some(Self::XE),
            "Cs" => Some(Self::CS),
            "Ba" => Some(Self::BA),
            "La" => Some(Self::LA),
            "Ce" => Some(Self::CE),
            "Pr" => Some(Self::PR),
            "Nd" => Some(Self::ND),
            "Pm" => Some(Self::PM),
            "Sm" => Some(Self::SM),
            "Eu" => Some(Self::EU),
            "Gd" => Some(Self::GD),
            "Tb" => Some(Self::TB),
            "Dy" => Some(Self::DY),
            "Ho" => Some(Self::HO),
            "Er" => Some(Self::ER),
            "Tm" => Some(Self::TM),
            "Yb" => Some(Self::YB),
            "Lu" => Some(Self::LU),
            "Hf" => Some(Self::HF),
            "Ta" => Some(Self::TA),
            "W" => Some(Self::W),
            "Re" => Some(Self::RE),
            "Os" => Some(Self::OS),
            "Ir" => Some(Self::IR),
            "Pt" => Some(Self::PT),
            "Au" => Some(Self::AU),
            "Hg" => Some(Self::HG),
            "Tl" => Some(Self::TL),
            "Pb" => Some(Self::PB),
            "Bi" => Some(Self::BI),
            "Po" => Some(Self::PO),
            "At" => Some(Self::AT),
            "Rn" => Some(Self::RN),
            "Fr" => Some(Self::FR),
            "Ra" => Some(Self::RA),
            "Ac" => Some(Self::AC),
            "Th" => Some(Self::TH),
            "Pa" => Some(Self::PA),
            "U" => Some(Self::U),
            "Np" => Some(Self::NP),
            "Pu" => Some(Self::PU),
            "Am" => Some(Self::AM),
            "Cm" => Some(Self::CM),
            "Bk" => Some(Self::BK),
            "Cf" => Some(Self::CF),
            "Es" => Some(Self::ES),
            "Fm" => Some(Self::FM),
            "Md" => Some(Self::MD),
            "No" => Some(Self::NO),
            "Lr" => Some(Self::LR),
            "Rf" => Some(Self::RF),
            "Db" => Some(Self::DB),
            "Sg" => Some(Self::SG),
            "Bh" => Some(Self::BH),
            "Hs" => Some(Self::HS),
            "Mt" => Some(Self::MT),
            "Ds" => Some(Self::DS),
            "Rg" => Some(Self::RG),
            "Cn" => Some(Self::CN),
            "Nh" => Some(Self::NH),
            "Fl" => Some(Self::FL),
            "Mc" => Some(Self::MC),
            "Lv" => Some(Self::LV),
            "Ts" => Some(Self::TS),
            "Og" => Some(Self::OG),
            _ => None,
        }
    }

    /// Return the IUPAC element symbol (title-case).
    #[inline]
    pub fn symbol(self) -> &'static str {
        SYMBOLS[(self.0 as usize) - 1]
    }

    /// Return the atomic number (1–118).
    #[inline]
    pub fn atomic_number(self) -> u8 {
        self.0
    }

    /// Returns true if this element is in the OpenSMILES organic subset
    /// (B, C, N, O, P, S, F, Cl, Br, I) — these may carry implicit H without brackets.
    #[inline]
    pub fn is_organic_subset(self) -> bool {
        matches!(self.0, 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53)
    }

    /// van der Waals radius in Å (Bondi 1964 + Alvarez 2013).
    /// Unknown/synthetic elements return 1.70 (carbon fallback).
    #[inline]
    pub fn vdw_radius(self) -> f32 {
        VDW_RADII[(self.0 as usize) - 1]
    }

    /// Single-bond covalent radius in Å (Alvarez 2008).
    /// Unknown elements return 0.77 (sp3 carbon fallback).
    #[inline]
    pub fn covalent_radius(self) -> f32 {
        COVALENT_RADII[(self.0 as usize) - 1]
    }

    /// Monoisotopic (exact) mass of the most abundant isotope (u).
    /// Used as CIP rule 4 tiebreaker: higher atomic mass has higher priority.
    pub fn atomic_mass(self) -> f64 {
        ATOMIC_MASSES[(self.0 as usize) - 1]
    }

    /// Normal valence list (ascending). Empty = undefined (transition metals, etc.).
    /// Used for computing implicit H counts in organic-subset atoms.
    pub fn normal_valences(self) -> &'static [u8] {
        match self.0 {
            1 => &[1],           // H
            5 => &[3],           // B
            6 => &[4],           // C
            7 => &[3, 5],        // N
            8 => &[2],           // O
            9 => &[1],           // F
            14 => &[4],          // Si
            15 => &[3, 5],       // P
            16 => &[2, 4, 6],    // S
            17 => &[1, 3, 5, 7], // Cl
            33 => &[3, 5],       // As
            34 => &[2, 4, 6],    // Se
            35 => &[1, 3, 5, 7], // Br
            // Te: source-verified against RDKit atomic_data.cpp @ pinned commit
            // 8afba32ec539dcb2369bc84549d802aca3f7eb39 ("52  Te  ...  2  4  6" in the
            // trailing "list of allowed valences" columns) AND cross-checked against
            // installed RDKit 2026.03.3: GetPeriodicTable().GetValenceList(52) == [2, 4, 6],
            // GetDefaultValence(52) == 2. Matches S/Se's list by coincidence, not analogy.
            52 => &[2, 4, 6],    // Te
            53 => &[1, 3, 5, 7], // I
            _ => &[],
        }
    }
}

impl core::fmt::Display for Element {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.symbol())
    }
}

/// Element symbol table indexed by (atomic_number - 1), covering elements 1–118.
static SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

/// van der Waals radii (Å) indexed by (atomic_number - 1).
/// Sources: Bondi J.Phys.Chem. 1964; Alvarez Dalton Trans. 2013, 10.1039/c3dt50599e.
/// Elements without published values use 1.70 (C fallback).
#[rustfmt::skip]
static VDW_RADII: [f32; 118] = [
    // H     He    Li    Be    B     C     N     O     F     Ne
    1.20, 1.40, 1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
    // Na    Mg    Al    Si    P     S     Cl    Ar    K     Ca
    2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88, 2.75, 2.31,
    // Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
    2.11, 2.00, 1.92, 1.93, 1.97, 1.96, 2.00, 1.63, 1.40, 1.39,
    // Ga    Ge    As    Se    Br    Kr    Rb    Sr    Y     Zr
    1.87, 2.11, 1.85, 1.90, 1.85, 2.02, 3.03, 2.49, 2.34, 2.23,
    // Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd    In    Sn
    2.18, 2.17, 2.16, 2.13, 2.10, 1.63, 1.72, 1.58, 1.93, 2.17,
    // Sb    Te    I     Xe    Cs    Ba    La    Ce    Pr    Nd
    2.06, 2.06, 1.98, 2.16, 3.43, 2.68, 2.57, 2.56, 2.54, 2.51,
    // Pm    Sm    Eu    Gd    Tb    Dy    Ho    Er    Tm    Yb
    2.50, 2.48, 2.45, 2.43, 2.42, 2.40, 2.38, 2.37, 2.36, 2.35,
    // Lu    Hf    Ta    W     Re    Os    Ir    Pt    Au    Hg
    2.34, 2.23, 2.22, 2.18, 2.16, 2.16, 2.13, 1.75, 1.66, 1.55,
    // Tl    Pb    Bi    Po    At    Rn    Fr    Ra    Ac    Th
    1.96, 2.02, 2.07, 1.97, 2.02, 2.20, 3.48, 2.83, 2.64, 2.54,
    // Pa    U     Np    Pu    Am    Cm    Bk    Cf    Es    Fm
    2.52, 2.41, 2.39, 2.43, 2.44, 2.45, 2.44, 2.45, 2.45, 1.70,
    // Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
    1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70,
    // Rg    Cn    Nh    Fl    Mc    Lv    Ts    Og
    1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70,
];

/// Single-bond covalent radii (Å) indexed by (atomic_number - 1).
/// Source: Alvarez et al. Dalton Trans. 2008, 2832–2838, 10.1039/b801115j.
/// Elements 104–118 use computational estimates; others default to 0.77.
#[rustfmt::skip]
static COVALENT_RADII: [f32; 118] = [
    // H     He    Li    Be    B     C     N     O     F     Ne
    0.31, 0.28, 1.28, 0.96, 0.84, 0.77, 0.71, 0.66, 0.57, 0.58,
    // Na    Mg    Al    Si    P     S     Cl    Ar    K     Ca
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76,
    // Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
    1.70, 1.60, 1.53, 1.39, 1.50, 1.52, 1.38, 1.24, 1.32, 1.22,
    // Ga    Ge    As    Se    Br    Kr    Rb    Sr    Y     Zr
    1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75,
    // Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd    In    Sn
    1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39,
    // Sb    Te    I     Xe    Cs    Ba    La    Ce    Pr    Nd
    1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
    // Pm    Sm    Eu    Gd    Tb    Dy    Ho    Er    Tm    Yb
    1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87,
    // Lu    Hf    Ta    W     Re    Os    Ir    Pt    Au    Hg
    1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
    // Tl    Pb    Bi    Po    At    Rn    Fr    Ra    Ac    Th
    1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06,
    // Pa    U     Np    Pu    Am    Cm    Bk    Cf    Es    Fm
    2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.68, 1.68, 1.65, 1.67,
    // Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
    1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28,
    // Rg    Cn    Nh    Fl    Mc    Lv    Ts    Og
    1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57,
];

/// Monoisotopic atomic masses (u) indexed by (atomic_number - 1).
/// Source: NIST, exact masses of most abundant isotopes, used for CIP rule 4 (atomic mass tiebreaker).
#[rustfmt::skip]
static ATOMIC_MASSES: [f64; 118] = [
    // H        He       Li       Be       B        C        N        O        F        Ne
    1.007825, 4.002603, 7.016003, 9.012182, 11.009305, 12.000000, 14.003074, 15.994915, 18.998403, 19.992440,
    // Na       Mg       Al       Si       P        S        Cl       Ar       K        Ca
    22.989770, 23.985042, 26.981538, 27.976927, 30.973762, 31.972071, 34.968853, 39.962383, 38.963707, 39.962591,
    // Sc       Ti       V        Cr       Mn       Fe       Co       Ni       Cu       Zn
    44.955910, 47.947947, 50.943963, 51.940512, 54.938045, 55.934840, 58.933195, 57.935343, 62.929598, 63.929145,
    // Ga       Ge       As       Se       Br       Kr       Rb       Sr       Y        Zr
    68.925581, 73.921218, 74.921596, 79.916521, 78.918337, 83.911507, 84.911789, 87.905614, 88.905848, 89.904704,
    // Nb       Mo       Tc       Ru       Rh       Pd       Ag       Cd       In       Sn
    92.906378, 97.905408, 98.906255, 101.904435, 102.905504, 105.903483, 106.905093, 113.903358, 114.903878, 119.902197,
    // Sb       Te       I        Xe       Cs       Ba       La       Ce       Pr       Nd
    120.903818, 129.906224, 126.904477, 131.904154, 132.905447, 137.905241, 138.906349, 139.905434, 140.907648, 141.907720,
    // Pm       Sm       Eu       Gd       Tb       Dy       Ho       Er       Tm       Yb
    145.000000, 151.919729, 152.921227, 157.924103, 158.925343, 163.929171, 164.930319, 165.930292, 168.934211, 173.938858,
    // Lu       Hf       Ta       W        Re       Os       Ir       Pt       Au       Hg
    174.940768, 179.946549, 180.947885, 183.950931, 186.955751, 191.961479, 192.960837, 194.964791, 196.966569, 201.972326,
    // Tl       Pb       Bi       Po       At       Rn       Fr       Ra       Ac       Th
    202.972329, 207.982019, 208.980384, 209.000000, 210.000000, 222.000000, 223.000000, 226.000000, 227.000000, 232.038054,
    // Pa       U        Np       Pu       Am       Cm       Bk       Cf       Es       Fm
    231.035884, 238.028915, 237.000000, 244.000000, 243.000000, 247.000000, 247.000000, 251.000000, 252.000000, 257.000000,
    // Md       No       Lr       Rf       Db       Sg       Bh       Hs       Mt       Ds
    258.000000, 259.000000, 266.000000, 267.000000, 268.000000, 269.000000, 270.000000, 270.000000, 278.000000, 281.000000,
    // Rg       Cn       Nh       Fl       Mc       Lv       Ts       Og
    280.000000, 285.000000, 286.000000, 289.000000, 290.000000, 293.000000, 294.000000, 294.000000,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_symbol() {
        for n in 1u8..=118 {
            let elem = Element::from_atomic_number(n).unwrap();
            let sym = elem.symbol();
            let back =
                Element::from_symbol(sym).unwrap_or_else(|| panic!("no elem for symbol {sym}"));
            assert_eq!(elem, back, "roundtrip failed for atomic number {n}");
        }
    }

    #[test]
    fn test_organic_subset() {
        for sym in &["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"] {
            assert!(
                Element::from_symbol(sym).unwrap().is_organic_subset(),
                "{sym} should be in organic subset"
            );
        }
        assert!(!Element::H.is_organic_subset());
        assert!(!Element::FE.is_organic_subset());
    }

    #[test]
    fn test_vdw_radii() {
        assert!((Element::H.vdw_radius() - 1.20).abs() < 1e-4);
        assert!((Element::C.vdw_radius() - 1.70).abs() < 1e-4);
        assert!((Element::N.vdw_radius() - 1.55).abs() < 1e-4);
        assert!((Element::O.vdw_radius() - 1.52).abs() < 1e-4);
        assert!((Element::S.vdw_radius() - 1.80).abs() < 1e-4);
        assert!((Element::CL.vdw_radius() - 1.75).abs() < 1e-4);
        assert!((Element::BR.vdw_radius() - 1.85).abs() < 1e-4);
        assert!((Element::I.vdw_radius() - 1.98).abs() < 1e-4);
        // Heavy elements fall back to 1.70
        assert!((Element::OG.vdw_radius() - 1.70).abs() < 1e-4);
    }

    #[test]
    fn test_covalent_radii() {
        assert!((Element::H.covalent_radius() - 0.31).abs() < 1e-4);
        assert!((Element::C.covalent_radius() - 0.77).abs() < 1e-4);
        assert!((Element::N.covalent_radius() - 0.71).abs() < 1e-4);
        assert!((Element::O.covalent_radius() - 0.66).abs() < 1e-4);
        assert!((Element::CL.covalent_radius() - 1.02).abs() < 1e-4);
    }

    #[test]
    fn test_valences() {
        assert_eq!(Element::C.normal_valences(), &[4]);
        assert_eq!(Element::N.normal_valences(), &[3, 5]);
        assert_eq!(Element::S.normal_valences(), &[2, 4, 6]);
    }

    #[test]
    fn test_te_valence_source_verified() {
        // Source-verified (not assumed from S/Se by analogy):
        //  - RDKit atomic_data.cpp @ pinned commit 8afba32ec539dcb2369bc84549d802aca3f7eb39,
        //    row "52  Te  5  1.38  1.378  2.1  127.6  6  130  129.9062244  2  4  6"
        //    (trailing columns = "list of allowed valences" per the file's own header comment).
        //  - Cross-checked against installed RDKit 2026.03.3:
        //    GetPeriodicTable().GetValenceList(52) == [2, 4, 6], GetDefaultValence(52) == 2.
        // It happens to equal S/Se's list, but that's a verified coincidence, not the basis
        // for this entry.
        assert_eq!(Element::TE.normal_valences(), &[2, 4, 6]);
    }

    #[test]
    fn test_s_se_valences_unchanged_regression_guard() {
        // Adding the Te arm must not disturb the pre-existing S/Se entries.
        assert_eq!(Element::S.normal_valences(), &[2, 4, 6]);
        assert_eq!(Element::SE.normal_valences(), &[2, 4, 6]);
    }

    #[test]
    fn test_te_not_in_organic_subset() {
        // Te is outside the OpenSMILES organic subset (only B,C,N,O,P,S,F,Cl,Br,I are),
        // so it always requires bracket notation and implicit_hcount() short-circuits
        // to 0 for it regardless of normal_valences(). This must stay true.
        assert!(!Element::TE.is_organic_subset());
    }
}
