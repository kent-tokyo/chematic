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
            "H"  => Some(Self::H),
            "He" => Some(Self::HE),
            "Li" => Some(Self::LI),
            "Be" => Some(Self::BE),
            "B"  => Some(Self::B),
            "C"  => Some(Self::C),
            "N"  => Some(Self::N),
            "O"  => Some(Self::O),
            "F"  => Some(Self::F),
            "Ne" => Some(Self::NE),
            "Na" => Some(Self::NA),
            "Mg" => Some(Self::MG),
            "Al" => Some(Self::AL),
            "Si" => Some(Self::SI),
            "P"  => Some(Self::P),
            "S"  => Some(Self::S),
            "Cl" => Some(Self::CL),
            "Ar" => Some(Self::AR),
            "K"  => Some(Self::K),
            "Ca" => Some(Self::CA),
            "Sc" => Some(Self::SC),
            "Ti" => Some(Self::TI),
            "V"  => Some(Self::V),
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
            "Y"  => Some(Self::Y),
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
            "I"  => Some(Self::I),
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
            "W"  => Some(Self::W),
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
            "U"  => Some(Self::U),
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

    /// Normal valence list (ascending). Empty = undefined (transition metals, etc.).
    /// Used for computing implicit H counts in organic-subset atoms.
    pub fn normal_valences(self) -> &'static [u8] {
        match self.0 {
            1  => &[1],           // H
            5  => &[3],           // B
            6  => &[4],           // C
            7  => &[3, 5],        // N
            8  => &[2],           // O
            9  => &[1],           // F
            14 => &[4],           // Si
            15 => &[3, 5],        // P
            16 => &[2, 4, 6],     // S
            17 => &[1, 3, 5, 7],  // Cl
            33 => &[3, 5],        // As
            34 => &[2, 4, 6],     // Se
            35 => &[1, 3, 5, 7],  // Br
            53 => &[1, 3, 5, 7],  // I
            _  => &[],
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
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_symbol() {
        for n in 1u8..=118 {
            let elem = Element::from_atomic_number(n).unwrap();
            let sym = elem.symbol();
            let back = Element::from_symbol(sym).unwrap_or_else(|| panic!("no elem for symbol {sym}"));
            assert_eq!(elem, back, "roundtrip failed for atomic number {n}");
        }
    }

    #[test]
    fn test_organic_subset() {
        for sym in &["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"] {
            assert!(Element::from_symbol(sym).unwrap().is_organic_subset(), "{sym} should be in organic subset");
        }
        assert!(!Element::H.is_organic_subset());
        assert!(!Element::FE.is_organic_subset());
    }

    #[test]
    fn test_valences() {
        assert_eq!(Element::C.normal_valences(), &[4]);
        assert_eq!(Element::N.normal_valences(), &[3, 5]);
        assert_eq!(Element::S.normal_valences(), &[2, 4, 6]);
    }
}
