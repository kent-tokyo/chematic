use chematic_core::atomic_weight::{
    average_atomic_weight as avg_mass, rdkit_average_atomic_weight as rdkit_avg_mass,
    rdkit_isotope_mass,
};

// --- Element Detection Helpers ---
// Consolidate atomic number matching to eliminate 50+ hardcoded checks throughout the file.

#[inline]
fn is_carbon(an: u8) -> bool {
    an == 6
}

#[inline]
fn is_nitrogen(an: u8) -> bool {
    an == 7
}

#[inline]
fn is_oxygen(an: u8) -> bool {
    an == 8
}

/// Mass of an explicitly labelled isotope: RDKit's `GetMassForIsotope`,
/// with its mass-number fallback for unknown labels.
fn exact_mass_isotope(element: Element, isotope: u16) -> f64 {
    rdkit_isotope_mass(element, isotope)
}

/// Monoisotopic (most-abundant-isotope) mass table (Da), indexed by atomic number minus one: RDKit 2026.03's `GetMostCommonIsotopeMass`
/// for every element, so exact mass is bit-identical to `ExactMolWt`.
/// (The table used to mix more precise nuclide masses for H, N, O with
/// rounded 4-5 decimal values for halogens, Si and metals: up to 4.7e-5 Da
/// off per Cl.) Generated with RDKit 2026.03.1; do not edit by hand.
static MONO_MASS_TABLE: [f64; 118] = [
    1.007825032, // 1 H (1H)
    4.002603254, // 2 He (4He)
    7.01600455,  // 3 Li (7Li)
    9.0121822,   // 4 Be (9Be)
    11.0093054,  // 5 B (11B)
    12.0,        // 6 C (12C)
    14.003074,   // 7 N (14N)
    15.99491462, // 8 O (16O)
    18.99840322, // 9 F (19F)
    19.99244018, // 10 Ne (20Ne)
    22.98976928, // 11 Na (23Na)
    23.9850417,  // 12 Mg (24Mg)
    26.98153863, // 13 Al (27Al)
    27.97692653, // 14 Si (28Si)
    30.97376163, // 15 P (31P)
    31.972071,   // 16 S (32S)
    34.96885268, // 17 Cl (35Cl)
    39.96238312, // 18 Ar (40Ar)
    38.96370668, // 19 K (39K)
    39.96259098, // 20 Ca (40Ca)
    44.9559119,  // 21 Sc (45Sc)
    47.9479463,  // 22 Ti (48Ti)
    50.9439595,  // 23 V (51V)
    51.9405075,  // 24 Cr (52Cr)
    54.9380451,  // 25 Mn (55Mn)
    55.9349375,  // 26 Fe (56Fe)
    58.933195,   // 27 Co (59Co)
    57.9353429,  // 28 Ni (58Ni)
    62.9295975,  // 29 Cu (63Cu)
    63.9291422,  // 30 Zn (64Zn)
    68.9255736,  // 31 Ga (69Ga)
    73.9211778,  // 32 Ge (74Ge)
    74.9215965,  // 33 As (75As)
    79.9165213,  // 34 Se (80Se)
    78.9183371,  // 35 Br (79Br)
    83.911507,   // 36 Kr (84Kr)
    84.91178974, // 37 Rb (85Rb)
    87.9056121,  // 38 Sr (88Sr)
    88.9058483,  // 39 Y (89Y)
    89.9047044,  // 40 Zr (90Zr)
    92.9063781,  // 41 Nb (93Nb)
    97.9054082,  // 42 Mo (98Mo)
    96.906365,   // 43 Tc (97Tc)
    101.9043493, // 44 Ru (102Ru)
    102.905504,  // 45 Rh (103Rh)
    105.903486,  // 46 Pd (106Pd)
    106.905097,  // 47 Ag (107Ag)
    113.9033585, // 48 Cd (114Cd)
    114.903878,  // 49 In (115In)
    119.9021947, // 50 Sn (120Sn)
    120.9038157, // 51 Sb (121Sb)
    129.9062244, // 52 Te (130Te)
    126.904473,  // 53 I (127I)
    131.9041535, // 54 Xe (132Xe)
    132.9054519, // 55 Cs (133Cs)
    137.9052472, // 56 Ba (138Ba)
    138.9063533, // 57 La (139La)
    139.9054387, // 58 Ce (140Ce)
    140.9076528, // 59 Pr (141Pr)
    141.9077233, // 60 Nd (142Nd)
    144.912749,  // 61 Pm (145Pm)
    151.9197324, // 62 Sm (152Sm)
    152.9212303, // 63 Eu (153Eu)
    157.9241039, // 64 Gd (158Gd)
    158.9253468, // 65 Tb (159Tb)
    163.9291748, // 66 Dy (164Dy)
    164.9303221, // 67 Ho (165Ho)
    165.9302931, // 68 Er (166Er)
    168.9342133, // 69 Tm (169Tm)
    173.9388621, // 70 Yb (174Yb)
    174.9407718, // 71 Lu (175Lu)
    179.94655,   // 72 Hf (180Hf)
    180.9479958, // 73 Ta (181Ta)
    183.9509312, // 74 W (184W)
    186.9557531, // 75 Re (187Re)
    191.9614807, // 76 Os (192Os)
    192.9629264, // 77 Ir (193Ir)
    194.9647911, // 78 Pt (195Pt)
    196.9665687, // 79 Au (197Au)
    201.970643,  // 80 Hg (202Hg)
    204.9744275, // 81 Tl (205Tl)
    207.9766521, // 82 Pb (208Pb)
    208.9803987, // 83 Bi (209Bi)
    208.9824304, // 84 Po (209Po)
    209.987148,  // 85 At (210At)
    222.0175706, // 86 Rn (222Rn)
    223.0197359, // 87 Fr (223Fr)
    226.0254026, // 88 Ra (226Ra)
    227.0277521, // 89 Ac (227Ac)
    232.0380553, // 90 Th (232Th)
    231.035884,  // 91 Pa (231Pa)
    238.0507882, // 92 U (238U)
    236.04657,   // 93 Np (236Np)
    238.0495599, // 94 Pu (238Pu)
    241.0568291, // 95 Am (241Am)
    243.0613891, // 96 Cm (243Cm)
    247.070307,  // 97 Bk (247Bk)
    249.0748535, // 98 Cf (249Cf)
    252.08298,   // 99 Es (252Es)
    257.095105,  // 100 Fm (257Fm)
    258.098431,  // 101 Md (258Md)
    259.10103,   // 102 No (259No)
    262.10963,   // 103 Lr (262Lr)
    267.12153,   // 104 Rf (267Rf)
    268.12545,   // 105 Db (268Db)
    271.13347,   // 106 Sg (271Sg)
    270.13362,   // 107 Bh (270Bh)
    269.13406,   // 108 Hs (269Hs)
    278.15481,   // 109 Mt (278Mt)
    281.16206,   // 110 Ds (281Ds)
    281.16537,   // 111 Rg (281Rg)
    285.17411,   // 112 Cn (285Cn)
    284.17873,   // 113 Nh (284Nh)
    289.19042,   // 114 Fl (289Fl)
    288.19274,   // 115 Mc (288Mc)
    293.20449,   // 116 Lv (293Lv)
    292.20746,   // 117 Ts (292Ts)
    294.21392,   // 118 Og (294Og)
];

fn mono_mass(element: Element) -> f64 {
    let an = element.atomic_number();
    MONO_MASS_TABLE
        .get(an as usize - 1)
        .copied()
        .unwrap_or(an as f64)
}

// ---------------------------------------------------------------------------
// 1. Molecular weight
// ---------------------------------------------------------------------------

/// Compute the average molecular weight (Da).
///
/// Sums the average atomic mass of all heavy atoms plus each atom's implicit
/// hydrogen contribution (1.008 Da per H).
pub fn molecular_weight(mol: &Molecule) -> f64 {
    let mut mw = 0.0f64;
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        mw += avg_mass(atom.element);
        let h = implicit_hcount(mol, idx);
        mw += h as f64 * 1.008;
    }
    mw
}

/// Compute molecular weight with the pinned RDKit average-mass convention.
///
/// This is deliberately separate from [`molecular_weight`].  It matches the
/// RDKit periodic-table values for the supported unlabelled-atom profile while
/// retaining chematic's native defaults for existing callers. Explicit
/// isotope-labelled atoms use RDKit's nuclide masses where present and its
/// mass-number fallback otherwise.
pub fn rdkit_molecular_weight(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, rdkit_molecular_weight_impl)
}

fn rdkit_molecular_weight_impl(mol: &Molecule) -> f64 {
    let mut mw = 0.0f64;
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            // RDKit gives a dummy atom zero mass but still includes bracket
            // hydrogens: [*H] and [*H2] weigh 1.008 and 2.016 Da.
            mw += f64::from(atom.hydrogen_count.unwrap_or(0)) * 1.008;
            continue;
        }
        mw += match atom.isotope {
            // In SMILES, isotope zero means unspecified. RDKit therefore
            // uses the element's average atomic weight for `[0CH4]`, just as
            // it does for `C`, rather than treating zero as a nuclide mass.
            None | Some(0) => rdkit_avg_mass(atom.element),
            Some(isotope) => rdkit_isotope_mass(atom.element, isotope),
        };
        mw += implicit_hcount(mol, idx) as f64 * 1.008;
    }
    mw
}

// ---------------------------------------------------------------------------
// 2. Exact mass (monoisotopic)
// ---------------------------------------------------------------------------

/// Compute the monoisotopic (exact) mass (Da).
///
/// Uses the most-abundant isotope for each element, or the explicit nuclide
/// mass when the atom carries a known isotope label. Unknown isotope labels
/// use RDKit's mass-number fallback. Implicit hydrogens use the ¹H monoisotopic
/// mass. Formal charge adjusts the neutral-atom sum by the electron rest mass,
/// matching RDKit's `ExactMolWt` convention.
pub fn exact_mass(mol: &Molecule) -> f64 {
    // RDKit's getExactMolWt, operation for operation: each atom's mass, then
    // its charge as electrons; the hydrogens are added once at the end.
    const ELECTRON_MASS: f64 = 0.000_548_579_91;
    let mut mass = 0.0f64;
    let mut hydrogens = 0u32;
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        mass += match atom.isotope {
            Some(iso) => exact_mass_isotope(atom.element, iso),
            None => mono_mass(atom.element),
        };
        mass -= f64::from(atom.charge) * ELECTRON_MASS;
        hydrogens += u32::from(implicit_hcount(mol, idx));
    }
    mass + f64::from(hydrogens) * mono_mass(Element::H)
}

// ---------------------------------------------------------------------------
// 3. Heavy atom count
// ---------------------------------------------------------------------------

/// Count non-hydrogen heavy atoms.
///
/// Hydrogen atoms are normally implicit in chematic, but some molecules may
/// carry explicit H atoms in the graph (e.g. from bracket notation `[H]`).
/// Those are excluded from the heavy-atom count.
pub fn heavy_atom_count(mol: &Molecule) -> usize {
    // RDKit's `CalcNumHeavyAtoms`: atomic number above 1, so dummy atoms
    // (atomic number 0) are not heavy atoms.
    mol.atoms()
        .filter(|(_, atom)| !atom.wildcard && atom.element != Element::H)
        .count()
}

/// Number of hydrogens attached to one atom for descriptor atom typing.
///
/// `implicit_hcount` already includes bracket H counts (for example `[NH]`),
/// but an explicit graph atom such as `[2H]O` is a separate neighbour and
/// therefore must be added here.  Descriptor definitions treat protium and
/// isotopic hydrogen identically for donor/acceptor and TPSA atom types.
fn descriptor_attached_hcount(mol: &Molecule, idx: AtomIdx) -> u8 {
    let explicit_graph_h = mol
        .neighbors(idx)
        .filter(|(neighbor, _)| mol.atom(*neighbor).element.atomic_number() == 1)
        .count() as u8;
    implicit_hcount(mol, idx).saturating_add(explicit_graph_h)
}

// ---------------------------------------------------------------------------
// 4. Hydrogen bond donor count
// ---------------------------------------------------------------------------

/// RDKit's `CalcNumHBD` pattern (Lipinski.cpp, version 2.0.1).
const RDKIT_HBD_SMARTS: &str = "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]";
/// RDKit's `CalcNumHBA` pattern (Lipinski.cpp, version 2.0.2).
const RDKIT_HBA_SMARTS: &str = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),\
$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]";

/// Unique matches of an RDKit descriptor pattern, as RDKit's
/// `countMatches` counts them.
fn count_rdkit_pattern(
    cell: &'static std::sync::OnceLock<chematic_smarts::QueryMolecule>,
    smarts: &str,
    mol: &Molecule,
) -> usize {
    let query = cell.get_or_init(|| {
        chematic_smarts::parse_smarts(smarts).expect("RDKit descriptor SMARTS is valid")
    });
    chematic_smarts::find_match_atom_sets_perceived(
        query,
        mol,
        &chematic_smarts::MatchConfig::default(),
    )
    .len()
}

/// Count hydrogen-bond donors: RDKit's `CalcNumHBD`, matched with its own
/// SMARTS (N-H on trivalent or protonated N, neutral O-H and S-H, aromatic
/// N-H). The former hand-written rule counted O-H/S-H on charged or
/// metal-bound atoms (7 exposed-10k rows).
pub fn hbd_count(mol: &Molecule) -> usize {
    chematic_smiles::with_rdkit_model_molecule(mol, hbd_count_impl)
}

fn hbd_count_impl(mol: &Molecule) -> usize {
    static QUERY: std::sync::OnceLock<chematic_smarts::QueryMolecule> = std::sync::OnceLock::new();
    count_rdkit_pattern(&QUERY, RDKIT_HBD_SMARTS, mol)
}

// ---------------------------------------------------------------------------
// 5. Hydrogen bond acceptor count (Ertl / RDKit-aligned definition)
// ---------------------------------------------------------------------------

/// Count hydrogen-bond acceptors: RDKit 2026.03's `CalcNumHBA`, matched
/// with its own SMARTS. The former hand-written rules disagreed on 277
/// exposed-10k rows (tetrazoles, enol/amidine tautomers, hypervalent
/// centres).
pub fn hba_count(mol: &Molecule) -> usize {
    static QUERY: std::sync::OnceLock<chematic_smarts::QueryMolecule> = std::sync::OnceLock::new();
    count_rdkit_pattern(&QUERY, RDKIT_HBA_SMARTS, mol)
}

/// Count hydrogen-bond acceptors using the RDKit 2026.03 profile; the same
/// value as [`hba_count`].
pub fn rdkit_hba_count(mol: &Molecule) -> usize {
    chematic_smiles::with_rdkit_model_molecule(mol, rdkit_hba_count_impl)
}

fn rdkit_hba_count_impl(mol: &Molecule) -> usize {
    hba_count(mol)
}

// ---------------------------------------------------------------------------
// 6. Rotatable bond count
// ---------------------------------------------------------------------------

/// RDKit's `NumRotatableBondsOptions::Strict` pattern (Lipinski.cpp),
/// matched with chematic's SMARTS engine. The former hand-written rules
/// disagreed with RDKit on 97 of the exposed 10k rows (oxime ethers,
/// N-acyl hydroxylamines, imine-aryl bonds, disulfides).
///
/// RDKit 2026.03 also excludes `[CH3]` on both ends: a methyl written with
/// explicit isotopic hydrogens (`[2H]C([2H])([2H])N`) has degree 4 but is
/// still not a rotor.
// Test degree and methyl exclusions before recursive multi-atom patterns.
// Every term is still the same conjunction, with the same endpoint order.
const STRICT_ROTATABLE_SMARTS: &str = "[!D1&!$([CH3])&!$(*#*)&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)\
&!$(C([CH3])([CH3])[CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])\
&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]\
-,:;!@[!D1&!$([CH3])&!$(*#*)&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])]";

fn strict_rotatable_query() -> &'static chematic_smarts::QueryMolecule {
    static QUERY: std::sync::OnceLock<chematic_smarts::QueryMolecule> = std::sync::OnceLock::new();
    QUERY.get_or_init(|| {
        chematic_smarts::parse_smarts(STRICT_ROTATABLE_SMARTS)
            .expect("the strict rotatable-bond SMARTS is valid")
    })
}

/// The atom-index pairs of every rotatable bond (RDKit strict definition),
/// sorted, one per bond. Used by [`rotatable_bond_count`] and by
/// `chematic-3d`'s torsion motif extraction so the two crates never define
/// "rotatable" differently.
pub fn rotatable_bond_atom_pairs(mol: &Molecule) -> Vec<(AtomIdx, AtomIdx)> {
    let mut pairs: Vec<(AtomIdx, AtomIdx)> = chematic_smarts::find_match_atom_sets_perceived(
        strict_rotatable_query(),
        mol,
        &chematic_smarts::MatchConfig::default(),
    )
    .into_iter()
    .filter_map(|set| match set.as_slice() {
        &[a, b] => Some((AtomIdx(a.min(b) as u32), AtomIdx(a.max(b) as u32))),
        _ => None,
    })
    .filter(|&(a, b)| {
        // RDKit removes ordinary explicit H atoms on SMILES input. Keep
        // isotope-labelled H, but do not let a removable H make a terminal
        // OH/NH group appear to have two neighbours in the strict SMARTS.
        [a, b].into_iter().all(|end| {
            mol.neighbors(end)
                .filter(|&(nb, _)| !crate::hydrogen::is_removable_explicit_h(mol, nb))
                .count()
                > 1
        })
    })
    .collect();
    pairs.sort_unstable_by_key(|&(a, b)| (a.0, b.0));
    pairs.dedup();
    pairs
}

/// Count rotatable bonds: RDKit's `CalcNumRotatableBonds` with the
/// (default) strict definition.
pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    rotatable_bond_atom_pairs(mol).len()
}
