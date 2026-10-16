// ---------------------------------------------------------------------------
// Element Counts — specific element frequencies
// ---------------------------------------------------------------------------

fn count_element(mol: &Molecule, atomic_num: u8) -> usize {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() == atomic_num)
        .count()
}

/// Count carbons (C atoms, including aromatic).
pub fn num_carbons(mol: &Molecule) -> usize {
    count_element(mol, 6)
}

/// Count nitrogens (N atoms, including aromatic).
pub fn num_nitrogens(mol: &Molecule) -> usize {
    count_element(mol, 7)
}

/// Count oxygens (O atoms).
pub fn num_oxygens(mol: &Molecule) -> usize {
    count_element(mol, 8)
}

/// Count fluorines (F atoms).
pub fn num_fluorines(mol: &Molecule) -> usize {
    count_element(mol, 9)
}

/// Count chlorines (Cl atoms).
pub fn num_chlorines(mol: &Molecule) -> usize {
    count_element(mol, 17)
}

/// Count bromines (Br atoms).
pub fn num_bromines(mol: &Molecule) -> usize {
    count_element(mol, 35)
}

/// Count iodines (I atoms).
pub fn num_iodines(mol: &Molecule) -> usize {
    count_element(mol, 53)
}

/// Count sulfurs (S atoms).
pub fn num_sulfurs(mol: &Molecule) -> usize {
    count_element(mol, 16)
}

/// Count phosphorus (P atoms).
pub fn num_phosphorus(mol: &Molecule) -> usize {
    count_element(mol, 15)
}

/// Total hydrogen count (explicit + implicit).
///
/// Counts each real hydrogen atom node (including isotopes/ions) once, plus
/// attached bracket or valence-inferred hydrogens. `implicit_hcount` already
/// returns stored bracket H counts; adding them again would double-count.
/// Wildcards contribute only their explicitly stored attached H count.
pub fn num_hydrogens(mol: &Molecule) -> usize {
    mol.atoms()
        .map(|(idx, atom)| {
            if atom.wildcard {
                usize::from(atom.hydrogen_count.unwrap_or(0))
            } else {
                usize::from(atom.element == Element::H) + usize::from(implicit_hcount(mol, idx))
            }
        })
        .sum()
}

// ---------------------------------------------------------------------------
// Functional Group Bond Counts (granular classification)
// ---------------------------------------------------------------------------

/// Number of amide bonds, RDKit's `CalcNumAmideBonds` (unique matches of
/// its `C(=[O;!R])N` pattern).
pub fn num_amide_bonds(mol: &Molecule) -> usize {
    static QUERY: std::sync::OnceLock<chematic_smarts::QueryMolecule> = std::sync::OnceLock::new();
    count_rdkit_pattern(&QUERY, "C(=[O;!R])N", mol)
}

/// Count ester C(=O)-O bonds in the molecule.
///
/// Identifies ester linkages: carbonyl C bonded to O (via C-O-) where the
/// oxygen is bonded to a carbon (R-O-C=O, not H-O-C=O which is carboxylic acid).
pub fn num_ester_bonds(mol: &Molecule) -> usize {
    let mut count = 0;
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        // Check if this carbon is part of a carbonyl (C=O)
        let has_carbonyl_o = mol.neighbors(idx).any(|(nb, bid)| {
            mol.atom(nb).element.atomic_number() == 8 && mol.bond(bid).order == BondOrder::Double
        });

        if !has_carbonyl_o {
            continue;
        }

        // Check if this carbon is bonded to oxygen (via single bond)
        for (o_idx, bid) in mol.neighbors(idx) {
            let is_oxygen = mol.atom(o_idx).element.atomic_number() == 8;
            let is_single = matches!(
                mol.bond(bid).order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            );
            if !is_oxygen || !is_single {
                continue;
            }

            // Found C(=O)-O. Check if the O is bonded to a carbon (ester) not just H (acid)
            let o_bonded_to_carbon = mol
                .neighbors(o_idx)
                .any(|(nb, _)| nb != idx && mol.atom(nb).element.atomic_number() == 6);

            if o_bonded_to_carbon {
                count += 1;
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Molecular Formula Generation
// ---------------------------------------------------------------------------

/// Generate molecular formula in Hill notation (C first, H second, then alphabetical).
///
/// Example: "C6H12O2" for acetic acid derivative, "H2O" for water (no carbon).
/// As RDKit's `CalcMolFormula`: a net charge is appended (`"C2H3O2-"`,
/// `"O4S-2"`, `"H4N+"`) and dummy atoms are written as `*`.
pub fn calc_mol_formula(mol: &Molecule) -> String {
    chematic_smiles::with_rdkit_model_molecule(mol, calc_mol_formula_impl)
}

fn calc_mol_formula_impl(mol: &Molecule) -> String {
    use std::collections::BTreeMap;

    // Count atoms by element
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();

    let mut charge = 0i32;
    for (_, atom) in mol.atoms() {
        charge += i32::from(atom.charge);
        if atom.wildcard {
            // RDKit's `getMolFormula` keys dummy atoms by the symbol `*`,
            // sorted with the other non-C, non-H symbols (byte order).
            *counts.entry("*".to_string()).or_insert(0) += 1;
            continue;
        }
        let symbol = atom.element.symbol().to_string();
        *counts.entry(symbol).or_insert(0) += 1;
    }

    // Count total implicit hydrogens
    // RDKit's `getTotalNumHs` also counts a dummy atom's hydrogens (`[*H]`).
    let total_h: usize = mol
        .atoms()
        .map(|(idx, a)| {
            if a.wildcard {
                a.hydrogen_count.unwrap_or(0) as usize
            } else {
                implicit_hcount(mol, idx) as usize
            }
        })
        .sum();

    if total_h > 0 {
        *counts.entry("H".to_string()).or_insert(0) += total_h;
    }

    // Build formula in Hill notation: C first, H second, rest alphabetical
    let mut formula = String::new();

    // Carbon
    if let Some(&c_count) = counts.get("C") {
        formula.push('C');
        if c_count > 1 {
            formula.push_str(&c_count.to_string());
        }
    }

    // Hydrogen (always include if present)
    if let Some(&h_count) = counts.get("H") {
        formula.push('H');
        if h_count > 1 {
            formula.push_str(&h_count.to_string());
        }
    }

    // Rest in alphabetical order (excluding C, H)
    for (symbol, &count) in counts.iter() {
        if symbol != "C" && symbol != "H" {
            formula.push_str(symbol);
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }
    }

    // If no atoms at all, return empty
    if formula.is_empty() {
        formula.push_str("H0"); // or empty string, depending on convention
    }

    // Net charge as RDKit's CalcMolFormula writes it: `+`, `-`, `+2`, `-3`.
    match charge {
        0 => {}
        1 => formula.push('+'),
        -1 => formula.push('-'),
        q if q > 0 => formula.push_str(&format!("+{q}")),
        q => formula.push_str(&format!("-{}", -q)),
    }

    formula
}

// ---------------------------------------------------------------------------
// Carbon hybridisation types (Mordred CarbonTypes)
// ---------------------------------------------------------------------------

/// Carbon atom counts grouped by hybridisation and heavy-atom degree.
///
/// Mirrors the Mordred `CarbonTypes` descriptor family: C*x*SP*y* where
/// *x* = number of non-hydrogen (heavy-atom) bonds and *y* = hybridisation
/// (1 = sp, 2 = sp², 3 = sp³).
///
/// Uses `hybridization_per_atom()` internally (1=sp, 2=sp², 3=sp³).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CarbonTypes {
    /// sp with 1 heavy-atom bond (terminal alkyne carbon)
    pub c1sp1: u32,
    /// sp with 2 heavy-atom bonds (internal alkyne carbon)
    pub c2sp1: u32,
    /// sp² with 1 heavy-atom bond (terminal alkene =CH₂)
    pub c1sp2: u32,
    /// sp² with 2 heavy-atom bonds (alkene =CH–)
    pub c2sp2: u32,
    /// sp² with 3 heavy-atom bonds (aromatic C or trisubstituted alkene)
    pub c3sp2: u32,
    /// sp³ with 1 heavy-atom bond (methyl group –CH₃)
    pub c1sp3: u32,
    /// sp³ with 2 heavy-atom bonds (methylene –CH₂–)
    pub c2sp3: u32,
    /// sp³ with 3 heavy-atom bonds (methine –CH<)
    pub c3sp3: u32,
}

/// Count carbon atoms by hybridisation × heavy-atom degree (Mordred CarbonTypes).
pub fn carbon_types(mol: &Molecule) -> CarbonTypes {
    let hyb = hybridization_per_atom(mol);
    let mut ct = CarbonTypes::default();
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue; // only carbon
        }
        let h = hyb[idx.0 as usize];
        let deg: u32 = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .count() as u32;
        match (h, deg) {
            (1, 1) => ct.c1sp1 += 1,
            (1, 2) => ct.c2sp1 += 1,
            (2, 1) => ct.c1sp2 += 1,
            (2, 2) => ct.c2sp2 += 1,
            (2, 3) => ct.c3sp2 += 1,
            (3, 1) => ct.c1sp3 += 1,
            (3, 2) => ct.c2sp3 += 1,
            (3, 3) => ct.c3sp3 += 1,
            _ => {}
        }
    }
    ct
}
