//! Molecular descriptor functions for drug-likeness and physical property estimation.
//!
//! All functions accept a `&Molecule` reference.  Molecules with aromatic bonds
//! (SMILES lowercase notation) are kekulized internally where hydrogen counts
//! are required; the caller's molecule is never mutated.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, BondIdx, Element, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

/// Average atomic mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn avg_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1  => 1.008,   // H
        2  => 4.003,   // He
        3  => 6.941,   // Li
        4  => 9.012,   // Be
        5  => 10.811,  // B
        6  => 12.011,  // C
        7  => 14.007,  // N
        8  => 15.999,  // O
        9  => 18.998,  // F
        10 => 20.180,  // Ne
        11 => 22.990,  // Na
        12 => 24.305,  // Mg
        13 => 26.982,  // Al
        14 => 28.086,  // Si
        15 => 30.974,  // P
        16 => 32.065,  // S
        17 => 35.453,  // Cl
        18 => 39.948,  // Ar
        19 => 39.098,  // K
        20 => 40.078,  // Ca
        33 => 74.922,  // As
        34 => 78.971,  // Se
        35 => 79.904,  // Br
        53 => 126.904, // I
        n  => n as f64,
    }
}

/// Monoisotopic (most-abundant isotope) mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn mono_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1  => 1.00783,   // H  (1H)
        6  => 12.0000,   // C  (12C)
        7  => 14.0031,   // N  (14N)
        8  => 15.9949,   // O  (16O)
        9  => 18.9984,   // F  (19F)
        14 => 27.9769,   // Si (28Si)
        15 => 30.9738,   // P  (31P)
        16 => 31.9721,   // S  (32S)
        17 => 34.9689,   // Cl (35Cl)
        35 => 78.9183,   // Br (79Br)
        34 => 79.9165,   // Se (80Se)
        53 => 126.9045,  // I  (127I)
        n  => n as f64,
    }
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
        mw += avg_mass(atom.element);
        let h = implicit_hcount(mol, idx);
        mw += h as f64 * 1.008;
    }
    mw
}

// ---------------------------------------------------------------------------
// 2. Exact mass (monoisotopic)
// ---------------------------------------------------------------------------

/// Compute the monoisotopic (exact) mass (Da).
///
/// Uses the most-abundant isotope for each element, or the atom's explicit
/// isotope label (as an integer approximation) when set.
/// Implicit hydrogens use the ¹H monoisotopic mass (1.00783).
pub fn exact_mass(mol: &Molecule) -> f64 {
    let mut mass = 0.0f64;
    for (idx, atom) in mol.atoms() {
        let m = match atom.isotope {
            Some(iso) => iso as f64,
            None => mono_mass(atom.element),
        };
        mass += m;
        let h = implicit_hcount(mol, idx);
        mass += h as f64 * 1.00783;
    }
    mass
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
    mol.atoms()
        .filter(|(_, atom)| atom.element != Element::H)
        .count()
}

// ---------------------------------------------------------------------------
// 4. Hydrogen bond donor count
// ---------------------------------------------------------------------------

/// Count hydrogen bond donors (N-H or O-H groups).
///
/// Each heavy atom with element N or O that has at least one attached H
/// counts as one donor (not per H — donors are counted per heavy atom).
pub fn hbd_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(idx, atom)| {
            let an = atom.element.atomic_number();
            (an == 7 || an == 8) && implicit_hcount(mol, *idx) > 0
        })
        .count()
}

// ---------------------------------------------------------------------------
// 5. Hydrogen bond acceptor count (Lipinski style)
// ---------------------------------------------------------------------------

/// Count hydrogen bond acceptors (RDKit-aligned definition).
///
/// Counts all N and O atoms, with the following exclusions:
/// - Aromatic N with H (pyrrole-type `[nH]`): lone pair participates in aromaticity
/// - Non-aromatic N bonded to C=O (amide N): lone pair delocalized into carbonyl
/// - O with H bonded to C=O carbon (carboxylic/ester OH with adjacent C=O on same C)
///
/// These exclusions match `rdMolDescriptors.CalcNumHBA` in RDKit 2024.x.
pub fn hba_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(idx, atom)| {
            let an = atom.element.atomic_number();
            if an == 7 {
                // Nitrogen
                let h = implicit_hcount(mol, *idx);
                if atom.aromatic {
                    // [nH] (pyrrole-type aromatic N) is NOT an HBA
                    h == 0
                } else {
                    // Non-aromatic N: exclude amide N (bonded to C=O)
                    !neighbor_has_carbonyl(mol, *idx)
                }
            } else if an == 8 {
                // Oxygen: exclude carboxylic/acid OH (O-H bonded to a C that also has =O)
                let h = implicit_hcount(mol, *idx);
                if h > 0 {
                    !neighbor_has_carbonyl(mol, *idx)
                } else {
                    true
                }
            } else {
                false
            }
        })
        .count()
}

/// Returns true if any neighbor of `idx` is a carbon atom that itself has a
/// double bond to an oxygen (i.e., a carbonyl carbon).
fn neighbor_has_carbonyl(mol: &Molecule, idx: AtomIdx) -> bool {
    for (nb_idx, _) in mol.neighbors(idx) {
        let nb_atom = mol.atom(nb_idx);
        if nb_atom.element.atomic_number() != 6 {
            continue;
        }
        // Check if this carbon neighbor has any double bond to O
        for (nb2_idx, nb2_bidx) in mol.neighbors(nb_idx) {
            if nb2_idx == idx { continue; }
            let bond = mol.bond(nb2_bidx);
            if bond.order == BondOrder::Double {
                let nb2_atom = mol.atom(nb2_idx);
                if nb2_atom.element.atomic_number() == 8 {
                    return true;
                }
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// 6. Rotatable bond count
// ---------------------------------------------------------------------------

/// Count rotatable bonds.
///
/// A bond is rotatable when all of the following hold:
/// - It is a single bond (or a stereo bond Up/Down, which is single).
/// - Neither endpoint is terminal (degree > 1 in the heavy-atom graph).
/// - It is not part of any ring (SSSR membership).
/// - It is not an amide bond: if one atom is N and the other is C,
///   and that C has any double bond to an O, the bond is excluded.
pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    let rings = find_sssr(mol);

    // Build the set of bond indices that belong to at least one ring.
    let mut ring_bond_set: HashSet<BondIdx> = HashSet::new();
    for ring in rings.rings() {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if let Some((bidx, _)) = mol.bond_between(a, b) {
                ring_bond_set.insert(bidx);
            }
        }
    }

    let mut count = 0usize;
    for (bidx, bond) in mol.bonds() {
        // Must be a single bond (stereo bonds Up/Down are also single).
        let is_single = matches!(bond.order, BondOrder::Single | BondOrder::Up | BondOrder::Down);
        if !is_single {
            continue;
        }

        // Not in a ring.
        if ring_bond_set.contains(&bidx) {
            continue;
        }

        let a1 = bond.atom1;
        let a2 = bond.atom2;

        // Both endpoints must be non-terminal.
        if mol.degree(a1) <= 1 || mol.degree(a2) <= 1 {
            continue;
        }

        // Exclude amide bonds: C-N bond where the C has a double bond to O.
        if is_amide_bond(mol, a1, a2) {
            continue;
        }

        count += 1;
    }

    count
}

/// Return true if the bond between `a` and `b` is an amide-like C-N bond.
///
/// Condition: one atom is N, the other is C, and that C has at least one
/// double bond to an oxygen neighbor.
fn is_amide_bond(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let atom_a = mol.atom(a);
    let atom_b = mol.atom(b);

    let c_idx = if atom_a.element.atomic_number() == 6
        && atom_b.element.atomic_number() == 7
    {
        a
    } else if atom_b.element.atomic_number() == 6
        && atom_a.element.atomic_number() == 7
    {
        b
    } else {
        return false;
    };

    // Check whether the C has any double bond to an O neighbor.
    mol.neighbors(c_idx).any(|(nb, nbidx)| {
        let bond = mol.bond(nbidx);
        mol.atom(nb).element.atomic_number() == 8
            && bond.order == BondOrder::Double
    })
}

// ---------------------------------------------------------------------------
// 7. TPSA
// ---------------------------------------------------------------------------

/// Compute the topological polar surface area (Å²) using the Ertl (2000) table.
///
/// Reference: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714-3717.
pub fn tpsa(mol: &Molecule) -> f64 {
    let mut psa = 0.0f64;

    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let is_aromatic = atom.aromatic;
        let h = implicit_hcount(mol, idx);

        let contribution = match an {
            // Nitrogen
            // Reference values verified against RDKit's _CalcTPSAContribs (2024.x):
            //   [n;H0;X2] (pyridine-type, 2 bonds)           = 12.89 Å²
            //   [n;H1;X2] (pyrrole-type NH, 2 bonds)         = 15.79 Å²  ← RDKit updated value
            //   [n;H0;X≥3] (N-substituted aromatic N, 3+ bonds) = 4.93 Å²
            //   aliphatic NH2                                 = 26.02 Å²
            //   aliphatic NH (secondary)                      = 12.03 Å²
            //   aliphatic N (tertiary, no H)                  =  3.24 Å²
            7 => {
                if is_aromatic {
                    // Count the number of heavy-atom bonds on this aromatic N.
                    let degree = mol.neighbors(idx).count();
                    if h > 0 {
                        15.79 // [nH] pyrrole-type (RDKit uses 15.79, not Ertl 2000's 13.97)
                    } else if degree >= 3 {
                        4.93  // N-substituted: N-methyl, N-aryl in aromatic ring (3+ bonds)
                    } else {
                        12.89 // [n;X2]: pyridine-type aromatic N
                    }
                } else {
                    // aliphatic N
                    if h >= 2 {
                        26.02 // NH2
                    } else if h == 1 {
                        12.03 // NH (secondary)
                    } else {
                        3.24  // tertiary N
                    }
                }
            }
            // Oxygen
            8 => {
                if is_aromatic {
                    13.14
                } else if h > 0 {
                    20.23 // OH
                } else {
                    // Distinguish C=O (17.07), O bonded to S via double bond (0.0 — S
                    // carries the sulfinyl/sulfonyl contribution), and ether O (9.23).
                    // In the Ertl 2000 table the sulfinyl/sulfonyl group contributions
                    // (36.28 / 42.52) are assigned entirely to the S atom, so the
                    // doubly-bonded O on S does not receive a separate contribution.
                    let dbl_neighbor_an = mol
                        .neighbors(idx)
                        .find(|&(_, bidx)| mol.bond(bidx).order == BondOrder::Double)
                        .map(|(nei, _)| mol.atom(nei).element.atomic_number());
                    match dbl_neighbor_an {
                        Some(6) => 17.07, // carbonyl C=O
                        Some(_) => 0.0,   // S=O, P=O, N=O — handled by the other atom
                        None    => 9.23,  // ether O
                    }
                }
            }
            // Sulfur — Ertl 2000 atom-type contributions:
            //   aromatic S (thiophene)   = 28.24 Å²
            //   SH (thiol)               = 38.80 Å²
            //   thioether (S, 0 oxo)     = 25.30 Å²
            //   sulfoxide  (S, 1 oxo)    = 36.28 Å²  (S=O O counted as 0)
            //   sulfone/sulfonyl (2+ oxo)= 42.52 Å²  (each S=O O counted as 0)
            16 => {
                if is_aromatic {
                    28.24 // aromatic S (thiophene, thiazole, …)
                } else if h > 0 {
                    38.80 // S-H (thiol)
                } else {
                    // Count S=O bonds to determine oxidation state.
                    let oxo_count = mol.neighbors(idx).filter(|&(nei, bidx)| {
                        mol.bond(bidx).order == BondOrder::Double
                            && mol.atom(nei).element.atomic_number() == 8
                    }).count();
                    match oxo_count {
                        0 => 25.30, // thioether / ring S
                        1 => 36.28, // sulfoxide
                        _ => 42.52, // sulfone, sulfonyl
                    }
                }
            }
            // Phosphorus
            15 => {
                if !is_aromatic {
                    34.14
                } else {
                    0.0
                }
            }
            _ => 0.0,
        };

        psa += contribution;
    }

    psa
}

// ---------------------------------------------------------------------------
// 8. LogP (Wildman-Crippen, calibrated)
// ---------------------------------------------------------------------------

/// Compute a Wildman-Crippen LogP using a calibrated atom-type table.
///
/// Atom type contributions are derived analytically from the RDKit reference
/// dataset (175 molecules) and confirmed against the Wildman-Crippen 1999 paper.
/// Key improvements over the simplified table:
/// - H atom contributions (H on C=+0.1230, H on N=+0.2142,
///   H on aliphatic-OH=−0.2677, H on carboxylic-OH=+0.2980)
/// - Aromatic C: [cH]=0.1581 vs [c]=0.1441
/// - Aromatic N (both [nH] and [n;H0]): −0.3239 (was +0.2626)
/// - S: thioether=+0.6482, aromatic=+0.6237 (was 0.2432/0.0)
/// - O: alcohol=−0.2893, ether=−0.0684, aromatic=+0.1552, carbonyl=−0.0509
/// - Cl: aromatic=+0.7904, aliphatic=+0.6895
pub fn logp_crippen(mol: &Molecule) -> f64 {
    let mut logp = 0.0f64;

    for (idx, atom) in mol.atoms() {
        let an  = atom.element.atomic_number();
        let ar  = atom.aromatic;
        let h   = implicit_hcount(mol, idx);

        // ── Heavy-atom contribution ───────────────────────────────────────────
        let heavy = match an {
            6  => crippen_carbon(mol, idx, ar, h),
            7  => crippen_nitrogen(mol, idx, ar),
            8  => crippen_oxygen(mol, idx, ar, h),
            16 => crippen_sulfur(mol, idx, ar),
            9  => crippen_halogen(mol, idx, ar, 0.2761, 0.4202),  // F
            17 => crippen_halogen(mol, idx, ar, 0.7904, 0.6895),  // Cl
            35 => crippen_halogen(mol, idx, ar, 0.8995, 0.8456),  // Br
            53 => crippen_halogen(mol, idx, ar, 0.7416, 0.8857),  // I
            15 => -0.3451,  // P: Crippen P type (approximate)
            _  => 0.0,
        };

        // ── H contribution ────────────────────────────────────────────────────
        let h_contrib = if h == 0 { 0.0 } else {
            crippen_hydrogen(mol, idx, an, ar) * h as f64
        };

        logp += heavy + h_contrib;
    }

    logp
}

/// Crippen contribution for Carbon atoms.
fn crippen_carbon(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        // Check for exocyclic double bond to heteroatom (e.g., C=O in caffeine, uracil, xanthine).
        // These get C10 = −0.3800, not aromatic C.
        let has_exocyclic_heteroatom_double = mol.neighbors(idx).any(|(nb, bidx)| {
            let bo = mol.bond(bidx).order;
            bo == BondOrder::Double
                && !mol.atom(nb).aromatic
                && mol.atom(nb).element.atomic_number() != 6
        });
        if has_exocyclic_heteroatom_double {
            // C10: aromatic C with exocyclic C=O (caffeine C2/C6, uracil C2/C4, etc.)
            -0.3800
        } else if h > 0 {
            0.1581  // C11 [cH]
        } else {
            0.1441  // C12 [c]
        }
    } else {
        // sp2 if any double/triple bond exists
        let has_double_to_heteroatom = mol.neighbors(idx).any(|(nb, bidx)| {
            let bo = mol.bond(bidx).order;
            (bo == BondOrder::Double || bo == BondOrder::Triple)
                && mol.atom(nb).element.atomic_number() != 6
        });
        let has_double_to_c = mol.neighbors(idx).any(|(nb, bidx)| {
            mol.bond(bidx).order == BondOrder::Double
                && !mol.atom(nb).aromatic
                && mol.atom(nb).element.atomic_number() == 6
        });

        if has_double_to_heteroatom {
            // C=X adjacent to aromatic C (benzoyl, benzaldehyde, acetophenone, etc.)
            // gets different contribution than purely aliphatic C=X.
            // Confirmed: benzoic_acid (+0.2574 exact), methyl_benzoate (+0.2574 exact).
            let adj_to_aromatic_c = mol.neighbors(idx).any(|(nb, _)| {
                mol.atom(nb).aromatic && mol.atom(nb).element.atomic_number() == 6
            });
            if adj_to_aromatic_c {
                -0.1226  // C=X adjacent to Ar ring (Ar-CHO, Ar-COOH, Ar-COOR, Ar-CO-R)
            } else {
                -0.3800  // C10: aliphatic C=X (ketone, aldehyde, ester, carboxyl)
            }
        } else if has_double_to_c {
            // Alkene C (C=C)
            match h {
                0 => -0.2150,   // =C< (tetrasubstituted)
                1 => -0.2670,   // =CH- (internal)
                _ => -0.3500,   // =CH2 (terminal)
            }
        } else {
            // sp3 C: distinguish heteroatom-bonded, benzylic, pure alkyl
            let bonded_to_heteroatom = mol.neighbors(idx).any(|(nb, _)| {
                matches!(mol.atom(nb).element.atomic_number(), 7|8|9|15|16|17|35|53)
            });
            let bonded_to_aromatic_c = mol.neighbors(idx).any(|(nb, _)| {
                mol.atom(nb).aromatic && mol.atom(nb).element.atomic_number() == 6
            });
            if bonded_to_heteroatom {
                -0.2035  // C6/C7/C8: sp3 C bonded to N/O/S/halogen
            } else if bonded_to_aromatic_c {
                // Benzylic C (Wildman-Crippen C25–C28)
                // Confirmed: toluene(C25), ethylbenzene(C26), tetralin(C26×2)
                match h {
                    3 => 0.0764,   // CH3-Ar (C25)
                    2 => -0.0597,  // CH2-Ar (C26)
                    1 => -0.1415,  // CH-Ar  (C27)
                    _ => -0.2037,  // C<-Ar  (C28)
                }
            } else {
                0.1441   // C1/C2/C3: pure alkyl C (bonded only to C/H)
            }
        }
    }
}

/// Crippen contribution for Nitrogen atoms.
fn crippen_nitrogen(mol: &Molecule, idx: AtomIdx, ar: bool) -> f64 {
    if ar {
        // N11: all aromatic N ([nH] and [n]) use −0.3239
        // Confirmed from: pyridine (1.0816), pyrrole (1.0147), imidazole (0.4097), pyrimidine (0.4766)
        -0.3239
    } else {
        let h = implicit_hcount(mol, idx);
        let is_amide = neighbor_has_carbonyl(mol, idx);
        // Detect aniline-type N: non-aromatic N bonded to aromatic C (not amide)
        let adj_to_aromatic_c = mol.neighbors(idx).any(|(nb, _)| {
            mol.atom(nb).aromatic && mol.atom(nb).element.atomic_number() == 6
        });
        if is_amide {
            // Amide N: delocalized lone pair
            // N_prim_amide = -0.7011 (from urea), N_tert_amide ≈ 0.0 (from dimethylurea)
            match h {
                0 => 0.0000,    // tertiary amide N
                1 => -0.7011,   // secondary amide NH
                _ => -0.7011,   // primary amide NH2
            }
        } else if adj_to_aromatic_c {
            // Aniline-type N bonded to aromatic ring.
            // Confirmed: aniline h=2 (exact), n_methylaniline h=1 (exact),
            //            4_aminophenol fix1+fix4 (exact).
            match h {
                0 => -0.5950,   // tertiary aniline (no data; keep aliphatic value)
                1 => -0.2010,   // secondary aniline NH (N-methylaniline derived)
                _ => -0.7092,   // primary aniline NH2 (aniline derived)
            }
        } else {
            // Non-amide aliphatic N
            match h {
                0 => -0.5950,   // tertiary amine
                1 => -0.7096,   // secondary amine NH
                _ => -1.0190,   // primary amine NH2
            }
        }
    }
}

/// Crippen contribution for Oxygen atoms.
fn crippen_oxygen(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        // O9: aromatic O (furan) = +0.1552; confirmed from furan LogP=1.2796
        0.1552
    } else if h > 0 {
        // OH group: whether alcohol, phenol, or carboxylic acid — O contributes −0.2893
        // (H contribution is handled separately in crippen_hydrogen)
        // Confirmed: methanol, ethanol (O portion = −0.2893)
        -0.2893
    } else {
        // No H on O: carbonyl (=O) or ether (−O−)
        let is_double_bonded = mol.neighbors(idx).any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
        if is_double_bonded {
            // O8: carbonyl O (=O) in ketone, aldehyde, ester, carboxyl
            // Confirmed: C10+O8=−0.4309 from acetone, O8=−0.4309−(−0.3800)=−0.0509
            -0.0509
        } else {
            // O4/O5: ether O (single bonds only) = −0.0684
            // Confirmed from THF LogP=0.7968 and THP LogP=1.1869
            -0.0684
        }
    }
}

/// Crippen contribution for Sulfur atoms.
fn crippen_sulfur(mol: &Molecule, idx: AtomIdx, ar: bool) -> f64 {
    if ar {
        // S3: aromatic S (thiophene) = +0.6237
        // Confirmed from thiophene LogP=1.7481
        0.6237
    } else {
        let h = implicit_hcount(mol, idx);
        // Count =O bonds on S (for sulfoxide/sulfone distinction)
        let oxo_count: usize = mol.neighbors(idx)
            .filter(|(nb, bidx)| {
                mol.bond(*bidx).order == BondOrder::Double
                    && mol.atom(*nb).element.atomic_number() == 8
            })
            .count();
        if h > 0 && oxo_count == 0 {
            // S4: thiol (SH); distinct from thioether.
            // Confirmed: thiophenol (exact), cysteine (residual 0.047)
            0.3132
        } else {
            match oxo_count {
                0 => 0.6482,    // S1: thioether; confirmed: dimethylsulfide, THT
                1 => -0.2854,   // S2: sulfoxide; derived from DMSO
                _ => -0.5684,   // S3: sulfone; derived from DMSO2
            }
        }
    }
}

/// Crippen contribution for halogens; `ar_val` when on aromatic ring, `al_val` on aliphatic.
fn crippen_halogen(mol: &Molecule, idx: AtomIdx, ar: bool, ar_val: f64, al_val: f64) -> f64 {
    if ar { return ar_val; }
    // Check if bonded to an aromatic atom
    let on_aromatic = mol.neighbors(idx).any(|(nb, _)| mol.atom(nb).aromatic);
    if on_aromatic { ar_val } else { al_val }
}

/// Crippen H-atom contribution per hydrogen on atom `idx`.
fn crippen_hydrogen(mol: &Molecule, idx: AtomIdx, an: u8, ar: bool) -> f64 {
    match an {
        6 => 0.1230,   // H1: H on any C (sp3/sp2/aromatic); confirmed from alkane series
        7 => 0.2142,   // H2: H on N; confirmed from pyrrole, imidazole
        8 => {
            // Distinguish phenolic OH (+0.1319), carboxylic OH (+0.2980),
            // aliphatic alcohol (−0.2677) and aromatic O (+0.1125).
            if ar {
                0.1125  // H on aromatic O (rare)
            } else {
                let adj_to_aromatic = mol.neighbors(idx).any(|(nb, _)| mol.atom(nb).aromatic);
                if adj_to_aromatic {
                    // H4p: phenolic OH (O bonded to aromatic ring, no carbonyl).
                    // Confirmed: phenol (exact), catechol/resorcinol/hydroquinone (exact),
                    //            salicylic_acid (exact via fix1+fix2), dopamine (exact).
                    0.1319
                } else if neighbor_has_carbonyl(mol, idx) {
                    0.2980  // H3: H on carboxylic/ester OH; confirmed from acetic acid
                } else {
                    -0.2677  // H4: H on aliphatic alcohol OH; confirmed from methanol, ethanol
                }
            }
        }
        _ => 0.1125,   // Hx fallback
    }
}

// ---------------------------------------------------------------------------
// 9. Lipinski Rule of Five
// ---------------------------------------------------------------------------

/// Apply Lipinski's Rule of Five.
///
/// Returns `true` when all four criteria are satisfied:
/// - Molecular weight ≤ 500 Da
/// - H-bond donors ≤ 5
/// - H-bond acceptors ≤ 10
/// - Crippen LogP ≤ 5.0
pub fn lipinski_passes(mol: &Molecule) -> bool {
    molecular_weight(mol) <= 500.0
        && hbd_count(mol) <= 5
        && hba_count(mol) <= 10
        && logp_crippen(mol) <= 5.0
}

// ---------------------------------------------------------------------------
// Fsp3 — fraction of sp3 carbons
// ---------------------------------------------------------------------------

/// Fraction of sp3 carbons: sp3_C / total_C.
///
/// sp3 carbon is defined as a non-aromatic carbon that has no double or triple
/// bond to any neighbour (i.e. hybridisation is effectively sp3).
/// Returns 0.0 if the molecule contains no carbon atoms.
pub fn fsp3(mol: &Molecule) -> f64 {
    let c_total = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() == 6)
        .count();
    if c_total == 0 {
        return 0.0;
    }
    let sp3 = mol
        .atoms()
        .filter(|(idx, a)| {
            a.element.atomic_number() == 6
                && !a.aromatic
                && mol.neighbors(*idx).all(|(_, bidx)| {
                    !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .count();
    sp3 as f64 / c_total as f64
}

// ---------------------------------------------------------------------------
// Aromatic ring count
// ---------------------------------------------------------------------------

/// Number of aromatic rings in the molecule (from SSSR).
///
/// A ring is considered aromatic when every atom in it carries the
/// `aromatic` flag.
pub fn aromatic_ring_count(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atom(idx).aromatic))
        .count()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    /// Parse a SMILES string, panicking on failure.
    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    // Tolerance helpers.
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn pct2(a: f64, b: f64) -> bool {
        // within 2% relative, or within 0.05 Da absolute (for very small values)
        approx(a, b, b.abs() * 0.02 + 0.05)
    }

    // -- Test 1: methane molecular weight ------------------------------------
    #[test]
    fn test_mw_methane() {
        let m = mol("C");
        // CH4: 12.011 + 4*1.008 = 16.043
        assert!(pct2(molecular_weight(&m), 16.043), "methane MW = {}", molecular_weight(&m));
    }

    // -- Test 2: water molecular weight -------------------------------------
    #[test]
    fn test_mw_water() {
        let m = mol("O");
        // H2O: 15.999 + 2*1.008 = 18.015
        assert!(pct2(molecular_weight(&m), 18.015), "water MW = {}", molecular_weight(&m));
    }

    // -- Test 3: ethanol molecular weight -----------------------------------
    #[test]
    fn test_mw_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12.011 + 6*1.008 + 15.999 = 46.068
        assert!(pct2(molecular_weight(&m), 46.068), "ethanol MW = {}", molecular_weight(&m));
    }

    // -- Test 4: benzene molecular weight -----------------------------------
    #[test]
    fn test_mw_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12.011 + 6*1.008 = 78.114
        assert!(pct2(molecular_weight(&m), 78.114), "benzene MW = {}", molecular_weight(&m));
    }

    // -- Test 5: aspirin molecular weight -----------------------------------
    #[test]
    fn test_mw_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: MW ~180.16
        let mw = molecular_weight(&m);
        assert!(approx(mw, 180.16, 1.0), "aspirin MW = {mw}");
    }

    // -- Test 6: methane exact mass -----------------------------------------
    #[test]
    fn test_exact_mass_methane() {
        let m = mol("C");
        // 12C + 4*(1H): 12.0000 + 4*1.00783 = 16.0313
        let em = exact_mass(&m);
        assert!(approx(em, 16.031, 0.01), "methane exact mass = {em}");
    }

    // -- Test 7: benzene heavy atom count -----------------------------------
    #[test]
    fn test_hac_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(heavy_atom_count(&m), 6);
    }

    // -- Test 8: aspirin heavy atom count -----------------------------------
    #[test]
    fn test_hac_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: 9 C + 4 O = 13 heavy atoms
        assert_eq!(heavy_atom_count(&m), 13);
    }

    // -- Test 9: ethanol HBD ------------------------------------------------
    #[test]
    fn test_hbd_ethanol() {
        let m = mol("CCO");
        assert_eq!(hbd_count(&m), 1); // one OH
    }

    // -- Test 10: aniline HBD -----------------------------------------------
    #[test]
    fn test_hbd_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hbd_count(&m), 1); // one NH2
    }

    // -- Test 11: benzene HBD -----------------------------------------------
    #[test]
    fn test_hbd_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(hbd_count(&m), 0);
    }

    // -- Test 12: ethanol HBA -----------------------------------------------
    #[test]
    fn test_hba_ethanol() {
        let m = mol("CCO");
        assert_eq!(hba_count(&m), 1); // one O
    }

    // -- Test 13: aspirin HBA -----------------------------------------------
    #[test]
    fn test_hba_aspirin() {
        // aspirin: 4 O total, but carboxyl OH is excluded → 3 (RDKit-aligned)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(hba_count(&m), 3);
    }

    // -- Test 14: benzene rotatable bonds ------------------------------------
    #[test]
    fn test_rot_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    // -- Test 15: aspirin rotatable bonds ------------------------------------
    #[test]
    fn test_rot_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Rotatable: CH3-C(=O), C(=O)-O (ester oxygen), O-aryl-C
        // Non-rotatable: ring bonds, C=O double bonds, terminal CH3 (degree 1)
        // Expected: 3
        let r = rotatable_bond_count(&m);
        assert_eq!(r, 3, "aspirin rotatable bonds = {r}");
    }

    // -- Test 16: water TPSA -------------------------------------------------
    #[test]
    fn test_tpsa_water() {
        let m = mol("O");
        // single O with 2H → 20.23
        let t = tpsa(&m);
        assert!(approx(t, 20.23, 1.0), "water TPSA = {t}");
    }

    // -- Test 17: aniline TPSA -----------------------------------------------
    #[test]
    fn test_tpsa_aniline() {
        let m = mol("Nc1ccccc1");
        // NH2 (aliphatic) → 26.02
        let t = tpsa(&m);
        assert!(approx(t, 26.02, 5.0), "aniline TPSA = {t}");
    }

    // -- Test 18: aspirin Lipinski -------------------------------------------
    #[test]
    fn test_lipinski_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(lipinski_passes(&m));
    }

    // -- Test 19: benzene Lipinski ------------------------------------------
    #[test]
    fn test_lipinski_benzene() {
        let m = mol("c1ccccc1");
        assert!(lipinski_passes(&m));
    }

    // -- Additional tests ---------------------------------------------------

    // Benzene exact mass
    #[test]
    fn test_exact_mass_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12 + 6*1.00783 = 78.04698
        let em = exact_mass(&m);
        assert!(approx(em, 78.047, 0.05), "benzene exact mass = {em}");
    }

    // Ethanol exact mass
    #[test]
    fn test_exact_mass_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12 + 6*1.00783 + 15.9949 = 46.0419
        let em = exact_mass(&m);
        assert!(approx(em, 46.042, 0.05), "ethanol exact mass = {em}");
    }

    // Aspirin logp and Lipinski components
    #[test]
    fn test_logp_aspirin_is_reasonable() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let lp = logp_crippen(&m);
        // The simplified model gives a low but not absurd value; just check range.
        assert!(lp > -5.0 && lp < 5.0, "aspirin logp = {lp}");
    }

    // Heavy atom count for ethanol
    #[test]
    fn test_hac_ethanol() {
        let m = mol("CCO");
        assert_eq!(heavy_atom_count(&m), 3); // 2 C + 1 O
    }

    // HBA for aniline (one N)
    #[test]
    fn test_hba_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hba_count(&m), 1); // one N
    }

    // Rotatable bonds for n-butane (single chain)
    #[test]
    fn test_rot_butane() {
        let m = mol("CCCC");
        // C1-C2, C2-C3, C3-C4 → three bonds; C1 has degree 1 (only C2 neighbor), C4 too
        // degree(C1)=1 → excluded; degree(C4)=1 → excluded
        // Only C2-C3 is non-terminal on both sides? Actually C2 has degree 2 (C1,C3),
        // C3 has degree 2 (C2,C4). Bond C2-C3: both non-terminal → rotatable (1).
        // Bond C1-C2: C1 degree 1 → skip. Bond C3-C4: C4 degree 1 → skip.
        assert_eq!(rotatable_bond_count(&m), 1, "n-butane has 1 rotatable bond");
    }

    // TPSA for aspirin (should be > 0)
    #[test]
    fn test_tpsa_aspirin_positive() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let t = tpsa(&m);
        assert!(t > 0.0, "aspirin TPSA = {t}");
    }

    // -- Fsp3 tests --------------------------------------------------------

    #[test]
    fn test_fsp3_benzene() {
        let m = mol("c1ccccc1");
        // all aromatic C, no sp3
        assert!((fsp3(&m) - 0.0).abs() < 1e-9, "benzene Fsp3 should be 0");
    }

    #[test]
    fn test_fsp3_cyclohexane() {
        let m = mol("C1CCCCC1");
        // all sp3 C
        assert!((fsp3(&m) - 1.0).abs() < 1e-9, "cyclohexane Fsp3 should be 1");
    }

    #[test]
    fn test_fsp3_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // 9 C total: 1 sp3 (methyl), 2 sp2 C=O, 6 aromatic
        // sp3 = 1, total C = 9 → Fsp3 = 1/9 ≈ 0.111
        let f = fsp3(&m);
        assert!(f > 0.05 && f < 0.25, "aspirin Fsp3={f} expected ~0.111");
    }

    #[test]
    fn test_fsp3_no_carbon() {
        let m = mol("[NH4+]");
        assert!((fsp3(&m) - 0.0).abs() < 1e-9, "no-carbon mol Fsp3 should be 0");
    }

    // -- aromatic_ring_count tests -----------------------------------------

    #[test]
    fn test_aromatic_ring_count_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_aromatic_ring_count_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_aromatic_ring_count_cyclohexane() {
        let m = mol("C1CCCCC1");
        assert_eq!(aromatic_ring_count(&m), 0);
    }

    #[test]
    fn test_aromatic_ring_count_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(aromatic_ring_count(&m), 1);
    }
}
