use std::collections::HashMap;
use std::os::raw::c_char;

use chematic_chem::tetrahedral_stereo_neighbors;
use chematic_core::{AtomIdx, BondOrder, CipCode, Molecule, apply_kekule, kekulize};

use super::ffi::{InchiAtom, InchiStereo0D, MAXVAL};

/// Error from molecule-to-InChI-input conversion.
#[derive(Debug)]
pub enum ConvertError {
    /// Aromatic bonds could not be Kekulized.
    KekulizationFailed(String),
}

/// Convert a `Molecule` into the atom + stereo lists required by the IUPAC InChI C API.
///
/// Aromatic bonds are Kekulized first. Tetrahedral stereo is derived from CIP R/S codes:
/// neighbors are sorted by decreasing CIP priority, then CIP R → EVEN parity, S → ODD.
/// Chiral centres with bracket-H (e.g. `[C@H]`) get an explicit H atom appended to the
/// atom list so the Stereo0D neighbor array can reference it by index.
pub fn mol_to_inchi_atoms(
    mol: &Molecule,
) -> Result<(Vec<InchiAtom>, Vec<InchiStereo0D>), ConvertError> {
    // Kekulize aromatic bonds.
    let kekulized_mol;
    let mol = if mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        match kekulize(mol) {
            Ok(kekule) => {
                kekulized_mol = apply_kekule(mol, &kekule);
                &kekulized_mol
            }
            Err(e) => return Err(ConvertError::KekulizationFailed(e.to_string())),
        }
    } else {
        mol
    };

    // --- Phase 1: map heavy atoms to 0-based InChI indices -------------------
    let mut heavy_order: Vec<AtomIdx> = Vec::new();
    let mut inchi_idx: HashMap<AtomIdx, i16> = HashMap::new();

    for (aidx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 1 {
            let idx = heavy_order.len() as i16;
            inchi_idx.insert(aidx, idx);
            heavy_order.push(aidx);
        }
    }

    // --- Phase 2: gather stereo info & find centres that need explicit H -----
    // tetrahedral_stereo_neighbors returns CIP code + 4 neighbors sorted by
    // DECREASING CIP priority. AtomIdx(u32::MAX) is the virtual-H sentinel.
    //
    // stereo_data: (center_aidx, CipCode, [p4, p3, p2, p1])
    // stereo_h: center_aidx → InChI index of the explicit H we will add
    let mut stereo_data: Vec<(AtomIdx, CipCode, [AtomIdx; 4])> = Vec::new();
    let mut stereo_h: Vec<(AtomIdx, i16)> = Vec::new(); // ordered by heavy_order

    for &aidx in &heavy_order {
        let Some((code, sorted_nbrs)) = tetrahedral_stereo_neighbors(mol, aidx) else {
            continue;
        };
        let needs_explicit_h = sorted_nbrs.iter().any(|&nb| nb.0 == u32::MAX);
        if needs_explicit_h {
            // Index = number of heavy atoms + number of H atoms added so far.
            let h_idx = (heavy_order.len() + stereo_h.len()) as i16;
            stereo_h.push((aidx, h_idx));
        }
        stereo_data.push((aidx, code, sorted_nbrs));
    }

    // Build lookup: center_aidx → explicit-H InChI index
    let h_idx_for: HashMap<AtomIdx, i16> = stereo_h.iter().copied().collect();

    // --- Phase 3: build InchiAtom list (heavy atoms) -------------------------
    let mut c_atoms: Vec<InchiAtom> = Vec::with_capacity(heavy_order.len() + stereo_h.len());

    for &aidx in &heavy_order {
        let atom = mol.atom(aidx);
        let mut ca = InchiAtom::default();

        for (i, b) in atom.element.symbol().bytes().enumerate().take(5) {
            ca.elname[i] = b as c_char;
        }

        // Half-adjacency bonds (only list j when j < self).
        let self_ni = *inchi_idx.get(&aidx).unwrap();
        let mut nb = 0i16;
        for (neigh_idx, bond_idx) in mol.neighbors(aidx) {
            if mol.atom(neigh_idx).element.atomic_number() == 1 {
                continue;
            }
            let Some(&ni) = inchi_idx.get(&neigh_idx) else {
                continue;
            };
            if ni >= self_ni || nb as usize >= MAXVAL {
                continue;
            }
            ca.neighbor[nb as usize] = ni;
            ca.bond_type[nb as usize] = match mol.bond(bond_idx).order {
                BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                _ => 1,
            };
            nb += 1;
        }
        ca.num_bonds = nb;

        // H count: total implicit + explicit-in-graph H.
        let explicit_h: u8 = mol
            .neighbors(aidx)
            .filter(|(ni, _)| mol.atom(*ni).element.atomic_number() == 1)
            .count() as u8;
        let implicit_h = mol.implicit_hydrogen_count(aidx);
        let total_h = implicit_h + explicit_h;
        // If we're adding an explicit H atom for this stereo centre, subtract one
        // so the library doesn't double-count.
        let h_for_inchi = if h_idx_for.contains_key(&aidx) {
            total_h.saturating_sub(1)
        } else {
            total_h
        };
        ca.num_iso_h[0] = h_for_inchi as i8;

        if let Some(mass) = atom.isotope {
            ca.isotopic_mass = mass as i16;
        }
        ca.charge = atom.charge;

        c_atoms.push(ca);
    }

    // --- Phase 4: add explicit H atoms in the same order as stereo_h ---------
    // Iterate in stereo_h order (which matches heavy_order) to guarantee indices.
    for &(center_aidx, h_inchi_idx) in &stereo_h {
        debug_assert_eq!(h_inchi_idx as usize, c_atoms.len());
        let center_inchi = *inchi_idx.get(&center_aidx).unwrap();
        let mut h_atom = InchiAtom::default();
        h_atom.elname[0] = b'H' as c_char;
        // List bond H→center (half-adjacency: H has higher index, so lists the lower).
        h_atom.num_bonds = 1;
        h_atom.neighbor[0] = center_inchi;
        h_atom.bond_type[0] = 1;
        c_atoms.push(h_atom);
    }

    // --- Phase 5: build Stereo0D descriptors (tetrahedral) -------------------
    //
    // CIP R (CW from lowest priority = neighbor[3]) → EVEN (2)
    // CIP S (CCW from lowest priority) → ODD (1)
    let mut stereo: Vec<InchiStereo0D> = Vec::new();

    for (center_aidx, code, sorted_nbrs) in &stereo_data {
        let Some(&center_ni) = inchi_idx.get(center_aidx) else {
            continue;
        };

        let mut ok = true;
        let mut arr = [0i16; 4];
        for (i, &nb) in sorted_nbrs.iter().enumerate() {
            arr[i] = if nb.0 == u32::MAX {
                match h_idx_for.get(center_aidx) {
                    Some(&h) => h,
                    None => {
                        ok = false;
                        break;
                    }
                }
            } else {
                match inchi_idx.get(&nb) {
                    Some(&ni) => ni,
                    None => {
                        ok = false;
                        break;
                    }
                }
            };
        }
        if !ok {
            continue;
        }

        let parity = if *code == CipCode::R { 2i8 } else { 1i8 };

        stereo.push(InchiStereo0D {
            neighbor: arr,
            central_atom: center_ni,
            stereo_type: 2,
            parity,
        });
    }

    // --- Phase 6: build Stereo0D descriptors (E/Z double bond) -------------------
    //
    // For each double bond, look for Up/Down stereo bonds on non-H substituents.
    // Same direction (both up or both down) → Z (zusammen) → ODD (1).
    // Opposite directions → E (entgegen) → EVEN (2).
    // central_atom = -1 (NO_ATOM), stereo_type = 1 (DoubleBond).
    //
    // is_up(alkene_end, sub) mirrors substituent_is_up in chematic-chem/cip.rs:
    //   Up bond: atom1 == alkene_end → true
    //   Down bond: atom1 == sub → true (i.e. atom1 != alkene_end)
    //
    // A ring bond adjacent to an exocyclic double bond (e.g. a mancude ring
    // flanking an imine) can carry its real `/`/`\` direction in
    // `Molecule::bond_direction`'s side channel while the bond's own `order`
    // is Aromatic pre-kekulization, or Single/Double post-kekulization
    // (`apply_kekule` preserves the stash verbatim on the same bond index --
    // it only ever updates `order`, never resolves the stash into it). Read
    // that stash first so this doesn't depend on which of the two ring
    // bonds happened to carry the literal marker. (No CIP-priority ranking
    // is needed here, unlike chematic-chem's own E/Z label: InChI's Stereo0D
    // format encodes the parity relative to whichever specific substituent
    // is fed in, so any one determinate substituent per end is sufficient
    // -- see `inchi_api.h`'s 0D stereo notes.)
    let find_stereo_sub = |alkene_end: AtomIdx, other: AtomIdx| -> Option<(i16, bool)> {
        for (nb, _) in mol.neighbors(alkene_end) {
            if nb == other {
                continue;
            }
            if mol.atom(nb).element.atomic_number() == 1 {
                continue;
            }
            let Some((bond_idx, nb_bond)) = mol.bond_between(alkene_end, nb) else {
                continue;
            };
            let effective_order = mol.bond_direction(bond_idx).unwrap_or(nb_bond.order);
            let is_up = match effective_order {
                BondOrder::Up => Some(nb_bond.atom1 == alkene_end),
                BondOrder::Down => Some(nb_bond.atom1 == nb),
                _ => None,
            };
            if let (Some(up), Some(&ni)) = (is_up, inchi_idx.get(&nb)) {
                return Some((ni, up));
            }
        }
        None
    };

    for (_, bond) in mol.bonds() {
        if bond.order != BondOrder::Double {
            continue;
        }
        let a = bond.atom1;
        let b = bond.atom2;
        let (Some(&a_ni), Some(&b_ni)) = (inchi_idx.get(&a), inchi_idx.get(&b)) else {
            continue;
        };

        let Some((x_ni, x_up)) = find_stereo_sub(a, b) else {
            continue;
        };
        let Some((y_ni, y_up)) = find_stereo_sub(b, a) else {
            continue;
        };

        let parity: i8 = if x_up == y_up { 1 } else { 2 };

        stereo.push(InchiStereo0D {
            neighbor: [x_ni, a_ni, b_ni, y_ni],
            central_atom: -1,
            stereo_type: 1,
            parity,
        });
    }

    Ok((c_atoms, stereo))
}
