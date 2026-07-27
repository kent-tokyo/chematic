use std::collections::HashMap;
use std::os::raw::c_char;

use chematic_chem::tetrahedral_stereo_neighbors;
use chematic_core::{AtomIdx, BondOrder, Chirality, CipCode, Molecule, apply_kekule, kekulize};

use super::ffi::{InchiAtom, InchiStereo0D, MAXVAL};

/// Error from molecule-to-InChI-input conversion.
#[derive(Debug)]
pub enum ConvertError {
    /// Aromatic bonds could not be Kekulized.
    KekulizationFailed(String),
}

/// Where a manufactured "H" InChI atom (added only so a tetrahedral stereocentre's
/// hydrogen substituent has an index `Stereo0D` can reference -- `heavy_order`/
/// `inchi_idx` never include H, so it has no index of its own otherwise) traces
/// back to the source molecule.
///
/// * `Sentinel` -- bracket-H notation (`[C@H]`, `[C@@H]`, ...): the H is not a
///   graph atom at all, just an implicit count on the bracket atom, so it is
///   always isotopically ordinary.
/// * `Explicit(AtomIdx)` -- a real graph H atom (`[H]`, `[2H]`, `[3H]`, ...)
///   that IS present in `mol`; its isotope (if any) must be forwarded onto the
///   manufactured atom so the `/i` layer is correct.
#[derive(Clone, Copy)]
enum StereoHSource {
    Sentinel,
    Explicit(AtomIdx),
}

/// Isotope of an explicit hydrogen atom -> the `num_iso_h[]` bucket it tallies
/// into (`inchi_api.h`: `[0]` ordinary, `[1]` explicit 1H/protium, `[2]` 2H/D,
/// `[3]` 3H/T). Anything other than 1/2/3 is not a real H isotope; fall back to
/// ordinary rather than panicking or silently mis-bucketing.
fn h_isotope_bucket(isotope: Option<u16>) -> usize {
    match isotope {
        Some(1) => 1,
        Some(2) => 2,
        Some(3) => 3,
        _ => 0,
    }
}

/// True if `mol` contains a tetrahedral stereocentre with 2 or more H-like
/// substituents (the bracket-H sentinel, or an explicit graph H/D/T atom).
///
/// This is exactly the shape [`mol_to_inchi_atoms`] cannot represent (see the
/// `two_h_like_substituents_on_one_centre_drops_stereo_not_corrupts_it` /
/// `bracket_h_plus_explicit_isotope_h_drops_stereo_not_corrupts_it` tests
/// below): its single-manufactured-atom-per-centre mechanism registers at
/// most one H-like substituent, so a second one silently drops the whole
/// stereo descriptor rather than corrupting it. A non-empty
/// [`super::standard_inchi`] string is therefore NOT sufficient evidence
/// that a comparison based on it is trustworthy for such a centre --
/// callers that need a verified identity comparison
/// (`crate::dedup::IdentityPolicy`) must check this FIRST and fail closed,
/// never trust the string just because generation "succeeded".
///
/// Deliberately structural, not a SMILES allowlist: any atom where
/// [`tetrahedral_stereo_neighbors`] recognizes a genuine stereocentre (which,
/// for two identically-untagged hydrogens, it already wouldn't -- CIP rule 2
/// only distinguishes them when their isotopes actually differ) AND 2+ of
/// its 4 ranked substituents are H-like is flagged, regardless of which
/// SMILES produced it.
pub(crate) fn has_unrepresentable_multi_h_stereocenter(mol: &Molecule) -> bool {
    mol.atoms().any(|(aidx, _)| {
        let Some((_, sorted_nbrs)) = tetrahedral_stereo_neighbors(mol, aidx) else {
            return false;
        };
        let h_like_count = sorted_nbrs
            .iter()
            .filter(|&&nb| nb.0 == u32::MAX || mol.atom(nb).element.atomic_number() == 1)
            .count();
        h_like_count >= 2
    })
}

/// True if `mol` contains an atom whose SMILES explicitly specified
/// tetrahedral chirality (`@`/`@@`) but the legacy CIP-based neighbor
/// ranking this crate's native conversion depends on
/// (`chematic_chem::tetrahedral_stereo_neighbors`) could not resolve a rank
/// for it -- i.e. a genuine tie/failure in that ranking, not an atom that
/// was never claimed to be a stereocentre in the first place.
///
/// This is the confirmed root cause of a real false `VerifiedDuplicate`
/// found via live corpus verification (two real diastereomers, differing
/// at exactly the two ring atoms where this predicate fires, that
/// `standard_inchi` collapsed to one string with `?` -- undefined parity --
/// at both positions, identically for both inputs). Diagnosed directly
/// (not assumed): for the affected molecule, `atom.chirality` was correctly
/// parsed as `Clockwise`/`CounterClockwise` for all 4 specified centres, but
/// `tetrahedral_stereo_neighbors` returned `None` for exactly the 2 that
/// produced `?` in the output -- confirming the failure is in CIP
/// substituent ranking, not a downstream Stereo0D-index bug (that would be
/// [`has_unrepresentable_multi_h_stereocenter`]'s territory) nor a
/// chirality-parsing loss.
///
/// Deliberately narrower than "any molecule with 2+ ring stereocentres": it
/// targets exactly the failure mode above (specified stereo the legacy
/// engine could not rank), not a broad heuristic that would also flag
/// ordinary, correctly-resolved multi-stereocentre molecules.
///
/// Matches `CounterClockwise`/`Clockwise` explicitly rather than
/// `!= Chirality::None` -- this predicate is specifically about
/// *tetrahedral* stereo the legacy CIP ranking failed to resolve, so it must
/// not silently widen its net if a future non-tetrahedral `Chirality`
/// variant (e.g. an eventual axial/allenic chirality) is added.
pub(crate) fn has_unresolved_specified_tetrahedral_stereo(mol: &Molecule) -> bool {
    !unresolved_specified_tetrahedral_stereo_atoms(mol).is_empty()
}

/// Same predicate as [`has_unresolved_specified_tetrahedral_stereo`], but
/// returns the actual atom indices instead of a bare bool.
///
/// Added for issue #161: `crate::dedup`'s accurate-CIP preflight needs to
/// know exactly *which* centres the legacy engine failed to rank (to re-check
/// each one individually via `CipMode::Accurate`), not just whether any
/// exist. Ascending order (same order as [`Molecule::atoms`]).
pub(crate) fn unresolved_specified_tetrahedral_stereo_atoms(mol: &Molecule) -> Vec<AtomIdx> {
    mol.atoms()
        .filter(|(idx, atom)| {
            matches!(
                atom.chirality,
                Chirality::CounterClockwise | Chirality::Clockwise
            ) && tetrahedral_stereo_neighbors(mol, *idx).is_none()
        })
        .map(|(idx, _)| idx)
        .collect()
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
    // DECREASING CIP priority. AtomIdx(u32::MAX) is the virtual-H sentinel
    // (bracket-H, e.g. `[C@H]`); a real graph H atom (e.g. `[C@](Br)(Cl)(F)[H]`
    // or an isotope thereof) shows up here as its own real AtomIdx instead,
    // since chematic-chem's stereo_neighbors reads it straight out of the
    // molecule graph. Either way, H is never in `heavy_order`/`inchi_idx`
    // (Phase 1 deliberately excludes it), so neither can be referenced by
    // Stereo0D without a manufactured stand-in atom (Phase 4).
    //
    // stereo_data: (center_aidx, CipCode, [p4, p3, p2, p1])
    // stereo_h: (center_aidx, InChI index of the manufactured H atom, its source)
    let mut stereo_data: Vec<(AtomIdx, CipCode, [AtomIdx; 4])> = Vec::new();
    let mut stereo_h: Vec<(AtomIdx, i16, StereoHSource)> = Vec::new(); // ordered by heavy_order

    for &aidx in &heavy_order {
        let Some((code, sorted_nbrs)) = tetrahedral_stereo_neighbors(mol, aidx) else {
            continue;
        };
        // Usually at most one substituent is hydrogen (bracket-H or a real
        // graph H atom), since a tetrahedral stereocentre needs 4 CIP-distinct
        // groups and two *ordinary* hydrogens would tie. But two isotopically
        // distinct hydrogens (e.g. a centre bearing both `[2H]` and `[3H]`, or
        // bracket-H plus an explicit `[2H]`) DO rank differently under CIP
        // rule 2 (mass) and are a legitimate stereocentre -- confirmed
        // reachable and RDKit-stereogenic (see `two_h_like_substituents_*`
        // tests below). This mechanism only has ONE manufactured-atom slot
        // per centre, so it cannot represent two independently-indexed H
        // substituents; register at most one (`h_subs.len() == 1`) and leave
        // multi-H centres unregistered so Phase 5's existing "no mapping ->
        // ok=false -> drop the descriptor" path applies -- safe (same
        // behaviour as any other unsupported case in this file) rather than
        // silently emitting a corrupted Stereo0D with a duplicated index.
        let mut h_subs = sorted_nbrs
            .iter()
            .copied()
            .filter(|&nb| nb.0 == u32::MAX || mol.atom(nb).element.atomic_number() == 1);
        let first = h_subs.next();
        let is_single = first.is_some() && h_subs.next().is_none();
        if is_single {
            let nb = first.unwrap();
            // Index = number of heavy atoms + number of H atoms added so far.
            let h_idx = (heavy_order.len() + stereo_h.len()) as i16;
            let source = if nb.0 == u32::MAX {
                StereoHSource::Sentinel
            } else {
                StereoHSource::Explicit(nb)
            };
            stereo_h.push((aidx, h_idx, source));
        }
        stereo_data.push((aidx, code, sorted_nbrs));
    }

    // Build lookups: center_aidx → manufactured-H InChI index, and
    // center_aidx → which num_iso_h[] bucket that manufactured H "steals"
    // from (so it isn't ALSO tallied via the count-based H-layer below).
    let h_idx_for: HashMap<AtomIdx, i16> = stereo_h.iter().map(|&(c, i, _)| (c, i)).collect();
    let h_bucket_for: HashMap<AtomIdx, usize> = stereo_h
        .iter()
        .map(|&(c, _, source)| {
            let bucket = match source {
                StereoHSource::Sentinel => 0,
                StereoHSource::Explicit(h_aidx) => h_isotope_bucket(mol.atom(h_aidx).isotope),
            };
            (c, bucket)
        })
        .collect();

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

        // H count: implicit + explicit-in-graph H, bucketed by isotope.
        // num_iso_h[0]=ordinary (no isotope tag), [1]=explicit 1H(protium),
        // [2]=explicit 2H(D), [3]=explicit 3H(T) -- see inchi_api.h. Implicit
        // H is never isotope-tagged (isotopes must be written explicitly in
        // the source graph), so it always tallies into bucket 0.
        let mut iso_h = [0u8; 4];
        for (nb, _) in mol.neighbors(aidx) {
            if mol.atom(nb).element.atomic_number() != 1 {
                continue;
            }
            let bucket = h_isotope_bucket(mol.atom(nb).isotope);
            iso_h[bucket] = iso_h[bucket].saturating_add(1);
        }
        iso_h[0] = iso_h[0].saturating_add(mol.implicit_hydrogen_count(aidx));
        // If this centre's stereo H-substituent is represented by a manufactured
        // InChI atom instead (Phase 4), don't ALSO count it here.
        if let Some(&bucket) = h_bucket_for.get(&aidx) {
            iso_h[bucket] = iso_h[bucket].saturating_sub(1);
        }
        ca.num_iso_h = [
            iso_h[0] as i8,
            iso_h[1] as i8,
            iso_h[2] as i8,
            iso_h[3] as i8,
        ];

        if let Some(mass) = atom.isotope {
            ca.isotopic_mass = mass as i16;
        }
        ca.charge = atom.charge;

        c_atoms.push(ca);
    }

    // --- Phase 4: add explicit H atoms in the same order as stereo_h ---------
    // Iterate in stereo_h order (which matches heavy_order) to guarantee indices.
    for &(center_aidx, h_inchi_idx, source) in &stereo_h {
        debug_assert_eq!(h_inchi_idx as usize, c_atoms.len());
        let center_inchi = *inchi_idx.get(&center_aidx).unwrap();
        let mut h_atom = InchiAtom::default();
        h_atom.elname[0] = b'H' as c_char;
        // List bond H→center (half-adjacency: H has higher index, so lists the lower).
        h_atom.num_bonds = 1;
        h_atom.neighbor[0] = center_inchi;
        h_atom.bond_type[0] = 1;
        // A real graph H atom (not the bracket-H sentinel) may itself carry an
        // isotope (`[2H]`, `[3H]`) -- forward it onto the manufactured atom,
        // same as the heavy-atom isotope path below (`ca.isotopic_mass`).
        if let StereoHSource::Explicit(h_aidx) = source
            && let Some(mass) = mol.atom(h_aidx).isotope
        {
            h_atom.isotopic_mass = mass as i16;
        }
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
            // Sentinel (bracket-H) or a real graph H atom: neither has an
            // entry in `inchi_idx` (H is never a heavy atom), so both route
            // through the manufactured stand-in atom from Phase 4 instead.
            let is_h = nb.0 == u32::MAX || mol.atom(nb).element.atomic_number() == 1;
            arr[i] = if is_h {
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

#[cfg(test)]
mod tests {
    use super::mol_to_inchi_atoms;
    use chematic_smiles::parse;

    /// Structural check (not just final-string matching): a real explicit
    /// graph H atom that is a stereocentre's substituent must NOT be folded
    /// into the heavy-atom numbering. The only correct way to give it a
    /// Stereo0D-referenceable index is via a manufactured extra atom appended
    /// AFTER all heavy atoms (Phase 4) -- so `c_atoms.len()` must be exactly
    /// `heavy_atom_count + 1`, and that manufactured atom must be the last
    /// entry, elname "H".
    #[test]
    fn explicit_h_stereo_substituent_is_manufactured_not_folded_into_heavy_atoms() {
        // [C@](Br)(Cl)(F)[H]: 4 heavy atoms (C, Br, Cl, F) + 1 manufactured H.
        let mol = parse("[C@](Br)(Cl)(F)[H]").unwrap();
        let (atoms, stereo) = mol_to_inchi_atoms(&mol).unwrap();
        assert_eq!(
            atoms.len(),
            5,
            "expected 4 heavy atoms + 1 manufactured H atom, got {}",
            atoms.len()
        );
        let last = atoms.last().unwrap();
        assert_eq!(last.elname[0] as u8 as char, 'H');
        assert_eq!(
            last.num_bonds, 1,
            "manufactured H atom must have exactly 1 bond (back to its centre)"
        );
        // Stereo0D must reference the manufactured atom (index 4), not any
        // heavy-atom index or the sentinel value itself.
        assert_eq!(stereo.len(), 1);
        assert!(
            stereo[0].neighbor.contains(&4),
            "Stereo0D neighbor list must reference the manufactured H atom's index (4): {:?}",
            stereo[0].neighbor
        );
    }

    /// Same structural check for an isotopically-labeled explicit H substituent:
    /// the manufactured atom must carry the isotope (D=2), not the centre atom's
    /// own isotopic_mass field, and must still not be folded into heavy atoms.
    #[test]
    fn explicit_deuterium_stereo_substituent_isotope_on_manufactured_atom() {
        let mol = parse("[C@](Br)(Cl)(F)[2H]").unwrap();
        let (atoms, _stereo) = mol_to_inchi_atoms(&mol).unwrap();
        assert_eq!(atoms.len(), 5);
        let centre = &atoms[0];
        assert_eq!(
            centre.isotopic_mass, 0,
            "the heavy stereocentre atom itself must not carry the H's isotope"
        );
        assert_eq!(
            centre.num_iso_h,
            [0, 0, 0, 0],
            "the D substituent must be represented via the manufactured atom, not double-counted in num_iso_h"
        );
        let manufactured = atoms.last().unwrap();
        assert_eq!(manufactured.elname[0] as u8 as char, 'H');
        assert_eq!(
            manufactured.isotopic_mass, 2,
            "manufactured atom must carry the D isotope mass"
        );
    }

    /// Two independent stereocentres, each with its own explicit-graph-H
    /// substituent (D and T respectively): two manufactured atoms must be
    /// appended, with distinct indices, each carrying its own isotope, and
    /// each Stereo0D descriptor must reference its OWN manufactured atom
    /// (not the other centre's).
    #[test]
    fn two_stereocenters_each_get_their_own_manufactured_atom() {
        let mol = parse("[C@](Br)(Cl)([2H])[C@@]([3H])(F)I").unwrap();
        let (atoms, stereo) = mol_to_inchi_atoms(&mol).unwrap();
        // 6 heavy atoms (C,Br,Cl,C,F,I) + 2 manufactured H atoms.
        assert_eq!(
            atoms.len(),
            8,
            "expected 6 heavy atoms + 2 manufactured H atoms"
        );
        let manufactured = &atoms[6..8];
        let masses: Vec<i16> = manufactured.iter().map(|a| a.isotopic_mass).collect();
        assert_eq!(
            masses,
            vec![2, 3],
            "manufactured atoms must carry D then T in stereo_h/heavy_order traversal order"
        );
        assert_eq!(stereo.len(), 2);
        // Each descriptor must reference exactly one manufactured index (6 or 7),
        // and the two descriptors must reference DIFFERENT indices.
        let referenced: Vec<i16> = stereo
            .iter()
            .map(|s| {
                let hits: Vec<i16> = s
                    .neighbor
                    .iter()
                    .copied()
                    .filter(|&n| n == 6 || n == 7)
                    .collect();
                assert_eq!(
                    hits.len(),
                    1,
                    "each stereocentre must reference exactly one manufactured atom: {:?}",
                    s.neighbor
                );
                hits[0]
            })
            .collect();
        assert_ne!(
            referenced[0], referenced[1],
            "the two centres must reference DIFFERENT manufactured atoms, not cross-talk"
        );
    }

    /// Non-stereo explicit isotopic H (bug 1 in isolation, no stereocentre
    /// involved at all): methane-d4 must tally 4 deuteriums into num_iso_h[2],
    /// not num_iso_h[0], and must not manufacture any extra atom (no stereo).
    #[test]
    fn methane_d4_isotope_tally_no_stereo() {
        let mol = parse("[2H]C([2H])([2H])[2H]").unwrap();
        let (atoms, stereo) = mol_to_inchi_atoms(&mol).unwrap();
        assert_eq!(atoms.len(), 1, "no stereocentre => no manufactured atom");
        assert!(stereo.is_empty());
        assert_eq!(
            atoms[0].num_iso_h,
            [0, 0, 4, 0],
            "4 explicit D atoms must tally into bucket 2 (D), not bucket 0 (ordinary)"
        );
    }

    /// A stereocentre with TWO isotopically-distinct H-like substituents
    /// (here D and T, which DO rank differently under CIP rule 2 and so form
    /// a real stereocentre -- confirmed reachable and RDKit-stereogenic) has
    /// no single manufactured-atom slot that can hold both. Must NOT emit a
    /// Stereo0D with a duplicated neighbor index (that would be silently
    /// wrong); must instead safely drop the descriptor while keeping the
    /// isotope counts on the heavy atom itself fully correct.
    #[test]
    fn two_h_like_substituents_on_one_centre_drops_stereo_not_corrupts_it() {
        let mol = parse("[C@](Br)(F)([2H])[3H]").unwrap();
        let (atoms, stereo) = mol_to_inchi_atoms(&mol).unwrap();
        assert_eq!(
            atoms.len(),
            3,
            "no manufactured atom for either H: 3 heavy atoms only (C, Br, F)"
        );
        assert!(
            stereo.is_empty(),
            "ambiguous double-H-slot centre must be dropped, not emitted with a duplicated index"
        );
        // The D and T themselves must still be correctly isotope-tallied on
        // the heavy carbon -- only the stereo *descriptor* is unsupported.
        assert_eq!(
            atoms[0].num_iso_h,
            [0, 0, 1, 1],
            "D and T must still tally correctly even though stereo is dropped"
        );
    }

    /// Same guard, but for bracket-H (sentinel) + an explicit isotopic H on
    /// the same centre -- also two H-like substituents, also must drop
    /// cleanly rather than duplicate an index.
    #[test]
    fn bracket_h_plus_explicit_isotope_h_drops_stereo_not_corrupts_it() {
        let mol = parse("[C@H](Br)([2H])F").unwrap();
        let (atoms, stereo) = mol_to_inchi_atoms(&mol).unwrap();
        assert_eq!(
            atoms.len(),
            3,
            "no manufactured atom: 3 heavy atoms (C, Br, F)"
        );
        assert!(stereo.is_empty(), "must be dropped, not corrupted");
        assert_eq!(
            atoms[0].num_iso_h,
            [1, 0, 1, 0],
            "bracket-H implicit count (bucket 0 = 1) + explicit D (bucket 2 = 1)"
        );
    }
}
