//! The parts of RDKit's `MolFromSmiles` sanitization that rewrite the graph
//! and change what RDKit-compatible perception sees: `MolOps::cleanUp`'s
//! `halogenCleanup` and `MolOps::cleanUpOrganometallics`.
//!
//! Applied only to the RDKit-parity view
//! ([`crate::apply_aromaticity_rdkit_parity_shared`]); the caller's molecule
//! is never changed.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// RDKit `halogenCleanup` precondition: a neutral Cl/Br/I of explicit
/// valence 3, 5 or 7 bonded only to oxygen, with at least one `=O` (without
/// one the rewrite is a no-op).
fn halogen_cleanup_applies(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    if atom.wildcard || atom.charge != 0 || !matches!(atom.element.atomic_number(), 17 | 35 | 53) {
        return false;
    }
    let mut ev = 0.0f64;
    let mut has_double = false;
    for (nb, b) in mol.neighbors(idx) {
        if mol.atom(nb).element.atomic_number() != 8 || mol.atom(nb).wildcard {
            return false;
        }
        let order = mol.bond(b).order;
        has_double |= order == BondOrder::Double;
        ev += match order {
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Quadruple => 4.0,
            BondOrder::Aromatic => 1.5,
            BondOrder::Zero => 0.0,
            _ => 1.0,
        };
    }
    ev += f64::from(atom.hydrogen_count.unwrap_or(0));
    let ev = (ev + 0.1) as i32;
    has_double && matches!(ev, 3 | 5 | 7)
}

/// Whether `idx` is a metal in RDKit's `QueryOps::isMetal` sense.
fn is_metal_atom(mol: &Molecule, idx: AtomIdx) -> bool {
    let a = mol.atom(idx);
    !a.wildcard
        && !matches!(
            a.element.atomic_number(),
            1 | 2 | 5..=10 | 14..=18 | 33..=36 | 52..=54 | 85 | 86
        )
}

/// Cheap linear pre-check: could [`rdkit_sanitize_cleanup`] change `mol`?
/// (A metal atom, or a halogen oxoacid group.)
pub fn rdkit_sanitize_cleanup_may_apply(mol: &Molecule) -> bool {
    may_need_cleanup(mol)
}

fn may_need_cleanup(mol: &Molecule) -> bool {
    mol.atoms()
        .any(|(idx, _)| is_metal_atom(mol, idx) || halogen_cleanup_applies(mol, idx))
}

/// `mol` with RDKit's sanitization rewrites applied, or `None` when they
/// change nothing:
///
/// * `halogenCleanup`: `X(=O)(=O)(=O)O` → `[X+3]([O-])([O-])([O-])O`
///   (perchlorate; likewise chlorate and chlorite, Br and I).
/// * `cleanUpOrganometallics`: the single bond from a hypervalent non-metal
///   to a metal becomes dative (ferrocenyl `[C-]`→Fe), chosen as
///   [`crate::rdkit_organometallic_dative_bonds`] ports it. chematic stores
///   a dative bond donor → acceptor, so only a bond already written from
///   the non-metal is converted.
pub fn rdkit_sanitize_cleanup(mol: &Molecule) -> Option<Molecule> {
    if !may_need_cleanup(mol) {
        return None;
    }
    let mut out: Option<Molecule> = None;
    for (idx, _) in mol.atoms() {
        if !halogen_cleanup_applies(mol, idx) {
            continue;
        }
        let m = out.get_or_insert_with(|| mol.clone());
        let mut charge = 0i8;
        for (nb, b) in mol.neighbors(idx) {
            if mol.bond(b).order == BondOrder::Double {
                m.set_bond_order(b, BondOrder::Single);
                m.set_charge(nb, -1);
                charge += 1;
            }
        }
        m.set_charge(idx, charge);
    }
    let base = out.as_ref().unwrap_or(mol);
    if base.atoms().any(|(idx, _)| is_metal_atom(base, idx)) {
        let dative = crate::rdkit_sssr_order::organometallic_dative_bonds(base);
        let convert: Vec<BondIdx> = dative
            .iter()
            .enumerate()
            .filter(|&(i, &d)| d && !is_metal_atom(base, base.bond(BondIdx(i as u32)).atom1))
            .map(|(i, _)| BondIdx(i as u32))
            .collect();
        if !convert.is_empty() {
            let m = out.get_or_insert_with(|| mol.clone());
            for b in convert {
                m.set_bond_order(b, BondOrder::Dative);
            }
        }
    }
    out
}

/// Whether [`rdkit_sanitize_cleanup`] would change `mol`.
pub fn rdkit_sanitize_cleanup_needed(mol: &Molecule) -> bool {
    may_need_cleanup(mol) && rdkit_sanitize_cleanup(mol).is_some()
}
