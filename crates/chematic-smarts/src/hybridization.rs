//! RDKit's hybridization model for the SMARTS `^n` primitive.
//!
//! Follows RDKit 2026.03.6 `ConjugHybrid.cpp`: the number of orbitals is the
//! total degree plus lone pairs, and an atom with four orbitals, fewer than
//! four neighbours and a conjugated bond is SP2 (amide N, ester O,
//! carboxylate O⁻). Conjugation follows `setConjugation`: aromatic bonds,
//! plus `markConjAtomBonds` around each 2- or 3-coordinate candidate atom.
//! Plain terminal H atoms report no hybridization, as for H atoms added by
//! RDKit's `AddHs`.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

/// Outer-shell electrons of a main-group element (`None` otherwise).
fn outer_electrons(z: u8) -> Option<i32> {
    let z = i32::from(z);
    Some(match z {
        1 => 1,
        2 => 2,
        3..=10 => z - 2,
        11..=18 => z - 10,
        19 | 37 | 55 | 87 => 1,
        20 | 38 | 56 | 88 => 2,
        31..=36 => z - 28,
        49..=54 => z - 46,
        81..=86 => z - 78,
        // Transition metals, as RDKit's periodic table counts them
        // (group 3–11 → 3–11, group 12 → 2).
        21..=29 => z - 18,
        39..=47 => z - 36,
        72..=79 => z - 68,
        57 => 3,
        30 | 48 | 80 => 2,
        _ => return None,
    })
}

/// RDKit 2026.03.6's allowed-valence list (`GetValenceList`); `None` for
/// elements it lists as `[-1]` (any valence: transition metals and the rest).
fn rdkit_valence_list(z: u8) -> Option<&'static [i32]> {
    Some(match z {
        1 | 9 | 17 | 35 => &[1],
        2 | 10 | 18 | 36 | 86 => &[0],
        3 | 11 | 19 | 37 => &[1, -1],
        55 | 87 => &[1],
        4 => &[2],
        12 | 20 | 38 | 56 | 88 => &[2, -1],
        5 | 13 | 31 | 49 => &[3],
        6 | 14 | 32 => &[4],
        7 => &[3],
        8 => &[2],
        15 | 33 | 51 | 83 => &[3, 5],
        16 | 34 | 52 | 84 => &[2, 4, 6],
        50 | 82 => &[2, 4],
        53 | 85 => &[1, 3, 5],
        54 => &[0, 2, 4, 6],
        _ => return None,
    })
}

/// RDKit `assignRadicals`: radical electrons of an atom with no implicit H
/// (a bracket atom), from its total valence.
fn rdkit_radicals(mol: &Molecule, idx: AtomIdx, nouter: i32, total_valence: i32) -> i32 {
    let atom = mol.atom(idx);
    if atom.hydrogen_count.is_none() {
        return 0;
    }
    let z = atom.element.atomic_number();
    let chg = i32::from(atom.charge);
    let Some(valens) = rdkit_valence_list(z) else {
        if mol.degree(idx) > 0 {
            return 0;
        }
        return (nouter - chg).max(0) % 2;
    };
    let base = if z <= 2 { 2 } else { 8 };
    let mut radicals = base - nouter - total_valence + chg;
    if radicals < 0 {
        radicals = 0;
        if valens.len() > 1
            && let Some(&v) = valens.iter().find(|&&v| v - total_valence + chg >= 0)
        {
            radicals = v - total_valence + chg;
        }
    }
    let early = nouter - total_valence - chg;
    if early >= 0 {
        radicals = radicals.min(early);
    }
    radicals
}

/// RDKit's default (first) valence for the elements conjugation looks at.
fn default_valence(mol: &Molecule, idx: AtomIdx) -> i32 {
    let element = mol.atom(idx).element;
    element
        .normal_valences()
        .first()
        .map_or(-1, |&v| i32::from(v))
}

/// Per-atom facts the model reads repeatedly, computed once per molecule
/// (H counts, degrees) or on first use (conjugation candidacy, bond
/// conjugation): a `^n` query asks for every candidate atom.
struct Facts<'a> {
    mol: &'a Molecule,
    lazy: &'a LazyFacts,
}

/// Per-atom and per-bond facts of one molecule, each computed on first use
/// and memoized on the molecule with the codes (a `^n` query asks for the
/// atoms that pass its other primitives, not every atom). Atomics so the
/// memo is shareable; every value is a pure function of the molecule.
pub(crate) struct LazyFacts {
    /// Implicit-H count + 1 (0 = not yet known).
    hs: Vec<std::sync::atomic::AtomicU8>,
    /// 0 = not yet known, 1 = false, 2 = true.
    candidate: Vec<std::sync::atomic::AtomicU8>,
    conjugated: Vec<std::sync::atomic::AtomicU8>,
    /// Hybridization code + 1, 255 = unspecified, 0 = not yet known.
    codes: Vec<std::sync::atomic::AtomicU8>,
    /// Total valence + 1 (0 = not yet known).
    valence: Vec<std::sync::atomic::AtomicU8>,
}

impl LazyFacts {
    fn new(mol: &Molecule) -> Self {
        let zeros = |n: usize| {
            (0..n)
                .map(|_| std::sync::atomic::AtomicU8::new(0))
                .collect()
        };
        Self {
            hs: zeros(mol.atom_count()),
            candidate: zeros(mol.atom_count()),
            conjugated: zeros(mol.bond_count()),
            codes: zeros(mol.atom_count()),
            valence: zeros(mol.atom_count()),
        }
    }
}

fn memo(cell: &std::sync::atomic::AtomicU8, compute: impl FnOnce() -> bool) -> bool {
    use std::sync::atomic::Ordering::Relaxed;
    match cell.load(Relaxed) {
        0 => {
            let value = compute();
            cell.store(if value { 2 } else { 1 }, Relaxed);
            value
        }
        known => known == 2,
    }
}

fn total_hs(f: &Facts, idx: AtomIdx) -> i32 {
    use std::sync::atomic::Ordering::Relaxed;
    let cell = &f.lazy.hs[idx.0 as usize];
    match cell.load(Relaxed) {
        0 => {
            let hs = implicit_hcount(f.mol, idx);
            // Counts above 254 are not cached (never seen; kept exact).
            if hs < 255 {
                cell.store(hs + 1, Relaxed);
            }
            i32::from(hs)
        }
        known => i32::from(known - 1),
    }
}

/// `match_vf2::total_valence`, memoized.
fn total_valence(f: &Facts, idx: AtomIdx) -> i32 {
    use std::sync::atomic::Ordering::Relaxed;
    let cell = &f.lazy.valence[idx.0 as usize];
    match cell.load(Relaxed) {
        0 => {
            let v = crate::match_vf2::total_valence(f.mol, idx);
            if v < 255 {
                cell.store(v + 1, Relaxed);
            }
            i32::from(v)
        }
        known => i32::from(known - 1),
    }
}

fn heavy_degree(f: &Facts, idx: AtomIdx) -> i32 {
    f.mol.neighbors(idx).count() as i32
}

/// RDKit `countAtomElec`: electrons available to a pi system (`-1`: none).
fn count_atom_elec(f: &Facts, idx: AtomIdx) -> i32 {
    let mol = f.mol;
    let atom = mol.atom(idx);
    let dv = default_valence(mol, idx);
    if dv <= 1 {
        return -1;
    }
    let degree = heavy_degree(f, idx) + total_hs(f, idx);
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = outer_electrons(atom.element.atomic_number()) else {
        return -1;
    };
    let nlp = (nouter - dv - i32::from(atom.charge)).max(0);
    let mut res = (dv - degree) + nlp;
    if res > 1 {
        let explicit_valence = total_valence(f, idx) - total_hs(f, idx);
        if explicit_valence - heavy_degree(f, idx) > 1 {
            res = 1;
        }
    }
    res
}

/// RDKit `isAtomConjugCand`.
fn conjugation_candidate(f: &Facts, idx: AtomIdx) -> bool {
    memo(&f.lazy.candidate[idx.0 as usize], || {
        conjugation_candidate_uncached(f, idx)
    })
}

fn conjugation_candidate_uncached(f: &Facts, idx: AtomIdx) -> bool {
    let z = f.mol.atom(idx).element.atomic_number();
    let nouter = outer_electrons(z).unwrap_or(0);
    // A group-16 atom beyond the first row only with one substituent
    // counting H (`C=C[S-]` conjugates, the thiol in `C=CS` does not).
    (z <= 10 || (nouter != 5 && nouter != 6) || (nouter == 6 && substituents(f, idx) < 2))
        && count_atom_elec(f, idx) > 0
}

fn valence_contrib(order: BondOrder) -> f32 {
    match order {
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Quadruple => 4.0,
        _ => 1.0,
    }
}

fn substituents(f: &Facts, idx: AtomIdx) -> i32 {
    heavy_degree(f, idx) + total_hs(f, idx)
}

/// Whether `markConjAtomBonds` run at `center` marks `bond`.
fn marked_at(f: &Facts, center: AtomIdx, bond: BondIdx) -> bool {
    let mol = f.mol;
    if !conjugation_candidate(f, center) {
        return false;
    }
    let sbo = substituents(f, center);
    if !(2..=3).contains(&sbo) {
        return false;
    }
    let other_ok = |b: BondIdx| {
        let e = mol.bond(b);
        let other = if e.atom1 == center { e.atom2 } else { e.atom1 };
        substituents(f, other) <= 3 && conjugation_candidate(f, other)
    };
    // A multiple bond counts only toward a candidate atom: a double bond to
    // `[S+]`/`[Se+]` with two substituents conjugates nothing (RDKit 2026.03.6
    // marks no bond in `CN(C)C(=[S+]C)C`; the metal dithiocarbamates of the
    // exposed 10k have sp3 amine N).
    let multiple = |b: BondIdx| {
        let e = mol.bond(b);
        let other = if e.atom1 == center { e.atom2 } else { e.atom1 };
        valence_contrib(e.order) >= 1.5 && conjugation_candidate(f, other)
    };
    let others = || {
        mol.neighbors(center)
            .map(|(_, b)| b)
            .filter(move |&b| b != bond)
    };
    // `bond` as the multiple bond, paired with any other candidate bond …
    (multiple(bond) && others().any(other_ok))
        // … or as the partner of another multiple bond.
        || (others().any(multiple) && other_ok(bond))
}

fn bond_is_conjugated(f: &Facts, bond: BondIdx) -> bool {
    memo(&f.lazy.conjugated[bond.0 as usize], || {
        let e = f.mol.bond(bond);
        e.order == BondOrder::Aromatic || marked_at(f, e.atom1, bond) || marked_at(f, e.atom2, bond)
    })
}

/// Whether `bond` is conjugated in RDKit's model (`Bond::getIsConjugated`
/// after `MolOps::setConjugation`), on the molecule as given (pass RDKit's
/// aromaticity view). Memoized with [`rdkit_hybridization`]'s facts.
pub fn rdkit_bond_is_conjugated(mol: &Molecule, bond: BondIdx) -> bool {
    let lazy = lazy_facts(mol);
    let f = Facts { mol, lazy: &lazy };
    bond_is_conjugated(&f, bond)
}

/// RDKit hybridization as the SMARTS `^n` code (0 = S, 1 = SP, 2 = SP2,
/// 3 = SP3, 4 = SP3D, 5 = SP3D2); `None` = unspecified.
///
/// Follows RDKit's `ConjugHybrid` model: orbital count = total degree plus
/// lone pairs, and a four-orbital atom with fewer than four neighbours and
/// a conjugated bond (amide N, aryl ether O, enamine N) is SP2.
///
/// The codes of every atom are computed together on the first call and
/// memoized on `mol` (a `^n` query reads them once per candidate atom).
pub fn rdkit_hybridization(mol: &Molecule, idx: AtomIdx) -> Option<u8> {
    hybridization_from(mol, &lazy_facts(mol), idx)
}

/// The memo of [`rdkit_hybridization`] for `mol`.
pub(crate) fn lazy_facts(mol: &Molecule) -> std::sync::Arc<LazyFacts> {
    mol.derived(chematic_core::DerivedSlot::RdkitHybridization, || {
        LazyFacts::new(mol)
    })
}

/// [`rdkit_hybridization`] with its memo already looked up.
pub(crate) fn hybridization_from(mol: &Molecule, lazy: &LazyFacts, idx: AtomIdx) -> Option<u8> {
    use std::sync::atomic::Ordering::Relaxed;
    let cell = lazy.codes.get(idx.0 as usize)?;
    match cell.load(Relaxed) {
        0 => {
            let facts = Facts { mol, lazy };
            let code = rdkit_hybridization_uncached(&facts, idx);
            cell.store(code.map_or(255, |c| c + 1), Relaxed);
            code
        }
        255 => None,
        known => Some(known - 1),
    }
}

fn rdkit_hybridization_uncached(f: &Facts, idx: AtomIdx) -> Option<u8> {
    let mol = f.mol;
    let atom = mol.atom(idx);
    let z = atom.element.atomic_number();
    if atom.wildcard || z == 0 {
        return None;
    }
    // A plain terminal H is what RDKit's `AddHs` adds, with no hybridization;
    // an H RDKit keeps from SMILES (`[2H]`, `[H+]`, `[H][H]`) is S.
    if z == 1 {
        let plain = atom.isotope.is_none() && atom.charge == 0 && mol.degree(idx) == 1 && {
            let (nb, _) = mol.neighbors(idx).next().unwrap();
            mol.atom(nb).element.atomic_number() != 1
        };
        return (!plain).then_some(0);
    }
    let total_degree = heavy_degree(f, idx) + total_hs(f, idx);
    // An aromatic atom with two or three connections has three orbitals or
    // a conjugated lone pair: SP2 in RDKit's model. Answering here skips the
    // Kekulé valence lookup the general rule needs for aromatic atoms.
    if atom.aromatic && (2..=3).contains(&total_degree) {
        return Some(2);
    }
    // RDKit's clean-up reads a neutral Cl/Br/I of valence 3, 5 or 7 bonded
    // only to O as `[X+k]([O-])…`: its `=O` atoms are sp3 `[O-]` (perchlorate).
    if z == 8
        && atom.charge == 0
        && total_degree == 1
        && let Some((halogen, bond)) = mol.neighbors(idx).next()
        && mol.bond(bond).order == BondOrder::Double
        && matches!(mol.atom(halogen).element.atomic_number(), 17 | 35 | 53)
        && mol.atom(halogen).charge == 0
        && matches!(crate::match_vf2::total_valence(mol, halogen), 3 | 5 | 7)
        && mol
            .neighbors(halogen)
            .all(|(nb, _)| mol.atom(nb).element.atomic_number() == 8)
    {
        return Some(3);
    }
    let norbs = match outer_electrons(z) {
        Some(nouter) if z < 89 => {
            let total_valence = total_valence(f, idx);
            // Not clamped: RDKit lets a negative lone-pair count lower the
            // orbital count (`[Zn++]` with four bonds is SP).
            let chg = i32::from(atom.charge);
            let free = nouter - (total_valence + chg);
            if total_valence + nouter - chg < 8 {
                // Below an octet RDKit counts radical electrons as orbitals
                // (`[Mg]` with four bonds to O+ has two: SP3).
                let radicals = rdkit_radicals(mol, idx, nouter, total_valence);
                total_degree + (free - radicals) / 2 + radicals
            } else {
                total_degree + free / 2
            }
        }
        _ => total_degree,
    };
    match norbs {
        0 | 1 => Some(0),
        2 => Some(1),
        3 => Some(2),
        // Conjugation is only looked up for an atom with fewer than four
        // neighbours (it decided nothing for an sp3 CH3/CH2, and was most
        // of a `[C^3]` match).
        4 => Some(
            if total_degree < 4 && mol.neighbors(idx).any(|(_, b)| bond_is_conjugated(f, b)) {
                2
            } else {
                3
            },
        ),
        5 => Some(4),
        6 => Some(5),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hybridization_matches_rdkit() {
        // (SMILES, per-atom RDKit 2026.03.6 `^n` code; `-` = unspecified).
        for (smiles, want) in [
            ("CC(=O)[O-]", "3222"),
            ("CC(=O)NC", "32223" /* amide N is SP2 */),
            ("COc1ccccc1", "32222222"),
            ("CC(S)=N", "3232"),
            ("C=C[S-]", "222"),
            ("C=CCl", "223"),
            (
                "CN(C)C(=[S+]C)C",
                "3332233", /* no conjugation through [S+] */
            ),
            (
                "OCl(=O)(=O)=O",
                "33333", /* clean-up: [Cl+3]([O-])... */
            ),
            ("CC(=O)N=S1OCCO1", "322233333"),
            ("CC#N", "311"),
            ("C[N+](C)(C)C", "33333"),
            ("[2H]C", "03"),
            // Below an octet radical electrons count as orbitals (RDKit
            // `assignRadicals` on bracket atoms).
            (
                "CC1=[O+][Mg]2([O+]=C(C)C1)[O+]=C(C)CC(=[O+]2)C",
                "322322332233223",
            ),
            ("C[CH2]", "33"),
            ("[CH2]", "3"),
            ("C[CH]C", "333"),
            ("C[O]", "33"),
            ("C[N]C", "333"),
            ("[Mg](C)C", "133"),
            ("C[B]C", "323"),
            ("[Na]C", "03"),
            ("C[Sn](C)(C)C", "33333"),
            ("[Li]", "0"),
            ("C[Al](C)C", "3233"),
            ("[Ca]", "1"),
            ("[SiH2]", "3"),
            ("CC[Ge]C", "3333"),
        ] {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let got: String = (0..mol.atom_count())
                .map(|i| match rdkit_hybridization(&mol, AtomIdx(i as u32)) {
                    Some(h) => char::from(b'0' + h),
                    None => '-',
                })
                .collect();
            assert_eq!(got, want, "{smiles}");
        }
    }
}
