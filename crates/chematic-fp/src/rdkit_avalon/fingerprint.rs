// Ported from the Avalon Cheminformatics Toolkit, Copyright 2001-2011 Novartis
// Pharma AG, BSD-3-Clause license; see THIRD_PARTY_NOTICES.md.
//! `CountFingerprintPatterns` of the Avalon toolkit (`ssmatch.c`), ported
//! operation for operation for `as_query == FALSE` and `exclude_atom == 0`
//! (what `SetFingerprintBits`, and so RDKit's `GetAvalonFP`, runs), with the
//! helpers it calls: `SetPathBitsRec`, `SetPathLengthFlags`,
//! `SpecialNeighboursRec`, `SetFeatureBits`, `ComputeImplicitH` /
//! `ImplicitHydrogens` (`utilities.c`) and the hash functions of
//! `hashcode.c`.

use super::rings::{
    Neighbourhood, perceive_aromatic_bonds, perceive_dy_aromaticity, ring_state,
    set_ring_size_flags, setup_neighbourhood,
};
use super::{AROMATIC, AvalonMolecule, DOUBLE, SINGLE, TRIPLE};

const ANY_BOND: i32 = 8;

pub(crate) const USE_RING_PATTERN: u32 = 0x000001;
pub(crate) const USE_RING_PATH: u32 = 0x000002;
pub(crate) const USE_ATOM_SYMBOL_PATH: u32 = 0x000004;
pub(crate) const USE_ATOM_CLASS_PATH: u32 = 0x000008;
pub(crate) const USE_ATOM_COUNT: u32 = 0x000010;
pub(crate) const USE_AUGMENTED_ATOM: u32 = 0x000020;
pub(crate) const USE_HCOUNT_PATH: u32 = 0x000040;
pub(crate) const USE_HCOUNT_CLASS_PATH: u32 = 0x000080;
pub(crate) const USE_HCOUNT_PAIR: u32 = 0x000100;
pub(crate) const USE_BOND_PATH: u32 = 0x000200;
pub(crate) const USE_AUGMENTED_BOND: u32 = 0x000400;
pub(crate) const USE_RING_SIZE_COUNTS: u32 = 0x000800;
pub(crate) const USE_DEGREE_PATH: u32 = 0x001000;
pub(crate) const USE_CLASS_SPIDERS: u32 = 0x002000;
pub(crate) const USE_FEATURE_PAIRS: u32 = 0x004000;
pub(crate) const USE_SCAFFOLD_IDS: u32 = 0x100000;
pub(crate) const USE_SCAFFOLD_COLORS: u32 = 0x200000;
pub(crate) const USE_SCAFFOLD_LINKS: u32 = 0x400000;
pub(crate) const USE_SHORTCUT_LABELS: u32 = 0x800000;
pub(crate) const USE_NON_SSS_BITS: u32 = 0xF00000;

const ATOM_SYMBOL_PATH_SEED: u64 = 17;
const ATOM_CLASS_PATH_SEED: u64 = 23;
const RING_PATTERN_SEED: u64 = 11;
const RING_PATH_SEED: u64 = 13;
const ATOM_COUNT_SEED: u64 = 31;
const AUGMENTED_ATOM_SEED: u64 = 37;
const HCOUNT_PATH_SEED: u64 = 41;
const HCOUNT_CLASS_PATH_SEED: u64 = 43;
const HCOUNT_PAIR_SEED: u64 = 47;
const BOND_PATH_SEED: u64 = 53;
const AUGMENTED_BOND_SEED: u64 = 61;
const RING_SIZE_SEED: u64 = 67;
const DEGREE_PATH_SEED: u64 = 71;
const CLASS_SPIDER_SEED: u64 = 79;
const RING_CLOSURE_SEED: i64 = 101;
const NON_SSS_SEED: u64 = 179;

const PROCESS_RING_CLOSURES: u32 = 0x0001;
const PROCESS_CHAINS: u32 = 0x0002;
const FORCED_HETERO_END: u32 = 0x0004;
const IGNORE_PATH_SYMBOL: u32 = 0x0008;
const IGNORE_TERM_SYMBOL: u32 = 0x0010;
const FORCED_RING_PATH: u32 = 0x0020;
const STOP_AT_HEAVY_ATOM: u32 = 0x0040;

const ANY_COLOR: i32 = 113;
const CSP3: i32 = 19;
const HETERO: i32 = 23;
const GENERIC: i32 = -1;
const SPECIAL_RING: u32 = 0xFC & !(1 << 6);
const MAX_SPIDER: usize = 7;

const HETERO_FLAG: i32 = 0x0100;
const RING_SUBST_FLAG: i32 = 0x0200;
const QUART_FLAG: i32 = 0x0400;
const CSP3_FLAG: i32 = 0x0800;
const RS_SPECIAL_FLAG: i32 = 0x1000;
const TYPE_MASK: i32 = 0x00FF;
const C_FLAG: i32 = 0x0001;
const O_FLAG: i32 = 0x0002;
const N_FLAG: i32 = 0x0003;
const S_FLAG: i32 = 0x0004;
const P_FLAG: i32 = 0x0005;
const X_FLAG: i32 = 0x0006;

const NCOUNT_HASH: usize = 128;
const NCOUNT_SEED_HASH: usize = 128 * 128;

/// `next_hash` (one-at-a-time step with 64-bit intermediates).
#[inline]
pub(crate) fn next_hash(hash: u64, data: u64) -> u64 {
    let mut h = hash.wrapping_add(data);
    h = h.wrapping_add(h << 10);
    h ^= h >> 6;
    h
}

/// `hash_position`.
#[inline]
pub(crate) fn hash_position(hash: u64, nslots: usize) -> usize {
    let mut h = hash;
    h = h.wrapping_add(h << 3);
    h ^= h >> 11;
    h = h.wrapping_add(h << 15);
    (h % nslots as u64) as usize
}

/// `hash_string`.
fn hash_string(s: &str) -> u64 {
    let mut h = 1001u64;
    for b in s.bytes() {
        // `*str` is a (signed) `char` converted to `uint64_t`.
        h = next_hash(h, b as i8 as i64 as u64);
    }
    h
}

/// `NEXT_SEED(seed, increment)` with a C `int`-typed increment: negative
/// values are sign-extended to `uint64_t`.
#[inline]
fn ns(seed: u64, inc: i64) -> u64 {
    next_hash(seed, inc as u64)
}

/// `periodic_table` / `StringToInt`.
pub(crate) fn periodic_number(symbol: &str) -> i32 {
    const TABLE: [&str; 94] = [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu",
    ];
    match symbol {
        "D" | "T" => 1,
        "X" => 120,
        "Q" => 121,
        "M" => 122,
        "R" => 123,
        "A" => 124,
        s => TABLE
            .iter()
            .position(|&t| t == s)
            .map_or(0, |p| p as i32 + 1),
    }
}

fn symbol_in(symbol: &str, list: &[&str]) -> bool {
    list.contains(&symbol)
}

/// One row of `valence_table`: symbol, from, to, step, lone pairs, pair
/// deficit.
type ValenceRow = (&'static str, i32, i32, i32, i32, i32);

const VALENCE_TABLE: &[ValenceRow] = &[
    ("C", 4, 4, 1, 0, 0),
    ("H", -1, 1, 2, 0, 0),
    ("N", 3, 5, 2, 1, 0),
    ("O", 2, 2, 1, 2, 0),
    ("Cl", 1, 7, 2, 0, 0),
    ("P", 3, 5, 2, 1, 0),
    ("S", 2, 6, 2, 2, 0),
    ("F", 1, 1, 1, 1, 0),
    ("H", -1, 1, 2, 0, 0),
    ("Li", -1, 1, 2, 0, 0),
    ("Na", -1, 1, 2, 0, 0),
    ("K", -1, 1, 2, 0, 0),
    ("Rb", -1, 1, 2, 0, 0),
    ("Cs", -1, 1, 2, 0, 0),
    ("Be", -2, 2, 2, 0, 0),
    ("Mg", -2, 2, 2, 0, 0),
    ("B", 3, 3, 1, 0, 1),
    ("Al", -3, 3, 2, 0, 1),
    ("Ga", -3, 3, 2, 0, 1),
    ("In", -3, 3, 2, 0, 1),
    ("Tl", -3, 3, 2, 0, 0),
    ("Si", 4, 4, 1, 0, 0),
    ("As", 3, 5, 2, 0, 0),
    ("Sb", 3, 5, 2, 0, 0),
    ("Bi", 3, 5, 2, 0, 0),
    ("Se", 2, 6, 2, 0, 0),
    ("Te", 2, 6, 2, 0, 0),
    ("La", -3, 3, 2, 0, 1),
    ("Ce", -3, 3, 2, 0, 1),
    ("Pr", -3, 3, 2, 0, 1),
    ("Nd", -3, 3, 2, 0, 1),
    ("Pm", -3, 3, 2, 0, 1),
    ("Sm", -3, 3, 2, 0, 1),
    ("Eu", -3, 3, 2, 0, 1),
    ("Gd", -3, 3, 2, 0, 1),
    ("Tb", -3, 3, 2, 0, 1),
    ("Dy", -3, 3, 2, 0, 1),
    ("Ho", -3, 3, 2, 0, 1),
    ("Er", -3, 3, 2, 0, 1),
    ("Tm", -3, 3, 2, 0, 1),
    ("Yb", -3, 3, 2, 0, 1),
    ("Lu", -3, 3, 2, 0, 1),
    ("Br", 1, 7, 2, 0, 0),
    ("I", 1, 7, 2, 0, 0),
];

/// `ImplicitHydrogens`.
fn implicit_hydrogens(
    symbol: &str,
    nsingle: i32,
    naromatic: i32,
    ndouble: i32,
    ntriple: i32,
    radical: bool,
    charge: i32,
) -> i32 {
    let mut bond_electrons = nsingle + 2 * ndouble + 3 * ntriple;
    if radical {
        bond_electrons += 1;
    }
    bond_electrons += match naromatic {
        0 => 0,
        1 => 2,
        2 => 3,
        3 => 4,
        n => n + 1,
    };
    for &(sym, from, to, step, lone_pairs, pair_deficit) in VALENCE_TABLE {
        if sym != symbol {
            continue;
        }
        let mut val = from;
        if charge == 0 {
            while val <= to {
                let h = val - bond_electrons;
                if h >= 0 {
                    return h;
                }
                val += step;
            }
        } else if charge > 0 {
            while val <= to {
                let h = val - bond_electrons + charge;
                if h >= 0 {
                    return if lone_pairs > 0 { h } else { 0 };
                }
                val += step;
            }
        } else {
            while val <= to {
                let mut h = val - bond_electrons - charge;
                if h >= 0 {
                    if pair_deficit < h {
                        h = pair_deficit;
                    }
                    return if lone_pairs > 0 { 0 } else { h };
                }
                val += step;
            }
        }
    }
    0
}

/// `ComputeImplicitH` into a zero-initialised `h_count[1..=n_atoms]`.
fn compute_implicit_h(mol: &AvalonMolecule, bond_types: &[i32], h_count: &mut [i32]) {
    let n = mol.atoms.len();
    let mut single = vec![0i32; n + 1];
    let mut aromatic = vec![0i32; n + 1];
    let mut double = vec![0i32; n + 1];
    let mut triple = vec![0i32; n + 1];
    for (i, b) in mol.bonds.iter().enumerate() {
        let (a0, a1) = (b.atoms[0] as usize, b.atoms[1] as usize);
        let arr = match bond_types[i] {
            DOUBLE => &mut double,
            TRIPLE => &mut triple,
            AROMATIC => &mut aromatic,
            _ => &mut single,
        };
        arr[a0] += 1;
        arr[a1] += 1;
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        if h_count[i + 1] == 0 {
            let h = implicit_hydrogens(
                &a.symbol,
                single[i + 1],
                aromatic[i + 1],
                double[i + 1],
                triple[i + 1],
                a.radical == 2,
                a.charge,
            );
            h_count[i + 1] = h.max(0);
        }
    }
}

/// The colour arrays and neighbourhood `SetPathBitsRec` walks.
struct Walk<'a> {
    nbp: &'a [Neighbourhood],
    atom_color: &'a [i32],
    bond_color: &'a [i32],
    bond_rsize: &'a [u32],
}

impl Walk<'_> {
    /// `SetPathBitsRec`.
    #[allow(clippy::too_many_arguments)]
    fn path_bits(
        &self,
        fp: &mut [i32],
        mut seed: u64,
        touched: &mut [i32],
        nbonds: i32,
        minbonds: i32,
        maxbonds: i32,
        sprout: usize,
        last: Option<usize>,
        flags: u32,
    ) -> i32 {
        let ncounts = fp.len();
        let mut result = 0;
        if nbonds > maxbonds {
            return result;
        }
        let nb = &self.nbp[sprout];
        for i in 0..nb.atoms.len() {
            let ai = nb.atoms[i];
            if Some(ai) == last {
                continue;
            }
            let bi = nb.bonds[i];
            if flags & FORCED_RING_PATH != 0 && self.bond_rsize[bi] == 0 {
                continue;
            }
            let acolor_ap = self.atom_color[ai];
            if acolor_ap >= 18 && acolor_ap != ANY_COLOR && flags & STOP_AT_HEAVY_ATOM != 0 {
                continue;
            }
            if touched[ai] > 0 {
                // ring closure
                if flags & PROCESS_RING_CLOSURES == 0 {
                    continue;
                }
                if touched[ai] > 1 {
                    continue;
                }
                let bcolor = self.bond_color[bi];
                if bcolor == 0 {
                    continue;
                }
                let acolor = self.atom_color[ai];
                if acolor == 0 && flags & IGNORE_PATH_SYMBOL == 0 {
                    continue;
                }
                let old_seed = seed;
                seed = ns(seed, i64::from(bcolor) * 16);
                if flags & IGNORE_PATH_SYMBOL == 0 {
                    seed = ns(seed, i64::from(acolor));
                }
                seed = ns(seed, i64::from(touched[ai]));
                fp[hash_position(seed, ncounts)] += 1;
                seed = ns(seed, RING_CLOSURE_SEED * i64::from(nbonds - touched[ai]));
                fp[hash_position(seed, ncounts)] += 1;
                result += 1;
                seed = old_seed;
            } else {
                let bcolor = self.bond_color[bi];
                if bcolor == 0 {
                    continue;
                }
                let acolor = self.atom_color[ai];
                if acolor == 0 && flags & IGNORE_PATH_SYMBOL == 0 {
                    continue;
                }
                touched[ai] = nbonds + 1;
                let old_seed = seed;
                seed = ns(seed, i64::from(bcolor) * 16);
                if flags & IGNORE_PATH_SYMBOL == 0 {
                    seed = ns(seed, i64::from(acolor));
                }
                if nbonds >= minbonds
                    && (acolor > 0 || flags & IGNORE_TERM_SYMBOL != 0)
                    && flags & PROCESS_CHAINS != 0
                    && (flags & FORCED_HETERO_END == 0
                        || flags & IGNORE_TERM_SYMBOL != 0
                        || acolor != 6)
                    && acolor > 1
                {
                    if flags & IGNORE_TERM_SYMBOL != 0 {
                        fp[hash_position(seed, ncounts)] += 1;
                    } else {
                        fp[hash_position(ns(seed, 17 * i64::from(acolor)), ncounts)] += 1;
                    }
                    result += 1;
                }
                if nbonds < maxbonds {
                    result += self.path_bits(
                        fp,
                        seed,
                        touched,
                        nbonds + 1,
                        minbonds,
                        maxbonds,
                        ai,
                        Some(sprout),
                        flags,
                    );
                }
                touched[ai] = 0;
                seed = old_seed;
            }
        }
        result
    }
}

/// `SetPathLengthFlags`.
#[allow(clippy::too_many_arguments)]
fn path_length_flags(
    atom_color: &[i32],
    nbp: &[Neighbourhood],
    touched: &mut [i32],
    start: usize,
    path_length: i32,
    current: usize,
    max_size: i32,
    lm: &mut [Vec<i32>],
) {
    for i in 0..nbp[current].atoms.len() {
        if path_length + 1 > max_size {
            continue;
        }
        let ai = nbp[current].atoms[i];
        if touched[ai] != 0 {
            continue;
        }
        if atom_color[ai] == 0 {
            continue;
        }
        touched[ai] = 1;
        lm[start][ai] |= 1 << (path_length + 1);
        path_length_flags(
            atom_color,
            nbp,
            touched,
            start,
            path_length + 1,
            ai,
            max_size,
            lm,
        );
        touched[ai] = 0;
    }
}

/// `SpecialNeighboursRec`.
#[allow(clippy::too_many_arguments)]
fn special_neighbours(
    atom_color: &[i32],
    nbp: &[Neighbourhood],
    touched: &mut [i32],
    path_length: usize,
    current: usize,
    max_size: usize,
    csp3: &mut [i32],
    hetero: &mut [i32],
) {
    for i in 0..nbp[current].atoms.len() {
        if path_length + 1 > max_size {
            continue;
        }
        let ai = nbp[current].atoms[i];
        if touched[ai] != 0 {
            continue;
        }
        let c = atom_color[ai];
        if c <= 0 {
            continue;
        }
        if c == CSP3 {
            csp3[path_length] += 1;
        } else if c == HETERO {
            hetero[path_length] += 1;
        }
        if c == HETERO && path_length > 1 {
            continue;
        }
        touched[ai] = 1;
        special_neighbours(
            atom_color,
            nbp,
            touched,
            path_length + 1,
            ai,
            max_size,
            csp3,
            hetero,
        );
        touched[ai] = 0;
    }
}

/// `SetFeatureBits` with `use_counts == FALSE`, `use_atom_types == TRUE`.
fn feature_bits(
    atom_color: &[i32],
    fp: &mut [i32],
    start_flags: i32,
    end_flags: i32,
    path_min: i32,
    path_max: i32,
    lm: &[Vec<i32>],
    start_seed: u64,
) {
    let ncounts = fp.len();
    let nslots = (ncounts * 4) as u64;
    let mut counts = vec![0i32; ncounts * 4];
    let n = atom_color.len();
    for i in 0..n {
        let coli = atom_color[i];
        if coli & start_flags == 0 {
            continue;
        }
        if coli & TYPE_MASK == 0 {
            continue;
        }
        let seed_i = ns(start_seed, i64::from(coli & TYPE_MASK));
        for j in 0..n {
            let colj = atom_color[j];
            if colj & end_flags == 0 {
                continue;
            }
            if colj & TYPE_MASK == 0 {
                continue;
            }
            let seed = ns(seed_i, i64::from(colj & TYPE_MASK));
            for k in path_min..=path_max {
                if (1 << k) & lm[i][j] != 0 {
                    let slot = ((k as u64 * 19).wrapping_add(seed) % nslots) as usize;
                    counts[slot] += 1;
                }
            }
        }
    }
    for (i, &c) in counts.iter().enumerate() {
        if c > 0 {
            fp[hash_position(i as u64, ncounts)] += 1;
        }
    }
}

/// The extended-connectivity relaxation shared by the scaffold bits.
fn relax_extcon(
    extcon: &mut [i32],
    extcon2: &mut [i32],
    atom_status: &[i32],
    bond_status: &[i32],
    nbp: &[Neighbourhood],
) {
    let n = extcon.len();
    for _ in 0..32 {
        extcon2.iter_mut().for_each(|e| *e = 0);
        for j in 0..n {
            if atom_status[j] <= 0 {
                continue;
            }
            let e2 = atom_status[j]
                .wrapping_mul(3)
                .wrapping_add(extcon[j].wrapping_mul(0xF));
            let mut sum: u64 = 0;
            for jj in 0..nbp[j].atoms.len() {
                if bond_status[nbp[j].bonds[jj]] <= 0 {
                    continue;
                }
                sum = sum.wrapping_add(extcon[nbp[j].atoms[jj]] as i64 as u64);
            }
            // `prod` starts at 0 and only gets multiplied: always 0.
            let v = ((e2 as i64 as u64).wrapping_add(sum.wrapping_mul(191)) << 8)
                .wrapping_add(((e2 & 0xFF0000) >> 16) as i64 as u64);
            extcon2[j] = (v as u32 as i32) & 0xFFFFFF;
        }
        extcon.copy_from_slice(extcon2);
    }
    loop {
        let mut changed = false;
        for j in 0..n {
            if atom_status[j] <= 0 {
                continue;
            }
            for jj in 0..nbp[j].atoms.len() {
                if bond_status[nbp[j].bonds[jj]] <= 0 {
                    continue;
                }
                let a = nbp[j].atoms[jj];
                if extcon[j] < extcon[a] {
                    changed = true;
                    extcon[a] = extcon[j];
                }
            }
        }
        if !changed {
            break;
        }
    }
}

/// `CountFingerprintPatterns(mp, fp_counts, ncounts, which_bits,
/// as_query = FALSE, fpflags, exclude_atom = 0, NULL)`; `dy_aromaticity`
/// is the `USE_DY_AROMATICITY` flag. Adds to `fp` (`ncounts = fp.len()`).
pub(crate) fn count_fingerprint_patterns(
    mol: &AvalonMolecule,
    fp: &mut [i32],
    which_bits: u32,
    dy_aromaticity: bool,
) {
    let n = mol.atoms.len();
    let nb = mol.bonds.len();
    let ncounts = fp.len();
    macro_rules! add_bit {
        ($seed:expr) => {
            fp[hash_position($seed, ncounts)] += 1
        };
    }
    macro_rules! add_bit_count {
        ($seed:expr, $count:expr) => {
            fp[hash_position($seed, ncounts)] += $count
        };
    }

    let nbp = setup_neighbourhood(mol);
    let mut touched = vec![0i32; n];
    let mut bond_type: Vec<i32> = mol.bonds.iter().map(|b| b.bond_type).collect();
    let atom_sym = |i: usize| mol.atoms[i].symbol.as_str();
    let a0 = |b: usize| (mol.bonds[b].atoms[0] - 1) as usize;
    let a1 = |b: usize| (mol.bonds[b].atoms[1] - 1) as usize;

    // hydrogen counts (implicit, from the input bond types, plus explicit H)
    let mut h_count = vec![0i32; n + 1];
    compute_implicit_h(mol, &bond_type, &mut h_count);
    for b in 0..nb {
        if atom_sym(a0(b)) == "H" {
            h_count[a1(b) + 1] += 1;
        } else if atom_sym(a1(b)) == "H" {
            h_count[a0(b) + 1] += 1;
        }
    }

    if dy_aromaticity {
        perceive_dy_aromaticity(mol, &nbp, &mut bond_type);
    } else {
        perceive_aromatic_bonds(mol, &mut bond_type);
    }
    let mut degree = vec![0i32; n];
    let mut cdegree = vec![0i32; n];
    let mut nspecial = vec![0i32; n];
    let mut unsaturated = vec![false; n];
    let (atom_status, bond_status) = ring_state(mol);
    let (atom_rsize, bond_rsize) = set_ring_size_flags(mol, 14, &nbp);

    let mut color = vec![0i32; n];
    let mut bcolor = vec![0i32; nb];
    let mut nrare_atoms = 0;
    for i in 0..n {
        let sym = atom_sym(i);
        let mut c = periodic_number(sym);
        if c <= 1 {
            c = 0;
        }
        if c > 115 {
            c = -1;
        }
        if sym == "A" {
            c = -1;
        }
        let is_rare = if sym == "L" {
            // Symbol lists never come from an RDKit molecule; an "L" atom
            // without a list leaves the previous value in C.
            c = -1;
            false
        } else {
            c > 0 && !symbol_in(sym, &["C", "H", "O", "N", "S", "P", "Cl", "F"])
        };
        if is_rare {
            nrare_atoms += 1;
        }
        if which_bits & USE_SHORTCUT_LABELS != 0 && sym == "R" && !mol.atoms[i].atext.is_empty() {
            c = ((0xFFFF00 & hash_string(&mol.atoms[i].atext)) | 119) as i32;
        }
        color[i] = c;
    }

    let mut ndouble = 0;
    let mut naromatic = 0;
    let mut nfusionb = 0;
    for b in 0..nb {
        let (x, y) = (a0(b), a1(b));
        bcolor[b] = match bond_type[b] {
            SINGLE => 1,
            DOUBLE => {
                ndouble += 1;
                2
            }
            TRIPLE => 3,
            AROMATIC => {
                naromatic += 1;
                4
            }
            _ => 0,
        };
        if bcolor[b] > 1 {
            unsaturated[x] = true;
            unsaturated[y] = true;
        }
        if color[x] != 0 && color[y] != 0 {
            degree[x] += 1;
            degree[y] += 1;
        }
        if bond_type[b] == DOUBLE {
            nspecial[x] += 1;
            nspecial[y] += 1;
        } else if bond_type[b] == TRIPLE {
            nspecial[x] += 2;
            nspecial[y] += 2;
        }
        if bond_type[b] != SINGLE {
            continue;
        }
        if color[x] != 0 && color[y] == 6 {
            cdegree[x] += 1;
        }
        if color[y] != 0 && color[x] == 6 {
            cdegree[y] += 1;
        }
        if atom_status[x] > 2 && atom_status[y] > 2 {
            nfusionb += 1;
        }
    }
    for c in color.iter_mut() {
        if *c < 0 {
            *c = 0;
        }
    }

    if which_bits & USE_ATOM_COUNT != 0 {
        let mut type_count_hash = [0i32; NCOUNT_HASH];
        // Upstream defines NCOUNT_SEED_HASH as the unparenthesized 128*128.
        // Thus seed%NCOUNT_SEED_HASH expands to (seed%128)*128. Preserve
        // that slot mapping for RDKit-compatible fingerprints above 2048 bits.
        let mut type_count_seed_hash = vec![0i32; NCOUNT_SEED_HASH];
        let (mut nringch2, mut nfusionch, mut nspiro) = (0, 0, 0);
        for i in 0..n {
            let c = color[i];
            if c == 0 {
                continue;
            }
            if c == 6 {
                if atom_status[i] > 0 && h_count[i + 1] >= 2 {
                    nringch2 += 1;
                }
                if atom_status[i] > 2 && h_count[i + 1] >= 1 {
                    nfusionch += 1;
                }
                if atom_status[i] > 3 {
                    nspiro += 1;
                }
                continue;
            }
            let mut hash: i32 = 0;
            let mut seed = ns(317 * ATOM_COUNT_SEED, 507);
            seed = ns(seed, i64::from(c) + 17);
            type_count_seed_hash[((seed % 128) * 128) as usize] += 1;
            if (c == 7 || c == 8) && h_count[i + 1] <= 0 {
                continue;
            }
            hash = hash.wrapping_mul(7).wrapping_add(c + 13);
            type_count_hash[(hash as usize) % NCOUNT_HASH] += 1;
            seed = ns(seed, i64::from(c) + 13);
            type_count_seed_hash[((seed % 128) * 128) as usize] += 1;
            if c != 7 && c != 8 {
                hash = hash.wrapping_mul(7).wrapping_add(c + 2 * 13);
                type_count_hash[(hash as usize) % NCOUNT_HASH] += 1;
                seed = ns(seed, i64::from(c) + 2 * 13);
                type_count_seed_hash[((seed % 128) * 128) as usize] += 1;
            }
        }
        if ncounts <= 2048 {
            for (i, &cnt) in type_count_hash.iter().enumerate() {
                let i = i as i64;
                if cnt > 0 {
                    add_bit_count!(ns((i * 19) as u64, 3), cnt);
                }
                if cnt > 1 {
                    add_bit_count!(ns((i * 23) as u64, 5), cnt);
                }
                if cnt > 2 {
                    add_bit_count!(ns((i * 29) as u64, 7), cnt);
                }
                if cnt > 4 {
                    add_bit_count!(ns((i * 31) as u64, 11), cnt);
                }
                if cnt > 8 {
                    add_bit_count!(ns((i * 37) as u64, 13), cnt);
                }
                if cnt > 16 {
                    add_bit_count!(ns((i * 41) as u64, 17), cnt);
                }
            }
        } else {
            for (i, &cnt) in type_count_seed_hash.iter().enumerate() {
                let i = i as i64;
                if cnt > 0 {
                    add_bit_count!(ns((i * 19) as u64, 3), cnt);
                }
                if cnt > 1 {
                    add_bit_count!(ns((i * 23) as u64, 5), cnt);
                }
                if cnt > 2 {
                    add_bit_count!(ns((i * 29) as u64, 7), cnt);
                }
                if cnt > 3 {
                    add_bit_count!(ns((i * 29) as u64, 71), cnt);
                }
                if cnt > 4 {
                    add_bit_count!(ns((i * 31) as u64, 11), cnt);
                }
                if cnt > 6 {
                    add_bit_count!(ns((i * 31) as u64, 113), cnt);
                }
                if cnt > 8 {
                    add_bit_count!(ns((i * 37) as u64, 13), cnt);
                }
                if cnt > 12 {
                    add_bit_count!(ns((i * 37) as u64, 133), cnt);
                }
                if cnt > 16 {
                    add_bit_count!(ns((i * 41) as u64, 17), cnt);
                }
            }
        }
        for (threshold, inc) in [(10, 341), (14, 441), (18, 541), (22, 641), (26, 741)] {
            if naromatic >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), naromatic);
            }
        }
        // The C code adds the aromatic-bond count for the double-bond bits.
        for (threshold, inc) in [(3, 371), (5, 471), (8, 571)] {
            if ndouble >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), naromatic);
            }
        }
        for (threshold, inc) in [(6, 411), (12, 511), (22, 611)] {
            if nringch2 >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), nringch2);
            }
        }
        for (threshold, inc) in [(2, 421), (4, 521)] {
            if nfusionch >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), nfusionch);
            }
        }
        for (threshold, inc) in [(1, 322), (1, 422), (2, 522), (2, 622)] {
            if nspiro >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), nspiro);
            }
        }
        for (threshold, inc) in [(1, 425), (2, 525), (3, 625), (5, 825)] {
            if nfusionb >= threshold {
                add_bit_count!(ns(ATOM_COUNT_SEED, inc), nfusionb);
            }
        }
        if nrare_atoms > 0 {
            let seed = ns(317 * ATOM_COUNT_SEED, 5);
            add_bit_count!(seed, nrare_atoms);
            add_bit_count!(ns(seed, 101), nrare_atoms);
            add_bit_count!(ns(seed, 103), nrare_atoms);
            if nrare_atoms > 1 {
                add_bit_count!(ns(seed, 203), nrare_atoms);
            }
            if nrare_atoms > 3 {
                add_bit_count!(ns(seed, 211), nrare_atoms);
            }
        }
    }

    if which_bits & USE_ATOM_SYMBOL_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = ATOM_SYMBOL_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            if c != 6 && c != 7 && c != 8 {
                add_bit!(seed);
            }
            if c != 6 {
                walk.path_bits(
                    fp,
                    seed.wrapping_add(12347),
                    &mut touched,
                    1,
                    1,
                    1,
                    i,
                    None,
                    FORCED_HETERO_END | PROCESS_CHAINS,
                );
            }
            if c >= 10
                || (c == 7 && degree[i] > 0)
                || (c == 8 && degree[i] > 1)
                || (c == 6 && degree[i] > 2 && atom_status[i] > 0)
            {
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    1,
                    2,
                    i,
                    None,
                    STOP_AT_HEAVY_ATOM | FORCED_HETERO_END | PROCESS_CHAINS,
                );
            }
            if c > 6 {
                walk.path_bits(
                    fp,
                    ((217 * c) as i64 as u64).wrapping_add(seed),
                    &mut touched,
                    1,
                    3,
                    4,
                    i,
                    None,
                    IGNORE_PATH_SYMBOL
                        | PROCESS_RING_CLOSURES
                        | STOP_AT_HEAVY_ATOM
                        | PROCESS_CHAINS,
                );
                if c > 10 && c <= 18 {
                    walk.path_bits(
                        fp,
                        17u64.wrapping_add(seed),
                        &mut touched,
                        1,
                        5,
                        7,
                        i,
                        None,
                        FORCED_HETERO_END
                            | IGNORE_PATH_SYMBOL
                            | STOP_AT_HEAVY_ATOM
                            | PROCESS_CHAINS,
                    );
                }
            }
            if (c == 7 || c == 8) && cdegree[i] > 2 {
                seed = ns(seed, i64::from(c) * 23);
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    2,
                    6,
                    i,
                    None,
                    STOP_AT_HEAVY_ATOM | IGNORE_TERM_SYMBOL | PROCESS_CHAINS,
                );
            }
            seed = old_seed;
            if degree[i] >= 4 && c != 5 && c < 18 {
                seed = ns(seed, i64::from(c));
                walk.path_bits(
                    fp,
                    ns(seed, 3 * 107),
                    &mut touched,
                    1,
                    2,
                    4,
                    i,
                    None,
                    STOP_AT_HEAVY_ATOM | IGNORE_TERM_SYMBOL | PROCESS_CHAINS,
                );
            }
            seed = old_seed;
            if atom_status[i] >= 4 && c != 5 && c < 18 {
                seed = ns(seed, i64::from(c) + 55);
                walk.path_bits(
                    fp,
                    ns(seed, 3 * 109),
                    &mut touched,
                    1,
                    2,
                    4,
                    i,
                    None,
                    STOP_AT_HEAVY_ATOM
                        | IGNORE_TERM_SYMBOL
                        | IGNORE_PATH_SYMBOL
                        | PROCESS_RING_CLOSURES
                        | PROCESS_CHAINS,
                );
            }
            // paths starting with a double or triple bond
            for j in 0..nbp[i].atoms.len() {
                let b = nbp[i].bonds[j];
                let ai = nbp[i].atoms[j];
                if bcolor[b] == 0 {
                    continue;
                }
                if bond_type[b] != DOUBLE && bond_type[b] != TRIPLE {
                    continue;
                }
                let bond_seed = ns(ns(old_seed, i64::from(c)), i64::from(bcolor[b]) * 613);
                touched[ai] = 1;
                walk.path_bits(
                    fp,
                    bond_seed,
                    &mut touched,
                    2,
                    5,
                    5,
                    ai,
                    Some(i),
                    STOP_AT_HEAVY_ATOM | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
                );
                touched[ai] = 0;
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    if which_bits & USE_AUGMENTED_ATOM != 0 {
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            if nspecial[i] >= 2 {
                add_bit!(ns(i64::from(c) as u64 * AUGMENTED_ATOM_SEED, 101));
            }
            let nbi = &nbp[i];
            let deg = nbi.atoms.len();
            if (h_count[i + 1] > 0 || c != 6) && degree[i] >= 2 {
                for i1 in 0..deg {
                    let (x1, y1) = (nbi.atoms[i1], nbi.bonds[i1]);
                    if color[x1] == 0 || bcolor[y1] == 0 {
                        continue;
                    }
                    for i2 in i1 + 1..deg {
                        let (x2, y2) = (nbi.atoms[i2], nbi.bonds[i2]);
                        let mut seed = ns(AUGMENTED_ATOM_SEED, 97);
                        seed = ns(seed, i64::from(c));
                        if color[x2] == 0 || bcolor[y2] == 0 {
                            continue;
                        }
                        if color[x1] == 6 && color[x2] == 6 {
                            continue;
                        }
                        let sum = (color[x1] * bcolor[y1] + color[x2] * bcolor[y2]) as i64 as u64;
                        let mut prod: u64 = 1;
                        prod = prod.wrapping_mul(color[x1] as i64 as u64) & 0xFFF;
                        prod = prod.wrapping_mul(color[x2] as i64 as u64) & 0xFFF;
                        seed = next_hash(seed, sum);
                        seed = next_hash(seed, prod);
                        add_bit!(seed);
                        let mut nmulti = 0;
                        if bcolor[y1] >= 2 {
                            nmulti += 1;
                        }
                        if bcolor[y2] >= 2 {
                            nmulti += 1;
                        }
                        if c != 6 && nmulti > 0 {
                            add_bit!(ns(seed, i64::from(nmulti + 173 * degree[i])));
                        }
                    }
                }
            }
            if degree[i] <= 2 {
                continue;
            }
            for i1 in 0..deg {
                for i2 in i1 + 1..deg {
                    for i3 in i2 + 1..deg {
                        let mut seed = AUGMENTED_ATOM_SEED;
                        seed = ns(seed, i64::from(c));
                        let xs = [nbi.atoms[i1], nbi.atoms[i2], nbi.atoms[i3]];
                        let ys = [nbi.bonds[i1], nbi.bonds[i2], nbi.bonds[i3]];
                        if xs.iter().any(|&x| color[x] == 0) || ys.iter().any(|&y| bcolor[y] == 0) {
                            continue;
                        }
                        let nqtmp = xs.iter().filter(|&&x| color[x] != 6).count();
                        let nmulti = ys
                            .iter()
                            .filter(|&&y| bcolor[y] == 2 || bcolor[y] == 3)
                            .count();
                        let mut sum: u64 = 0;
                        let mut prod: u64 = 1;
                        for k in 0..3 {
                            sum = sum.wrapping_add((color[xs[k]] * bcolor[ys[k]]) as i64 as u64);
                        }
                        for &x in &xs {
                            prod = prod.wrapping_mul(color[x] as i64 as u64) & 0xFFF;
                        }
                        seed = next_hash(seed, sum);
                        seed = next_hash(seed, prod);
                        if (nqtmp > 2 || nmulti >= 2) && c == 6 {
                            add_bit!(seed);
                            add_bit!(ns(seed, 73));
                        }
                        seed = next_hash(seed, sum.wrapping_mul(prod) & 0xFFF);
                        add_bit!(seed);
                        if nmulti >= 2 || c > 6 {
                            add_bit!(ns(seed, 53));
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_AUGMENTED_BOND != 0 {
        for b in 0..nb {
            let (x, y) = (a0(b), a1(b));
            if degree[x] <= 2 || degree[y] <= 2 {
                continue;
            }
            let side = |a: usize, k1: usize, k2: usize| -> Option<(u64, u64)> {
                let nba = &nbp[a];
                if nba.bonds[k1] == b || nba.bonds[k2] == b {
                    return None;
                }
                let (u1, v1) = (nba.atoms[k1], nba.bonds[k1]);
                let (u2, v2) = (nba.atoms[k2], nba.bonds[k2]);
                if color[u1] == 0 || color[u2] == 0 || bcolor[v1] == 0 || bcolor[v2] == 0 {
                    return None;
                }
                let sum =
                    ((color[u1] * bcolor[v1] + color[u2] * bcolor[v2]) as i64 as u64) & 0x0FFF;
                let prod = (((color[u1] + bcolor[v1]) as i64 as u64)
                    .wrapping_mul((color[u2] + bcolor[v2]) as i64 as u64))
                    & 0x0FFF;
                Some((sum, prod))
            };
            for i1 in 0..nbp[x].atoms.len() {
                for i2 in i1 + 1..nbp[x].atoms.len() {
                    let Some((sumi, prodi)) = side(x, i1, i2) else {
                        continue;
                    };
                    for j1 in 0..nbp[y].atoms.len() {
                        for j2 in j1 + 1..nbp[y].atoms.len() {
                            let Some((sumj, prodj)) = side(y, j1, j2) else {
                                continue;
                            };
                            let mut seed = AUGMENTED_BOND_SEED;
                            seed = ns(seed, i64::from(bcolor[b]));
                            seed = next_hash(seed, prodi + prodj);
                            seed = next_hash(seed, sumi * sumj);
                            add_bit!(seed);
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_HCOUNT_PAIR != 0 {
        for b in 0..nb {
            let (x, y) = (a0(b), a1(b));
            let (hx, hy) = (h_count[x + 1], h_count[y + 1]);
            if hx == 0 && hy == 0 {
                continue;
            }
            if color[x] == 6 && color[y] == 6 && bcolor[b] == 1 {
                continue;
            }
            if color[x] == 0 || color[y] == 0 {
                continue;
            }
            for j1 in 0..=hx {
                for j2 in 0..=hy {
                    if j1 + j2 == 0 {
                        continue;
                    }
                    if !unsaturated[x] && !unsaturated[y] && j1 * j2 == 0 && bcolor[b] <= 1 {
                        continue;
                    }
                    if j1 + j2 > 3 {
                        continue;
                    }
                    let mut seed = ns(HCOUNT_PAIR_SEED, i64::from(7 * (j1 + 1) * (j2 + 1)));
                    seed = ns(seed, i64::from(53 * (j1 + j2)));
                    seed = ns(seed, i64::from(bcolor[b]));
                    seed = ns(seed, i64::from(color[x] + color[y]));
                    seed = ns(seed, i64::from(color[x] * color[y]));
                    add_bit!(seed);
                    if bcolor[b] > 1 {
                        seed = ns(seed, 83);
                        add_bit!(seed);
                        if j1 + j2 == 1 && (color[x] != 6 || color[y] != 6) && bcolor[b] <= 3 {
                            seed = ns(seed, 91);
                            add_bit!(seed);
                            seed = ns(seed, 97);
                            add_bit!(seed);
                            seed = ns(seed, 103);
                            add_bit!(seed);
                        }
                    }
                }
            }
        }
    }

    if which_bits & USE_HCOUNT_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = HCOUNT_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            let h = h_count[i + 1];
            if h == 0 {
                continue;
            }
            if c == 6 && h < 3 && degree[i] < 3 {
                continue;
            }
            touched[i] = 1;
            let mut old_seed = seed;
            seed = ns(seed, i64::from(c));
            if c != 6 {
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    2,
                    5,
                    i,
                    None,
                    IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
                );
            } else {
                if degree[i] > 2 {
                    walk.path_bits(
                        fp,
                        ns(seed, 103),
                        &mut touched,
                        1,
                        2,
                        5,
                        i,
                        None,
                        IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
                    );
                }
                if h >= 3 {
                    walk.path_bits(
                        fp,
                        ns(seed, 1105),
                        &mut touched,
                        1,
                        3,
                        6,
                        i,
                        None,
                        IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
                    );
                }
            }
            touched[i] = 0;
            seed = old_seed;
            if c == 6 {
                continue;
            }
            if h > 1 {
                touched[i] = 1;
                old_seed = seed;
                seed = ns(seed, 113);
                seed = ns(seed, i64::from(c));
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    1,
                    5,
                    i,
                    None,
                    IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
                );
                touched[i] = 0;
                seed = old_seed;
            }
            for j in 1..h {
                seed = ns(seed, i64::from(61 * j));
                seed = ns(seed, i64::from(c));
                add_bit!(seed);
            }
            seed = old_seed;
        }
    }

    // ring paths: only ring atoms and bonds touching ring atoms
    for i in 0..n {
        if atom_status[i] <= 0 {
            color[i] = 0;
        }
    }
    for b in 0..nb {
        if bond_status[b] <= 0 && atom_status[a0(b)] == 0 && atom_status[a1(b)] == 0 {
            bcolor[b] = 0;
        }
    }

    if which_bits & USE_RING_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = RING_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            if c > 5 && c < 10 && atom_status[i] > 2 {
                seed = ns(seed, 61);
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    3,
                    8,
                    i,
                    None,
                    STOP_AT_HEAVY_ATOM | PROCESS_RING_CLOSURES,
                );
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    // all non-hydrogen atoms alike, bond types kept
    for i in 0..n {
        color[i] = if periodic_number(atom_sym(i)) == 1 {
            0
        } else {
            ANY_COLOR
        };
    }
    for b in 0..nb {
        bcolor[b] = match bond_type[b] {
            SINGLE => 1,
            DOUBLE => 2,
            TRIPLE => 3,
            AROMATIC => 4,
            _ => 0,
        };
    }

    if which_bits & USE_BOND_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = BOND_PATH_SEED;
        for i in 0..n {
            if color[i] <= 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            if degree[i] > 2 && atom_status[i] > 1 {
                walk.path_bits(fp, seed, &mut touched, 1, 4, 4, i, None, PROCESS_CHAINS);
                if atom_rsize[i] & SPECIAL_RING != 0 {
                    walk.path_bits(
                        fp,
                        ns(seed, 217),
                        &mut touched,
                        1,
                        5,
                        5,
                        i,
                        None,
                        IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
                    );
                }
            }
            seed = old_seed;
            seed = ns(seed, 11);
            walk.path_bits(
                fp,
                seed,
                &mut touched,
                1,
                4,
                6,
                i,
                None,
                FORCED_RING_PATH | IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
            );
            for j in 0..nbp[i].atoms.len() {
                let b = nbp[i].bonds[j];
                let ai = nbp[i].atoms[j];
                if bcolor[b] == 0 {
                    continue;
                }
                if bond_type[b] != DOUBLE && bond_type[b] != TRIPLE {
                    continue;
                }
                if atom_status[i] <= 0 && bond_type[b] != TRIPLE {
                    continue;
                }
                seed = old_seed;
                seed = ns(seed, i64::from(bcolor[b]) * 413);
                touched[ai] = 1;
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    2,
                    4,
                    5,
                    ai,
                    Some(i),
                    IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
                );
                touched[ai] = 0;
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    // atom classes: carbon / hetero
    for i in 0..n {
        let sym = atom_sym(i);
        let mut c = match sym {
            // symbol lists never come from an RDKit molecule
            "L" => 0,
            "H" | "D" | "T" | "A" => 0,
            "C" => 6,
            "N" | "O" | "S" | "Q" => 8,
            s => {
                let t = periodic_number(s);
                if 1 < t && t < 115 { 8 } else { 0 }
            }
        };
        if c > 115 {
            c = 0;
        }
        color[i] = c;
    }

    if which_bits & USE_HCOUNT_CLASS_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = HCOUNT_CLASS_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            let h = h_count[i + 1];
            if h == 0 {
                continue;
            }
            if c == 6 && h < 2 {
                continue;
            }
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            touched[i] = 1;
            if c == 6 {
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    2,
                    4,
                    i,
                    None,
                    FORCED_HETERO_END | PROCESS_CHAINS,
                );
            } else {
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    2,
                    5,
                    i,
                    None,
                    IGNORE_PATH_SYMBOL | FORCED_HETERO_END | PROCESS_CHAINS,
                );
            }
            seed = old_seed;
            touched[i] = 0;
        }
    }

    if which_bits & USE_ATOM_CLASS_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = ns(ATOM_CLASS_PATH_SEED, 117);
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            if c == 6 && degree[i] < 3 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            if c != 6 {
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    3,
                    4,
                    i,
                    None,
                    FORCED_RING_PATH | PROCESS_RING_CLOSURES | PROCESS_CHAINS,
                );
            }
            touched[i] = 0;
            seed = old_seed;
        }
    }

    // a single bond class
    for b in 0..nb {
        bcolor[b] = if (SINGLE..=ANY_BOND).contains(&bond_type[b]) {
            5
        } else {
            0
        };
    }

    if which_bits & USE_ATOM_CLASS_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = ATOM_CLASS_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 || c == 6 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            seed = ns(seed, 23 + i64::from(c) * 19);
            walk.path_bits(
                fp,
                seed,
                &mut touched,
                1,
                3,
                9,
                i,
                None,
                IGNORE_PATH_SYMBOL | PROCESS_RING_CLOSURES,
            );
            seed = old_seed;
            touched[i] = 0;
        }
        let mut qq_count = 0;
        let mut qc_count = 0;
        for b in 0..nb {
            if bond_status[b] == 0 || bcolor[b] == 0 {
                continue;
            }
            let (c1, c2) = (color[a0(b)], color[a1(b)]);
            if c1 == 0 || c2 == 0 || (c1 == 6 && c2 == 6) {
                continue;
            }
            let base = if c1 != 6 && c2 != 6 {
                qq_count += 1;
                ATOM_CLASS_PATH_SEED * 17
            } else {
                qc_count += 1;
                ATOM_CLASS_PATH_SEED * 19
            };
            for j in 3..9 {
                if bond_rsize[b] & (1 << j) != 0 {
                    add_bit!(ns(base, j * 8));
                }
            }
        }
        let mut seed = 2 * ATOM_CLASS_PATH_SEED + 3;
        let mut i: i32 = 1;
        while i <= qq_count {
            seed = ns(seed, i64::from(i) * 153);
            add_bit!(seed);
            seed = ns(seed, 53);
            if i <= 1 {
                add_bit!(seed);
            }
            i = (1.0 + f64::from(i) * 1.5) as i32;
        }
        seed = 3 * ATOM_CLASS_PATH_SEED + 5;
        for i in 1..=qc_count.min(2) {
            seed = ns(seed, i64::from(i) * 157);
            add_bit!(seed);
        }
        let mut i: i32 = 3;
        while i <= qc_count {
            seed = ns(seed, i64::from(i) * 157);
            add_bit!(seed);
            i = (f64::from(i) * 1.8) as i32;
        }
    }

    // ring patterns: bonds with at least one ring atom
    for b in 0..nb {
        bcolor[b] = 5;
        let (x, y) = (a0(b), a1(b));
        if atom_status[x] == 0 && atom_status[y] == 0 {
            bcolor[b] = 0;
        } else {
            for a in [x, y] {
                if color[a] != 0 {
                    if color[a] != 6 {
                        color[a] = 8;
                    }
                    if atom_rsize[a] == 0 {
                        color[a] = 0;
                    }
                }
            }
        }
    }

    if which_bits & USE_RING_PATTERN != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = RING_PATTERN_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            if c == 6 && atom_rsize[i] & SPECIAL_RING == 0 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            walk.path_bits(fp, seed, &mut touched, 1, 3, 3, i, None, PROCESS_CHAINS);
            seed = old_seed;
            touched[i] = 0;
        }
        // The colour changes made here for the (disabled) complete-ring
        // classes are overwritten below before anything reads them.
    }

    if which_bits & USE_RING_SIZE_COUNTS != 0 {
        for j in 3i64..10 {
            let nrbonds = (0..nb).filter(|&b| bond_rsize[b] & (1 << j) != 0).count() as i64;
            let mut seed = ns(RING_SIZE_SEED, j * 13);
            let mut i: i64 = 1;
            while i < 100 {
                if nrbonds >= j * i {
                    seed = ns(seed, i);
                    if (j != 6 && j != 5) || nrbonds > j * i {
                        add_bit_count!(seed, nrbonds as i32);
                    }
                    if j != 6 && j != 5 {
                        seed = ns(seed, 17);
                        add_bit_count!(seed, nrbonds as i32);
                    }
                } else {
                    break;
                }
                i *= 2;
            }
        }
        let mut rscounts = [[0i32; 15]; 15];
        for b in 0..nb {
            let (rx, ry) = (atom_rsize[a0(b)], atom_rsize[a1(b)]);
            for j in 3..15 {
                for k in 3..15 {
                    if j == k {
                        continue;
                    }
                    if rx & (1 << j) != 0 && ry & (1 << k) != 0 {
                        rscounts[j][k] += 1;
                        rscounts[k][j] += 1;
                    }
                }
            }
        }
        for j in 3usize..9 {
            for k in j + 1..9 {
                if rscounts[j][k] == 0 {
                    continue;
                }
                let (jj, kk) = (j as i64, k as i64);
                let mut seed = 2 * RING_SIZE_SEED;
                seed = ns(seed, (jj + kk) * 11);
                seed = ns(seed, jj * kk * 13);
                seed = ns(seed, 19);
                add_bit!(seed);
                if rscounts[j][k] == 1 || j == 6 || k == 6 {
                    continue;
                }
                let mut seed = 2 * RING_SIZE_SEED + 23;
                seed = ns(seed, (jj + kk) * 11);
                seed = ns(seed, jj * kk * 13);
                seed = ns(seed, 19);
                add_bit!(seed);
            }
        }
    }

    // degree colours
    for i in 0..n {
        let sym = atom_sym(i);
        let mut c = periodic_number(sym);
        if c <= 1 {
            c = 0;
        }
        if c > 115 {
            c = -1;
        }
        if sym == "A" {
            c = -1;
        }
        color[i] = if c > 1 { c + 32 * degree[i] } else { 0 };
    }
    for b in 0..nb {
        bcolor[b] = if (SINGLE..=ANY_BOND).contains(&bond_type[b]) {
            5
        } else {
            0
        };
    }

    if which_bits & USE_DEGREE_PATH != 0 {
        let walk = Walk {
            nbp: &nbp,
            atom_color: &color,
            bond_color: &bcolor,
            bond_rsize: &bond_rsize,
        };
        let mut seed = DEGREE_PATH_SEED;
        for i in 0..n {
            let c = color[i];
            if c <= 0 {
                continue;
            }
            if degree[i] <= 3 && atom_status[i] <= 2 && degree[i] != 1 {
                continue;
            }
            if degree[i] == 1 && atom_sym(i) != "C" {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(c));
            walk.path_bits(
                fp,
                seed,
                &mut touched,
                1,
                2,
                4,
                i,
                None,
                IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
            );
            if atom_status[i] > 2 && h_count[i + 1] >= 1 {
                seed = ns(seed, 219);
                walk.path_bits(
                    fp,
                    seed,
                    &mut touched,
                    1,
                    2,
                    5,
                    i,
                    None,
                    IGNORE_PATH_SYMBOL | PROCESS_CHAINS,
                );
            }
            seed = old_seed;
            touched[i] = 0;
        }
        for i in 0..n {
            let t = periodic_number(atom_sym(i));
            if 1 >= t || t >= 115 || t == 6 || t < 10 {
                continue;
            }
            touched[i] = 1;
            let old_seed = seed;
            seed = ns(seed, i64::from(t));
            walk.path_bits(fp, seed, &mut touched, 1, 2, 2, i, None, PROCESS_CHAINS);
            seed = old_seed;
            touched[i] = 0;
        }
    }

    let need_lengths = which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS) != 0;
    let mut lm: Vec<Vec<i32>> = Vec::new();
    if need_lengths {
        lm = vec![vec![0i32; n]; n];
        touched.iter_mut().for_each(|t| *t = 0);
        for i in 0..n {
            touched[i] = 1;
            path_length_flags(&color, &nbp, &mut touched, i, 0, i, 12, &mut lm);
            touched[i] = 0;
        }
    }

    if which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS) != 0 {
        for i in 0..n {
            let sym = atom_sym(i);
            color[i] = match sym {
                "H" | "D" | "T" => 0,
                "Q" => HETERO,
                "A" | "L" => GENERIC,
                "C" => {
                    if cdegree[i] >= 3 {
                        CSP3
                    } else {
                        6
                    }
                }
                s => {
                    let c = periodic_number(s);
                    if c > 1 && c < 115 { HETERO } else { 0 }
                }
            };
        }
        for i in 0..n {
            if degree[i] < 3 {
                continue;
            }
            let c = color[i];
            if c != CSP3 && c != 6 {
                continue;
            }
            touched[i] = 1;
            let mut csp3 = [0i32; MAX_SPIDER + 1];
            let mut hetero = [0i32; MAX_SPIDER + 1];
            if which_bits & USE_CLASS_SPIDERS != 0 {
                special_neighbours(
                    &color,
                    &nbp,
                    &mut touched,
                    1,
                    i,
                    MAX_SPIDER,
                    &mut csp3,
                    &mut hetero,
                );
            }
            touched[i] = 0;
            let centre = |second: i64| -> u64 {
                let seed = CLASS_SPIDER_SEED;
                if c == HETERO {
                    ns(ns(seed, i64::from(HETERO) * 8), second)
                } else {
                    ns(ns(seed, 6 * 8), second)
                }
            };
            if which_bits & USE_CLASS_SPIDERS != 0 {
                for j in 1..=MAX_SPIDER {
                    if csp3[j] == 0 {
                        continue;
                    }
                    let seed = centre(i64::from(CSP3) * 11);
                    for j1 in 1..=MAX_SPIDER {
                        if hetero[j1] <= 0 {
                            continue;
                        }
                        for j2 in j1..=MAX_SPIDER {
                            let mut tmp2 = hetero[j2];
                            if j2 == j1 {
                                tmp2 -= 1;
                            }
                            if tmp2 <= 0 {
                                continue;
                            }
                            let s = ns(seed, j as i64);
                            let s = ns(s, (j1 + j2) as i64);
                            add_bit!(s);
                        }
                    }
                }
            }
            if c != CSP3 {
                continue;
            }
            if which_bits & USE_CLASS_SPIDERS != 0 {
                for j in 1..=MAX_SPIDER {
                    if hetero[j] == 0 {
                        continue;
                    }
                    let seed = centre(i64::from(HETERO) * 11);
                    for j1 in j..=MAX_SPIDER {
                        let mut tmp1 = hetero[j1];
                        if j1 == j {
                            tmp1 -= 1;
                        }
                        if tmp1 <= 0 {
                            continue;
                        }
                        for j2 in j1..=MAX_SPIDER {
                            let mut tmp2 = hetero[j2];
                            if j2 == j {
                                tmp2 -= 1;
                            }
                            if j2 == j1 {
                                tmp2 -= 1;
                            }
                            if tmp2 <= 0 {
                                continue;
                            }
                            let s = ns(seed, (j * j1 * j2) as i64);
                            add_bit!(s);
                            if degree[i] > 3 {
                                add_bit!(ns(s, 501));
                            }
                        }
                    }
                }
            }
        }

        if which_bits & USE_FEATURE_PAIRS != 0 {
            for i in 0..n {
                let sym = atom_sym(i);
                let mut flags = match sym {
                    "C" => C_FLAG,
                    "O" => O_FLAG,
                    "N" => N_FLAG,
                    "S" => S_FLAG,
                    "P" => P_FLAG,
                    "F" | "Cl" | "Br" | "I" | "At" => X_FLAG,
                    _ => 0,
                };
                if color[i] == HETERO {
                    flags |= HETERO_FLAG;
                }
                if cdegree[i] >= 3 {
                    flags |= CSP3_FLAG;
                }
                if degree[i] >= 4 {
                    flags |= QUART_FLAG;
                }
                if atom_status[i] > 0 && degree[i] >= 3 {
                    flags |= RING_SUBST_FLAG;
                    if atom_rsize[i] & SPECIAL_RING != 0 {
                        flags |= RS_SPECIAL_FLAG;
                    }
                }
                color[i] = flags;
            }
            let _ = RS_SPECIAL_FLAG;
            feature_bits(
                &color,
                fp,
                RING_SUBST_FLAG,
                RING_SUBST_FLAG,
                5,
                7,
                &lm,
                2237,
            );
            feature_bits(&color, fp, QUART_FLAG, HETERO_FLAG, 1, 8, &lm, 4237);
            feature_bits(&color, fp, QUART_FLAG, RING_SUBST_FLAG, 1, 6, &lm, 5237);
            feature_bits(&color, fp, X_FLAG, CSP3_FLAG, 1, 1, &lm, 15237);
        }
    }

    if which_bits & USE_NON_SSS_BITS != 0 {
        let seed = NON_SSS_SEED;
        let mut extcon = vec![0i32; n];
        let mut extcon2 = vec![0i32; n];
        for j in 0..n {
            if atom_status[j] > 0 {
                extcon[j] = atom_rsize[j] as i32;
            }
        }
        relax_extcon(&mut extcon, &mut extcon2, &atom_status, &bond_status, &nbp);
        for j in 0..n {
            if atom_status[j] <= 0 {
                continue;
            }
            if which_bits & USE_SCAFFOLD_IDS != 0 {
                // RDKit's compiled Avalon hashes the full positive product.
                // Truncating it to signed i32 changes non-power-of-two sizes.
                add_bit!(ns(seed, i64::from(extcon[j]) * 10013));
            }
            if which_bits & USE_SCAFFOLD_LINKS != 0 {
                if degree[j] <= atom_status[j] {
                    continue;
                }
                for jj in 0..n {
                    if atom_status[jj] <= 0 || j == jj || degree[jj] <= atom_status[jj] {
                        continue;
                    }
                    let Some(k) = (0..12).find(|&k| lm[j][jj] == 1 << k) else {
                        continue;
                    };
                    if k >= 3 {
                        continue;
                    }
                    let k = i64::from(k);
                    let v1 = i64::from(extcon[j]) * 1013 + i64::from(extcon[jj]) * 2003;
                    add_bit!(ns(ns(seed, 3 * k), v1));
                    let v2 = i64::from(extcon2[j]) * 2013 + i64::from(extcon[jj]) * 1003;
                    add_bit!(ns(ns(seed, 5 * k), v2));
                    let v3 = i64::from(extcon[j]) * 2013 + i64::from(extcon2[jj]) * 1003;
                    add_bit!(ns(ns(seed, 5 * k), v3));
                    let v4 = i64::from(extcon2[j]) * 3013 + i64::from(extcon2[jj]) * 3003;
                    add_bit!(ns(ns(seed, 7 * k), v4));
                }
            }
        }
        if which_bits & USE_SCAFFOLD_COLORS != 0 {
            for j in 0..n {
                if atom_status[j] <= 0 {
                    continue;
                }
                let tmp1 = match atom_sym(j) {
                    "C" => 101,
                    "O" => 301,
                    "N" => 401,
                    "S" => 601,
                    "P" => 701,
                    "F" | "Cl" | "Br" | "I" | "At" => 901,
                    _ => 0,
                };
                extcon[j] = atom_rsize[j] as i32 + tmp1;
            }
            relax_extcon(&mut extcon, &mut extcon2, &atom_status, &bond_status, &nbp);
            for &e in &extcon {
                add_bit!(ns(ns(seed, 4), i64::from(e) * 3013));
            }
        }
    }
}
