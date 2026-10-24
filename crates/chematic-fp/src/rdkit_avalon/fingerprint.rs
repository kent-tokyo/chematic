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

include!("fingerprint/count_patterns.rs");
