//! Non-tetrahedral stereo classes (`@SP`, `@TB`, `@OH`) and their
//! permutation algebra, as RDKit 2026.03.1 defines it
//! (`GraphMol/NontetrahedralStereo.cpp`): how a permutation number changes
//! when two ligand positions swap, which ligand sits across from which, and
//! the axial positions of a trigonal bipyramid.
//!
//! A permutation number refers to an ordered ligand list; position `i` is
//! the `i`-th ligand in that order (implicit or missing ligands take their
//! own positions).

/// A non-tetrahedral stereo class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NonTetrahedralClass {
    /// `@SP1`..`@SP3`: four ligands in a plane.
    SquarePlanar,
    /// `@TB1`..`@TB20`: five ligands, two axial.
    TrigonalBipyramidal,
    /// `@OH1`..`@OH30`: six ligands.
    Octahedral,
}

impl NonTetrahedralClass {
    /// The number of ligand positions (`Chirality::getMaxNbors`).
    pub fn max_nbors(self) -> usize {
        match self {
            Self::SquarePlanar => 4,
            Self::TrigonalBipyramidal => 5,
            Self::Octahedral => 6,
        }
    }

    /// The largest permutation number.
    pub fn max_permutation(self) -> u32 {
        match self {
            Self::SquarePlanar => 3,
            Self::TrigonalBipyramidal => 20,
            Self::Octahedral => 30,
        }
    }

    /// The SMILES class token (`SP`, `TB`, `OH`).
    pub fn token(self) -> &'static str {
        match self {
            Self::SquarePlanar => "SP",
            Self::TrigonalBipyramidal => "TB",
            Self::Octahedral => "OH",
        }
    }
}

const SWAP_SP: [[u8; 6]; 4] = [
    [0, 0, 0, 0, 0, 0],
    [3, 1, 2, 2, 1, 3],
    [2, 3, 1, 1, 3, 2],
    [1, 2, 3, 3, 2, 1],
];

const SWAP_TB: [[u8; 10]; 21] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 20, 17, 2, 2, 2, 7, 2, 6, 3],
    [11, 15, 18, 1, 1, 1, 8, 1, 5, 4],
    [10, 19, 4, 18, 4, 8, 4, 5, 4, 1],
    [12, 16, 3, 17, 3, 7, 3, 6, 3, 2],
    [13, 6, 16, 20, 7, 6, 6, 3, 2, 6],
    [14, 5, 19, 15, 8, 5, 5, 4, 1, 5],
    [8, 14, 10, 11, 5, 4, 1, 8, 8, 8],
    [7, 13, 12, 9, 6, 3, 2, 7, 7, 7],
    [1, 11, 11, 8, 15, 18, 11, 11, 14, 10],
    [3, 12, 7, 12, 16, 12, 17, 13, 12, 9],
    [2, 9, 9, 7, 20, 17, 9, 9, 13, 12],
    [4, 10, 8, 10, 19, 10, 18, 14, 10, 11],
    [5, 8, 14, 14, 14, 19, 15, 10, 11, 14],
    [6, 7, 13, 13, 13, 16, 20, 12, 9, 13],
    [20, 2, 20, 6, 9, 20, 13, 17, 20, 16],
    [19, 4, 5, 19, 10, 14, 19, 19, 18, 15],
    [18, 18, 1, 4, 18, 11, 10, 15, 19, 18],
    [17, 17, 2, 3, 17, 9, 12, 20, 16, 17],
    [16, 3, 6, 16, 12, 13, 16, 16, 17, 20],
    [15, 1, 15, 5, 11, 15, 14, 18, 15, 19],
];

const SWAP_OH: [[u8; 15]; 31] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3],
    [7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16],
    [18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1],
    [15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5],
    [14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4],
    [16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7],
    [2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6],
    [23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9],
    [24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8],
    [20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11],
    [19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10],
    [22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13],
    [21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12],
    [5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15],
    [4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14],
    [6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2],
    [1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18],
    [3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17],
    [11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20],
    [10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19],
    [13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22],
    [12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21],
    [8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24],
    [9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23],
    [30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26],
    [29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25],
    [28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28],
    [27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27],
    [26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30],
    [25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29],
];

const ACROSS_SP: [[u8; 4]; 4] = [[4, 4, 4, 4], [2, 3, 0, 1], [1, 0, 3, 2], [3, 2, 1, 0]];

const ACROSS_TB: [[u8; 5]; 21] = [
    [5, 5, 5, 5, 5],
    [4, 5, 5, 5, 0],
    [4, 5, 5, 5, 0],
    [3, 5, 5, 0, 5],
    [3, 5, 5, 0, 5],
    [2, 5, 0, 5, 5],
    [2, 5, 0, 5, 5],
    [1, 0, 5, 5, 5],
    [1, 0, 5, 5, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 4, 5, 5, 1],
    [5, 3, 5, 1, 5],
    [5, 2, 1, 5, 5],
    [5, 2, 1, 5, 5],
    [5, 5, 4, 5, 2],
    [5, 5, 3, 2, 5],
    [5, 5, 5, 4, 3],
    [5, 5, 5, 4, 3],
    [5, 5, 3, 2, 5],
    [5, 5, 4, 5, 2],
];

const ACROSS_OH: [[u8; 6]; 31] = [
    [6, 6, 6, 6, 6, 6],
    [5, 3, 4, 1, 2, 0],
    [5, 3, 4, 1, 2, 0],
    [4, 3, 5, 1, 0, 2],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [3, 4, 5, 0, 1, 2],
    [3, 5, 4, 0, 2, 1],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [5, 2, 1, 4, 3, 0],
    [4, 2, 1, 5, 0, 3],
    [3, 2, 1, 0, 5, 4],
    [3, 2, 1, 0, 5, 4],
    [5, 4, 3, 2, 1, 0],
    [4, 5, 3, 2, 0, 1],
    [4, 3, 5, 1, 0, 2],
    [3, 5, 4, 0, 2, 1],
    [3, 4, 5, 0, 1, 2],
    [2, 4, 0, 5, 1, 3],
    [2, 5, 0, 4, 3, 1],
    [2, 3, 0, 1, 5, 4],
    [2, 3, 0, 1, 5, 4],
    [2, 5, 0, 4, 3, 1],
    [2, 4, 0, 5, 1, 3],
    [1, 0, 4, 5, 2, 3],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 3, 2, 5, 4],
    [1, 0, 5, 4, 3, 2],
    [1, 0, 4, 5, 2, 3],
];

const AXIAL_TB: [[u8; 2]; 21] = [
    [5, 5],
    [0, 4],
    [0, 4],
    [0, 3],
    [0, 3],
    [0, 2],
    [0, 2],
    [0, 1],
    [0, 1],
    [1, 4],
    [1, 4],
    [1, 3],
    [1, 3],
    [1, 2],
    [1, 2],
    [2, 4],
    [2, 3],
    [3, 4],
    [3, 4],
    [2, 3],
    [2, 4],
];

const INVERT_TB: [u8; 21] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 11, 12, 9, 10, 14, 13, 20, 19, 18, 17, 16, 15,
];

const INVERT_OH: [u8; 31] = [
    0, 2, 1, 16, 14, 15, 18, 17, 10, 11, 8, 9, 13, 12, 4, 5, 3, 7, 6, 24, 23, 22, 21, 20, 19, 30,
    29, 28, 27, 26, 25,
];

/// `Atom::invertChirality` on a non-tetrahedral tag: the mirror-image
/// permutation (square-planar tags have none and stay as they are).
pub fn invert(class: NonTetrahedralClass, perm: u32) -> u32 {
    match class {
        NonTetrahedralClass::SquarePlanar => perm,
        NonTetrahedralClass::TrigonalBipyramidal => {
            INVERT_TB.get(perm as usize).map_or(0, |&p| u32::from(p))
        }
        NonTetrahedralClass::Octahedral => {
            INVERT_OH.get(perm as usize).map_or(0, |&p| u32::from(p))
        }
    }
}

/// The permutation after swapping ligand positions `x` and `y`
/// (`swap_squareplanar` / `swap_trigonalbipyramidal` / `swap_octahedral`);
/// 0 where RDKit's tables give none.
pub fn swap(class: NonTetrahedralClass, perm: u32, x: usize, y: usize) -> u32 {
    if x == y {
        return perm;
    }
    let (lo, hi) = if x < y { (x, y) } else { (y, x) };
    let n = class.max_nbors();
    if hi > n - 1 {
        return 0;
    }
    let offset = |i: usize| -> usize { (0..i).map(|k| n - 1 - k).sum() };
    let idx = offset(lo) + (hi - 1 - lo);
    let p = perm as usize;
    match class {
        NonTetrahedralClass::SquarePlanar if p < 4 => u32::from(SWAP_SP[p][idx]),
        NonTetrahedralClass::TrigonalBipyramidal if p < 21 => u32::from(SWAP_TB[p][idx]),
        NonTetrahedralClass::Octahedral if p < 31 => u32::from(SWAP_OH[p][idx]),
        _ => 0,
    }
}

/// The position across from ligand position `pos` (`*_across` tables),
/// `None` when there is none.
pub fn across(class: NonTetrahedralClass, perm: u32, pos: usize) -> Option<usize> {
    let p = perm as usize;
    let v = match class {
        NonTetrahedralClass::SquarePlanar if (1..=3).contains(&p) => {
            *ACROSS_SP[p].get(pos)? as usize
        }
        NonTetrahedralClass::TrigonalBipyramidal if (1..=20).contains(&p) => {
            *ACROSS_TB[p].get(pos)? as usize
        }
        NonTetrahedralClass::Octahedral if (1..=30).contains(&p) => {
            *ACROSS_OH[p].get(pos)? as usize
        }
        _ => return None,
    };
    (v < class.max_nbors()).then_some(v)
}

/// The two axial ligand positions of a trigonal bipyramid
/// (`trigonalbipyramidal_axial`), `None` outside `TB1`..`TB20`.
pub fn tb_axial(perm: u32) -> Option<[usize; 2]> {
    (1..=20)
        .contains(&perm)
        .then(|| AXIAL_TB[perm as usize].map(usize::from))
}

/// `Chirality::getChiralPermutation`'s permutation walk: `perm` refers to
/// the ligand order `reference` (ids, with `None` for implicit ligands); the
/// result refers to `probe`. Both hold the same ids. 0 when `perm` is 0 or
/// the walk hits a missing table entry.
pub fn permute(
    class: NonTetrahedralClass,
    mut perm: u32,
    reference: &[Option<usize>],
    probe: &[Option<usize>],
) -> Option<u32> {
    if reference.len() != probe.len() {
        return None;
    }
    let mut nbr = reference.to_vec();
    for i in 0..probe.len().saturating_sub(1) {
        let pval = probe[i];
        if nbr[i] == pval {
            continue;
        }
        let tgt = (i..nbr.len()).find(|&k| nbr[k] == pval)?;
        perm = swap(class, perm, i, tgt);
        nbr.swap(tgt, i);
    }
    Some(perm)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn swap_offsets_match_rdkit() {
        // RDKit's offsets: SP {0, 2, 3}, TB {0, 3, 5, 6}, OH {0, 4, 7, 9, 10}.
        assert_eq!(swap(NonTetrahedralClass::SquarePlanar, 1, 2, 3), 3);
        assert_eq!(swap(NonTetrahedralClass::SquarePlanar, 1, 1, 2), 2);
        assert_eq!(swap(NonTetrahedralClass::TrigonalBipyramidal, 1, 3, 4), 3);
        assert_eq!(swap(NonTetrahedralClass::TrigonalBipyramidal, 1, 2, 4), 6);
        assert_eq!(swap(NonTetrahedralClass::Octahedral, 1, 4, 5), 3);
        assert_eq!(swap(NonTetrahedralClass::Octahedral, 1, 3, 5), 7);
        assert_eq!(swap(NonTetrahedralClass::Octahedral, 1, 0, 5), 2);
    }

    #[test]
    fn swaps_are_involutions() {
        for class in [
            NonTetrahedralClass::SquarePlanar,
            NonTetrahedralClass::TrigonalBipyramidal,
            NonTetrahedralClass::Octahedral,
        ] {
            let n = class.max_nbors();
            for p in 1..=class.max_permutation() {
                for x in 0..n {
                    for y in 0..n {
                        assert_eq!(swap(class, swap(class, p, x, y), x, y), p);
                    }
                }
            }
        }
    }
}
