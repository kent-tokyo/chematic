//! RDKit-compatible "Layered fingerprint" port (`Chem.LayeredFingerprint`).
//!
//! Ported from RDKit's `Code/GraphMol/Fingerprints/Fingerprints.cpp`'s
//! `LayeredFingerprintMol`, default parameters (`layerFlags=0xFFFFFFFF,
//! minPath=1, maxPath=7, fpSize=2048`). Upstream itself documents this
//! fingerprint as "experimental — the API or results may change from release
//! to release"; this port targets its behavior as of the source read for this
//! implementation.
//!
//! # Algorithm
//!
//! Reuses the same branched-subgraph enumeration as this crate's
//! `rdkit_rdk_fp` port of `Chem.RDKFingerprint` (RDKit's own
//! `findAllSubgraphsOfLengthsMtoN`/`recurseWalkRange` — a "root at minimum
//! bond index" backtracking scheme) — duplicated here rather than shared
//! across crates/PRs, so this fingerprint stays independently mergeable. For
//! each subgraph (1..=7 bonds), up to 6 independent "layers" are computed,
//! each packing a different feature into a small bitfield per bond, then
//! folded separately:
//!
//! 1. **Topology**: `bondNbrs % 8` (bits 0-2) plus the two endpoint atoms'
//!    *path-local* degrees (bits 3-5, 6-8; larger degree first).
//! 2. **Bond order**: same as layer 1, but a bond hash (bits 0-2) is
//!    prepended — aromatic bonds and `SINGLE` bonds share the *same* code
//!    (`SINGLE`'s raw value), while `DOUBLE`/`TRIPLE`/`QUADRUPLE` etc. get
//!    their own raw `Bond::BondType` values. This intentionally differs from
//!    `rdkit_rdk_fp`'s own bond-type mapping (where aromatic gets a distinct
//!    code, `12`, from single) — confirmed against RDKit's own source, not
//!    assumed from the other fingerprint's convention.
//! 3. **Atom types**: `(atomic_number % 128)` for each endpoint (larger first,
//!    ties broken by larger path-local degree first) plus both degrees and
//!    the bond-neighbor count.
//! 4. **Ring membership**: emits a single `1` per ring bond (skipped
//!    entirely, not even zero, for non-ring bonds).
//! 5. **Ring size**: the *smallest* SSSR ring size containing the bond,
//!    `% 8` (`0` if the bond is in no ring).
//! 6. **Aromaticity**: each endpoint's atom-aromaticity flag (smaller-first
//!    canonical order) plus the bond-neighbor count.
//!
//! Each layer's collected per-bond values are sorted, then a trailing
//! distinct-atom count and the layer's own 1-indexed number are appended
//! (matching `rdkit_rdk_fp`'s own `hash_range` fold, reused verbatim here),
//! and the fold's `% fpSize` result sets one bit. Unlike `rdkit_rdk_fp`, there
//! is no `numBitsPerFeature`/PRNG step here — `LayeredFingerprintMol` sets
//! exactly one bit per (subgraph, non-empty layer) pair, and (unlike
//! `getEnvironments`) there is no single-bond shortcut: even a 1-bond
//! subgraph's single layer value goes through the full sort-then-fold.
//!
//! Ring membership/size use [`chematic_perception::find_sssr`] rather than
//! RDKit's own ring perception; this is the same pre-existing
//! `chematic-perception` dependency gap already documented for the other
//! fingerprints in this series and can disagree with RDKit's choice for some
//! fused/bridged ring systems.
//!
//! # Verified
//!
//! See this module's own test suite and the crate's fingerprint-parity
//! corpora for current bit-exactness figures.

use crate::bitvec::BitVec2048;
use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::aromaticity::{
    AromaticityAlgorithm, AromaticityModel, assign_aromaticity_ex,
};
use chematic_perception::find_sssr;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;

const FP_SIZE: u32 = 2048;
const MIN_PATH: usize = 1;
const MAX_PATH: usize = 7;
const NUM_LAYERS: usize = 6;

fn hash_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e3779b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// Layer 2's own bond-type code: aromatic and `SINGLE` bonds share one code;
/// other bond orders keep their own raw `Bond::BondType` value. Deliberately
/// different from `rdkit_rdk_fp`'s own bond-type mapping (which gives
/// aromatic its own distinct code) — confirmed from RDKit's own source for
/// this specific fingerprint, not copied from the other one.
fn layer2_bond_type_code(order: BondOrder) -> u32 {
    match order {
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        _ => 1,
    }
}

fn is_atom_aromatic(mol: &Molecule, aromaticity: Option<&AromaticityModel>, idx: AtomIdx) -> bool {
    if mol.atom(idx).aromatic {
        return true;
    }
    match aromaticity {
        Some(model) => model.is_atom_aromatic(idx),
        None => false,
    }
}

fn is_bond_aromatic(mol: &Molecule, aromaticity: Option<&AromaticityModel>, bidx: BondIdx) -> bool {
    let bond = mol.bond(bidx);
    bond.order == BondOrder::Aromatic
        || matches!(aromaticity, Some(model) if model.is_bond_aromatic(bidx))
}

/// Unlike `rdkit_rdk_fp`'s own port of `Chem.RDKFingerprint` (`useHs=true`),
/// RDKit's `LayeredFingerprintMol` enumerates with `useHs=false`: any bond
/// touching an explicit hydrogen atom (e.g. an isotope-labeled `[2H]`/`[3H]`,
/// which chematic — like RDKit — represents as a real graph atom, not folded
/// into an implicit H count) never gets an adjacency entry at all, so it can
/// never be a root or a candidate. Confirmed against a live oracle: a naive
/// `useHs=true` port silently included such bonds and produced extra bits.
fn build_bond_adjacency(mol: &Molecule) -> FxHashMap<u32, Vec<u32>> {
    let is_h = |idx: AtomIdx| mol.atom(idx).element.atomic_number() == 1;
    let keep = |bidx: BondIdx| {
        let bond = mol.bond(bidx);
        !is_h(bond.atom1) && !is_h(bond.atom2)
    };
    let mut nbrs: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for (bidx, _) in mol.bonds() {
        if keep(bidx) {
            nbrs.entry(bidx.0).or_default();
        }
    }
    for (aidx, _) in mol.atoms() {
        if is_h(aidx) {
            continue;
        }
        let incident: Vec<u32> = mol
            .neighbors(aidx)
            .map(|(_, b)| b)
            .filter(|&b| keep(b))
            .map(|b| b.0)
            .collect();
        for &b1 in &incident {
            for &b2 in &incident {
                if b1 != b2 {
                    nbrs.get_mut(&b1).unwrap().push(b2);
                }
            }
        }
    }
    nbrs
}

fn find_all_subgraphs(
    mol: &Molecule,
    min_len: usize,
    max_len: usize,
) -> BTreeMap<usize, Vec<Vec<u32>>> {
    let nbrs = build_bond_adjacency(mol);
    let n_bonds = mol.bond_count();
    let mut forbidden = vec![false; n_bonds];
    let mut res: BTreeMap<usize, Vec<Vec<u32>>> =
        (min_len..=max_len).map(|l| (l, Vec::new())).collect();

    let mut roots: Vec<u32> = nbrs.keys().copied().collect();
    roots.sort_unstable();
    for root in roots {
        if forbidden[root as usize] {
            continue;
        }
        forbidden[root as usize] = true;
        let cands = nbrs.get(&root).cloned().unwrap_or_default();
        let spath = vec![root];
        recurse_walk_range(
            &nbrs,
            spath,
            cands,
            min_len,
            max_len,
            forbidden.clone(),
            &mut res,
        );
    }
    res
}

fn recurse_walk_range(
    nbrs: &FxHashMap<u32, Vec<u32>>,
    spath: Vec<u32>,
    mut cands: Vec<u32>,
    min_len: usize,
    max_len: usize,
    mut forbidden: Vec<bool>,
    res: &mut BTreeMap<usize, Vec<Vec<u32>>>,
) {
    let nsize = spath.len();
    if nsize >= min_len && nsize <= max_len {
        res.get_mut(&nsize).unwrap().push(spath.clone());
    }
    if nsize >= max_len {
        return;
    }
    while let Some(next) = cands.pop() {
        if !forbidden[next as usize] {
            forbidden[next as usize] = true;
            let mut tstack = cands.clone();
            if let Some(next_nbrs) = nbrs.get(&next) {
                for &b in next_nbrs {
                    if !forbidden[b as usize] {
                        tstack.push(b);
                    }
                }
            }
            let mut tpath = spath.clone();
            tpath.push(next);
            recurse_walk_range(
                nbrs,
                tpath,
                tstack,
                min_len,
                max_len,
                forbidden.clone(),
                res,
            );
        }
    }
}

struct RingInfo {
    bond_in_ring: Vec<bool>,
    bond_min_ring_size: Vec<u32>,
}

fn compute_ring_info(mol: &Molecule) -> RingInfo {
    let n_bonds = mol.bond_count();
    let mut bond_in_ring = vec![false; n_bonds];
    let mut bond_min_ring_size = vec![u32::MAX; n_bonds];

    let ring_set = find_sssr(mol);
    for ring in ring_set.rings() {
        let size = ring.len() as u32;
        for i in 0..ring.len() {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % ring.len()];
            if let Some((_, bidx)) = mol.neighbors(a1).find(|&(nbr, _)| nbr == a2) {
                bond_in_ring[bidx.0 as usize] = true;
                let slot = &mut bond_min_ring_size[bidx.0 as usize];
                if size < *slot {
                    *slot = size;
                }
            }
        }
    }

    for v in &mut bond_min_ring_size {
        if *v == u32::MAX {
            *v = 0;
        }
    }
    RingInfo {
        bond_in_ring,
        bond_min_ring_size,
    }
}

/// Compute the RDKit-compatible "Layered fingerprint" (`Chem.LayeredFingerprint`).
///
/// See the module doc comment for the full algorithm and its verification
/// status.
pub fn rdkit_layered_fp(mol: &Molecule) -> BitVec2048 {
    let mut fp = BitVec2048::new();
    if mol.atom_count() == 0 {
        return fp;
    }

    let has_literal_aromatic_bond = mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic);
    let aromaticity = if has_literal_aromatic_bond {
        None
    } else {
        Some(assign_aromaticity_ex(mol, AromaticityAlgorithm::RdkitLike))
    };

    let anums: Vec<u32> = (0..mol.atom_count())
        .map(|i| mol.atom(AtomIdx(i as u32)).element.atomic_number() as u32)
        .collect();
    let atom_aromatic: Vec<bool> = (0..mol.atom_count())
        .map(|i| is_atom_aromatic(mol, aromaticity.as_ref(), AtomIdx(i as u32)))
        .collect();
    let ring_info = compute_ring_info(mol);

    let paths = find_all_subgraphs(mol, MIN_PATH, MAX_PATH);
    for plist in paths.values() {
        for path in plist {
            let mut atom_degrees: FxHashMap<u32, u32> = FxHashMap::default();
            let mut atoms_in_path: std::collections::BTreeSet<u32> =
                std::collections::BTreeSet::new();
            for &bidx in path {
                let bond = mol.bond(BondIdx(bidx));
                *atom_degrees.entry(bond.atom1.0).or_insert(0) += 1;
                *atom_degrees.entry(bond.atom2.0).or_insert(0) += 1;
                atoms_in_path.insert(bond.atom1.0);
                atoms_in_path.insert(bond.atom2.0);
            }

            let mut bond_nbrs = vec![0u32; path.len()];
            for i in 0..path.len() {
                let bi = mol.bond(BondIdx(path[i]));
                for j in (i + 1)..path.len() {
                    let bj = mol.bond(BondIdx(path[j]));
                    if bi.atom1 == bj.atom1
                        || bi.atom1 == bj.atom2
                        || bi.atom2 == bj.atom1
                        || bi.atom2 == bj.atom2
                    {
                        bond_nbrs[i] += 1;
                        bond_nbrs[j] += 1;
                    }
                }
            }

            let mut hash_layers: [Vec<u32>; NUM_LAYERS] = Default::default();

            for (i, &bidx) in path.iter().enumerate() {
                let bond = mol.bond(BondIdx(bidx));
                let a1 = bond.atom1.0;
                let a2 = bond.atom2.0;
                let deg1 = atom_degrees[&a1];
                let deg2 = atom_degrees[&a2];

                // Layer 1: straight topology.
                {
                    let (mut d1, mut d2) = (deg1, deg2);
                    if d1 < d2 {
                        std::mem::swap(&mut d1, &mut d2);
                    }
                    let mut h = bond_nbrs[i] % 8;
                    h |= (d1 % 8) << 3;
                    h |= (d2 % 8) << 6;
                    hash_layers[0].push(h);
                }

                // Layer 2: bond order.
                {
                    let is_aromatic = is_bond_aromatic(mol, aromaticity.as_ref(), BondIdx(bidx));
                    let bond_hash = if is_aromatic {
                        1
                    } else {
                        layer2_bond_type_code(bond.order)
                    };
                    let (mut d1, mut d2) = (deg1, deg2);
                    if d1 < d2 {
                        std::mem::swap(&mut d1, &mut d2);
                    }
                    let mut h = bond_hash % 8;
                    h |= (bond_nbrs[i] % 8) << 3;
                    h |= (d1 % 8) << 6;
                    h |= (d2 % 8) << 9;
                    hash_layers[1].push(h);
                }

                // Layer 3: atom types.
                {
                    let mut a1h = anums[a1 as usize] % 128;
                    let mut a2h = anums[a2 as usize] % 128;
                    let (mut d1, mut d2) = (deg1, deg2);
                    if a1h < a2h {
                        std::mem::swap(&mut a1h, &mut a2h);
                        std::mem::swap(&mut d1, &mut d2);
                    } else if a1h == a2h && d1 < d2 {
                        std::mem::swap(&mut d1, &mut d2);
                    }
                    let mut h = a1h;
                    h |= a2h << 7;
                    h |= (d1 % 8) << 14;
                    h |= (d2 % 8) << 17;
                    h |= (bond_nbrs[i] % 8) << 20;
                    hash_layers[2].push(h);
                }

                // Layer 4: ring membership.
                if ring_info.bond_in_ring[bidx as usize] {
                    hash_layers[3].push(1);
                }

                // Layer 5: ring size.
                {
                    let h = ring_info.bond_min_ring_size[bidx as usize] % 8;
                    hash_layers[4].push(h);
                }

                // Layer 6: aromaticity.
                {
                    let mut a1h = atom_aromatic[a1 as usize];
                    let mut a2h = atom_aromatic[a2 as usize];
                    if !a1h && a2h {
                        std::mem::swap(&mut a1h, &mut a2h);
                    }
                    let mut h = a1h as u32;
                    h |= (a2h as u32) << 1;
                    h |= (bond_nbrs[i] % 8) << 5;
                    hash_layers[5].push(h);
                }
            }

            for (l, layer) in hash_layers.iter_mut().enumerate() {
                if layer.is_empty() {
                    continue;
                }
                layer.sort_unstable();
                layer.push(atoms_in_path.len() as u32);
                layer.push((l + 1) as u32);
                let mut seed = 0u32;
                for &v in layer.iter() {
                    seed = hash_combine(seed, v);
                }
                fp.set((seed % FP_SIZE) as usize);
            }
        }
    }

    fp
}

/// Tanimoto similarity between two RDKit-compatible Layered fingerprints.
pub fn tanimoto_rdkit_layered(a: &BitVec2048, b: &BitVec2048) -> f64 {
    a.tanimoto(b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    #[test]
    fn deterministic() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(rdkit_layered_fp(&m), rdkit_layered_fp(&m));
    }

    #[test]
    fn different_molecules_differ() {
        let a = rdkit_layered_fp(&mol("CCO"));
        let b = rdkit_layered_fp(&mol("c1ccccc1"));
        assert_ne!(a, b);
    }

    #[test]
    fn single_atom_has_no_bits() {
        let fp = rdkit_layered_fp(&mol("C"));
        assert_eq!(fp.popcount(), 0);
    }

    #[test]
    fn matches_rdkit_oracle_on_hand_picked_molecules() {
        // (smiles, expected on-bit count from a live `Chem.LayeredFingerprint`
        // oracle run at minPath=1, maxPath=7, fpSize=2048).
        let cases: &[(&str, u32)] = &[
            ("CC", 5),
            ("C=C", 5),
            ("C#C", 5),
            ("CCC", 10),
            ("CC(C)C", 15),
            ("c1ccccc1", 36),
            ("C1CC1", 18),
            ("CCO", 11),
            ("CC(=O)O", 19),
            ("c1ccc2ccccc2c1", 104),
            ("CC(=O)Oc1ccccc1C(=O)O", 335),
        ];
        for &(smiles, expected) in cases {
            let fp = rdkit_layered_fp(&mol(smiles));
            assert_eq!(fp.popcount(), expected, "popcount mismatch for {smiles:?}");
        }
    }

    #[test]
    fn explicit_isotope_labeled_hydrogen_is_excluded_from_enumeration() {
        // Regression pin for a real bug this port had: LayeredFingerprintMol
        // enumerates with useHs=false, unlike rdkit_rdk_fp's useHs=true, so a
        // bond touching an explicit H atom (isotope-labeled deuterium/tritium
        // here, since chematic -- like RDKit -- represents those as real
        // graph atoms) must never appear in the enumeration. A naive
        // useHs=true port produced 56 bits here instead of the correct 19.
        let fp = rdkit_layered_fp(&mol("[2H]C([2H])([2H])NC=O"));
        assert_eq!(fp.popcount(), 19);
    }
}
