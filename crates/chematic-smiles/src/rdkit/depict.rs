//! RDKit's native 2D depiction (`RDDepict::compute2DCoords`, RDKit
//! 2026.03.1 `Code/GraphMol/Depictor`: `RDDepictor.cpp`, `EmbeddedFrag.cpp`,
//! `DepictUtils.cpp`, and `Geometry/Transform2D.cpp`), as
//! `rdDepictor.Compute2DCoords(mol)` runs it from Python with its defaults:
//! `canonOrient=True`, an empty coordinate map, no random sampling (so
//! collisions are removed by bond flips), no ring templates, no CoordGen.
//!
//! The port follows the C++ operation for operation, including the order
//! of floating-point operations, the `std::map` iteration order of the
//! embedded atoms, the way `operator[]` inserts default atoms, and the
//! aliasing of references the C++ keeps into atoms it moves. `sin`, `cos`
//! and `acos` are the platform's libm functions, as RDKit calls them. Recorded
//! Linux comparisons are bit-identical; other platforms may differ in the
//! last floating-point bits while remaining within the documented tolerance.

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::RdkitSmilesError;
use super::mol::{BondStereo, BondType, ChiralTag, Hybridization, Mol};

mod merge_helpers;

const BOND_LEN: f64 = 1.5;
const COLLISION_THRES: f64 = 0.70;
const BOND_THRES: f64 = 0.50;
const ANGLE_OPEN: f64 = 0.1222;
const MAX_COLL_ITERS: u32 = 15;
const HETEROATOM_COLL_SCALE: f64 = 1.3;
const NUM_BONDS_FLIPS: u32 = 3;
const NEIGH_RADIUS: f64 = 2.5;
const ZERO_TOLERANCE: f64 = 1.0e-16;
const LOCAL_INF: f64 = 1.0e8;

/// An invariant the C++ checks (and throws on) failed.
#[derive(Debug)]
struct DepictError(String);

type DResult<T> = Result<T, DepictError>;

fn derr<T>(what: &str) -> DResult<T> {
    Err(DepictError(what.to_string()))
}

include!("depict/geometry.rs");

// ---------------------------------------------------------------------------
// The molecule as the depictor sees it

struct Ctx<'a> {
    mol: &'a Mol,
    /// `_CIPRank` per atom, if set.
    cip: Option<&'a [u32]>,
    n: usize,
    /// `MolOps::getDistanceMat` (topological), filled on first use.
    dmat: Vec<f64>,
}

impl Ctx<'_> {
    fn degree(&self, a: usize) -> usize {
        self.mol.degree(a)
    }

    fn anum(&self, a: usize) -> u32 {
        self.mol.atoms[a].anum
    }

    fn num_bond_rings(&self, b: usize) -> usize {
        self.mol.num_bond_rings(b)
    }

    /// `getAtomDepictRank`.
    fn depict_rank(&self, a: usize) -> i64 {
        let anum = self.anum(a) as i64;
        let anum = if anum == 1 { 1000 } else { anum };
        100 * anum + self.degree(a) as i64
    }

    /// `rankAtomsByRank(mol, atoms, ascending=true)`.
    fn rank_atoms(&self, atoms: &[usize]) -> Vec<usize> {
        let mut rank_aid: Vec<(i32, i32)> = atoms
            .iter()
            .map(|&aid| {
                let rank: u32 = match self.cip {
                    Some(c) => c[aid],
                    None => (aid as u32)
                        .wrapping_add((self.n as u32).wrapping_mul(self.depict_rank(aid) as u32)),
                };
                (rank as i32, aid as i32)
            })
            .collect();
        rank_aid.sort(); // stable; (rank, aid) pairs
        rank_aid.into_iter().map(|(_, a)| a as usize).collect()
    }

    /// `setNbrOrder`.
    fn set_nbr_order(&self, aid: usize, nbrs: &[usize]) -> DResult<Vec<usize>> {
        let mut r: i64 = -1;
        for nb in self.mol.nbrs(aid) {
            if !nbrs.contains(&nb) {
                r = nb as i64;
            }
        }
        let mut thold: Vec<usize> = nbrs.to_vec();
        if r >= 0 {
            thold.push(r as usize);
        }
        if thold.len() <= 3 {
            return derr("setNbrOrder: too few neighbors");
        }
        let mut thold = self.rank_atoms(&thold);
        let ln = thold.len();
        thold.swap(ln - 3, ln - 2);
        let mut res = Vec::with_capacity(thold.len());
        let pos = if r >= 0 {
            thold.iter().position(|&x| x == r as usize)
        } else {
            None
        };
        match pos {
            Some(p) => {
                res.extend_from_slice(&thold[p + 1..]);
                res.extend_from_slice(&thold[..p]);
            }
            None => res.extend_from_slice(&thold),
        }
        Ok(res)
    }

    fn ensure_dmat(&mut self) {
        if !self.dmat.is_empty() || self.n == 0 {
            return;
        }
        let n = self.n;
        let mut d = vec![LOCAL_INF; n * n];
        let mut queue = std::collections::VecDeque::new();
        for s in 0..n {
            d[s * n + s] = 0.0;
            queue.clear();
            queue.push_back(s);
            while let Some(c) = queue.pop_front() {
                let dc = d[s * n + c];
                for nb in self.mol.nbrs(c) {
                    if d[s * n + nb] == LOCAL_INF {
                        d[s * n + nb] = dc + 1.0;
                        queue.push_back(nb);
                    }
                }
            }
        }
        self.dmat = d;
    }

    /// `MolOps::getShortestPath`.
    fn shortest_path(&self, aid1: usize, aid2: usize) -> Vec<usize> {
        let n = self.n;
        let mut pred = vec![-1i64; n];
        pred[aid1] = -2;
        pred[aid2] = -3;
        let mut q = std::collections::VecDeque::new();
        q.push_back(aid1);
        let mut done = false;
        while !done && !q.is_empty() {
            let cur = *q.front().unwrap();
            for nb in self.mol.nbrs(cur) {
                if done {
                    break;
                }
                match pred[nb] {
                    -1 => {
                        pred[nb] = cur as i64;
                        q.push_back(nb);
                    }
                    -3 => {
                        pred[nb] = cur as i64;
                        done = true;
                    }
                    _ => {}
                }
            }
            q.pop_front();
        }
        let mut res = std::collections::VecDeque::new();
        if done {
            let mut prev = aid2 as i64;
            res.push_back(aid2);
            loop {
                prev = pred[prev as usize];
                if prev != aid1 as i64 {
                    res.push_front(prev as usize);
                } else {
                    break;
                }
            }
            res.push_front(aid1);
        }
        res.into_iter().collect()
    }

    /// `getRotatableBonds(mol, aid1, aid2)`.
    fn rotatable_bonds(&self, aid1: usize, aid2: usize) -> Vec<usize> {
        let mut path = self.shortest_path(aid1, aid2);
        let mut res = Vec::new();
        if path.len() >= 4 {
            path.remove(0);
            path.pop();
            let mut pid = path[0];
            for &aid in &path {
                if aid == pid {
                    continue;
                }
                let bid = self.mol.bond_between(pid, aid).expect("path bond");
                if self.mol.bonds[bid].stereo <= BondStereo::Any && self.num_bond_rings(bid) == 0 {
                    res.push(bid);
                }
                pid = aid;
            }
        }
        res
    }

    /// `_recurseAtomOneSide`.
    fn recurse_atom_one_side(&self, end_aid: usize, beg_aid: usize, flip: &mut Vec<usize>) {
        flip.push(end_aid);
        let nbrs: Vec<usize> = self.mol.nbrs(end_aid).collect();
        for nb in nbrs {
            if nb != beg_aid && !flip.contains(&nb) {
                self.recurse_atom_one_side(nb, beg_aid, flip);
            }
        }
    }

    /// `_findDeg1Neighbor`.
    fn deg1_neighbor(&self, aid: usize) -> DResult<usize> {
        if self.degree(aid) != 1 {
            return derr("_findDeg1Neighbor: degree is not 1");
        }
        Ok(self.mol.nbrs(aid).next().unwrap())
    }

    /// `_findClosestNeighbor`.
    fn closest_neighbor(&self, aid1: usize, aid2: usize) -> usize {
        let mut res = 0;
        let mut mdist = 1.0e8;
        let naid = aid1 * self.n;
        for nb in self.mol.nbrs(aid2) {
            let d = self.dmat[naid + nb];
            if d < mdist {
                mdist = d;
                res = nb;
            }
        }
        res
    }

    /// `_findClosestPair`.
    fn closest_pair(&self, beg1: usize, end1: usize, beg2: usize, end2: usize) -> (usize, usize) {
        let na = self.n;
        let d1 = self.dmat[beg1 * na + beg2];
        let d2 = self.dmat[beg1 * na + end2];
        let d3 = self.dmat[end1 * na + beg2];
        let d4 = self.dmat[end1 * na + end2];
        // std::min(a, b, comp) keeps `a` unless comp(b, a).
        let mut min = (d1, (beg1, beg2));
        if d2 < min.0 {
            min = (d2, (beg1, end2));
        }
        if d3 < min.0 {
            min = (d3, (end1, beg2));
        }
        if d4 < min.0 {
            min = (d4, (end1, end2));
        }
        min.1
    }

    /// `_recurseDegTwoRingAtoms`.
    fn recurse_deg_two_ring_atoms(
        &self,
        aid: usize,
        rpath: &mut Vec<usize>,
        nbr_map: &mut BTreeMap<usize, Vec<usize>>,
    ) {
        let mut nbrs = Vec::new();
        for &b in &self.mol.atom_bonds[aid] {
            if self.num_bond_rings(b) > 0 {
                nbrs.push(self.mol.bonds[b].other(aid));
            }
        }
        if nbrs.len() != 2 {
            return;
        }
        rpath.push(aid);
        nbr_map.insert(aid, nbrs.clone());
        for nb in nbrs {
            if !rpath.contains(&nb) {
                self.recurse_deg_two_ring_atoms(nb, rpath, nbr_map);
            }
        }
    }

    /// `_anyNonRingBonds`.
    fn any_non_ring_bonds(&self, aid: usize, path: &[usize]) -> DResult<u32> {
        let mut prev = aid;
        let mut n_open = 0;
        for &pi in path {
            let Some(b) = self.mol.bond_between(prev, pi) else {
                return derr("no bond found");
            };
            if self.num_bond_rings(b) == 0 {
                n_open += 1;
            }
            prev = pi;
        }
        Ok(n_open)
    }
}

include!("depict/embedded_frag.rs");

/// `EmbeddedAtom::operator=` (copies everything but `aid`, which the port
/// does not keep).
fn assign_eatom(dst: &mut EAtom, src: &EAtom) {
    *dst = src.clone();
}

fn cross_val(v1: P, v2: P) -> f64 {
    v1.x * v2.y - v2.x * v1.y
}

/// `RDKit::Union(rings, res, exclude)`.
fn union_rings(rings: &[Vec<usize>], exclude: Option<&[usize]>) -> Vec<usize> {
    let mut res: Vec<usize> = Vec::new();
    for (id, ring) in rings.iter().enumerate() {
        if let Some(ex) = exclude
            && ex.contains(&id)
        {
            continue;
        }
        for &a in ring {
            if !res.contains(&a) {
                res.push(a);
            }
        }
    }
    res
}

/// `pickFirstRingToEmbed`.
fn pick_first_ring_to_embed(ctx: &Ctx, fused: &[Vec<usize>]) -> usize {
    let mut res = 0usize;
    let mut max_size = 0usize;
    let mut minsubs = 100_000_000i64;
    for (cnt, ring) in fused.iter().enumerate() {
        let subs = ring.iter().filter(|&&a| ctx.degree(a) > 2).count() as i64;
        if subs < minsubs {
            res = cnt;
            minsubs = subs;
            max_size = ring.len();
        } else if subs == minsubs && ring.len() > max_size {
            res = cnt;
            max_size = ring.len();
        }
    }
    res
}

/// `findNextRingToEmbed`.
fn find_next_ring_to_embed(
    done_rings: &[usize],
    fused: &[Vec<usize>],
) -> DResult<(Vec<usize>, usize)> {
    let not_done: Vec<usize> = (0..fused.len())
        .filter(|i| !done_rings.contains(i))
        .collect();
    let done_atoms = union_rings(fused, Some(&not_done));
    let mut max_common = 0usize;
    let mut next_id: Option<usize> = None;
    let mut res: Vec<usize> = Vec::new();
    for (curr, ring) in fused.iter().enumerate() {
        if done_rings.contains(&curr) {
            continue;
        }
        let common: Vec<usize> = ring
            .iter()
            .filter(|a| done_atoms.contains(a))
            .copied()
            .collect();
        if common.len() == 2 {
            return Ok((common, curr));
        }
        if common.len() > max_common {
            max_common = common.len();
            next_id = Some(curr);
            res = common;
        }
    }
    let Some(next_id) = next_id else {
        return derr("findNextRingToEmbed: no ring found");
    };
    let mut cmn_lst = 0usize;
    let n_cmn = res.len();
    for i in 0..n_cmn {
        if res[i] == fused[next_id][i] {
            cmn_lst += 1;
        } else {
            break;
        }
    }
    if cmn_lst > 0 && cmn_lst < res.len() {
        let temp = res.clone();
        let n_mov = n_cmn - cmn_lst;
        res[..n_mov].copy_from_slice(&temp[cmn_lst..n_cmn]);
        res[n_mov..n_cmn].copy_from_slice(&temp[..cmn_lst]);
    }
    if res.is_empty() {
        return derr("findNextRingToEmbed: empty result");
    }
    Ok((res, next_id))
}

/// `mirrorTransRingAtoms`.
fn mirror_trans_ring_atoms(
    ctx: &Ctx,
    ring: &[usize],
    coords: &mut BTreeMap<usize, P>,
) -> DResult<()> {
    let n = ring.len();
    for i in 0..n {
        let atom1 = ring[i];
        let atom2 = ring[(i + 1) % n];
        let Some(b) = ctx.mol.bond_between(atom1, atom2) else {
            return derr("ring atoms not bonded");
        };
        let bond = &ctx.mol.bonds[b];
        if bond.bt != BondType::Double {
            continue;
        }
        let stype = bond.stereo;
        if stype <= BondStereo::Any {
            continue;
        }
        if bond.stereo_atoms.len() != 2 {
            continue;
        }
        let left_in = ring.contains(&bond.stereo_atoms[0]);
        let right_in = ring.contains(&bond.stereo_atoms[1]);
        let is_trans = if stype == BondStereo::E {
            left_in == right_in
        } else {
            left_in != right_in
        };
        if !is_trans {
            continue;
        }
        let left = ring[(i + n - 1) % n];
        let right = atom2;
        let last = *coords.entry(left).or_default();
        let r = *coords.entry(right).or_default();
        let interest = *coords.entry(atom1).or_default();
        let d = last.sub(r);
        let a = (d.x * d.x - d.y * d.y) / d.dot(d);
        let bb = 2.0 * d.x * d.y / d.dot(d);
        let x = a * (interest.x - r.x) + bb * (interest.y - r.y) + r.x;
        let y = bb * (interest.x - r.x) - a * (interest.y - r.y) + r.y;
        coords.insert(atom1, P::new(x, y));
    }
    Ok(())
}

/// `makeRingNeighborMap` + `pickFusedRings` + `embedFusedSystems`.
fn embed_fused_systems(ctx: &Ctx, arings: &[Vec<usize>], efrags: &mut Vec<EFrag>) -> DResult<()> {
    let nr = arings.len();
    let mut neigh: Vec<Vec<usize>> = vec![Vec::new(); nr];
    for i in 0..nr {
        for j in i + 1..nr {
            if arings[i].iter().any(|a| arings[j].contains(a)) {
                neigh[i].push(j);
                neigh[j].push(i);
            }
        }
    }
    let mut fus_done = vec![false; nr];
    let mut curr = 0usize;
    while curr < nr {
        let mut fused = Vec::new();
        pick_fused_rings(curr, &neigh, &mut fused, &mut fus_done);
        let frings: Vec<Vec<usize>> = fused.iter().map(|&r| arings[r].clone()).collect();
        let mut efrag = EFrag::from_fused_rings(ctx, &frings)?;
        efrag.setup_new_neighs(ctx)?;
        efrags.push(efrag);
        match (0..nr).find(|&r| !fus_done[r]) {
            Some(r) => curr = r,
            None => break,
        }
    }
    Ok(())
}

fn pick_fused_rings(curr: usize, neigh: &[Vec<usize>], res: &mut Vec<usize>, done: &mut [bool]) {
    done[curr] = true;
    res.push(curr);
    for &n in &neigh[curr] {
        if !done[n] {
            pick_fused_rings(n, neigh, res, done);
        }
    }
}

/// `_findLargestFrag`.
fn find_largest_frag(efrags: &[EFrag]) -> Option<usize> {
    let mut msiz = 0;
    let mut res = None;
    for (i, f) in efrags.iter().enumerate() {
        if !f.dead && !f.done && f.size() > msiz {
            msiz = f.size();
            res = Some(i);
        }
    }
    res
}

/// `_shiftCoords`.
fn shift_coords(efrags: &mut [EFrag]) {
    let alive: Vec<usize> = (0..efrags.len()).filter(|&i| !efrags[i].dead).collect();
    if alive.is_empty() {
        return;
    }
    for &i in &alive {
        efrags[i].compute_box();
    }
    let first = &efrags[alive[0]];
    let mut xmax = first.px;
    let xmin = first.nx;
    let mut ymax = first.py;
    let ymin = first.ny;
    for &i in &alive[1..] {
        let mut xshift = true;
        if xmax + xmin > ymax + ymin {
            xshift = false;
        }
        let f = &mut efrags[i];
        let (xn, xp, yn, yp) = (f.nx, f.px, f.ny, f.py);
        let shift;
        if xshift {
            shift = P::new(xmax + xn + 1.0, 0.0);
            xmax += xp + xn + 1.0;
        } else {
            shift = P::new(0.0, ymax + yn + 1.0);
            ymax += yp + yn + 1.0;
        }
        for a in f.ea.values_mut() {
            a.loc = a.loc.add(shift);
        }
    }
}

/// RDKit's truncated `ISQRT2` (not `FRAC_1_SQRT_2`: the template
/// coordinates use these digits).
#[allow(clippy::approx_constant)]
const ISQRT2: f64 = 0.707107;
const SQRT3_2: f64 = 0.866025;

/// `embedNontetrahedralStereo`: a fixed template fragment per square-planar,
/// trigonal-bipyramidal or octahedral centre.
fn embed_nontetrahedral_stereo(
    ctx: &Ctx,
    atom_ranks: &[i64],
    efrags: &mut Vec<EFrag>,
) -> DResult<()> {
    let mol = ctx.mol;
    for a in 0..mol.atoms.len() {
        let tag = mol.atoms[a].chiral;
        if tag.nontet().is_none() {
            continue;
        }
        // getRankedAtomNeighbors (a stable sort for these short lists).
        let mut nbrs: Vec<usize> = mol.nbrs(a).collect();
        nbrs.sort_by_key(|&x| atom_ranks[x]);
        let mut map: BTreeMap<usize, P> = BTreeMap::new();
        map.insert(a, P::new(0.0, 0.0));
        match tag {
            ChiralTag::SquarePlanar => {
                let pts = [
                    P::new(ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
                    P::new(ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
                    P::new(-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
                    P::new(-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
                ];
                if nbrs.is_empty() {
                    return derr("square-planar centre without neighbours");
                }
                map.insert(nbrs[0], pts[0]);
                let mut q2_full = false;
                for &nbr in &nbrs[1..] {
                    let angle = mol.ideal_angle_between_ligands(a, nbrs[0], nbr);
                    if (angle - 180.0).abs() < 0.1 {
                        map.insert(nbr, pts[2]);
                    } else if !q2_full {
                        map.insert(nbr, pts[1]);
                        q2_full = true;
                    } else {
                        map.insert(nbr, pts[3]);
                    }
                }
            }
            ChiralTag::TrigonalBipyramidal => {
                let pts = [
                    P::new(0.0, BOND_LEN),
                    P::new(0.0, -BOND_LEN),
                    P::new(-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),
                    P::new(-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0),
                    P::new(BOND_LEN, 0.0),
                ];
                let axial1 = mol.tb_axial_atom(a, 1);
                let axial2 = mol.tb_axial_atom(a, -1);
                if let Some(x) = axial1 {
                    map.insert(x, pts[0]);
                }
                if let Some(x) = axial2 {
                    map.insert(x, pts[1]);
                }
                let mut which = 2;
                for &nbr in &nbrs {
                    if Some(nbr) != axial1 && Some(nbr) != axial2 {
                        // Past the template (no axial ligands known) RDKit
                        // reads beyond its static array, which holds zeros.
                        map.insert(nbr, pts.get(which).copied().unwrap_or_default());
                        which += 1;
                    }
                }
            }
            _ => {
                let pts = [
                    P::new(0.0, BOND_LEN),
                    P::new(0.0, -BOND_LEN),
                    P::new(SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),
                    P::new(SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0),
                    P::new(-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0),
                    P::new(-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),
                ];
                let mut axial1 = None;
                let mut axial2 = None;
                for i in 0..nbrs.len() {
                    let mut all90 = true;
                    for j in i + 1..nbrs.len() {
                        let ang = mol.ideal_angle_between_ligands(a, nbrs[i], nbrs[j]);
                        if (ang - 180.0).abs() < 0.1 {
                            axial1 = Some(nbrs[i]);
                            axial2 = Some(nbrs[j]);
                            all90 = false;
                            break;
                        } else if (ang - 90.0).abs() > 0.1 {
                            all90 = false;
                        }
                    }
                    if all90 {
                        axial1 = Some(nbrs[i]);
                    }
                    if axial1.is_some() {
                        break;
                    }
                }
                if let Some(x) = axial1 {
                    map.insert(x, pts[0]);
                }
                if let Some(x) = axial2 {
                    map.insert(x, pts[1]);
                }
                let mut ref1: Option<usize> = None;
                let mut ref2: Option<usize> = None;
                for &nbr in &nbrs {
                    if Some(nbr) == axial1 || Some(nbr) == axial2 {
                        continue;
                    }
                    if ref1.is_none() {
                        ref1 = Some(nbr);
                        map.insert(nbr, pts[2]);
                        ref2 = mol.chiral_across_atom(a, nbr);
                        if let Some(x) = ref2 {
                            map.insert(x, pts[4]);
                        }
                    } else {
                        if Some(nbr) == ref2 || Some(nbr) == ref1 {
                            continue;
                        }
                        map.insert(nbr, pts[3]);
                        if let Some(x) = mol.chiral_across_atom(a, nbr) {
                            map.insert(x, pts[5]);
                        }
                        break;
                    }
                }
            }
        }
        efrags.push(EFrag::from_coord_map(ctx, &map)?);
    }
    Ok(())
}

/// `compute2DCoords` with `rdDepictor.Compute2DCoords`' defaults on an
/// RDKit molecule as `MolFromSmiles` leaves it; `cip` holds the atoms'
/// `_CIPRank` values (empty when unset).
pub(crate) fn compute_2d_coords(mol: &Mol, cip: &[u32]) -> Result<Vec<[f64; 2]>, RdkitSmilesError> {
    compute(mol, cip, true)
        .map_err(|e| RdkitSmilesError::Unsupported(format!("depiction: {}", e.0)))
}

fn compute(mol_in: &Mol, cip: &[u32], canon_orient: bool) -> DResult<Vec<[f64; 2]>> {
    let n = mol_in.atoms.len();
    // symmetrizeSSSR(mol, arings, includeDativeBonds=true) on the copy.
    let mut mol = mol_in.clone();
    let bonds: Vec<(usize, usize, bool)> =
        mol.bonds.iter().map(|b| (b.begin, b.end, true)).collect();
    let Some(arings) = chematic_perception::rdkit_symmetrized_sssr(n, &bonds) else {
        return derr("ring perception falls back to RDKit's approximate ring finder");
    };
    mol.set_rings(arings.clone());
    let mut ctx = Ctx {
        mol: &mol,
        cip: if cip.len() == n && n > 0 {
            Some(cip)
        } else {
            None
        },
        n,
        dmat: Vec::new(),
    };
    let atom_ranks: Vec<i64> = (0..n).map(|a| ctx.depict_rank(a)).collect();

    let mut efrags: Vec<EFrag> = Vec::new();
    if !arings.is_empty() {
        embed_fused_systems(&ctx, &arings, &mut efrags)?;
    }
    embed_nontetrahedral_stereo(&ctx, &atom_ranks, &mut efrags)?;
    // embedCisTransSystems
    for (b, bond) in mol.bonds.iter().enumerate() {
        if bond.bt == BondType::Double
            && bond.stereo > BondStereo::Any
            && ctx.num_bond_rings(b) == 0
        {
            if bond.stereo_atoms.len() != 2 {
                continue;
            }
            let mut f = EFrag::from_double_bond(b, &ctx);
            f.setup_new_neighs(&ctx)?;
            efrags.push(f);
        }
    }
    // getNonEmbeddedAtoms
    let mut done = vec![false; n];
    for f in &efrags {
        for &k in f.ea.keys() {
            if k < n {
                done[k] = true;
            }
        }
    }
    let mut nratms: Vec<usize> = (0..n).filter(|&a| !done[a]).collect();
    let mut mri = find_largest_frag(&efrags);
    while mri.is_some() || !nratms.is_empty() {
        let i = match mri {
            Some(i) => i,
            None => {
                let mut mrank = i32::MAX as i64;
                let mut mnri = 0;
                for (p, &a) in nratms.iter().enumerate() {
                    let rank = atom_ranks[a] * n as i64 + a as i64;
                    if rank < mrank {
                        mrank = rank;
                        mnri = p;
                    }
                }
                let a = nratms.remove(mnri);
                efrags.push(EFrag::from_atom(&ctx, a)?);
                efrags.len() - 1
            }
        };
        efrags[i].done = true;
        let mut frag = std::mem::take(&mut efrags[i]);
        frag.done = true;
        let res = frag.expand(&ctx, &mut nratms, &mut efrags);
        efrags[i] = frag;
        res?;
        mri = find_largest_frag(&efrags);
    }

    ctx.ensure_dmat();
    for f in efrags.iter_mut().filter(|f| !f.dead) {
        f.remove_collisions_bond_flip(&ctx)?;
    }
    for f in efrags.iter_mut().filter(|f| !f.dead) {
        f.remove_collisions_open_angles(&ctx)?;
        f.remove_collisions_shorten_bonds(&ctx)?;
    }
    if canon_orient {
        for f in efrags.iter_mut().filter(|f| !f.dead) {
            f.canonicalize_orientation()?;
        }
    }
    shift_coords(&mut efrags);
    let mut out = vec![[0.0f64; 2]; n];
    for f in efrags.iter().filter(|f| !f.dead) {
        for (&k, a) in &f.ea {
            if k < n {
                out[k] = [a.loc.x, a.loc.y];
            } else {
                return derr("embedded atom index out of range");
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod fixed_fragment_boundary_tests {
    use super::*;

    #[test]
    fn geometry_helpers_cover_degenerate_and_reflected_inputs() {
        let mut zero = P::new(0.0, 0.0);
        assert!(zero.normalize().is_err());
        assert!(compute_normal(P::new(1.0, 1.0), P::new(1.0, 1.0)).is_err());
        assert!(angle_to_2d(P::new(0.0, 0.0), P::new(1.0, 0.0)).is_err());

        let identity = T2::align(
            P::new(0.0, 0.0),
            P::new(0.0, 0.0),
            P::new(1.0, 1.0),
            P::new(1.0, 1.0),
        );
        let mut point = P::new(2.0, 3.0);
        identity.apply(&mut point);
        assert_eq!(point, P::new(2.0, 3.0));

        let rotated = T2::rotation_about(P::new(1.0, 1.0), PI);
        rotated.apply(&mut point);
        assert!((point.x - 0.0).abs() < 1e-10);
        assert!((point.y + 1.0).abs() < 1e-10);

        let reflected = reflect_point(P::new(1.0, 2.0), P::new(0.0, 0.0), P::new(2.0, 0.0));
        assert!((reflected.x - 1.0).abs() < 1e-10);
        assert!((reflected.y + 2.0).abs() < 1e-10);
        assert_eq!(
            rotation_dir(P::default(), P::new(1.0, 0.0), P::new(0.0, 1.0), 0.0),
            -1
        );
        assert_eq!(
            rotation_dir(P::default(), P::new(1.0, 0.0), P::new(0.0, -1.0), 0.0),
            1
        );
        assert_eq!(
            compute_bisect_point(P::default(), PI + 0.1, P::new(1.0, 0.0), P::new(0.0, 1.0)),
            P::new(-0.5, -0.5)
        );
    }

    #[test]
    fn partial_fixed_fragments_keep_coordinates_and_place_remaining_ligands() {
        let input = crate::parse("C(F)(Cl)(Br)I").unwrap();
        let mut mol = super::super::parse::from_chematic(&input).unwrap();
        mol.set_rings(Vec::new());
        let ctx = Ctx {
            mol: &mol,
            cip: None,
            n: mol.atoms.len(),
            dmat: Vec::new(),
        };
        for placed in 0..=3 {
            for rotation in [0.0_f64, 0.37, -1.1] {
                let centre = P::new(3.0, -2.0);
                let mut map = BTreeMap::from([(0, centre)]);
                for i in 0..placed {
                    let angle = rotation + i as f64 * PI / 2.0;
                    map.insert(
                        i + 1,
                        centre.add(P::new(angle.cos(), angle.sin()).scale(BOND_LEN)),
                    );
                }
                let mut fragment = EFrag::from_coord_map(&ctx, &map).unwrap();
                if placed >= 2 {
                    // Two ligands store their occupied angle; three store
                    // the complementary sector of the widest pair.
                    let expected = if placed == 2 { PI / 2.0 } else { PI };
                    let actual = fragment.get(0).unwrap().angle;
                    // acos loses about sqrt(epsilon) angular precision at pi.
                    assert!(
                        (actual - expected).abs() < 1e-7,
                        "fixed={placed}, rotation={rotation}: {actual} != {expected}"
                    );
                }
                let mut remaining: Vec<_> = (0..ctx.n).filter(|a| !map.contains_key(a)).collect();
                fragment.expand(&ctx, &mut remaining, &mut []).unwrap();
                assert!(remaining.is_empty());
                assert!(fragment.attach.is_empty());
                assert_eq!(fragment.size(), ctx.n);
                for (&atom, &position) in &map {
                    assert_eq!(fragment.loc(atom).unwrap(), position);
                    assert!(fragment.get(atom).unwrap().fixed);
                }
                for atom in 1..ctx.n {
                    let position = fragment.loc(atom).unwrap();
                    assert!(position.x.is_finite() && position.y.is_finite());
                    assert!((position.sub(centre).length() - BOND_LEN).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn fixed_fragment_rejects_coincident_vectors_with_context() {
        let input = crate::parse("CCO").unwrap();
        let mut mol = super::super::parse::from_chematic(&input).unwrap();
        mol.set_rings(Vec::new());
        let ctx = Ctx {
            mol: &mol,
            cip: None,
            n: mol.atoms.len(),
            dmat: Vec::new(),
        };
        let map = BTreeMap::from([(0, P::new(1.0, 1.0)), (1, P::new(1.0, 1.0))]);
        let error = EFrag::from_coord_map(&ctx, &map).err().unwrap();
        assert!(error.0.contains("zero length vector"));
        assert_eq!(map.len(), 2);
        assert_eq!(map[&0], map[&1]);
    }
}
