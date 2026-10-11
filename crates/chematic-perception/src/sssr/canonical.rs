// ---------------------------------------------------------------------------
// Canonical tie-break (determinism, independent of atom-numbering)
// ---------------------------------------------------------------------------

/// A cheap, self-contained (no cross-crate dependency) approximation of
/// canonical atom ranking, used only to make candidate-cycle tie-breaking
/// deterministic: the same molecular graph always produces the same SSSR,
/// regardless of how its atoms happen to be numbered by the parser (SMILES
/// traversal order, SDF atom-block order, etc). This is a local
/// Weisfeiler-Leman-style refinement (seed on (element, ring degree, charge,
/// aromatic), then repeatedly fold in each atom's sorted ring-bond-aware
/// neighbor keys) —
/// it does not aim for full canonical-labeling discriminating power (ties
/// among genuinely symmetric atoms are expected and fine; the point is
/// input-order-independence, not maximal refinement).
fn canonical_atom_ranks(mol: &Molecule) -> Vec<u64> {
    // Same refinement as `canonical_atom_ranks_reference`, with the eligible
    // adjacency flattened once and scratch buffers reused across rounds.
    let n = mol.atom_count();
    let mut start = Vec::with_capacity(n + 1);
    let mut nbrs: Vec<(u32, u64)> = Vec::with_capacity(2 * mol.bond_count());
    for i in 0..n {
        start.push(nbrs.len());
        for (nb, bidx) in mol.neighbors(AtomIdx(i as u32)) {
            let order = mol.bond(bidx).order;
            if is_ring_eligible(order) {
                nbrs.push((nb.0, order.order_int() as u64));
            }
        }
    }
    start.push(nbrs.len());
    let mut keys: Vec<u64> = (0..n)
        .map(|i| {
            let atom = mol.atom(AtomIdx(i as u32));
            let z = atom.element.atomic_number() as u64;
            let degree = (start[i + 1] - start[i]) as u64;
            let charge = (atom.charge as i64 + 8) as u64;
            let aromatic = u64::from(atom.aromatic);
            (z << 24) | (degree << 16) | (charge << 8) | aromatic
        })
        .collect();
    let mut next = vec![0u64; n];
    let mut scratch: Vec<u64> = Vec::new();
    const ROUNDS: usize = 8;
    for _ in 0..ROUNDS {
        for i in 0..n {
            scratch.clear();
            scratch.extend(
                nbrs[start[i]..start[i + 1]]
                    .iter()
                    .map(|&(nb, ord)| keys[nb as usize].wrapping_mul(257).wrapping_add(ord)),
            );
            scratch.sort_unstable();
            let mut h = keys[i];
            for &nk in &scratch {
                h = h.wrapping_mul(1_000_003).wrapping_add(nk);
            }
            next[i] = h;
        }
        std::mem::swap(&mut keys, &mut next);
    }
    keys
}

/// [`canonical_atom_ranks`] restricted to what the SSSR tie-break reads: the
/// returned values are identical for every *ring* atom (an endpoint of a
/// `cyclic` bond); other entries are unspecified.
///
/// After `ROUNDS` refinement rounds a ring atom's key depends only on atoms
/// within `ROUNDS` bonds of it, so round `t` only needs the atoms within
/// `ROUNDS - t` bonds of some ring atom (a multi-source BFS, whose queue is
/// already ordered by that distance). Far acyclic tails are skipped.
fn canonical_ring_atom_ranks(mol: &Molecule, cyclic: &[bool]) -> Vec<u64> {
    const ROUNDS: usize = 8;
    let n = mol.atom_count();
    let mut start = Vec::with_capacity(n + 1);
    let mut nbrs: Vec<(u32, u64)> = Vec::with_capacity(2 * mol.bond_count());
    for i in 0..n {
        start.push(nbrs.len());
        for (nb, bidx) in mol.neighbors(AtomIdx(i as u32)) {
            let order = mol.bond(bidx).order;
            if is_ring_eligible(order) {
                nbrs.push((nb.0, order.order_int() as u64));
            }
        }
    }
    start.push(nbrs.len());

    const FAR: u32 = u32::MAX;
    let mut dist = vec![FAR; n];
    let mut queue: Vec<u32> = Vec::with_capacity(n);
    for (bidx, bond) in mol.bonds() {
        if cyclic[bidx.0 as usize] {
            for a in [bond.atom1.0, bond.atom2.0] {
                if dist[a as usize] == FAR {
                    dist[a as usize] = 0;
                    queue.push(a);
                }
            }
        }
    }
    let mut head = 0;
    while head < queue.len() {
        let cur = queue[head] as usize;
        head += 1;
        let d = dist[cur] + 1;
        if d as usize > ROUNDS {
            continue;
        }
        for &(nb, _) in &nbrs[start[cur]..start[cur + 1]] {
            if dist[nb as usize] == FAR {
                dist[nb as usize] = d;
                queue.push(nb);
            }
        }
    }
    // level_end[k] = number of queued atoms with distance <= k.
    let mut level_end = [0usize; ROUNDS + 1];
    for (k, end) in level_end.iter_mut().enumerate() {
        *end = queue.partition_point(|&a| dist[a as usize] as usize <= k);
    }

    let mut keys = vec![0u64; n];
    for &a in &queue {
        let i = a as usize;
        let atom = mol.atom(AtomIdx(a));
        let z = atom.element.atomic_number() as u64;
        let degree = (start[i + 1] - start[i]) as u64;
        let charge = (atom.charge as i64 + 8) as u64;
        let aromatic = u64::from(atom.aromatic);
        keys[i] = (z << 24) | (degree << 16) | (charge << 8) | aromatic;
    }
    let mut next = vec![0u64; n];
    let mut scratch: Vec<u64> = Vec::new();
    for round in 1..=ROUNDS {
        for &a in &queue[..level_end[ROUNDS - round]] {
            let i = a as usize;
            let row = &nbrs[start[i]..start[i + 1]];
            let nk =
                |&(nb, ord): &(u32, u64)| keys[nb as usize].wrapping_mul(257).wrapping_add(ord);
            let mut h = keys[i];
            match row {
                [] => {}
                [x] => h = h.wrapping_mul(1_000_003).wrapping_add(nk(x)),
                [x, y] => {
                    let (p, q) = (nk(x), nk(y));
                    let (lo, hi) = if p <= q { (p, q) } else { (q, p) };
                    h = h.wrapping_mul(1_000_003).wrapping_add(lo);
                    h = h.wrapping_mul(1_000_003).wrapping_add(hi);
                }
                _ => {
                    let mut buf = [0u64; 8];
                    let vals: &mut [u64] = if row.len() <= 8 {
                        for (slot, x) in buf.iter_mut().zip(row) {
                            *slot = nk(x);
                        }
                        &mut buf[..row.len()]
                    } else {
                        scratch.clear();
                        scratch.extend(row.iter().map(nk));
                        &mut scratch[..]
                    };
                    // Insertion sort: degrees are tiny.
                    for j in 1..vals.len() {
                        let v = vals[j];
                        let mut k = j;
                        while k > 0 && vals[k - 1] > v {
                            vals[k] = vals[k - 1];
                            k -= 1;
                        }
                        vals[k] = v;
                    }
                    for &v in vals.iter() {
                        h = h.wrapping_mul(1_000_003).wrapping_add(v);
                    }
                }
            }
            next[i] = h;
        }
        std::mem::swap(&mut keys, &mut next);
    }
    keys
}

/// Deterministic sort key for a candidate cycle: the sorted multiset of its
/// atoms' canonical ranks, combined order-independently — isomorphic cycles
/// (same molecule, different traversal) always collapse to the same key.
fn canonical_cycle_key(atom_seq: &[AtomIdx], ranks: &[u64]) -> u64 {
    let mut vals: Vec<u64> = atom_seq.iter().map(|a| ranks[a.0 as usize]).collect();
    vals.sort_unstable();
    let mut h: u64 = 0;
    for v in vals {
        h = h.wrapping_mul(1_000_003).wrapping_add(v);
    }
    h
}

/// Order-sensitive canonical key for a ring traversal.
///
/// The former multiset-only key is useful as a compact diagnostic, but it
/// cannot distinguish two different cycles whose atoms have the same refined
/// ranks.  That ambiguity is material for the bounded D2 representative
/// selector: selecting by raw candidate order can then choose a different
/// relevant-cycle family after atom relabeling.  Compare every rotation in
/// both directions and retain the lexicographically smallest sequence of
/// `(atom-rank, bond-order-to-next)` pairs.  Tied genuinely symmetric cycles
/// remain tied and are handled by the existing fail-closed basis checks.
fn canonical_cycle_order_key(
    mol: &Molecule,
    atom_seq: &[AtomIdx],
    ranks: &[u64],
) -> Vec<(u64, u8)> {
    // Output-identical to `canonical_cycle_order_key_reference`: the
    // lexicographically smallest of the 2n traversals, found by comparing
    // traversals lazily instead of materializing each one.
    if atom_seq.is_empty() {
        return Vec::new();
    }
    let n = atom_seq.len();
    let mut sequence: Vec<(u64, u8)> = Vec::with_capacity(n);
    for i in 0..n {
        let current = atom_seq[i];
        let next = atom_seq[(i + 1) % n];
        let order = mol
            .bond_between(current, next)
            .map(|(_, bond)| bond.order.order_int())
            .unwrap_or_default();
        sequence.push((ranks[current.0 as usize], order));
    }
    let elem = |reversed: bool, offset: usize, step: usize| -> (u64, u8) {
        if reversed {
            let index = (offset + n - (step % n)) % n;
            let previous = (index + n - 1) % n;
            (sequence[index].0, sequence[previous].1)
        } else {
            sequence[(offset + step) % n]
        }
    };
    let mut best = (false, 0usize);
    for reversed in [false, true] {
        for offset in 0..n {
            if (reversed, offset) == (false, 0) {
                continue;
            }
            let mut less = false;
            for step in 0..n {
                let c = elem(reversed, offset, step);
                let b = elem(best.0, best.1, step);
                if c != b {
                    less = c < b;
                    break;
                }
            }
            if less {
                best = (reversed, offset);
            }
        }
    }
    (0..n).map(|step| elem(best.0, best.1, step)).collect()
}

// ---------------------------------------------------------------------------
// GF(2) Gaussian elimination
// ---------------------------------------------------------------------------

/// Reduce `cycle` over GF(2) against the current `basis`.
///
/// Each basis entry maps a pivot bond (the minimum BondIdx in that row)
/// to the full row (sorted Vec<BondIdx>).
///
/// Returns the reduced cycle (empty if dependent on existing basis).
fn gf2_reduce(cycle: &[BondIdx], basis: &FxHashMap<BondIdx, Vec<BondIdx>>) -> Vec<BondIdx> {
    let mut current: Vec<BondIdx> = cycle.to_vec();
    while let Some(&pivot) = current.iter().min() {
        match basis.get(&pivot) {
            None => return current, // independent
            // XOR: symmetric difference of the two sorted sets.
            Some(basis_row) => current = sym_diff(&current, basis_row),
        }
    }
    current
}

/// Symmetric difference of two sorted slices (GF(2) addition / XOR for sets).
fn sym_diff(a: &[BondIdx], b: &[BondIdx]) -> Vec<BondIdx> {
    let mut result = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                result.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                result.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                // Both contain this element — XOR removes it.
                i += 1;
                j += 1;
            }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

