const DEFAULT_LOG_FREQ: f64 = -5.0;

#[inline]
fn lookup_score(hash: u64) -> f64 {
    let idx = FRAGMENT_SCORES.partition_point(|&(h, _)| h < hash);
    if idx < FRAGMENT_SCORES.len() && FRAGMENT_SCORES[idx].0 == hash {
        FRAGMENT_SCORES[idx].1 as f64 / 1000.0
    } else {
        DEFAULT_LOG_FREQ
    }
}

// ---------------------------------------------------------------------------
// Complexity corrections
// ---------------------------------------------------------------------------

fn complexity_penalty_with(mol: &Molecule, rb: &RingBundle) -> f64 {
    let n = mol.atom_count() as f64;
    if n == 0.0 {
        return 0.0;
    }

    let spiro = rb.num_spiro_atoms as f64;
    let bridge = rb.num_bridgehead_atoms as f64;
    let rings = rb.ring_count as f64;
    let stereo = num_stereocenters(mol) as f64;

    // count_macrorings inlined using rb.ring_count to avoid a second find_sssr call.
    let macro_count = if rb.ring_count > 0 && n / rings > 8.0 {
        1.0
    } else {
        0.0
    };

    let bond_n = mol.bond_count() as f64;
    let ring_bond_frac = if bond_n > 0.0 {
        ring_bond_count(mol) as f64 / bond_n
    } else {
        0.0
    };

    let penalty = spiro * 0.25
        + bridge * 0.35
        + macro_count * 0.30
        + stereo * 0.10
        + (rings - 1.0).max(0.0) * 0.05
        + ring_bond_frac * 0.50;

    let size_penalty = (n / 10.0 - 1.0).max(0.0) * 0.05;

    penalty + size_penalty
}

fn ring_bond_count(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    let mut in_ring = vec![false; mol.bond_count()];
    let mut visited = vec![false; n];
    let mut bond_idx_map: HashMap<(usize, usize), usize> = HashMap::new();
    for (bidx, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        bond_idx_map.insert((i.min(j), i.max(j)), bidx.0 as usize);
    }

    let adj: Vec<Vec<usize>> = (0..n)
        .map(|i| {
            mol.neighbors(AtomIdx(i as u32))
                .map(|(nb, _)| nb.0 as usize)
                .collect()
        })
        .collect();

    let mut stack: Vec<(usize, usize, usize)> = Vec::new();
    for root in 0..n {
        if visited[root] {
            continue;
        }
        stack.push((root, usize::MAX, 0));
        while let Some((cur, par, d)) = stack.pop() {
            if visited[cur] {
                continue;
            }
            visited[cur] = true;
            let _ = d;
            for &nb in &adj[cur] {
                if nb == par {
                    continue;
                }
                if visited[nb] {
                    let key = (cur.min(nb), cur.max(nb));
                    if let Some(&bi) = bond_idx_map.get(&key) {
                        in_ring[bi] = true;
                    }
                } else {
                    stack.push((nb, cur, 0));
                }
            }
        }
    }
    in_ring.iter().filter(|&&r| r).count()
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the Synthetic Accessibility Score for `mol`.
///
/// Returns a value in \[1, 10\]: 1 = trivially easy, 10 = very hard to synthesize.
///
/// Fragment environments use the same FNV-1a hash scheme as ECFP fingerprints
/// (`chematic_fp::morgan_fp_counts`), so the fragment table is directly
/// compatible with ECFP4 lookups.
pub fn sa_score(mol: &Molecule) -> f64 {
    sa_score_with_bundle(mol, &ring_bundle(mol))
}

/// SA Score with pre-computed ring values from [`ring_bundle`] — avoids 4 redundant `find_sssr` calls.
pub fn sa_score_with_bundle(mol: &Molecule, rb: &RingBundle) -> f64 {
    let n = mol.atom_count();
    if n == 0 {
        return 10.0;
    }

    let counts = morgan_fp_counts(mol, 2);
    let total_frags: u32 = counts.values().sum();
    if total_frags == 0 {
        return 10.0;
    }

    let total_score: f64 = counts
        .iter()
        .map(|(&hash, &count)| lookup_score(hash) * count as f64)
        .sum();
    let fscore = total_score / total_frags as f64;

    let cscore = complexity_penalty_with(mol, rb);

    let raw = fscore * 2.0 - cscore * 4.5;
    let sa = 1.0 + (-raw).max(0.0);
    sa.clamp(1.0, 10.0)
}

