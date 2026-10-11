fn build_bond_terms(mol: &Molecule, types: &[u8]) -> Vec<PreparedBond> {
    mol.bonds()
        .filter_map(|(_, bond)| {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            let bt = bond_type_for(types[i], types[j], bond.order);
            mmff94_bond_energy_resolved(bt, types[i], types[j]).map(|(params, _)| PreparedBond {
                i,
                j,
                params,
            })
        })
        .collect()
}

fn mol_neighbors(mol: &Molecule) -> Vec<Vec<usize>> {
    (0..mol.atom_count())
        .map(|j| {
            mol.neighbors(AtomIdx(j as u32))
                .map(|(nb, _)| nb.0 as usize)
                .collect()
        })
        .collect()
}

/// Angle-bend and stretch-bend parameter lists in RDKit's builder order and
/// orientation: centre `j` ascending, neighbour pairs `(i, k)` taken in the
/// order of `adj[j]` (RDKit's adjacency order). The empirical-rule angle
/// constants are not symmetric in floating point under `i <-> k`, so the
/// orientation must be RDKit's.
#[allow(clippy::type_complexity)]
pub(crate) fn rdkit_angle_parts(
    mol: &Molecule,
    adj: &[Vec<usize>],
) -> Result<
    (
        Vec<(usize, usize, usize, f64, f64, bool)>,
        Vec<(usize, usize, usize, f64, f64, f64, f64, f64)>,
    ),
    MinimizerError,
> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let rings = find_sssr(mol).rings().to_vec();
    let angles = build_angle_terms_with(&mmff_mol, &types, &rings, adj)
        .into_iter()
        .map(|a| (a.i, a.j, a.k, a.params.theta0, a.params.ka, a.linear))
        .collect();
    let stbn = build_stretch_bend_terms_with(&mmff_mol, &types, &rings, adj)
        .into_iter()
        .map(|t| {
            (
                t.i, t.j, t.k, t.r0_ij, t.r0_kj, t.theta0, t.kba_ij, t.kba_kj,
            )
        })
        .collect();
    Ok((angles, stbn))
}

fn build_angle_terms(mol: &Molecule, types: &[u8], rings: &[Vec<AtomIdx>]) -> Vec<PreparedAngle> {
    build_angle_terms_with(mol, types, rings, &mol_neighbors(mol))
}

fn build_angle_terms_with(
    mol: &Molecule,
    types: &[u8],
    rings: &[Vec<AtomIdx>],
    adj: &[Vec<usize>],
) -> Vec<PreparedAngle> {
    let mut terms = Vec::new();
    for j_idx in 0..mol.atom_count() {
        let neighbors = &adj[j_idx];
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let Some((bond_ij, _)) = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx])
                else {
                    continue;
                };
                let Some((bond_kj, _)) = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx])
                else {
                    continue;
                };
                let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                if let Some((params, _)) = mmff94_angle_energy_resolved(
                    at,
                    types[i],
                    types[j_idx],
                    types[k],
                    bond_ij.r0,
                    bond_kj.r0,
                    ring_size,
                ) {
                    terms.push(PreparedAngle {
                        i,
                        j: j_idx,
                        k,
                        params,
                        linear: central_type_is_linear(types[j_idx]),
                    });
                }
            }
        }
    }
    terms
}

fn build_torsion_terms(mol: &Molecule, types: &[u8]) -> Vec<PreparedTorsion> {
    let mut terms = Vec::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        let nbrs_j: Vec<usize> = mol
            .neighbors(bond.atom1)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        let nbrs_k: Vec<usize> = mol
            .neighbors(bond.atom2)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                // l == i closes a 3-membered ring: there is no dihedral and
                // RDKit adds no torsion term (#637).
                if l == j || l == i {
                    continue;
                }
                if let Some(params) = mmff94_torsion_term_params(mol, i, j, k, l, types) {
                    terms.push(PreparedTorsion { i, j, k, l, params });
                }
            }
        }
    }
    terms
}

fn build_oop_terms(mol: &Molecule, types: &[u8]) -> Vec<PreparedOop> {
    let mut terms = Vec::new();
    for j in 0..mol.atom_count() {
        if OOP_SP2_TYPES.binary_search(&types[j]).is_err() {
            continue;
        }
        let neighbors: Vec<usize> = mol
            .neighbors(AtomIdx(j as u32))
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        if neighbors.len() != 3 {
            continue;
        }
        if let Some(koop) = mmff94_oop(
            types[j],
            types[neighbors[0]],
            types[neighbors[1]],
            types[neighbors[2]],
        ) {
            // One Wilson angle per neighbour out of the plane of the other
            // two, all with the centre's koop (MMFF eq. 6; RDKit adds the
            // three contributions too).
            for [i, k, l] in oop_permutations(&neighbors) {
                terms.push(PreparedOop { i, j, k, l, koop });
            }
        }
    }
    terms
}

/// The three (i, k, l) choices of a trigonal centre's out-of-plane terms,
/// in RDKit's order: l = third, second, first neighbour.
fn oop_permutations(n: &[usize]) -> [[usize; 3]; 3] {
    [[n[0], n[1], n[2]], [n[0], n[2], n[1]], [n[1], n[2], n[0]]]
}

/// MMFF94 adds no stretch-bend term at a linear central atom (MMFFPROP
/// `linh`; RDKit skips these angles), e.g. alkyne or nitrile carbons (#637).
fn central_type_is_linear(atom_type: u8) -> bool {
    crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(atom_type)
        .is_some_and(|info| info.linear)
}

fn build_stretch_bend_terms(
    mol: &Molecule,
    types: &[u8],
    rings: &[Vec<AtomIdx>],
) -> Vec<PreparedStretchBend> {
    build_stretch_bend_terms_with(mol, types, rings, &mol_neighbors(mol))
}

fn build_stretch_bend_terms_with(
    mol: &Molecule,
    types: &[u8],
    rings: &[Vec<AtomIdx>],
    adj: &[Vec<usize>],
) -> Vec<PreparedStretchBend> {
    let mut terms = Vec::new();
    for j_idx in 0..mol.atom_count() {
        if central_type_is_linear(types[j_idx]) {
            continue;
        }
        let neighbors = &adj[j_idx];
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let Some((kba_ij, kba_kj)) = mmff94_stbn_oriented(
                    stretch_bend_type_for(at, types[i], types[k], bt_ij, bt_kj),
                    types[i],
                    types[j_idx],
                    types[k],
                    bt_ij,
                    bt_kj,
                    mol.atom(AtomIdx(i as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(j_idx as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(k as u32)).element.atomic_number(),
                ) else {
                    continue;
                };
                let Some((bond_ij, _)) = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx])
                else {
                    continue;
                };
                let Some((bond_kj, _)) = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx])
                else {
                    continue;
                };
                let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                let Some((angle, _)) = mmff94_angle_energy_resolved(
                    at,
                    types[i],
                    types[j_idx],
                    types[k],
                    bond_ij.r0,
                    bond_kj.r0,
                    ring_size,
                ) else {
                    continue;
                };
                terms.push(PreparedStretchBend {
                    i,
                    j: j_idx,
                    k,
                    kba_ij,
                    kba_kj,
                    r0_ij: bond_ij.r0,
                    r0_kj: bond_kj.r0,
                    theta0: angle.theta0,
                });
            }
        }
    }
    terms
}

fn prepared_bond_energy(coords: &[[f64; 3]], terms: &[PreparedBond]) -> f64 {
    const KB_CONV: f64 = MDYNE_A_TO_KCAL_MOL;
    const CS: f64 = 2.0;
    terms
        .iter()
        .map(|term| {
            let dr = dist(coords[term.i], coords[term.j]) - term.params.r0;
            let cubic = 1.0 - CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
            (KB_CONV * term.params.kb / 2.0) * dr * dr * cubic
        })
        .sum()
}

fn prepared_angle_energy(coords: &[[f64; 3]], terms: &[PreparedAngle]) -> f64 {
    terms
        .iter()
        .map(|term| {
            angle_bend_energy(
                cos_angle(coords[term.i], coords[term.j], coords[term.k]),
                term.params.theta0,
                term.params.ka,
                term.linear,
            )
        })
        .sum()
}

/// RDKit's `calcAngleBendEnergy`: `0.5 c2 ka Δθ² (1 + cb Δθ)` (Δθ in
/// degrees), or `143.9325 ka (1 + cos θ)` at a linear central atom.
fn angle_bend_energy(cos_theta: f64, theta0: f64, ka: f64, linear: bool) -> f64 {
    if linear {
        return MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cos_theta);
    }
    let dt = cos_theta.dacos() / DEG2RAD - theta0;
    0.5 * MMFF_ANGLE_CONV * ka * dt * dt * (1.0 - MMFF_ANGLE_CUBIC * dt)
}

fn prepared_stretch_bend_energy(coords: &[[f64; 3]], terms: &[PreparedStretchBend]) -> f64 {
    const CONV: f64 = MMFF_STBN_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    terms
        .iter()
        .map(|term| {
            let dr_ij = dist(coords[term.i], coords[term.j]) - term.r0_ij;
            let dr_kj = dist(coords[term.k], coords[term.j]) - term.r0_kj;
            let dtheta = cos_angle(coords[term.i], coords[term.j], coords[term.k]).dacos()
                * RAD_TO_DEG
                - term.theta0;
            CONV * (term.kba_ij * dr_ij + term.kba_kj * dr_kj) * dtheta
        })
        .sum()
}

fn prepared_torsion_energy(coords: &[[f64; 3]], terms: &[PreparedTorsion]) -> f64 {
    terms
        .iter()
        .map(|term| {
            let phi = dihedral(
                coords[term.i],
                coords[term.j],
                coords[term.k],
                coords[term.l],
            );
            0.5 * term.params.v1 * (1.0 + phi.dcos())
                + 0.5 * term.params.v2 * (1.0 - (2.0 * phi).dcos())
                + 0.5 * term.params.v3 * (1.0 + (3.0 * phi).dcos())
        })
        .sum()
}

fn prepared_torsion_gradient(coords: &[[f64; 3]], terms: &[PreparedTorsion]) -> Vec<[f64; 3]> {
    let mut gradient = vec![[0.0; 3]; coords.len()];
    for term in terms {
        let p1 = coords[term.i];
        let p2 = coords[term.j];
        let p3 = coords[term.k];
        let p4 = coords[term.l];
        let b1 = sub3(p1, p2);
        let b2 = sub3(p3, p2);
        let b3 = sub3(p4, p3);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let b2_sq = dot3(b2, b2);
        let b2_len = b2_sq.sqrt();
        let n1_sq = dot3(n1, n1);
        let n2_sq = dot3(n2, n2);
        if b2_len < 1e-12 || n1_sq < 1e-24 || n2_sq < 1e-24 {
            continue;
        }

        let phi = dihedral(p1, p2, p3, p4);
        let d_ed_phi = -0.5 * term.params.v1 * phi.dsin() + term.params.v2 * (2.0 * phi).dsin()
            - 1.5 * term.params.v3 * (3.0 * phi).dsin();
        let dphi = dual_dihedral(p1, p2, p3, p4);
        for atom in 0..4 {
            let target = &mut gradient[[term.i, term.j, term.k, term.l][atom]];
            for (axis, component) in target.iter_mut().enumerate() {
                *component += d_ed_phi * dphi.grad[atom * 3 + axis];
            }
        }
    }
    gradient
}

fn prepared_oop_energy(coords: &[[f64; 3]], terms: &[PreparedOop]) -> f64 {
    const CONV: f64 = MMFF_ANGLE_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    terms
        .iter()
        .map(|term| {
            let n = cross(
                sub3(coords[term.i], coords[term.j]),
                sub3(coords[term.k], coords[term.j]),
            );
            let l = sub3(coords[term.l], coords[term.j]);
            let n_len = dot3(n, n).sqrt();
            let l_len = dot3(l, l).sqrt();
            if n_len < 1e-12 || l_len < 1e-12 {
                return 0.0;
            }
            let chi = (dot3(n, l) / (n_len * l_len)).clamp(-1.0, 1.0).dasin();
            (CONV * term.koop / 2.0) * (chi * RAD_TO_DEG).powi(2)
        })
        .sum()
}

fn prepared_oop_gradient(coords: &[[f64; 3]], terms: &[PreparedOop]) -> Vec<[f64; 3]> {
    const CONV: f64 = MMFF_ANGLE_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut gradient = vec![[0.0; 3]; coords.len()];
    for term in terms {
        let n = cross(
            sub3(coords[term.i], coords[term.j]),
            sub3(coords[term.k], coords[term.j]),
        );
        let l = sub3(coords[term.l], coords[term.j]);
        if dot3(n, n).sqrt() < 1e-12 || dot3(l, l).sqrt() < 1e-12 {
            continue;
        }
        let nd = dual_cross(
            dual_sub(dual_point(coords[term.i], 0), dual_point(coords[term.j], 3)),
            dual_sub(dual_point(coords[term.k], 6), dual_point(coords[term.j], 3)),
        );
        let ld = dual_sub(dual_point(coords[term.l], 9), dual_point(coords[term.j], 3));
        let denominator = dual_mul(dual_sqrt(dual_dot(nd, nd)), dual_sqrt(dual_dot(ld, ld)));
        let sin_chi = dual_div(dual_dot(nd, ld), denominator);
        // asin(sin χ) has an unbounded Cartesian derivative as |sin χ|→1;
        // fail closed at the Wilson angle's ±90° branch point.
        if sin_chi.value.abs() >= 1.0 - 1e-10 {
            continue;
        }
        let chi = dual_asin_clamped(sin_chi);
        let d_ed_chi = CONV * term.koop * RAD_TO_DEG * RAD_TO_DEG * chi.value;
        for atom in 0..4 {
            let target = &mut gradient[[term.i, term.j, term.k, term.l][atom]];
            for (axis, component) in target.iter_mut().enumerate() {
                *component += d_ed_chi * chi.grad[atom * 3 + axis];
            }
        }
    }
    gradient
}

fn build_nonbonded_pairs(
    mol: &Molecule,
    types: &[u8],
    charges: &[f64],
) -> (VdwPairs, ElectrostaticPairs) {
    let n = mol.atom_count();
    let mut excluded = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
        for (neighbor, _) in mol.neighbors(bond.atom1) {
            let k = neighbor.0 as usize;
            excluded.insert((k.min(j), k.max(j)));
        }
        for (neighbor, _) in mol.neighbors(bond.atom2) {
            let k = neighbor.0 as usize;
            excluded.insert((i.min(k), i.max(k)));
        }
    }

    let mut one_four = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        for (neighbor_j, _) in mol.neighbors(bond.atom1) {
            let i = neighbor_j.0 as usize;
            if i == k {
                continue;
            }
            for (neighbor_k, _) in mol.neighbors(bond.atom2) {
                let l = neighbor_k.0 as usize;
                if l == j {
                    continue;
                }
                let pair = (i.min(l), i.max(l));
                if !excluded.contains(&pair) {
                    one_four.insert(pair);
                }
            }
        }
    }

    let mut vdw_pairs = Vec::new();
    let mut electrostatic_pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            if let Some((r_star, epsilon)) = mmff94_vdw_combined(types[i], types[j])
                && r_star > 0.0
                && epsilon > 0.0
            {
                vdw_pairs.push(PreparedVdwPair {
                    i,
                    j,
                    r_star,
                    epsilon,
                });
            }
            let scale = if one_four.contains(&(i, j)) {
                0.75
            } else {
                1.0
            };
            electrostatic_pairs.push(PreparedElectrostaticPair {
                i,
                j,
                scale_charge_product: scale * charges[i] * charges[j],
            });
        }
    }
    (vdw_pairs, electrostatic_pairs)
}

/// Rebuild the coordinate-dependent vdW candidate list for one evaluation.
///
/// The topology-prepared pair map remains immutable; only spatial bin
/// membership changes with coordinates. A two-bin halo is conservative for
/// the 10 Å cutoff even at floating-point cell boundaries. Non-finite input
/// deliberately falls back to the reference all-pair traversal so error
/// behavior is not changed by integer bin conversion.
fn prepared_vdw_neighbor_list<'a>(
    coords: &[[f64; 3]],
    pairs: &'a [PreparedVdwPair],
) -> Vec<&'a PreparedVdwPair> {
    const CUTOFF: f64 = 10.0;
    let valid = coords
        .iter()
        .all(|point| point.iter().all(|value| value.is_finite()));
    if !valid {
        return pairs.iter().collect();
    }
    let cell_size = CUTOFF * (1.0 + 1e-12);
    let mut bins: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (index, point) in coords.iter().enumerate() {
        let key = [
            (point[0] / cell_size).floor() as i64,
            (point[1] / cell_size).floor() as i64,
            (point[2] / cell_size).floor() as i64,
        ];
        bins.entry(key).or_default().push(index);
    }
    let mut pair_map: HashMap<(usize, usize), &PreparedVdwPair> = HashMap::new();
    for pair in pairs {
        pair_map.insert((pair.i, pair.j), pair);
    }
    let mut active = Vec::new();
    for (i, point) in coords.iter().enumerate() {
        let key = [
            (point[0] / cell_size).floor() as i64,
            (point[1] / cell_size).floor() as i64,
            (point[2] / cell_size).floor() as i64,
        ];
        for dx in -2..=2 {
            for dy in -2..=2 {
                for dz in -2..=2 {
                    let Some(indices) = bins.get(&[key[0] + dx, key[1] + dy, key[2] + dz]) else {
                        continue;
                    };
                    for &j in indices {
                        if j <= i {
                            continue;
                        }
                        if let Some(pair) = pair_map.get(&(i, j)) {
                            let r = dist(*point, coords[j]);
                            if r <= CUTOFF {
                                active.push(*pair);
                            }
                        }
                    }
                }
            }
        }
    }
    active
}

/// vdW energy over every prepared pair within RDKit's `nonBondedThresh`
/// (100 Å). The buffered 14-7 form is finite at any distance, coincident
/// atoms included, as in RDKit.
fn vdw_energy_pairs(coords: &[[f64; 3]], pairs: &[PreparedVdwPair]) -> f64 {
    pairs
        .iter()
        .map(|pair| {
            let r = dist(coords[pair.i], coords[pair.j]);
            if r <= NONBONDED_THRESH {
                mmff94_vdw_energy_value(r, pair.r_star, pair.epsilon)
            } else {
                0.0
            }
        })
        .sum()
}

/// vdW energy over the pairs within the bounded 10 Å cutoff.
fn vdw_energy_cutoff_pairs(coords: &[[f64; 3]], pairs: &[PreparedVdwPair]) -> f64 {
    let mut energy = 0.0;
    let mut add_pair = |pair: &PreparedVdwPair| {
        let r = dist(coords[pair.i], coords[pair.j]);
        if r > 10.0 {
            return;
        }
        if r > 0.01 {
            energy += mmff94_vdw_energy_value(r, pair.r_star, pair.epsilon);
        }
    };
    if coords.len() <= DIRECT_VDW_ATOM_THRESHOLD {
        for pair in pairs {
            add_pair(pair);
        }
    } else {
        for pair in prepared_vdw_neighbor_list(coords, pairs) {
            add_pair(pair);
        }
    }
    energy
}

fn electrostatic_energy_pairs(coords: &[[f64; 3]], pairs: &[PreparedElectrostaticPair]) -> f64 {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    pairs
        .iter()
        .map(|pair| {
            let r = dist(coords[pair.i], coords[pair.j]);
            if r <= NONBONDED_THRESH {
                pair.scale_charge_product * COULOMB / (r + DELTA)
            } else {
                0.0
            }
        })
        .sum()
}

fn prepared_electrostatic_neighbor_list<'a>(
    coords: &[[f64; 3]],
    pairs: &'a [PreparedElectrostaticPair],
) -> Vec<&'a PreparedElectrostaticPair> {
    const CUTOFF: f64 = 10.0;
    let valid = coords
        .iter()
        .all(|point| point.iter().all(|value| value.is_finite()));
    if !valid {
        return pairs.iter().collect();
    }
    let cell_size = CUTOFF * (1.0 + 1e-12);
    let mut bins: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (index, point) in coords.iter().enumerate() {
        let key = [
            (point[0] / cell_size).floor() as i64,
            (point[1] / cell_size).floor() as i64,
            (point[2] / cell_size).floor() as i64,
        ];
        bins.entry(key).or_default().push(index);
    }
    let mut pair_map: HashMap<(usize, usize), &PreparedElectrostaticPair> = HashMap::new();
    for pair in pairs {
        pair_map.insert((pair.i, pair.j), pair);
    }
    let mut active = Vec::new();
    for (i, point) in coords.iter().enumerate() {
        let key = [
            (point[0] / cell_size).floor() as i64,
            (point[1] / cell_size).floor() as i64,
            (point[2] / cell_size).floor() as i64,
        ];
        for dx in -2..=2 {
            for dy in -2..=2 {
                for dz in -2..=2 {
                    let Some(indices) = bins.get(&[key[0] + dx, key[1] + dy, key[2] + dz]) else {
                        continue;
                    };
                    for &j in indices {
                        if j <= i {
                            continue;
                        }
                        if let Some(pair) = pair_map.get(&(i, j))
                            && dist(*point, coords[j]) <= CUTOFF
                        {
                            active.push(*pair);
                        }
                    }
                }
            }
        }
    }
    active
}

fn electrostatic_energy_cutoff_pairs(
    coords: &[[f64; 3]],
    pairs: &[PreparedElectrostaticPair],
) -> f64 {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    prepared_electrostatic_neighbor_list(coords, pairs)
        .iter()
        .map(|pair| {
            pair.scale_charge_product * COULOMB / (dist(coords[pair.i], coords[pair.j]) + DELTA)
        })
        .sum()
}

fn prepared_cutoff_nonbonded_gradient(
    coords: &[[f64; 3]],
    vdw_pairs: &[PreparedVdwPair],
    electrostatic_pairs: &[PreparedElectrostaticPair],
) -> Vec<[f64; 3]> {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    let mut result = vec![[0.0; 3]; coords.len()];
    let add_pair =
        |gradient: &mut [[f64; 3]], i: usize, j: usize, radial: f64, delta: [f64; 3], r: f64| {
            if r <= 1e-12 {
                return;
            }
            let value = [
                radial * delta[0] / r,
                radial * delta[1] / r,
                radial * delta[2] / r,
            ];
            for axis in 0..3 {
                gradient[i][axis] += value[axis];
                gradient[j][axis] -= value[axis];
            }
        };
    for pair in prepared_vdw_neighbor_list(coords, vdw_pairs) {
        let delta = [
            coords[pair.i][0] - coords[pair.j][0],
            coords[pair.i][1] - coords[pair.j][1],
            coords[pair.i][2] - coords[pair.j][2],
        ];
        let r = dist(coords[pair.i], coords[pair.j]);
        if r > 0.01 && r <= 10.0 {
            let radial = mmff94_vdw_radial_derivative(r, pair.r_star, pair.epsilon);
            add_pair(&mut result, pair.i, pair.j, radial, delta, r);
        }
    }
    for pair in prepared_electrostatic_neighbor_list(coords, electrostatic_pairs) {
        let delta = [
            coords[pair.i][0] - coords[pair.j][0],
            coords[pair.i][1] - coords[pair.j][1],
            coords[pair.i][2] - coords[pair.j][2],
        ];
        let r = dist(coords[pair.i], coords[pair.j]);
        let radial = -pair.scale_charge_product * COULOMB / (r + DELTA).powi(2);
        add_pair(&mut result, pair.i, pair.j, radial, delta, r);
    }
    result
}

fn prepared_nonbonded_gradient(
    coords: &[[f64; 3]],
    vdw_pairs: &[PreparedVdwPair],
    electrostatic_pairs: &[PreparedElectrostaticPair],
) -> Vec<[f64; 3]> {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    let mut gradient = vec![[0.0; 3]; coords.len()];
    let add_pair =
        |gradient: &mut [[f64; 3]], i: usize, j: usize, radial: f64, delta: [f64; 3], r: f64| {
            if r <= 1e-12 {
                return;
            }
            let value = [
                radial * delta[0] / r,
                radial * delta[1] / r,
                radial * delta[2] / r,
            ];
            for axis in 0..3 {
                gradient[i][axis] += value[axis];
                gradient[j][axis] -= value[axis];
            }
        };
    let mut add_vdw = |pair: &PreparedVdwPair| {
        let delta = [
            coords[pair.i][0] - coords[pair.j][0],
            coords[pair.i][1] - coords[pair.j][1],
            coords[pair.i][2] - coords[pair.j][2],
        ];
        let r = dist(coords[pair.i], coords[pair.j]);
        if r <= NONBONDED_THRESH {
            let radial = mmff94_vdw_radial_derivative(r, pair.r_star, pair.epsilon);
            add_pair(&mut gradient, pair.i, pair.j, radial, delta, r);
        }
    };
    for pair in vdw_pairs {
        add_vdw(pair);
    }
    for pair in electrostatic_pairs {
        let delta = [
            coords[pair.i][0] - coords[pair.j][0],
            coords[pair.i][1] - coords[pair.j][1],
            coords[pair.i][2] - coords[pair.j][2],
        ];
        let r = dist(coords[pair.i], coords[pair.j]);
        if r > NONBONDED_THRESH {
            continue;
        }
        let radial = -pair.scale_charge_product * COULOMB / (r + DELTA).powi(2);
        add_pair(&mut gradient, pair.i, pair.j, radial, delta, r);
    }
    gradient
}

fn prepared_bond_angle_gradient(
    coords: &[[f64; 3]],
    bonds: &[PreparedBond],
    angles: &[PreparedAngle],
) -> Vec<[f64; 3]> {
    const KB_CONV: f64 = MDYNE_A_TO_KCAL_MOL;
    const CS: f64 = 2.0;
    const KA_CONV: f64 = MMFF_ANGLE_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut gradient = vec![[0.0; 3]; coords.len()];
    let add_pair =
        |gradient: &mut [[f64; 3]], i: usize, j: usize, radial: f64, delta: [f64; 3], r: f64| {
            if r <= 1e-12 {
                return;
            }
            for axis in 0..3 {
                let value = radial * delta[axis] / r;
                gradient[i][axis] += value;
                gradient[j][axis] -= value;
            }
        };
    for term in bonds {
        let delta = [
            coords[term.i][0] - coords[term.j][0],
            coords[term.i][1] - coords[term.j][1],
            coords[term.i][2] - coords[term.j][2],
        ];
        let r = dist(coords[term.i], coords[term.j]);
        let dr = r - term.params.r0;
        let prefactor = KB_CONV * term.params.kb / 2.0;
        let radial =
            prefactor * (2.0 * dr - 3.0 * CS * dr * dr + (7.0 / 3.0) * CS * CS * dr * dr * dr);
        add_pair(&mut gradient, term.i, term.j, radial, delta, r);
    }
    for term in angles {
        let u = [
            coords[term.i][0] - coords[term.j][0],
            coords[term.i][1] - coords[term.j][1],
            coords[term.i][2] - coords[term.j][2],
        ];
        let v = [
            coords[term.k][0] - coords[term.j][0],
            coords[term.k][1] - coords[term.j][1],
            coords[term.k][2] - coords[term.j][2],
        ];
        let u2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        let lu = u2.sqrt();
        let lv = v2.sqrt();
        if lu <= 1e-12 || lv <= 1e-12 {
            continue;
        }
        let cos_theta = cos_angle(coords[term.i], coords[term.j], coords[term.k]);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        if sin_theta <= 1e-12 {
            continue;
        }
        let d_e_d_cos = if term.linear {
            MDYNE_A_TO_KCAL_MOL * term.params.ka
        } else {
            let theta = cos_theta.dacos() * RAD_TO_DEG;
            let dt = theta - term.params.theta0;
            let prefactor = KA_CONV * term.params.ka / 2.0;
            let d_e_d_theta = prefactor * (2.0 * dt - 3.0 * MMFF_ANGLE_CUBIC * dt * dt);
            -RAD_TO_DEG * d_e_d_theta / sin_theta
        };
        let inv = 1.0 / (lu * lv);
        let dcos_di = [
            v[0] * inv - cos_theta * u[0] / u2,
            v[1] * inv - cos_theta * u[1] / u2,
            v[2] * inv - cos_theta * u[2] / u2,
        ];
        let dcos_dk = [
            u[0] * inv - cos_theta * v[0] / v2,
            u[1] * inv - cos_theta * v[1] / v2,
            u[2] * inv - cos_theta * v[2] / v2,
        ];
        for axis in 0..3 {
            let gi = d_e_d_cos * dcos_di[axis];
            let gk = d_e_d_cos * dcos_dk[axis];
            gradient[term.i][axis] += gi;
            gradient[term.k][axis] += gk;
            gradient[term.j][axis] -= gi + gk;
        }
    }
    gradient
}

fn prepared_stretch_bend_gradient(
    coords: &[[f64; 3]],
    terms: &[PreparedStretchBend],
) -> Vec<[f64; 3]> {
    const CONV: f64 = MMFF_STBN_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut gradient = vec![[0.0; 3]; coords.len()];
    for term in terms {
        let a = coords[term.i];
        let center = coords[term.j];
        let c = coords[term.k];
        let u = [a[0] - center[0], a[1] - center[1], a[2] - center[2]];
        let v = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
        let u2 = dot3(u, u);
        let v2 = dot3(v, v);
        let lu = u2.sqrt();
        let lv = v2.sqrt();
        if lu <= 1e-12 || lv <= 1e-12 {
            continue;
        }
        let r_ij = lu;
        let r_kj = lv;
        let cos_theta = dot3(u, v) / (lu * lv);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        if sin_theta <= 1e-12 {
            continue;
        }
        let dtheta = cos_theta.dacos() * RAD_TO_DEG - term.theta0;
        let stretch = term.kba_ij * (r_ij - term.r0_ij) + term.kba_kj * (r_kj - term.r0_kj);
        let d_e_d_r = CONV * dtheta;
        let d_e_d_theta = CONV * stretch;
        let d_e_d_cos = -RAD_TO_DEG * d_e_d_theta / sin_theta;
        let dcos_di = [
            v[0] / (lu * lv) - cos_theta * u[0] / u2,
            v[1] / (lu * lv) - cos_theta * u[1] / u2,
            v[2] / (lu * lv) - cos_theta * u[2] / u2,
        ];
        let dcos_dk = [
            u[0] / (lu * lv) - cos_theta * v[0] / v2,
            u[1] / (lu * lv) - cos_theta * v[1] / v2,
            u[2] / (lu * lv) - cos_theta * v[2] / v2,
        ];
        let mut gi = [0.0; 3];
        let mut gk = [0.0; 3];
        for axis in 0..3 {
            gi[axis] = d_e_d_r * term.kba_ij * u[axis] / r_ij + d_e_d_cos * dcos_di[axis];
            gk[axis] = d_e_d_r * term.kba_kj * v[axis] / r_kj + d_e_d_cos * dcos_dk[axis];
            gradient[term.i][axis] += gi[axis];
            gradient[term.k][axis] += gk[axis];
            gradient[term.j][axis] -= gi[axis] + gk[axis];
        }
    }
    gradient
}

