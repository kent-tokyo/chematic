// ─── Energy components ───────────────────────────────────────────────────────

fn total_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    charges: &[f64],
    rings: &[Vec<AtomIdx>],
) -> f64 {
    bond_energy(mol, coords, types)
        + angle_energy(mol, coords, types, rings)
        + stretch_bend_energy(mol, coords, types, rings)
        + torsion_energy(mol, coords, types)
        + oop_energy(mol, coords, types)
        + vdw_energy(mol, coords, types)
        + elec_energy(mol, coords, charges)
}

/// Stretch-bend coupling (Halgren MMFF.V eq. 4)
/// E_sb = 2.512076 × (kba_ijk × Δr_ij + kba_kji × Δr_kj) × Δθ   [kcal/mol, Δθ in degrees]
fn stretch_bend_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    rings: &[Vec<AtomIdx>],
) -> f64 {
    const CONV: f64 = MMFF_STBN_CONV; // md/Å → kcal/(mol·Å·deg)
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    const KB_CONV: f64 = MDYNE_A_TO_KCAL_MOL;
    const CS: f64 = 2.0;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() < 2 || central_type_is_linear(types[j_idx]) {
            continue;
        }
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let sbt = stretch_bend_type_for(at, types[i], types[k], bt_ij, bt_kj);
                if let Some((kba_ijk, kba_kji)) = mmff94_stbn_oriented(
                    sbt,
                    types[i],
                    types[j_idx],
                    types[k],
                    bt_ij,
                    bt_kj,
                    mol.atom(AtomIdx(i as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(j_idx as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(k as u32)).element.atomic_number(),
                ) {
                    let bond_ij = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx]);
                    let bond_kj = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx]);
                    // Δr_ij
                    let r_ij = dist(coords[i], coords[j_idx]);
                    let dr_ij = bond_ij.map(|(p, _)| r_ij - p.r0).unwrap_or(0.0);
                    // Δr_kj
                    let r_kj = dist(coords[k], coords[j_idx]);
                    let dr_kj = bond_kj.map(|(p, _)| r_kj - p.r0).unwrap_or(0.0);
                    // Δθ in degrees -- both flanking bonds must resolve to feed
                    // r0_ij/r0_jk into the angle empirical rule if it's needed
                    // (matches RDKit's real `getMMFFAngleBendParams`, which
                    // requires both `getMMFFBondStretchParams` calls to succeed
                    // before even attempting the empirical path).
                    let cos_t = cos_angle(coords[i], coords[j_idx], coords[k]);
                    if let (Some((rij, _)), Some((rkj, _))) = (bond_ij, bond_kj) {
                        let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                        if let Some((ap, _)) = mmff94_angle_energy_resolved(
                            at,
                            types[i],
                            types[j_idx],
                            types[k],
                            rij.r0,
                            rkj.r0,
                            ring_size,
                        ) {
                            let dtheta = cos_t.dacos() * RAD_TO_DEG - ap.theta0;
                            energy += CONV * (kba_ijk * dr_ij + kba_kji * dr_kj) * dtheta;
                        }
                    }
                    let _ = (KB_CONV, CS); // suppress warnings
                }
            }
        }
    }
    energy
}

/// sp2 atom types eligible for out-of-plane bending (Halgren 1996), used by
/// [`oop_energy`] below. `pub` so downstream coverage-checkers (e.g.
/// `chematic-3d`'s MMFF94 bridge) can mirror exactly which atoms this
/// module's own energy loop would evaluate, instead of hand-copying this
/// list and risking drift. Every centre type with an MMFFOOP row is listed,
/// as RDKit adds a term at any three-connected atom its table covers: NCN+
/// (55) and the N-oxide nitrogen NPOX (69) carry barriers; NR (8), S=O (17),
/// tricoordinate P (26) and SO2M (73) have zero rows.
pub const OOP_SP2_TYPES: &[u8] = &[
    2, 3, 8, 9, 10, 17, 26, 30, 37, 38, 39, 40, 41, 43, 45, 49, 54, 55, 56, 57, 58, 59, 63, 64, 65,
    66, 67, 69, 73, 76, 78, 79, 80, 81, 82,
];

/// Out-of-plane bending for trigonal sp2 centers (Halgren MMFF.VI eq. 6)
/// E_oop = (0.043844 × koop / 2) × χ²  (χ in degrees: Wilson angle of out-of-plane distortion)
fn oop_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    const CONV: f64 = MMFF_ANGLE_CONV;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        if OOP_SP2_TYPES.binary_search(&types[j_idx]).is_err() {
            continue;
        }
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() != 3 {
            continue; // OOP only for exactly 3 substituents (trigonal)
        }
        let Some(koop) = mmff94_oop(
            types[j_idx],
            types[neighbors[0]],
            types[neighbors[1]],
            types[neighbors[2]],
        ) else {
            continue;
        };
        for [i, k, l] in oop_permutations(&neighbors) {
            // Wilson out-of-plane angle: angle between j→l vector and plane (i,j,k)
            let pj = coords[j_idx];
            let pi = coords[i];
            let pk = coords[k];
            let pl = coords[l];
            let rji = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
            let rjk = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
            let rjl = [pl[0] - pj[0], pl[1] - pj[1], pl[2] - pj[2]];
            let n = cross(rji, rjk); // normal to ijk plane
            let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let l_len = (rjl[0] * rjl[0] + rjl[1] * rjl[1] + rjl[2] * rjl[2]).sqrt();
            if n_len < 1e-12 || l_len < 1e-12 {
                continue;
            }
            let sin_chi = dot3(n, rjl) / (n_len * l_len);
            let chi_deg = sin_chi.clamp(-1.0, 1.0).dasin() * RAD_TO_DEG;
            energy += (CONV * koop / 2.0) * chi_deg * chi_deg;
        }
    }
    energy
}

/// Bond stretching: cubic-corrected harmonic (Halgren MMFF.II eq. 1)
/// E = (143.9325 × kb / 2) × ΔR² × (1 − cs×ΔR + (7/12)×cs²×ΔR²)
fn bond_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    const KB_CONV: f64 = MDYNE_A_TO_KCAL_MOL;
    const CS: f64 = 2.0;
    let mut energy = 0.0;
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let bt = bond_type_for(types[i], types[j], bond.order);
        if let Some((p, _)) = mmff94_bond_energy_resolved(bt, types[i], types[j]) {
            let r = dist(coords[i], coords[j]);
            let dr = r - p.r0;
            let cubic = 1.0 - CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
            energy += (KB_CONV * p.kb / 2.0) * dr * dr * cubic;
        }
    }
    energy
}

/// Angle bending: cubic-corrected harmonic (Halgren MMFF.III eq. 2)
/// E = (0.043844 × ka / 2) × Δθ² × (1 − 0.007×Δθ)   [Δθ in degrees], with
/// RDKit's unrounded constants; `143.9325 ka (1 + cos θ)` at a linear centre
fn angle_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8], rings: &[Vec<AtomIdx>]) -> f64 {
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let r0_ij = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx]);
                let r0_kj = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx]);
                let (Some((bond_ij, _)), Some((bond_kj, _))) = (r0_ij, r0_kj) else {
                    continue;
                };
                let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                if let Some((p, _)) = mmff94_angle_energy_resolved(
                    at,
                    types[i],
                    types[j_idx],
                    types[k],
                    bond_ij.r0,
                    bond_kj.r0,
                    ring_size,
                ) {
                    let cos_t = cos_angle(coords[i], coords[j_idx], coords[k]);
                    energy += angle_bend_energy(
                        cos_t,
                        p.theta0,
                        p.ka,
                        central_type_is_linear(types[j_idx]),
                    );
                }
            }
        }
    }
    energy
}

/// Torsion: three-term Fourier (Halgren MMFF.IV)
/// E = (v1/2)(1+cosφ) + (v2/2)(1-cos2φ) + (v3/2)(1+cos3φ)
fn torsion_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    let mut energy = 0.0;
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
                if let Some(p) = mmff94_torsion_term_params(mol, i, j, k, l, types) {
                    let phi = dihedral(coords[i], coords[j], coords[k], coords[l]);
                    energy += 0.5 * p.v1 * (1.0 + phi.dcos())
                        + 0.5 * p.v2 * (1.0 - (2.0 * phi).dcos())
                        + 0.5 * p.v3 * (1.0 + (3.0 * phi).dcos());
                }
            }
        }
    }
    energy
}

/// Van der Waals: buffered 14-7 (Halgren MMFF.I eq. 8).
/// The potential implementation is shared with the prepared-pair path below.
fn vdw_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    let n = mol.atom_count();
    let mut excl = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excl.insert((i.min(j), i.max(j)));
        for (nb_i, _) in mol.neighbors(bond.atom1) {
            let ni = nb_i.0 as usize;
            excl.insert((ni.min(j), ni.max(j)));
        }
        for (nb_j, _) in mol.neighbors(bond.atom2) {
            let nj = nb_j.0 as usize;
            excl.insert((i.min(nj), i.max(nj)));
        }
    }
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excl.contains(&(i, j)) {
                continue;
            }
            let r = dist(coords[i], coords[j]);
            if let Some((r_star, eps)) = mmff94_vdw_combined(types[i], types[j])
                && r_star > 0.0
                && eps > 0.0
                && r <= NONBONDED_THRESH
            {
                energy += mmff94_vdw_energy_value(r, r_star, eps);
            }
        }
    }
    energy
}

/// Electrostatic: Coulomb with δ=0.05 Å buffer (Halgren MMFF.V eq. 14)
/// E = 332.0716 × q_i × q_j / (D × (r + δ))   [D=1.0]
fn elec_energy(mol: &Molecule, coords: &[[f64; 3]], charges: &[f64]) -> f64 {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    let n = mol.atom_count();
    let mut excl = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excl.insert((i.min(j), i.max(j)));
        for (nb_i, _) in mol.neighbors(bond.atom1) {
            excl.insert(((nb_i.0 as usize).min(j), (nb_i.0 as usize).max(j)));
        }
        for (nb_j, _) in mol.neighbors(bond.atom2) {
            excl.insert((i.min(nb_j.0 as usize), i.max(nb_j.0 as usize)));
        }
    }
    // 1-4 pairs: scale by 0.75 (MMFF94 convention)
    let mut one_four = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        for (nb_j, _) in mol.neighbors(bond.atom1) {
            let i = nb_j.0 as usize;
            if i == k {
                continue;
            }
            for (nb_k, _) in mol.neighbors(bond.atom2) {
                let l = nb_k.0 as usize;
                if l == j {
                    continue;
                }
                let key = (i.min(l), i.max(l));
                if !excl.contains(&key) {
                    one_four.insert(key);
                }
            }
        }
    }
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excl.contains(&(i, j)) {
                continue;
            }
            let r = dist(coords[i], coords[j]);
            if r > NONBONDED_THRESH {
                continue;
            }
            let scale = if one_four.contains(&(i, j)) {
                0.75
            } else {
                1.0
            };
            energy += scale * COULOMB * charges[i] * charges[j] / (r + DELTA);
        }
    }
    energy
}

// ─── MMFF94 buffered 14-7 helpers ───────────────────────────────────────────

/// Evaluate Halgren's buffered 14-7 van der Waals potential.
///
/// The second factor is buffered independently from the first one.  It is
/// tempting to write this as `t^7 * (t^7 - 2)`, but that is a different
/// potential and causes large errors for close aromatic contacts.
#[inline]
fn mmff94_vdw_energy_value(r: f64, r_star: f64, epsilon: f64) -> f64 {
    let first = (1.07 * r_star / (r + 0.07 * r_star)).powi(7);
    let r_star_7 = r_star.powi(7);
    let second = 1.12 * r_star_7 / (r.powi(7) + 0.12 * r_star_7) - 2.0;
    epsilon * first * second
}

/// Radial derivative of [`mmff94_vdw_energy_value`], dE/dr.
#[inline]
fn mmff94_vdw_radial_derivative(r: f64, r_star: f64, epsilon: f64) -> f64 {
    let buffered_denominator = r + 0.07 * r_star;
    let first = (1.07 * r_star / buffered_denominator).powi(7);
    let r_star_7 = r_star.powi(7);
    let second_denominator = r.powi(7) + 0.12 * r_star_7;
    let second = 1.12 * r_star_7 / second_denominator - 2.0;
    let d_first = -7.0 * first / buffered_denominator;
    let d_second = -7.0 * 1.12 * r_star_7 * r.powi(6) / second_denominator.powi(2);
    epsilon * (d_first * second + first * d_second)
}

// ─── Geometry helpers ─────────────────────────────────────────────────────────

#[inline]
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[inline]
fn cos_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot_val = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
    let na = (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]).sqrt();
    let nc = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]).sqrt();
    if na < 1e-12 || nc < 1e-12 {
        return 0.0;
    }
    (dot_val / (na * nc)).clamp(-1.0, 1.0)
}

/// Signed dihedral angle φ (radians) for the quartet i-j-k-l.
#[inline]
fn dihedral(i: [f64; 3], j: [f64; 3], k: [f64; 3], l: [f64; 3]) -> f64 {
    let b1 = [j[0] - i[0], j[1] - i[1], j[2] - i[2]];
    let b2 = [k[0] - j[0], k[1] - j[1], k[2] - j[2]];
    let b3 = [l[0] - k[0], l[1] - k[1], l[2] - k[2]];
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m1 = cross(n1, b2);
    let b2_len = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    if b2_len < 1e-12 {
        return 0.0;
    }
    let x = dot3(n1, n2);
    let y = dot3(m1, n2) / b2_len;
    y.datan2(x)
}

#[derive(Clone, Copy)]
struct Dual12 {
    value: f64,
    grad: [f64; 12],
}

fn dual_point(point: [f64; 3], offset: usize) -> [Dual12; 3] {
    std::array::from_fn(|axis| {
        let mut grad = [0.0; 12];
        grad[offset + axis] = 1.0;
        Dual12 {
            value: point[axis],
            grad,
        }
    })
}

fn dual_sub(a: [Dual12; 3], b: [Dual12; 3]) -> [Dual12; 3] {
    std::array::from_fn(|axis| Dual12 {
        value: a[axis].value - b[axis].value,
        grad: std::array::from_fn(|i| a[axis].grad[i] - b[axis].grad[i]),
    })
}

fn dual_cross(a: [Dual12; 3], b: [Dual12; 3]) -> [Dual12; 3] {
    [
        dual_sub_scalar(dual_mul(a[1], b[2]), dual_mul(a[2], b[1])),
        dual_sub_scalar(dual_mul(a[2], b[0]), dual_mul(a[0], b[2])),
        dual_sub_scalar(dual_mul(a[0], b[1]), dual_mul(a[1], b[0])),
    ]
}

fn dual_dot(a: [Dual12; 3], b: [Dual12; 3]) -> Dual12 {
    dual_add(
        dual_add(dual_mul(a[0], b[0]), dual_mul(a[1], b[1])),
        dual_mul(a[2], b[2]),
    )
}

fn dual_add(a: Dual12, b: Dual12) -> Dual12 {
    Dual12 {
        value: a.value + b.value,
        grad: std::array::from_fn(|i| a.grad[i] + b.grad[i]),
    }
}

fn dual_sub_scalar(a: Dual12, b: Dual12) -> Dual12 {
    Dual12 {
        value: a.value - b.value,
        grad: std::array::from_fn(|i| a.grad[i] - b.grad[i]),
    }
}

fn dual_mul(a: Dual12, b: Dual12) -> Dual12 {
    Dual12 {
        value: a.value * b.value,
        grad: std::array::from_fn(|i| a.grad[i] * b.value + a.value * b.grad[i]),
    }
}

fn dual_div(a: Dual12, b: Dual12) -> Dual12 {
    let inv = 1.0 / b.value;
    Dual12 {
        value: a.value * inv,
        grad: std::array::from_fn(|i| (a.grad[i] * b.value - a.value * b.grad[i]) * inv * inv),
    }
}

fn dual_sqrt(a: Dual12) -> Dual12 {
    let value = a.value.sqrt();
    Dual12 {
        value,
        grad: std::array::from_fn(|i| a.grad[i] / (2.0 * value)),
    }
}

fn dual_atan2(y: Dual12, x: Dual12) -> Dual12 {
    let denominator = x.value * x.value + y.value * y.value;
    Dual12 {
        value: y.value.datan2(x.value),
        grad: std::array::from_fn(|i| (x.value * y.grad[i] - y.value * x.grad[i]) / denominator),
    }
}

fn dual_asin_clamped(a: Dual12) -> Dual12 {
    let value = a.value.clamp(-1.0, 1.0).dasin();
    let denominator = (1.0 - a.value * a.value).max(1e-24).sqrt();
    Dual12 {
        value,
        grad: std::array::from_fn(|i| a.grad[i] / denominator),
    }
}

fn dual_dihedral(i: [f64; 3], j: [f64; 3], k: [f64; 3], l: [f64; 3]) -> Dual12 {
    let b1 = dual_sub(dual_point(i, 0), dual_point(j, 3));
    let b2 = dual_sub(dual_point(k, 6), dual_point(j, 3));
    let b3 = dual_sub(dual_point(l, 9), dual_point(k, 6));
    let n1 = dual_cross(b1, b2);
    let n2 = dual_cross(b2, b3);
    let x = dual_dot(n1, n2);
    let y = dual_div(
        dual_dot(dual_cross(n1, b2), n2),
        dual_sqrt(dual_dot(b2, b2)),
    );
    dual_atan2(y, x)
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

