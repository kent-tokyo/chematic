use super::*;

#[derive(Clone, Copy)]
enum Path14BoundsKind {
    Cis { extra_distance: f64 },
    Trans { extra_distance: f64 },
    Flexible,
    Custom { torsion: f64 },
}

fn path14_bounds(
    mol: &RdkitMolView,
    m: &BoundsMatrix,
    atoms: [usize; 4],
    angles: [f64; 2],
    kind: Path14BoundsKind,
) -> (f64, f64) {
    let [aid1, atm2, atm3, aid4] = atoms;
    let lower_bonds = [
        m.lower(aid1, atm2),
        m.lower(atm2, atm3),
        m.lower(atm3, aid4),
    ];
    let upper_bonds = [
        m.upper(aid1, atm2),
        m.upper(atm2, atm3),
        m.upper(atm3, aid4),
    ];
    let (angle_lower_12, angle_upper_12) =
        angle_range(angles[0], scaled_angle_tolerance(mol, [aid1, atm2, atm3]));
    let (angle_lower_23, angle_upper_23) =
        angle_range(angles[1], scaled_angle_tolerance(mol, [atm2, atm3, aid4]));

    let (mut lower, mut upper) = match kind {
        Path14BoundsKind::Cis { extra_distance } => (
            compute14_dist_cis(
                lower_bonds[0],
                lower_bonds[1],
                lower_bonds[2],
                angle_lower_12,
                angle_lower_23,
            ) + extra_distance,
            compute14_dist_cis(
                upper_bonds[0],
                upper_bonds[1],
                upper_bonds[2],
                angle_upper_12,
                angle_upper_23,
            ) + extra_distance,
        ),
        Path14BoundsKind::Trans { extra_distance } => (
            compute14_dist_trans(
                lower_bonds[0],
                lower_bonds[1],
                lower_bonds[2],
                angle_lower_12,
                angle_lower_23,
            ) + extra_distance,
            compute14_dist_trans(
                upper_bonds[0],
                upper_bonds[1],
                upper_bonds[2],
                angle_upper_12,
                angle_upper_23,
            ) + extra_distance,
        ),
        Path14BoundsKind::Flexible => (
            compute14_dist_cis(
                lower_bonds[0],
                lower_bonds[1],
                lower_bonds[2],
                angle_lower_12,
                angle_lower_23,
            ),
            compute14_dist_trans(
                upper_bonds[0],
                upper_bonds[1],
                upper_bonds[2],
                angle_upper_12,
                angle_upper_23,
            ),
        ),
        Path14BoundsKind::Custom { torsion } => (
            compute14_dist_3d(
                lower_bonds[0],
                lower_bonds[1],
                lower_bonds[2],
                angle_lower_12,
                angle_lower_23,
                torsion,
            ),
            compute14_dist_3d(
                upper_bonds[0],
                upper_bonds[1],
                upper_bonds[2],
                angle_upper_12,
                angle_upper_23,
                torsion,
            ),
        ),
    };
    if matches!(kind, Path14BoundsKind::Flexible) && upper < lower {
        std::mem::swap(&mut upper, &mut lower);
    }
    if (upper - lower).abs() < DIST12_DELTA {
        lower -= GEN_DIST_TOL;
        upper += GEN_DIST_TOL;
    }
    (lower, upper)
}

fn set_in_ring14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    ring_size: usize,
    collected: &mut HashMap<usize, Vec<Collected14Bounds>>,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let ahyb2 = mol.atoms[atm2].hybridization;
    let ahyb3 = mol.atoms[atm3].hybridization;
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    if cx.d(aid1, aid4) < 2.9 {
        return Ok(());
    }
    let angles = path_angles(acc, b1, b2, b3)?;
    let stype = atom_stereo(mol, b2, aid1, aid4);
    // RDKit 2026.09 no longer special-cases fused-ring membership here. A
    // path wholly in a small ring, or a path classified as sharing/different
    // ring bonds (`ring_size == 0`), is cis when both central atoms are sp2
    // unless explicit E stereo says otherwise.
    let prefer_cis = (ring_size <= 8 && ahyb2 == HYB_SP2 && ahyb3 == HYB_SP2 && stype != STEREO_E)
        || stype == STEREO_Z;
    let prefer_trans = stype == STEREO_E;
    let (ty, kind) = if prefer_cis {
        ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
        (
            Path14Type::Cis,
            Path14BoundsKind::Cis {
                extra_distance: 0.0,
            },
        )
    } else if prefer_trans {
        ComputedData::add_path(&mut acc.trans_paths, cx.nb, b1, b2, b3);
        (
            Path14Type::Trans,
            Path14BoundsKind::Trans {
                extra_distance: 0.0,
            },
        )
    } else {
        (Path14Type::Other, Path14BoundsKind::Flexible)
    };
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
        forced: false,
    });
    let (dl, du) = path14_bounds(mol, m, [aid1, atm2, atm3, aid4], [angles.0, angles.1], kind);
    acc.visited14[p] = true;
    collected.entry(p).or_default().push(Collected14Bounds {
        lower: dl,
        upper: du,
        aid1,
        aid4,
    });
    Ok(())
}

fn set_two_in_same_ring14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macrocycle: bool,
    prefer_trans: bool,
    collected: &mut HashMap<usize, Vec<Collected14Bounds>>,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    if cx.d(aid1, aid4) < 2.9 {
        return Ok(());
    }
    if mol.bond_between(aid1, atm3).is_some() || mol.bond_between(aid4, atm2).is_some() {
        return Ok(());
    }
    let angles = path_angles(acc, b1, b2, b3)?;
    let amide = |bb1: usize, bb3: usize, a1: usize, a2: usize, a3: usize, a4: usize| {
        let n = |a: usize| mol.atoms[a].atomic_num;
        n(a1) != 1
            && n(a3) == 6
            && mol.bonds[bb3].bond_type == DOUBLE
            && (n(a4) == 8 || n(a4) == 7)
            && mol.bonds[bb1].bond_type == SINGLE
            && (n(a2) == 8 || n(a2) == 7)
    };
    let stype = atom_stereo(mol, b2, aid1, aid4);
    let (ty, kind) = if macrocycle {
        if amide(b1, b3, aid1, atm2, atm3, aid4) || amide(b3, b1, aid4, atm3, atm2, aid1) {
            // RDKit's MACROCYCLE_TWO_IN_SAME_RING call intentionally leaves
            // `useMacrocycle14Config` false. Amide/ester paths therefore span
            // cis through trans instead of being pinned cis.
            (Path14Type::Other, Path14BoundsKind::Flexible)
        } else if stype == STEREO_Z {
            ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
            (
                Path14Type::Cis,
                Path14BoundsKind::Cis {
                    extra_distance: 0.0,
                },
            )
        } else if stype == STEREO_E {
            ComputedData::add_path(&mut acc.trans_paths, cx.nb, b1, b2, b3);
            (
                Path14Type::Trans,
                Path14BoundsKind::Trans {
                    extra_distance: 0.0,
                },
            )
        } else {
            (Path14Type::Other, Path14BoundsKind::Flexible)
        }
    } else if (prefer_trans
        && mol.atoms[atm2].hybridization == HYB_SP2
        && mol.atoms[atm3].hybridization == HYB_SP2
        && stype != STEREO_Z)
        || stype == STEREO_E
    {
        ComputedData::add_path(&mut acc.trans_paths, cx.nb, b1, b2, b3);
        (
            Path14Type::Trans,
            Path14BoundsKind::Trans {
                extra_distance: 0.0,
            },
        )
    } else if stype == STEREO_Z {
        ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
        (
            Path14Type::Cis,
            Path14BoundsKind::Cis {
                extra_distance: 0.0,
            },
        )
    } else {
        (Path14Type::Other, Path14BoundsKind::Flexible)
    };
    let (dl, du) = path14_bounds(mol, m, [aid1, atm2, atm3, aid4], [angles.0, angles.1], kind);
    collected.entry(p).or_default().push(Collected14Bounds {
        lower: dl,
        upper: du,
        aid1,
        aid4,
    });
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
        forced: false,
    });
    acc.visited14[p] = true;
    Ok(())
}

fn check_h2_nx3h1_ox2(mol: &RdkitMolView, a: usize) -> bool {
    let at = &mol.atoms[a];
    (at.atomic_num == 6 && at.total_num_hs_with_neighbors == 2)
        || (at.atomic_num == 8 && at.total_num_hs_with_neighbors == 0)
        || (at.atomic_num == 7 && mol.degree(a) == 3 && at.total_num_hs_with_neighbors == 1)
}

#[allow(dead_code)]
fn check_nh_ch_ch_nh(mol: &RdkitMolView, a1: usize, a2: usize, a3: usize, a4: usize) -> bool {
    mol.atoms[a1].atomic_num != 1
        && mol.atoms[a4].atomic_num != 1
        && check_h2_nx3h1_ox2(mol, a2)
        && check_h2_nx3h1_ox2(mol, a3)
}

fn check_amide_ester14(
    mol: &RdkitMolView,
    bnd1: usize,
    bnd3: usize,
    a2: usize,
    a3: usize,
    a4: usize,
) -> bool {
    let n = |a: usize| mol.atoms[a].atomic_num;
    n(a3) == 6
        && mol.bonds[bnd3].bond_type == DOUBLE
        && (n(a4) == 8 || n(a4) == 7)
        && mol.bonds[bnd1].bond_type == SINGLE
        && (n(a2) == 8 || (n(a2) == 7 && mol.atoms[a2].total_num_hs_with_neighbors == 1))
}

fn check_macrocycle_all_in_same_ring_amide_ester14(
    mol: &RdkitMolView,
    a1: usize,
    a2: usize,
    a3: usize,
    a4: usize,
) -> bool {
    let n = |a: usize| mol.atoms[a].atomic_num;
    if n(a3) != 6 {
        return false;
    }
    if (n(a2) == 7 || n(a2) == 8) && mol.degree(a2) == 3 && mol.degree(a3) == 3 {
        for nbr in mol.neighbors(a2) {
            if nbr != a1 && nbr != a3 {
                let b = mol.bond_between(a2, nbr).expect("bond");
                if (n(nbr) != 6 && n(nbr) != 1) || mol.bonds[b].bond_type != SINGLE {
                    return false;
                }
                break;
            }
        }
        for nbr in mol.neighbors(a3) {
            if nbr != a2 && nbr != a4 {
                let b = mol.bond_between(a3, nbr).expect("bond");
                if n(nbr) != 8 || mol.bonds[b].bond_type != DOUBLE {
                    return false;
                }
                break;
            }
        }
        return true;
    }
    false
}

fn is_carbonyl(mol: &RdkitMolView, a: usize) -> bool {
    if mol.atoms[a].atomic_num == 6 && mol.degree(a) > 2 {
        for nbr in mol.neighbors(a) {
            let n = mol.atoms[nbr].atomic_num;
            if (n == 8 || n == 7)
                && mol.bonds[mol.bond_between(a, nbr).expect("bond")].bond_type == DOUBLE
            {
                return true;
            }
        }
    }
    false
}

fn check_amide_ester15(mol: &RdkitMolView, bnd1: usize, bnd3: usize, a2: usize, a3: usize) -> bool {
    let n2 = mol.atoms[a2].atomic_num;
    (n2 == 8 || (n2 == 7 && mol.atoms[a2].total_num_hs_with_neighbors == 1))
        && mol.bonds[bnd1].bond_type == SINGLE
        && mol.atoms[a3].atomic_num == 6
        && mol.bonds[bnd3].bond_type == SINGLE
        && is_carbonyl(mol, a3)
}

struct Chain14Decision {
    kind: Option<Path14BoundsKind>,
    ty: Path14Type,
    forced: bool,
}

fn add_typed_path(acc: &mut ComputedData, nb: u64, bonds: [usize; 3], ty: Path14Type) {
    match ty {
        Path14Type::Cis => {
            ComputedData::add_path(&mut acc.cis_paths, nb, bonds[0], bonds[1], bonds[2])
        }
        Path14Type::Trans => {
            ComputedData::add_path(&mut acc.trans_paths, nb, bonds[0], bonds[1], bonds[2])
        }
        Path14Type::Other => {}
    }
}

fn pinned_chain14(
    kind: Path14BoundsKind,
    ty: Path14Type,
    acc: &mut ComputedData,
    nb: u64,
    bonds: [usize; 3],
) -> Chain14Decision {
    add_typed_path(acc, nb, bonds, ty);
    Chain14Decision {
        kind: Some(kind),
        ty,
        forced: false,
    }
}

fn forced_chain14(
    kind: Path14BoundsKind,
    ty: Path14Type,
    acc: &mut ComputedData,
    nb: u64,
    bonds: [usize; 3],
) -> Chain14Decision {
    let mut decision = pinned_chain14(kind, ty, acc, nb, bonds);
    decision.forced = true;
    decision
}

fn flexible_chain14() -> Chain14Decision {
    Chain14Decision {
        kind: Some(Path14BoundsKind::Flexible),
        ty: Path14Type::Other,
        forced: false,
    }
}

#[allow(clippy::too_many_arguments)]
fn double_bond_chain14(
    cx: &Ctx,
    atoms: [usize; 4],
    bonds: [usize; 3],
    acc: &mut ComputedData,
) -> Chain14Decision {
    if cx.mol.bonds[bonds[0]].bond_type == DOUBLE || cx.mol.bonds[bonds[2]].bond_type == DOUBLE {
        return pinned_chain14(
            Path14BoundsKind::Cis {
                extra_distance: 0.0,
            },
            Path14Type::Cis,
            acc,
            cx.nb,
            bonds,
        );
    }
    if cx.mol.bonds[bonds[1]].stereo <= STEREO_ANY {
        return flexible_chain14();
    }
    let ty = if atom_stereo(cx.mol, bonds[1], atoms[0], atoms[3]) == STEREO_Z {
        Path14Type::Cis
    } else {
        Path14Type::Trans
    };
    let kind = match ty {
        Path14Type::Cis => Path14BoundsKind::Cis {
            extra_distance: 0.0,
        },
        Path14Type::Trans => Path14BoundsKind::Trans {
            extra_distance: 0.0,
        },
        Path14Type::Other => unreachable!(),
    };
    pinned_chain14(kind, ty, acc, cx.nb, bonds)
}

fn secondary_amide_h(mol: &RdkitMolView, hydrogen: usize, nitrogen: usize) -> bool {
    mol.atoms[hydrogen].atomic_num == 1
        && mol.atoms[nitrogen].atomic_num == 7
        && mol.degree(nitrogen) == 3
        && mol.atoms[nitrogen].total_num_hs_with_neighbors == 1
}

fn pin_chain14_type(
    ty: Path14Type,
    acc: &mut ComputedData,
    nb: u64,
    bonds: [usize; 3],
) -> Chain14Decision {
    let kind = match ty {
        Path14Type::Cis => Path14BoundsKind::Cis {
            extra_distance: 0.0,
        },
        Path14Type::Trans => Path14BoundsKind::Trans {
            extra_distance: 0.0,
        },
        Path14Type::Other => Path14BoundsKind::Custom {
            torsion: M_PI / 2.0,
        },
    };
    forced_chain14(kind, ty, acc, nb, bonds)
}

fn amide14_chain_decision(
    cx: &Ctx,
    atoms: [usize; 4],
    bonds: [usize; 3],
    acc: &mut ComputedData,
    force_trans_amides: bool,
) -> Option<Chain14Decision> {
    let [aid1, atm2, atm3, aid4] = atoms;
    if !(check_amide_ester14(cx.mol, bonds[0], bonds[2], atm2, atm3, aid4)
        || check_amide_ester14(cx.mol, bonds[2], bonds[0], atm3, atm2, aid1))
    {
        return None;
    }
    if !force_trans_amides {
        return Some(flexible_chain14());
    }
    let ty = if secondary_amide_h(cx.mol, aid1, atm2) || secondary_amide_h(cx.mol, aid4, atm3) {
        Path14Type::Trans
    } else {
        Path14Type::Cis
    };
    Some(pin_chain14_type(ty, acc, cx.nb, bonds))
}

#[allow(clippy::too_many_arguments)]
fn amide15_chain_decision(
    cx: &Ctx,
    atoms: [usize; 4],
    bonds: [usize; 3],
    acc: &mut ComputedData,
    force_trans_amides: bool,
    macrocycle: bool,
) -> Option<Chain14Decision> {
    let [aid1, atm2, atm3, aid4] = atoms;
    if !(check_amide_ester15(cx.mol, bonds[0], bonds[2], atm2, atm3)
        || check_amide_ester15(cx.mol, bonds[2], bonds[0], atm3, atm2))
    {
        return None;
    }
    if macrocycle
        && cx.mol.atoms[atm2].atomic_num == 7
        && cx.mol.degree(atm2) == 3
        && cx.mol.atoms[aid1].atomic_num == 1
        && cx.mol.atoms[atm2].total_num_hs_with_neighbors == 1
    {
        return Some(forced_chain14(
            Path14BoundsKind::Cis {
                extra_distance: 0.0,
            },
            Path14Type::Cis,
            acc,
            cx.nb,
            bonds,
        ));
    }
    if !macrocycle && !force_trans_amides {
        return Some(flexible_chain14());
    }
    let ty = if !macrocycle
        && (secondary_amide_h(cx.mol, aid1, atm2) || secondary_amide_h(cx.mol, aid4, atm3))
    {
        Path14Type::Cis
    } else {
        Path14Type::Trans
    };
    Some(pin_chain14_type(ty, acc, cx.nb, bonds))
}

#[allow(clippy::too_many_arguments)]
fn single_bond_chain14(
    cx: &Ctx,
    atoms: [usize; 4],
    bonds: [usize; 3],
    acc: &mut ComputedData,
    force_trans_amides: bool,
    macrocycle: bool,
) -> Chain14Decision {
    let [aid1, atm2, atm3, aid4] = atoms;
    let atomic_num = |atom: usize| cx.mol.atoms[atom].atomic_num;
    if atomic_num(atm2) == 16
        && atomic_num(atm3) == 16
        && cx.mol.degree(atm2) == 2
        && cx.mol.degree(atm3) == 2
    {
        return pinned_chain14(
            Path14BoundsKind::Custom {
                torsion: M_PI / 2.0,
            },
            Path14Type::Other,
            acc,
            cx.nb,
            bonds,
        );
    }
    if macrocycle
        && (check_macrocycle_all_in_same_ring_amide_ester14(cx.mol, aid1, atm2, atm3, aid4)
            || check_macrocycle_all_in_same_ring_amide_ester14(cx.mol, aid4, atm3, atm2, aid1))
    {
        return forced_chain14(
            Path14BoundsKind::Trans {
                extra_distance: 0.1,
            },
            Path14Type::Trans,
            acc,
            cx.nb,
            bonds,
        );
    }
    if !macrocycle
        && let Some(decision) = amide14_chain_decision(cx, atoms, bonds, acc, force_trans_amides)
    {
        return decision;
    }
    if let Some(decision) =
        amide15_chain_decision(cx, atoms, bonds, acc, force_trans_amides, macrocycle)
    {
        return decision;
    }
    flexible_chain14()
}

#[allow(unused_assignments)]
/// `_setChain14Bounds` (`macrocycle == false`) and
/// `_setMacrocycleAllInSameRing14Bounds` (`macrocycle == true`).
fn set_chain14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    force_trans_amides: bool,
    macrocycle: bool,
    collected: &mut HashMap<usize, Vec<Collected14Bounds>>,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    let angles = path_angles(acc, b1, b2, b3)?;
    let bonds = [b1, b2, b3];
    let decision = match mol.bonds[b2].bond_type {
        DOUBLE => double_bond_chain14(cx, [aid1, atm2, atm3, aid4], bonds, acc),
        SINGLE => single_bond_chain14(
            cx,
            [aid1, atm2, atm3, aid4],
            bonds,
            acc,
            force_trans_amides,
            macrocycle,
        ),
        _ => flexible_chain14(),
    };
    if let Some(kind) = decision.kind {
        let (dl, du) = path14_bounds(mol, m, [aid1, atm2, atm3, aid4], [angles.0, angles.1], kind);
        collected.entry(p).or_default().push(Collected14Bounds {
            lower: dl,
            upper: du,
            aid1,
            aid4,
        });
        acc.paths14.push(Path14 {
            bid1: b1,
            bid2: b2,
            bid3: b3,
            ty: decision.ty,
            forced: decision.forced,
        });
        acc.visited14[p] = true;
    }
    Ok(())
}

/// `_record14Path`.
fn record14_path(cx: &Ctx, b1: usize, b2: usize, b3: usize, acc: &mut ComputedData) {
    let atm2 = acc.bond_adj.get(b1, b2) as usize;
    let atm3 = acc.bond_adj.get(b2, b3) as usize;
    let ty = if cx.mol.atoms[atm2].hybridization == HYB_SP2
        && cx.mol.atoms[atm3].hybridization == HYB_SP2
    {
        ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
        Path14Type::Cis
    } else {
        Path14Type::Other
    };
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
        forced: false,
    });
}

#[allow(clippy::too_many_arguments)]
fn set_ring_14_bounds(
    cx: &Ctx,
    ring: &[usize],
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macro14: bool,
    force_trans_amides: bool,
    macrocycle_bonds: &mut HashSet<usize>,
    ring_bond_pairs: &mut HashSet<u64>,
    cis_ring_bond_pairs: &mut HashSet<u64>,
    done_paths: &mut HashSet<u64>,
    collected: &mut HashMap<usize, Vec<Collected14Bounds>>,
) -> Result<(), BoundsInvariant> {
    if ring.len() < 3 {
        return Ok(());
    }
    let mut bid1 = ring[ring.len() - 1];
    for position in 0..ring.len() {
        let bid2 = ring[position];
        let bid3 = ring[(position + 1) % ring.len()];
        let (b1, b2, b3) = (bid1 as u64, bid2 as u64, bid3 as u64);
        ring_bond_pairs.extend([b1 * cx.nb + b2, b2 * cx.nb + b1]);
        done_paths.extend([
            b1 * cx.nb * cx.nb + b2 * cx.nb + b3,
            b3 * cx.nb * cx.nb + b2 * cx.nb + b1,
        ]);
        if ring.len() > 5 {
            if macro14 && ring.len() >= MIN_MACROCYCLE_RING_SIZE {
                set_chain14(
                    cx,
                    bid1,
                    bid2,
                    bid3,
                    acc,
                    m,
                    force_trans_amides,
                    true,
                    collected,
                )?;
                macrocycle_bonds.insert(bid2);
            } else {
                set_in_ring14(cx, bid1, bid2, bid3, acc, m, ring.len(), collected)?;
                if ring.len() <= 8 {
                    cis_ring_bond_pairs.extend([b1 * cx.nb + b2, b2 * cx.nb + b1]);
                }
            }
        } else {
            record14_path(cx, bid1, bid2, bid3, acc);
            cis_ring_bond_pairs.extend([b1 * cx.nb + b2, b2 * cx.nb + b1]);
        }
        bid1 = bid2;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn set_remaining_14_path(
    cx: &Ctx,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macro14: bool,
    force_trans_amides: bool,
    macrocycle_bonds: &HashSet<usize>,
    ring_bond_pairs: &HashSet<u64>,
    cis_ring_bond_pairs: &HashSet<u64>,
    collected: &mut HashMap<usize, Vec<Collected14Bounds>>,
) -> Result<(), BoundsInvariant> {
    let (b1, b2, b3) = (bid1 as u64, bid2 as u64, bid3 as u64);
    let pairs = [
        b1 * cx.nb + b2,
        b2 * cx.nb + b1,
        b2 * cx.nb + b3,
        b3 * cx.nb + b2,
    ];
    if pairs.iter().any(|pair| ring_bond_pairs.contains(pair)) {
        set_two_in_same_ring14(
            cx,
            bid1,
            bid2,
            bid3,
            acc,
            m,
            macro14 && macrocycle_bonds.contains(&bid2),
            pairs.iter().any(|pair| cis_ring_bond_pairs.contains(pair)),
            collected,
        )
    } else if (cx.mol.num_bond_rings(bid1) > 0 && cx.mol.num_bond_rings(bid2) > 0)
        || (cx.mol.num_bond_rings(bid2) > 0 && cx.mol.num_bond_rings(bid3) > 0)
        || cx.mol.num_bond_rings(bid2) > 0
    {
        set_in_ring14(cx, bid1, bid2, bid3, acc, m, 0, collected)
    } else {
        set_chain14(
            cx,
            bid1,
            bid2,
            bid3,
            acc,
            m,
            force_trans_amides,
            false,
            collected,
        )
    }
}

pub(super) fn set14_bounds(
    cx: &Ctx,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macro14: bool,
    force_trans_amides: bool,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let nb = cx.nb;
    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut cis_ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    let mut collected: HashMap<usize, Vec<Collected14Bounds>> = HashMap::new();
    let mut bond_rings = mol.bond_rings.clone();
    bond_rings.sort_by_key(|ring| std::cmp::Reverse(ring.len()));
    for bring in &bond_rings {
        set_ring_14_bounds(
            cx,
            bring,
            acc,
            m,
            macro14,
            force_trans_amides,
            &mut bid_is_macrocycle,
            &mut ring_bond_pairs,
            &mut cis_ring_bond_pairs,
            &mut done_paths,
            &mut collected,
        )?;
    }
    for bid2 in 0..mol.num_bonds() {
        let aid2 = mol.bonds[bid2].begin;
        let aid3 = mol.bonds[bid2].end;
        for &bid1 in &mol.atom_bonds[aid2] {
            if bid1 == bid2 {
                continue;
            }
            for &bid3 in &mol.atom_bonds[aid3] {
                if bid3 == bid2 {
                    continue;
                }
                let (b1, b2, b3) = (bid1 as u64, bid2 as u64, bid3 as u64);
                let id1 = b1 * nb * nb + b2 * nb + b3;
                let id2 = b3 * nb * nb + b2 * nb + b1;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }
                set_remaining_14_path(
                    cx,
                    bid1,
                    bid2,
                    bid3,
                    acc,
                    m,
                    macro14,
                    force_trans_amides,
                    &bid_is_macrocycle,
                    &ring_bond_pairs,
                    &cis_ring_bond_pairs,
                    &mut collected,
                )?;
            }
        }
    }
    for bounds in collected.into_values() {
        let merged = merge_14_bounds(bounds);
        check_and_set_bounds(merged.aid1, merged.aid4, merged.lower, merged.upper, m)?;
    }
    Ok(())
}
