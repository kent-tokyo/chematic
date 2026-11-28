//! Verified stereo repair passes (split out of `stereo_constraints` to keep file size in budget).

use super::*;

/// BFS component reachable from `starts`, never stepping into `blocked`. Returns
/// `None` if the search reaches any atom in `forbidden` -- meaning the starting
/// point(s) are not a clean graph bridge off `blocked` (some path loops back to the
/// rest of the molecule through a ring, not just through `blocked` itself).
pub(super) fn bridge_component(
    mol: &Molecule,
    starts: &[AtomIdx],
    blocked: AtomIdx,
    forbidden: &[AtomIdx],
) -> Option<HashSet<AtomIdx>> {
    let mut seen: HashSet<AtomIdx> = starts.iter().copied().collect();
    let mut queue: VecDeque<AtomIdx> = starts.iter().copied().collect();
    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if nb == blocked {
                continue;
            }
            if forbidden.contains(&nb) {
                return None;
            }
            if seen.insert(nb) {
                queue.push_back(nb);
            }
        }
    }
    Some(seen)
}

fn repair_tetrahedral_center(
    mol: &Molecule,
    coords: &mut Coords3D,
    idx: AtomIdx,
    variant: usize,
) -> Result<RepairedElement, RepairRejectionReason> {
    let order = tetrahedral_constraint_for(mol, idx)
        .expect("already validated by verify_stereo before repair was attempted")
        .order;
    let real_neighbors: Vec<u32> = order
        .iter()
        .copied()
        .filter(|&n| n != STEREO_H_SENTINEL)
        .collect();

    // Stable size order retains the historical smallest-component choice at
    // variant zero; other variants allow bounded, deterministic retries.
    let mut candidates: Vec<(u32, HashSet<AtomIdx>)> = Vec::new();
    for &cand in &real_neighbors {
        let others: Vec<AtomIdx> = real_neighbors
            .iter()
            .copied()
            .filter(|&n| n != cand)
            .map(AtomIdx)
            .collect();
        if let Some(component) = bridge_component(mol, &[AtomIdx(cand)], idx, &others) {
            candidates.push((cand, component));
        }
    }
    if candidates.is_empty() {
        return Err(RepairRejectionReason::NoBridgeEligibleSubstituent);
    }
    candidates.sort_by_key(|(_, component)| component.len());
    let selected = variant % candidates.len();
    let (moved, component) = candidates.swap_remove(selected);

    // Declared-order positions BEFORE any mutation (phantom computed from the
    // original, not-yet-moved real neighbor positions).
    let center = coords.get(idx);
    let positions4: Vec<Point3> = order
        .iter()
        .map(|&n| {
            if n == STEREO_H_SENTINEL {
                let real_positions: Vec<Point3> = real_neighbors
                    .iter()
                    .map(|&r| coords.get(AtomIdx(r)))
                    .collect();
                phantom_neighbor_position(center, &real_positions)
                    .expect("already validated by verify_stereo before repair was attempted")
            } else {
                coords.get(AtomIdx(n))
            }
        })
        .collect();

    let moved_slot = order
        .iter()
        .position(|&n| n == moved)
        .expect("moved came from order");
    let others: Vec<Point3> = (0..4)
        .filter(|&i| i != moved_slot)
        .map(|i| positions4[i])
        .collect();

    let normal = others[1]
        .sub(&others[0])
        .cross(&others[2].sub(&others[0]))
        .try_normalize()
        .ok_or(RepairRejectionReason::DegenerateReflectionPlane)?;

    let p_moved = positions4[moved_slot];
    let signed_dist = normal.dot(&p_moved.sub(&others[0]));
    let reflected = p_moved.sub(&normal.scale(2.0 * signed_dist));
    let direction = reflected
        .sub(&center)
        .try_normalize()
        .ok_or(RepairRejectionReason::DegenerateReflectionDirection)?;
    let bond_len = p_moved.sub(&center).norm();
    let new_p_moved = center.add(&direction.scale(bond_len));
    let displacement = new_p_moved.sub(&p_moved);
    let displacement_mag = displacement.norm();

    for &atom in &component {
        let old = coords.get(atom);
        coords.set(atom, old.add(&displacement));
    }

    match verify_tetrahedral_center(mol, coords, idx) {
        StereoStatus::Satisfied => Ok(RepairedElement {
            element: StereoElement::Tetrahedral(idx),
            atoms_moved: component.len(),
            max_displacement: displacement_mag,
        }),
        _ => Err(RepairRejectionReason::StillViolatedAfterRepair),
    }
}

fn rotate_180_around_axis(p: Point3, axis_point: Point3, axis_dir: Point3) -> Point3 {
    let v = p.sub(&axis_point);
    let parallel = axis_dir.scale(v.dot(&axis_dir));
    let perp = v.sub(&parallel);
    let new_v = parallel.sub(&perp);
    axis_point.add(&new_v)
}

fn repair_double_bond(
    mol: &Molecule,
    coords: &mut Coords3D,
    bond_idx: BondIdx,
    variant: usize,
) -> Result<RepairedElement, RepairRejectionReason> {
    let bond = mol.bond(bond_idx);
    let a1 = bond.atom1;
    let a2 = bond.atom2;

    let a1_subs: Vec<AtomIdx> = mol
        .neighbors(a1)
        .filter(|&(nb, _)| nb != a2)
        .map(|(nb, _)| nb)
        .collect();
    let a2_subs: Vec<AtomIdx> = mol
        .neighbors(a2)
        .filter(|&(nb, _)| nb != a1)
        .map(|(nb, _)| nb)
        .collect();

    let a1_component = if a1_subs.is_empty() {
        None
    } else {
        bridge_component(mol, &a1_subs, a1, &[a2])
    };
    let a2_component = if a2_subs.is_empty() {
        None
    } else {
        bridge_component(mol, &a2_subs, a2, &[a1])
    };

    let (pivot, other_end, component) = match (a1_component, a2_component) {
        (Some(c1), Some(c2)) => {
            if (c1.len() <= c2.len()) == variant.is_multiple_of(2) {
                (a1, a2, c1)
            } else {
                (a2, a1, c2)
            }
        }
        (Some(c1), None) => (a1, a2, c1),
        (None, Some(c2)) => (a2, a1, c2),
        (None, None) => return Err(RepairRejectionReason::NoBridgeEligibleSubstituent),
    };

    let axis_point = coords.get(pivot);
    let axis_dir = coords
        .get(other_end)
        .sub(&axis_point)
        .try_normalize()
        .ok_or(RepairRejectionReason::DegenerateReflectionDirection)?;

    let mut max_disp = 0.0_f64;
    for &atom in &component {
        let old = coords.get(atom);
        let new = rotate_180_around_axis(old, axis_point, axis_dir);
        let disp = new.sub(&old).norm();
        if disp > max_disp {
            max_disp = disp;
        }
        coords.set(atom, new);
    }

    match verify_double_bond(mol, coords, bond_idx) {
        Some(StereoStatus::Satisfied) => Ok(RepairedElement {
            element: StereoElement::DoubleBond(bond_idx),
            atoms_moved: component.len(),
            max_displacement: max_disp,
        }),
        _ => Err(RepairRejectionReason::StillViolatedAfterRepair),
    }
}

/// Verify `(mol, coords)`, then repair every `Violated` element found. Returns
/// `Ok(RepairOutcome)` only when every originally-`Violated` element is now
/// `Satisfied` **and** no originally-`Satisfied` element regressed; otherwise
/// `Err(RepairFailureReport)` naming exactly which element(s) failed and why.
/// `Unevaluable` elements are left untouched either way (nothing to compare a repair
/// against).
pub fn repair_stereo(
    mol: &Molecule,
    coords: &Coords3D,
) -> Result<RepairOutcome, RepairFailureReport> {
    repair_stereo_with_variant(mol, coords, 0)
}

/// Run the same verified repair with a deterministic alternative choice of
/// bridge-connected substituent at each stereo element. Variant zero retains
/// the default behavior; nonzero variants do not guarantee clash-free output.
pub fn repair_stereo_with_variant(
    mol: &Molecule,
    coords: &Coords3D,
    variant: usize,
) -> Result<RepairOutcome, RepairFailureReport> {
    let initial = verify_stereo(mol, coords);
    let mut working = coords.clone();
    let mut repaired = Vec::new();
    let mut failures: Vec<(StereoElement, RepairRejectionReason)> = Vec::new();

    let unevaluable = unevaluable_elements(&initial);
    apply_stereo_repairs(
        mol,
        &initial,
        &mut working,
        &mut repaired,
        &mut failures,
        variant,
    );

    // Unconditional final reconciliation pass -- runs whether or not the loop above
    // already produced failures, and checks BOTH directions:
    //   (a) an originally-Satisfied element disturbed by some OTHER element's repair.
    //   (b) an originally-Violated element that `repair_tetrahedral_center`/
    //       `repair_double_bond`'s OWN immediate re-verify confirmed fixed (already
    //       pushed to `repaired`), but which a STILL LATER, unrelated repair in this
    //       same call silently moved back into violation.
    // (b) is not hypothetical: found via this module's own validation harness on
    // real corpus data (cholesterol) -- two ring-fused stereocenters whose bridge-
    // safe substituents can still overlap enough that fixing one undoes an
    // already-confirmed fix on another, even though each individual repair's own
    // bridge check passed. Both directions are reported the same way
    // (`RepairCausedNewViolation`): from the caller's perspective, both mean "some
    // OTHER repair in this call broke something that was fine a moment ago." This
    // is a detect-and-fail-closed pass, not an iterate-to-a-fixed-point one -- a
    // caller can retry `repair_stereo` on the returned `partial_coords` if a fixed
    // point is reachable, but this function does not loop internally.
    let final_check = verify_stereo(mol, &working);

    let repaired = retain_successful_repairs(repaired, &final_check, &mut failures);
    record_stereo_regressions(&initial, &final_check, &mut failures);

    if !failures.is_empty() {
        return Err(RepairFailureReport {
            partial_coords: working,
            repaired,
            failures,
        });
    }

    Ok(RepairOutcome {
        coords: working,
        repaired,
        unevaluable,
    })
}

fn unevaluable_elements(verification: &StereoVerification) -> Vec<StereoElement> {
    verification
        .tetrahedral
        .iter()
        .filter(|result| matches!(result.status, StereoStatus::Unevaluable(_)))
        .map(|result| StereoElement::Tetrahedral(result.atom))
        .chain(
            verification
                .double_bond
                .iter()
                .filter(|result| matches!(result.status, StereoStatus::Unevaluable(_)))
                .map(|result| StereoElement::DoubleBond(result.bond)),
        )
        .collect()
}

fn apply_stereo_repairs(
    mol: &Molecule,
    initial: &StereoVerification,
    coords: &mut Coords3D,
    repaired: &mut Vec<RepairedElement>,
    failures: &mut Vec<(StereoElement, RepairRejectionReason)>,
    variant: usize,
) {
    for result in &initial.tetrahedral {
        if result.status.is_violated() {
            match repair_tetrahedral_center(mol, coords, result.atom, variant) {
                Ok(info) => repaired.push(info),
                Err(reason) => failures.push((StereoElement::Tetrahedral(result.atom), reason)),
            }
        }
    }
    for result in &initial.double_bond {
        if result.status.is_violated() {
            match repair_double_bond(mol, coords, result.bond, variant) {
                Ok(info) => repaired.push(info),
                Err(reason) => failures.push((StereoElement::DoubleBond(result.bond), reason)),
            }
        }
    }
}

fn retain_successful_repairs(
    repaired: Vec<RepairedElement>,
    final_check: &StereoVerification,
    failures: &mut Vec<(StereoElement, RepairRejectionReason)>,
) -> Vec<RepairedElement> {
    repaired
        .into_iter()
        .filter(|info| {
            let satisfied = match info.element {
                StereoElement::Tetrahedral(atom) => final_check
                    .tetrahedral
                    .iter()
                    .find(|result| result.atom == atom)
                    .is_some_and(|result| result.status.is_satisfied()),
                StereoElement::DoubleBond(bond) => final_check
                    .double_bond
                    .iter()
                    .find(|result| result.bond == bond)
                    .is_some_and(|result| result.status.is_satisfied()),
            };
            if !satisfied {
                failures.push((
                    info.element,
                    RepairRejectionReason::RepairCausedNewViolation,
                ));
            }
            satisfied
        })
        .collect()
}

pub(super) fn record_stereo_regressions(
    initial: &StereoVerification,
    final_check: &StereoVerification,
    failures: &mut Vec<(StereoElement, RepairRejectionReason)>,
) {
    for result in &initial.tetrahedral {
        let final_status = final_check
            .tetrahedral
            .iter()
            .find(|candidate| candidate.atom == result.atom)
            .map(|candidate| candidate.status);
        let regressed = stereo_status_regressed(result.status, final_status);
        if regressed {
            failures.push((
                StereoElement::Tetrahedral(result.atom),
                RepairRejectionReason::RepairCausedNewViolation,
            ));
        }
    }
    for result in &initial.double_bond {
        let final_status = final_check
            .double_bond
            .iter()
            .find(|candidate| candidate.bond == result.bond)
            .map(|candidate| candidate.status);
        let regressed = stereo_status_regressed(result.status, final_status);
        if regressed {
            failures.push((
                StereoElement::DoubleBond(result.bond),
                RepairRejectionReason::RepairCausedNewViolation,
            ));
        }
    }
}

pub(super) fn stereo_status_regressed(
    initial: StereoStatus,
    final_status: Option<StereoStatus>,
) -> bool {
    match initial {
        StereoStatus::Satisfied => !final_status.is_some_and(|s| s.is_satisfied()),
        StereoStatus::Unevaluable(_) => final_status.is_some_and(|s| s.is_violated()),
        StereoStatus::Violated => false, // Repaired elements are checked separately.
    }
}
