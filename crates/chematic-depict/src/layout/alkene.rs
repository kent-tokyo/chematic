use super::*;
use chematic_core::ez_markers::EzFact;

pub(super) fn direction_facts(mol: &Molecule) -> Vec<EzFact> {
    if mol.has_bond_directions()
        || mol.bonds().any(|(_, bond)| {
            matches!(
                bond.order,
                chematic_core::BondOrder::Up | chematic_core::BondOrder::Down
            )
        })
    {
        chematic_core::ez_markers::read_ez_facts(mol)
            .map(|(facts, _)| facts)
            .unwrap_or_default()
    } else {
        Vec::new()
    }
}

pub(super) fn match_alkene_geometry(
    mol: &Molecule,
    facts: &[&EzFact],
    atoms: &[AtomIdx],
    bonds: &[(AtomIdx, AtomIdx)],
    placed: &mut Vec<Option<Point>>,
) {
    let at = |coords: &[Option<Point>], a: AtomIdx| coords[a.0 as usize].unwrap();
    let matches = |fact: &EzFact, coords: &[Option<Point>]| {
        let [left, right] = fact.ends;
        let p = at(coords, left.atom);
        let q = at(coords, right.atom);
        let cross = |end: chematic_core::ez_markers::EzFactEnd| {
            let sub = mol.bond(end.reference).other(end.atom).unwrap();
            let s = at(coords, sub);
            (q.x - p.x) * (s.y - p.y) - (q.y - p.y) * (s.x - p.x)
        };
        let (a, b) = (cross(left), cross(right));
        a.abs() > 1e-6 && b.abs() > 1e-6 && ((a > 0.0) == (b > 0.0)) == (left.up == right.up)
    };

    for &fact in facts {
        if matches(fact, placed) {
            continue;
        }
        let [left, right] = fact.ends.map(|end| end.atom);
        let mut best = None;
        for (end, far) in [(left, right), (right, left)] {
            let (group, ring) = branch(mol, end, far, None, true);
            let (group, pivot) = if ring {
                // A simple ring can flip the alkene end about the line from
                // its other ring neighbour to the far end of the double bond.
                let mut pivots = mol
                    .neighbors(end)
                    .map(|(a, _)| a)
                    .filter(|&a| a != far && branch(mol, a, far, Some(end), false).1);
                let (Some(pivot), None) = (pivots.next(), pivots.next()) else {
                    continue;
                };
                let (group, reaches_far) = branch(mol, end, far, Some(pivot), true);
                if reaches_far {
                    continue;
                }
                (group, pivot)
            } else {
                (group, end)
            };

            let p = at(placed, pivot);
            let q = at(placed, far);
            let (dx, dy) = (q.x - p.x, q.y - p.y);
            let length_sq = dx * dx + dy * dy;
            if length_sq < 1e-6 {
                continue;
            }
            let mut trial = placed.clone();
            for &a in &group {
                let s = at(placed, a);
                let t = ((s.x - p.x) * dx + (s.y - p.y) * dy) / length_sq;
                trial[a.0 as usize] = Some(Point::new(
                    2.0 * (p.x + t * dx) - s.x,
                    2.0 * (p.y + t * dy) - s.y,
                ));
            }
            if !matches(fact, &trial)
                || facts
                    .iter()
                    .any(|p| matches(p, placed) && !matches(p, &trial))
            {
                continue;
            }
            let (clashes, crossings) = layout_defects(atoms, bonds, &Bonded(mol), &trial);
            let score = (clashes, crossings, group.len());
            if best.as_ref().is_none_or(|(old, _)| score < *old) {
                best = Some((score, trial));
            }
        }
        if let Some((_, trial)) = best {
            *placed = trial;
        }
    }
}

fn branch(
    mol: &Molecule,
    start: AtomIdx,
    far: AtomIdx,
    blocked: Option<AtomIdx>,
    skip_direct: bool,
) -> (Vec<AtomIdx>, bool) {
    let mut seen = HashSet::default();
    seen.insert(start);
    seen.extend(blocked);
    let mut stack = vec![start];
    let mut atoms = Vec::new();
    let mut reaches_far = false;
    while let Some(a) = stack.pop() {
        atoms.push(a);
        for (next, _) in mol.neighbors(a) {
            if next == far {
                reaches_far |= a != start || !skip_direct;
            } else if seen.insert(next) {
                stack.push(next);
            }
        }
    }
    (atoms, reaches_far)
}
