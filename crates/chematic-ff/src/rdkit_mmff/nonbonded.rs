//! Preparation of the fixed nonbonded pair list.
use super::{
    BondAdjacency, NonbondedTerm, P3, fragments, is_double_zero, topological_distances, vdw_params,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn mmff_nonbonded_terms(
    coords: &[[f64; 3]],
    adjacency: &BondAdjacency,
    types: &[u8],
    charges: &[f64],
    threshold: f64,
    ignore_interfrag: bool,
    terms: [bool; 7],
    exact_charges: bool,
) -> Vec<NonbondedTerm> {
    let n = adjacency.len();
    let topology = topological_distances(n, adjacency);
    let fragments = fragments(n, adjacency);
    let mut output = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if ignore_interfrag && fragments[i] != fragments[j] {
                continue;
            }
            let relation = topology[i * n + j];
            if matches!(relation, 1 | 2) {
                continue;
            }
            let delta = P3 {
                x: coords[i][0] - coords[j][0],
                y: coords[i][1] - coords[j][1],
                z: coords[i][2] - coords[j][2],
            };
            if delta.length() > threshold {
                continue;
            }
            let vdw = terms[5].then(|| vdw_params(types[i], types[j])).flatten();
            let charged = if exact_charges {
                charges[i] != 0.0 && charges[j] != 0.0
            } else {
                !is_double_zero(charges[i]) && !is_double_zero(charges[j])
            };
            let ele = (terms[6] && charged).then_some((
                charges[i] * charges[j],
                charges[i] * charges[j],
                relation == 3,
            ));
            if vdw.is_some() || ele.is_some() {
                output.push(NonbondedTerm { i, j, vdw, ele });
            }
        }
    }
    output
}
