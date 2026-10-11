use super::*;

pub(super) fn add_torsions(
    field: &mut RdkitUffField,
    mol: &Molecule,
    adjacency: &UffAdjacency,
    params: &[Option<AtomicParams>],
    atomic_numbers: &[usize],
    hybridization: &[Option<u8>],
) {
    let torsion_end = |atom: usize| {
        adjacency[atom].len() > 1
            && !adjacency[atom]
                .iter()
                .any(|&(_, bond)| mol.bond(BondIdx(bond as u32)).order == BondOrder::Triple)
    };
    let mut seen = std::collections::HashSet::new();
    for j in 0..adjacency.len() {
        if !torsion_end(j) {
            continue;
        }
        for &(k, central_bond) in &adjacency[j] {
            if !torsion_end(k)
                || !seen.insert((j.min(k), j.max(k)))
                || !matches!(hybridization[j], Some(2) | Some(3))
                || !matches!(hybridization[k], Some(2) | Some(3))
            {
                continue;
            }
            let (Some(j_params), Some(k_params)) = (&params[j], &params[k]) else {
                continue;
            };
            add_torsions_for_bond(
                field,
                mol,
                adjacency,
                atomic_numbers,
                hybridization,
                j_params,
                k_params,
                j,
                k,
                central_bond,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn add_torsions_for_bond(
    field: &mut RdkitUffField,
    mol: &Molecule,
    adjacency: &UffAdjacency,
    atomic_numbers: &[usize],
    hybridization: &[Option<u8>],
    j_params: &AtomicParams,
    k_params: &AtomicParams,
    j: usize,
    k: usize,
    central_bond: usize,
) {
    let start = field.torsions.len();
    let central_order = bond_as_double(mol.bond(BondIdx(central_bond as u32)).order);
    for &(i, first_bond) in &adjacency[j] {
        if first_bond == central_bond {
            continue;
        }
        for &(l, second_bond) in &adjacency[k] {
            if second_bond == central_bond || second_bond == first_bond || l == i {
                continue;
            }
            let end_sp2 = hybridization[i] == Some(2) || hybridization[l] == Some(2);
            let (force_constant, order, cos_term) = torsion_params(
                central_order,
                atomic_numbers[j],
                atomic_numbers[k],
                hybridization[j],
                hybridization[k],
                j_params,
                k_params,
                end_sp2,
            );
            field.torsions.push(TorsionTerm {
                i,
                j,
                k,
                l,
                force_constant,
                order,
                cos_term,
            });
        }
    }
    let count = (field.torsions.len() - start) as f64;
    for term in &mut field.torsions[start..] {
        term.force_constant /= count;
    }
}
