//! RDKit-exact graph descriptors from `rdkit.Chem.GraphDescriptors`:
//! [`balaban_j`], [`bertz_ct`] and [`ipc`].
//!
//! RDKit implements these in Python on top of numpy and a few C++ helpers
//! (`Chem.GetDistanceMatrix`, `rdInfoTheory.InfoEntropy`). Each function here
//! reproduces the reference floating-point operation order so the results are
//! bit-identical to RDKit 2026.03 (CPython 3.13, numpy 2.x):
//!
//! * `GetDistanceMatrix` is RDKit's double-buffered Floyd–Warshall
//!   (`Matrices.cpp`), whose path-length sums depend on the composition order.
//! * Python's builtin `sum` over numpy values is a naive left-to-right sum;
//!   `numpy.trace` uses numpy's pairwise summation.
//! * `InfoEntropy` is the C++ template in `InfoGainFuncs.h`.

use chematic_core::{BondOrder, Molecule};
use std::collections::HashMap;

/// RDKit's `LOCAL_INF` distance for unconnected atom pairs.
const LOCAL_INF: f64 = 1e8;

/// Whether a bond is aromatic in RDKit's perception (descriptor view).
fn is_aromatic(order: BondOrder) -> bool {
    order == BondOrder::Aromatic
}

/// RDKit `Bond::getBondTypeAsDouble` for a non-aromatic bond.
fn bond_type_as_double(order: BondOrder) -> f64 {
    match order {
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Zero => 0.0,
        BondOrder::Aromatic => 1.5,
        _ => 1.0,
    }
}

/// `float(bond.GetBondType())` (the `BondType` enum value) for a
/// non-aromatic bond, as `GraphDescriptors._LookUpBondOrder` uses it.
fn bond_type_enum_value(order: BondOrder) -> f64 {
    match order {
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Zero => 21.0,
        BondOrder::Dative => 17.0,
        BondOrder::Aromatic => 12.0,
        _ => 1.0,
    }
}

/// `Chem.GetDistanceMatrix(mol, useBO, useAtomWts=0)`, row-major `n × n`.
///
/// Port of `MolOps::getDistanceMat` + `FloydWarshall` (`Matrices.cpp`),
/// including its double buffering and `v1 <= v2` tie rule, so every summed
/// path length is composed exactly as RDKit composes it.
fn rdkit_distance_matrix(view: &Molecule, use_bo: bool) -> Vec<f64> {
    let n = view.atom_count();
    let mut last = vec![LOCAL_INF; n * n];
    for i in 0..n {
        last[i * n + i] = 0.0;
    }
    for (_, bond) in view.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let contrib = if use_bo {
            if !is_aromatic(bond.order) {
                1.0 / bond_type_as_double(bond.order)
            } else {
                2.0 / 3.0
            }
        } else {
            1.0
        };
        last[i * n + j] = contrib;
        last[j * n + i] = contrib;
    }
    let mut curr = vec![0.0; n * n];
    for k in 0..n {
        let ktab = k * n;
        for i in 0..n {
            let itab = i * n;
            let dik = last[itab + k];
            for j in 0..n {
                let v1 = last[itab + j];
                let v2 = dik + last[ktab + j];
                curr[itab + j] = if v1 <= v2 { v1 } else { v2 };
            }
        }
        std::mem::swap(&mut curr, &mut last);
    }
    last
}

/// `rdInfoTheory.InfoEntropy` on a float64 array (`InfoGainFuncs.h`).
fn info_entropy(values: &[f64]) -> f64 {
    let mut n_instances = 0.0f64;
    for &v in values {
        n_instances += v;
    }
    let mut accum = 0.0f64;
    if n_instances != 0.0 {
        for &v in values {
            let d = v / n_instances;
            if d != 0.0 {
                accum += -d * d.ln();
            }
        }
    }
    accum / std::f64::consts::LN_2
}

/// numpy's pairwise summation (`DOUBLE_pairwise_sum`), as `numpy.add.reduce`
/// (and so `numpy.trace`) applies it to a float64 vector.
fn numpy_pairwise_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut res = 0.0;
        for &x in a {
            res += x;
        }
        res
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - (n % 8) {
            for (k, rk) in r.iter_mut().enumerate() {
                *rk += a[i + k];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for &x in &a[i..] {
            res += x;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        numpy_pairwise_sum(&a[..n2]) + numpy_pairwise_sum(&a[n2..])
    }
}

/// RDKit `GraphDescriptors.BalabanJ` (Balaban, *Chem. Phys. Lett.* **89**,
/// 399–404, 1982), bit-identical to RDKit.
///
/// J = q / (μ + 1) · Σ_{bonded i ≤ j} 1/√(sᵢ·sⱼ), with sᵢ the column sums of
/// the bond-order-weighted distance matrix (edge weight 1/order, aromatic
/// 2/3; unconnected pairs count RDKit's `1e8`), q the bond count and
/// μ = q − n + 1.
///
/// Returns 0.0 for molecules larger than 1000 atoms (the O(n³) all-pairs
/// path computation is impractical at that scale).
pub fn balaban_j(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n > 1000 {
        return 0.0;
    }
    let view = crate::descriptors::descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    let d = rdkit_distance_matrix(view, true);
    // `sum(mat)`: rows added elementwise, left to right.
    let mut s = vec![0.0f64; n];
    for i in 0..n {
        for (j, sj) in s.iter_mut().enumerate() {
            *sj += d[i * n + j];
        }
    }
    let q = view.bond_count() as i64;
    let mu = q - n as i64 + 1;
    let mut adjacent = vec![false; n * n];
    for (_, bond) in view.bonds() {
        let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        adjacent[i * n + j] = true;
        adjacent[j * n + i] = true;
    }
    let mut sum = 0.0f64;
    for i in 0..n {
        let si = s[i];
        for j in i..n {
            if adjacent[i * n + j] {
                sum += 1.0 / (si * s[j]).sqrt();
            }
        }
    }
    if mu + 1 != 0 {
        q as f64 / (mu + 1) as f64 * sum
    } else {
        0.0
    }
}

/// Connection-dictionary key of [`bertz_ct`]: a bond (pair of symmetry
/// classes) or a two-bond path (outer classes plus hinge class).
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum ConnectionKey {
    Pair(usize, usize),
    Triple(usize, usize, usize),
}

/// RDKit `GraphDescriptors.BertzCT` (Bertz, *J. Am. Chem. Soc.* **103**,
/// 3599–3601, 1981) with the default `cutoff=100`, bit-identical to RDKit.
///
/// Atoms are grouped into symmetry classes by their sorted bond-order
/// weighted distances (formatted `%.4f`, first `cutoff` entries); bond and
/// two-bond-path "connections" between classes are counted with aromatic
/// bonds as order 1.5; the result is the sum of the class-weighted
/// information contents of connections and of atom types.
///
/// Returns 0.0 for fewer than two atoms (RDKit returns the integer 0), and
/// for molecules larger than 1000 atoms.
pub fn bertz_ct(mol: &Molecule) -> f64 {
    const CUTOFF: usize = 100;
    let n = mol.atom_count();
    if !(2..=1000).contains(&n) {
        return 0.0;
    }
    let view = crate::descriptors::descriptor_aromaticity(mol);
    let view: &Molecule = &view;

    // _CreateBondDictEtc
    let mut bond_order: HashMap<(usize, usize), f64> = HashMap::new();
    let mut neighbors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (_, bond) in view.bonds() {
        let (mut a1, mut a2) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if a1 > a2 {
            std::mem::swap(&mut a1, &mut a2);
        }
        let order = if is_aromatic(bond.order) {
            1.5
        } else {
            bond_type_enum_value(bond.order)
        };
        bond_order.insert((a1, a2), order);
        if !neighbors[a1].contains(&a2) {
            neighbors[a1].push(a2);
        }
        if !neighbors[a2].contains(&a1) {
            neighbors[a2].push(a1);
        }
    }
    for nl in &mut neighbors {
        nl.sort_unstable();
    }
    let lookup = |a: usize, b: usize| bond_order[&(a.min(b), a.max(b))];

    // _AssignSymmetryClasses on the "Balaban" (bond-order) distance matrix.
    let d = rdkit_distance_matrix(view, true);
    let mut keys_seen: HashMap<Vec<String>, usize> = HashMap::new();
    let mut sym = vec![0usize; n];
    for i in 0..n {
        let mut row: Vec<f64> = d[i * n..(i + 1) * n].to_vec();
        row.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let key: Vec<String> = row
            .iter()
            .take(CUTOFF)
            .map(|x| python_percent_4f(*x))
            .collect();
        let next = keys_seen.len();
        sym[i] = *keys_seen.entry(key).or_insert(next) + 1;
    }

    // Insertion-ordered dictionaries, as Python's dict.
    let mut atom_types: Vec<(u8, f64)> = Vec::new();
    let mut conn_index: HashMap<ConnectionKey, usize> = HashMap::new();
    let mut connections: Vec<f64> = Vec::new();
    let mut add_conn = |key: ConnectionKey, v: f64| match conn_index.get(&key) {
        Some(&k) => connections[k] += v,
        None => {
            conn_index.insert(key, connections.len());
            connections.push(v);
        }
    };
    for (atom_idx, (_, atom)) in view.atoms().enumerate() {
        let z = if atom.wildcard {
            0
        } else {
            atom.element.atomic_number()
        };
        match atom_types.iter_mut().find(|(k, _)| *k == z) {
            Some((_, c)) => *c += 1.0,
            None => atom_types.push((z, 1.0)),
        }
        let hinge = sym[atom_idx];
        let nbrs = &neighbors[atom_idx];
        for (i, &ni) in nbrs.iter().enumerate() {
            let ni_class = sym[ni];
            let bi = lookup(atom_idx, ni);
            if bi > 1.0 && ni > atom_idx {
                let num = bi * (bi - 1.0) / 2.0;
                add_conn(
                    ConnectionKey::Pair(hinge.min(ni_class), hinge.max(ni_class)),
                    num,
                );
            }
            for &nj in &nbrs[i + 1..] {
                let nj_class = sym[nj];
                let bj = lookup(atom_idx, nj);
                add_conn(
                    ConnectionKey::Triple(ni_class.min(nj_class), hinge, ni_class.max(nj_class)),
                    bi * bj,
                );
            }
        }
    }
    if connections.is_empty() {
        connections.push(1.0);
    }

    // _CalculateEntropies. The connection counts are multiples of 1/4, so
    // CPython's compensated float `sum` equals the plain sum here.
    let mut tot = 0.0f64;
    for &c in &connections {
        tot += c;
    }
    let connection_ie = tot * (info_entropy(&connections) + tot.ln() / std::f64::consts::LN_2);
    let type_counts: Vec<f64> = atom_types.iter().map(|&(_, c)| c).collect();
    let atom_type_ie = n as f64 * info_entropy(&type_counts);
    atom_type_ie + connection_ie
}

/// Python's `'%.4f' % x`.
fn python_percent_4f(x: f64) -> String {
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x.is_nan() {
        return "nan".into();
    }
    format!("{x:.4}")
}

/// RDKit `GraphDescriptors.Ipc` (Bonchev & Trinajstić, *J. Chem. Phys.*
/// **67**, 4517–4533, 1977), bit-identical to RDKit for molecules whose
/// characteristic polynomial stays in the range where the reference
/// `numpy.dot` (BLAS) result is order-independent.
///
/// The characteristic polynomial of the hydrogen-suppressed adjacency
/// matrix is computed by `Graphs.CharacteristicPolynomial`'s
/// Le Verrier–Faddeev–Frame recursion; Ipc = Σ|cᵢ| · H(|c|).
///
/// Returns 0.0 for molecules larger than 1000 atoms.
pub fn ipc(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n > 1000 {
        return 0.0;
    }
    let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (_, bond) in mol.bonds() {
        let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if i != j && !nbrs[i].contains(&j) {
            nbrs[i].push(j);
            nbrs[j].push(i);
        }
    }
    for nb in &mut nbrs {
        nb.sort_unstable();
    }
    let mut res = vec![0.0f64; n + 1];
    res[0] = 1.0;
    // An = A (0/1); the first trace is 0.
    let mut an = vec![0.0f64; n * n];
    for (i, nb) in nbrs.iter().enumerate() {
        for &j in nb {
            an[i * n + j] = 1.0;
        }
    }
    let mut diag = vec![0.0f64; n];
    let mut bn = vec![0.0f64; n * n];
    for step in 1..=n {
        for (i, di) in diag.iter_mut().enumerate() {
            *di = an[i * n + i];
        }
        let c = 1.0 / step as f64 * numpy_pairwise_sum(&diag);
        res[step] = c;
        if step == n {
            break;
        }
        // Bn = An - c·I
        bn.copy_from_slice(&an);
        for i in 0..n {
            bn[i * n + i] = an[i * n + i] - c;
        }
        // An = A · Bn, accumulating over k in index order (as the BLAS
        // kernel does for each output element).
        for (i, nb) in nbrs.iter().enumerate() {
            let row = &mut an[i * n..(i + 1) * n];
            row.fill(0.0);
            for &k in nb {
                let brow = &bn[k * n..(k + 1) * n];
                for (r, b) in row.iter_mut().zip(brow) {
                    *r += *b;
                }
            }
        }
    }
    for r in res.iter_mut().skip(1) {
        *r = -*r;
    }
    let cpoly: Vec<f64> = res.iter().map(|x| x.abs()).collect();
    let mut total = 0.0f64;
    for &x in &cpoly {
        total += x;
    }
    total * info_entropy(&cpoly)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mol(smi: &str) -> Molecule {
        chematic_smiles::parse(smi).unwrap()
    }

    // Reference values: RDKit 2026.03.1 `GraphDescriptors`.
    #[test]
    fn balaban_j_matches_rdkit() {
        assert_eq!(balaban_j(&mol("c1ccccc1")), 3.000000000000001);
        assert_eq!(balaban_j(&mol("C1CCCCC1")), 2.0);
        assert_eq!(balaban_j(&mol("CCC")), 1.6329931618554523);
        assert_eq!(balaban_j(&mol("C")), 0.0);
        assert_eq!(balaban_j(&mol("CC(=O)Oc1ccccc1C(=O)O")), 3.0435273546341013);
    }

    #[test]
    fn ipc_matches_rdkit() {
        assert_eq!(ipc(&mol("c1ccccc1")), 34.3994618804395);
        assert_eq!(ipc(&mol("CC")), 2.0);
        assert_eq!(ipc(&mol("C")), 0.0);
        assert_eq!(ipc(&mol("CC(=O)Oc1ccccc1C(=O)O")), 729.6807528797516);
    }

    #[test]
    fn bertz_ct_charged_amidine_matches_rdkit() {
        assert_eq!(bertz_ct(&mol("C[NH+]=C(N)c1ccccc1")), 226.27272407307768);
    }
}
