//! B5 Advanced: MMFF94 batch processing, ensemble calculations, and advanced electrostatics.
//!
//! Provides high-performance operations for:
//! - Batch molecular property calculations
//! - Molecular ensemble properties (multiple conformers)
//! - Full electrostatic scaling matrices
//! - Cached/vectorized computations

use chematic_core::{AtomIdx, Molecule};

use crate::mmff94::MMFF94Type;
use crate::mmff94_params::mmff94_electrostatic_scaling_1_4;

/// Electrostatic scaling matrix for a complete molecule.
/// Stores 1-4 and 1-3 scaling factors for all atom pairs.
#[derive(Clone, Debug)]
pub struct ElectrostaticMatrix {
    /// Scaling factors: matrix[i][j] contains dielectric scaling from atom i to j
    pub scaling: Vec<Vec<f64>>,
    /// 1-4 pair indicators: true if atoms are 1-4 pairs (through 3 bonds)
    pub is_1_4_pair: Vec<Vec<bool>>,
    /// 1-3 pair indicators: true if atoms are bonded through one atom
    pub is_1_3_pair: Vec<Vec<bool>>,
}

impl ElectrostaticMatrix {
    /// Create an electrostatic scaling matrix for a molecule.
    ///
    /// Computes all pairwise dielectric scaling factors based on:
    /// - Topological distance (1-2, 1-3, 1-4, farther)
    /// - Atom types and bonding patterns
    pub fn new(mol: &Molecule, mmff94_types: &[MMFF94Type]) -> Self {
        let n = mol.atom_count();
        let mut scaling = vec![vec![4.0; n]; n]; // Default dielectric
        let mut is_1_4_pair = vec![vec![false; n]; n];
        let mut is_1_3_pair = vec![vec![false; n]; n];

        // Build adjacency matrix
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (_, bond) in mol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            adj[i].push(j);
            adj[j].push(i);
        }

        // Compute topological distances
        let dist = topological_distance_matrix(mol);

        // Fill scaling matrix
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    scaling[i][j] = 0.0; // Self-interaction (not used)
                    continue;
                }

                let d = dist[i][j];

                if d == 1 {
                    // 1-2 pairs: typically excluded
                    scaling[i][j] = 0.0;
                } else if d == 2 {
                    // 1-3 pairs: excluded or reduced
                    is_1_3_pair[i][j] = true;
                    scaling[i][j] = 0.0; // Usually excluded
                } else if d == 3 {
                    // 1-4 pairs: scaled
                    is_1_4_pair[i][j] = true;
                    let type_i = mmff94_types[i];
                    let type_j = mmff94_types[j];
                    let params = mmff94_electrostatic_scaling_1_4(type_i, type_j);
                    scaling[i][j] = params.scale_dielectric;
                } else {
                    // Farther pairs: full electrostatic
                    scaling[i][j] = 4.0; // Default organic environment
                }
            }
        }

        ElectrostaticMatrix {
            scaling,
            is_1_4_pair,
            is_1_3_pair,
        }
    }

    /// Get scaling factor for atom pair (i, j)
    pub fn get_scaling(&self, i: usize, j: usize) -> f64 {
        if i < self.scaling.len() && j < self.scaling[i].len() {
            self.scaling[i][j]
        } else {
            4.0 // Default
        }
    }

    /// Check if atoms are 1-4 pairs
    pub fn is_1_4(&self, i: usize, j: usize) -> bool {
        if i < self.is_1_4_pair.len() && j < self.is_1_4_pair[i].len() {
            self.is_1_4_pair[i][j]
        } else {
            false
        }
    }
}

/// Batch properties for multiple molecules.
///
/// Enables efficient calculation of properties for molecular ensembles
/// (e.g., multiple conformers, related molecules).
#[derive(Clone, Debug)]
pub struct MMFF94BatchProperties {
    /// Per-molecule charges
    pub charges: Vec<Vec<f64>>,
    /// Per-molecule bond dipoles
    pub bond_dipoles: Vec<Vec<(usize, usize, f64)>>,
    /// Per-molecule electrostatic matrices
    pub electrostatic_matrices: Vec<ElectrostaticMatrix>,
}

impl MMFF94BatchProperties {
    /// Create empty batch properties
    pub fn new() -> Self {
        MMFF94BatchProperties {
            charges: Vec::new(),
            bond_dipoles: Vec::new(),
            electrostatic_matrices: Vec::new(),
        }
    }

    /// Get total size (number of molecules in batch)
    pub fn len(&self) -> usize {
        self.charges.len()
    }

    /// Check if batch is empty
    pub fn is_empty(&self) -> bool {
        self.charges.is_empty()
    }

    /// Add properties for a new molecule
    pub fn push(
        &mut self,
        charges: Vec<f64>,
        bond_dipoles: Vec<(usize, usize, f64)>,
        electrostatic_matrix: ElectrostaticMatrix,
    ) {
        self.charges.push(charges);
        self.bond_dipoles.push(bond_dipoles);
        self.electrostatic_matrices.push(electrostatic_matrix);
    }

    /// Calculate average charges across all molecules in batch
    pub fn average_charges(&self) -> Vec<f64> {
        if self.charges.is_empty() {
            return Vec::new();
        }

        let n_atoms = self.charges[0].len();
        let mut avg = vec![0.0; n_atoms];

        for charges in &self.charges {
            for (i, &q) in charges.iter().enumerate() {
                avg[i] += q;
            }
        }

        let n_mols = self.charges.len() as f64;
        for q in &mut avg {
            *q /= n_mols;
        }

        avg
    }

    /// Calculate charge variance across molecules
    pub fn charge_variance(&self) -> Vec<f64> {
        if self.charges.is_empty() {
            return Vec::new();
        }

        let n_atoms = self.charges[0].len();
        let avg = self.average_charges();
        let mut var = vec![0.0; n_atoms];

        for charges in &self.charges {
            for (i, &q) in charges.iter().enumerate() {
                let diff = q - avg[i];
                var[i] += diff * diff;
            }
        }

        let n_mols = self.charges.len() as f64;
        for v in &mut var {
            *v /= n_mols;
            *v = v.sqrt(); // Standard deviation
        }

        var
    }
}

impl Default for MMFF94BatchProperties {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute topological distance matrix (BFS-based).
fn topological_distance_matrix(mol: &Molecule) -> Vec<Vec<i32>> {
    let n = mol.atom_count();
    let mut dist = vec![vec![i32::MAX; n]; n];

    // Initialize
    for (i, row) in dist.iter_mut().enumerate().take(n) {
        row[i] = 0;
    }

    // BFS from each atom
    for start in 0..n {
        let mut queue = vec![start];
        let mut visited = vec![false; n];
        visited[start] = true;

        while !queue.is_empty() {
            let current = queue.remove(0);
            let atom_idx = AtomIdx(current as u32);

            for (neighbor, _) in mol.neighbors(atom_idx) {
                let neighbor_idx = neighbor.0 as usize;
                if !visited[neighbor_idx] {
                    visited[neighbor_idx] = true;
                    dist[start][neighbor_idx] = dist[start][current] + 1;
                    queue.push(neighbor_idx);
                }
            }
        }
    }

    dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff94::assign_mmff94_types;
    use chematic_smiles::parse;

    #[test]
    fn test_electrostatic_matrix_creation() {
        let mol = parse("CC").unwrap();
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let matrix = ElectrostaticMatrix::new(&mol, &mmff94_types);
        assert_eq!(matrix.scaling.len(), 2);
        assert_eq!(matrix.scaling[0].len(), 2);
    }

    #[test]
    fn test_electrostatic_matrix_1_2_exclusion() {
        let mol = parse("CC").unwrap();
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let matrix = ElectrostaticMatrix::new(&mol, &mmff94_types);
        // 1-2 pairs (bonded) should be excluded
        assert_eq!(matrix.get_scaling(0, 1), 0.0);
        assert_eq!(matrix.get_scaling(1, 0), 0.0);
    }

    #[test]
    fn test_batch_properties_empty() {
        let batch = MMFF94BatchProperties::new();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_batch_properties_average() {
        let mut batch = MMFF94BatchProperties::new();

        let mol = parse("C").unwrap();
        let mmff94_types = assign_mmff94_types(&mol).unwrap();
        let matrix = ElectrostaticMatrix::new(&mol, &mmff94_types);

        // Add two identical molecules with different charges
        batch.push(vec![0.1], vec![], matrix.clone());
        batch.push(vec![0.3], vec![], matrix);

        let avg = batch.average_charges();
        assert_eq!(avg.len(), 1);
        assert!((avg[0] - 0.2).abs() < 1e-6); // (0.1 + 0.3) / 2
    }

    #[test]
    fn test_batch_properties_variance() {
        let mut batch = MMFF94BatchProperties::new();

        let mol = parse("C").unwrap();
        let mmff94_types = assign_mmff94_types(&mol).unwrap();
        let matrix = ElectrostaticMatrix::new(&mol, &mmff94_types);

        batch.push(vec![0.0], vec![], matrix.clone());
        batch.push(vec![1.0], vec![], matrix);

        let var = batch.charge_variance();
        assert_eq!(var.len(), 1);
        // Variance of [0.0, 1.0] around mean 0.5
        // σ = sqrt(((0-0.5)² + (1-0.5)²) / 2) = sqrt(0.25) = 0.5
        assert!((var[0] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_topological_distance_ethane() {
        let mol = parse("CC").unwrap();
        let dist = topological_distance_matrix(&mol);

        assert_eq!(dist[0][0], 0);
        assert_eq!(dist[1][1], 0);
        assert_eq!(dist[0][1], 1);
        assert_eq!(dist[1][0], 1);
    }

    #[test]
    fn test_topological_distance_propane() {
        let mol = parse("CCC").unwrap();
        let dist = topological_distance_matrix(&mol);

        // C-C-C structure
        assert_eq!(dist[0][1], 1); // 1-2
        assert_eq!(dist[1][2], 1); // 1-2
        assert_eq!(dist[0][2], 2); // 1-3 (through one atom)
    }

    #[test]
    fn test_electrostatic_matrix_scaling_values() {
        let mol = parse("CCCC").unwrap();
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let matrix = ElectrostaticMatrix::new(&mol, &mmff94_types);

        // Check some expected scaling patterns
        for i in 0..4 {
            // Diagonal should be zero
            assert_eq!(matrix.get_scaling(i, i), 0.0);
        }

        // All other scaling should be non-negative
        for i in 0..4 {
            for j in 0..4 {
                assert!(matrix.get_scaling(i, j) >= 0.0);
            }
        }
    }
}
