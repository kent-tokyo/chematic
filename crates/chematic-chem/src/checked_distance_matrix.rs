//! Bounded native heavy-atom graph distances, motivated by RDKit #9647.
//!
//! This operates directly on the native graph without sanitizing or removing
//! atoms. Every atom except hydrogen (including wildcard atoms) is heavy.

use chematic_core::{AtomIdx, Molecule};

/// Hard input cap, checked before traversing the graph or allocating buffers.
pub const MAX_NATIVE_DISTANCE_ATOMS: usize = 1024;
/// Hard bond cap, including bonds incident to hydrogen atoms.
pub const MAX_NATIVE_DISTANCE_BONDS: usize = 16_384;
/// Hard dense-output cap (4 MiB of u32 payload, excluding row headers).
pub const MAX_NATIVE_DISTANCE_ENTRIES: usize = 1_048_576;

/// Caller budgets can tighten, but cannot raise, the hard limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NativeDistanceMatrixLimits {
    pub max_atoms: usize,
    pub max_bonds: usize,
    pub max_entries: usize,
}

impl Default for NativeDistanceMatrixLimits {
    fn default() -> Self {
        Self {
            max_atoms: MAX_NATIVE_DISTANCE_ATOMS,
            max_bonds: MAX_NATIVE_DISTANCE_BONDS,
            max_entries: MAX_NATIVE_DISTANCE_ENTRIES,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NativeDistanceMatrixError {
    ResourceLimit {
        resource: &'static str,
        requested: usize,
        limit: usize,
    },
    DimensionOverflow,
    Allocation {
        buffer: &'static str,
    },
}

impl std::fmt::Display for NativeDistanceMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "native distance matrix resource limit: {requested} {resource} exceeds {limit}"
            ),
            Self::DimensionOverflow => write!(f, "native distance matrix dimension overflow"),
            Self::Allocation { buffer } => {
                write!(f, "native distance matrix allocation failed: {buffer}")
            }
        }
    }
}

impl std::error::Error for NativeDistanceMatrixError {}

fn check_limit(
    resource: &'static str,
    requested: usize,
    limit: usize,
) -> Result<(), NativeDistanceMatrixError> {
    if requested > limit {
        return Err(NativeDistanceMatrixError::ResourceLimit {
            resource,
            requested,
            limit,
        });
    }
    Ok(())
}

fn reserved<T>(count: usize, buffer: &'static str) -> Result<Vec<T>, NativeDistanceMatrixError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| NativeDistanceMatrixError::Allocation { buffer })?;
    Ok(result)
}

/// Checked version of [`crate::topological_distance_matrix`].
///
/// Returns unweighted shortest-path lengths in the induced non-hydrogen
/// subgraph. Rows follow original atom order with H omitted; disconnected
/// entries are `u32::MAX`. Empty and hydrogen-only graphs return an empty matrix.
/// No SMILES provenance, sanitization, coordinates or derived caches are needed.
pub fn try_topological_distance_matrix(
    mol: &Molecule,
) -> Result<Vec<Vec<u32>>, NativeDistanceMatrixError> {
    try_topological_distance_matrix_with_limits(mol, &NativeDistanceMatrixLimits::default())
}

/// As [`try_topological_distance_matrix`], with tighter caller budgets.
///
/// All input atoms and bonds count towards the preflight caps, even when H
/// exclusion would make the output smaller. Dense dimensions use checked
/// multiplication and every native heap buffer uses fallible reservation.
/// BFS costs O(H * (V + E)); no partial matrix is returned on failure.
pub fn try_topological_distance_matrix_with_limits(
    mol: &Molecule,
    limits: &NativeDistanceMatrixLimits,
) -> Result<Vec<Vec<u32>>, NativeDistanceMatrixError> {
    check_limit(
        "atoms",
        mol.atom_count(),
        limits.max_atoms.min(MAX_NATIVE_DISTANCE_ATOMS),
    )?;
    check_limit(
        "bonds",
        mol.bond_count(),
        limits.max_bonds.min(MAX_NATIVE_DISTANCE_BONDS),
    )?;
    let n = mol
        .atoms()
        .filter(|(_, atom)| atom.element.atomic_number() != 1)
        .count();
    let entries = n
        .checked_mul(n)
        .ok_or(NativeDistanceMatrixError::DimensionOverflow)?;
    check_limit(
        "entries",
        entries,
        limits.max_entries.min(MAX_NATIVE_DISTANCE_ENTRIES),
    )?;
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut positions = reserved(mol.atom_count(), "atom positions")?;
    positions.resize(mol.atom_count(), usize::MAX);
    let mut heavy = reserved(n, "heavy atom indices")?;
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 1 {
            positions[idx.0 as usize] = heavy.len();
            heavy.push(idx);
        }
    }
    let mut matrix = reserved(n, "matrix rows")?;
    let mut queue = reserved(n, "BFS queue")?;
    for (start_pos, &start) in heavy.iter().enumerate() {
        let mut row = reserved(n, "distance row")?;
        row.resize(n, u32::MAX);
        fill_row(mol, &positions, start, start_pos, &mut row, &mut queue);
        matrix.push(row);
    }
    Ok(matrix)
}

fn fill_row(
    mol: &Molecule,
    positions: &[usize],
    start: AtomIdx,
    start_pos: usize,
    row: &mut [u32],
    queue: &mut Vec<AtomIdx>,
) {
    queue.clear();
    row[start_pos] = 0;
    queue.push(start);
    let mut head = 0;
    while head < queue.len() {
        let current = queue[head];
        head += 1;
        let next_distance = row[positions[current.0 as usize]] + 1;
        for (neighbor, _) in mol.neighbors(current) {
            let pos = positions[neighbor.0 as usize];
            if pos != usize::MAX && row[pos] == u32::MAX {
                row[pos] = next_distance;
                queue.push(neighbor);
            }
        }
    }
}
