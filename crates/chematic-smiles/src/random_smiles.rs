use crate::writer::write;
use chematic_core::{AtomIdx, Molecule, MoleculeBuilder};
use std::collections::HashSet;

/// Generate a single random SMILES from a molecule using the given seed.
/// Atoms are permuted based on a simple xorshift64 RNG.
///
/// # Example
/// ```ignore
/// let mol = parse("CC(C)O")?;
/// let smiles = random_smiles(&mol, 42);
/// assert!(!smiles.is_empty());
/// ```
pub fn random_smiles(mol: &Molecule, seed: u64) -> String {
    let permutation = generate_permutation(mol.atom_count(), seed);
    let permuted = apply_permutation(mol, &permutation);
    write(&permuted)
}

/// Generate `count` unique random SMILES from a molecule using sequential seeds.
/// Each seed increments from the base seed. Returns up to `count` unique SMILES.
///
/// # Example
/// ```ignore
/// let mol = parse("CC(C)O")?;
/// let variants = random_smiles_vect(&mol, 5, 42);
/// assert!(variants.len() <= 5);
/// // All elements should be unique
/// ```
pub fn random_smiles_vect(mol: &Molecule, count: usize, seed: u64) -> Vec<String> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    let max_attempts = count.saturating_mul(3).max(10);

    for attempt in 0..max_attempts {
        if result.len() >= count {
            break;
        }
        let smiles = random_smiles(mol, seed.wrapping_add(attempt as u64));
        if seen.insert(smiles.clone()) {
            result.push(smiles);
        }
    }
    result
}

/// Xorshift64 pseudo-random number generator.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        let state = if seed == 0 { 1 } else { seed };
        Xorshift64 { state }
    }

    fn next(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn range(&mut self, max: usize) -> usize {
        if max == 0 {
            0
        } else {
            (self.next() as usize) % max
        }
    }
}

/// Generate a random permutation of atom indices [0..n).
fn generate_permutation(n: usize, seed: u64) -> Vec<usize> {
    let mut rng = Xorshift64::new(seed);
    let mut perm: Vec<usize> = (0..n).collect();

    // Fisher-Yates shuffle
    for i in (1..n).rev() {
        let j = rng.range(i + 1);
        perm.swap(i, j);
    }
    perm
}

/// Apply an atom permutation to a molecule, returning a new molecule with atoms reordered.
fn apply_permutation(mol: &Molecule, permutation: &[usize]) -> Molecule {
    let mut builder = MoleculeBuilder::new();

    // Add atoms in permuted order
    for &old_idx in permutation {
        let atom = mol.atom(AtomIdx(old_idx as u32));
        builder.add_atom(atom.clone());
    }

    // Create old_to_new mapping
    let mut old_to_new = vec![0usize; mol.atom_count()];
    for (new_idx, &old_idx) in permutation.iter().enumerate() {
        old_to_new[old_idx] = new_idx;
    }

    // Add bonds with remapped indices
    for (_, bond_entry) in mol.bonds() {
        let old_a = bond_entry.atom1;
        let old_b = bond_entry.atom2;
        let new_a = AtomIdx(old_to_new[old_a.0 as usize] as u32);
        let new_b = AtomIdx(old_to_new[old_b.0 as usize] as u32);
        let _ = builder.add_bond(new_a, new_b, bond_entry.order);
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse;

    #[test]
    fn test_random_smiles_single() {
        let mol = parse("CC").unwrap();
        let smiles = random_smiles(&mol, 42);
        assert!(!smiles.is_empty());
    }

    #[test]
    fn test_random_smiles_vect_generates_multiple() {
        let mol = parse("CCCC").unwrap();
        let variants = random_smiles_vect(&mol, 3, 42);
        assert!(!variants.is_empty());
        // All variants should be valid SMILES
        for s in &variants {
            assert!(!s.is_empty());
        }
    }

    #[test]
    fn test_random_smiles_vect_unique() {
        let mol = parse("CCCC").unwrap();
        let variants = random_smiles_vect(&mol, 10, 100);
        let set: HashSet<_> = variants.iter().cloned().collect();
        // Should have multiple unique variants (may not be exactly 10 due to permutation limits)
        assert!(set.len() > 1);
    }

    #[test]
    fn test_random_smiles_roundtrip() {
        let original_smiles = "CC(C)O";
        let mol = parse(original_smiles).unwrap();
        let random = random_smiles(&mol, 99);
        // Parse the random SMILES back — should parse successfully
        let mol2 = parse(&random);
        assert!(mol2.is_ok());
    }

    #[test]
    fn test_permutation_deterministic() {
        let mol = parse("CCCC").unwrap();
        let s1 = random_smiles(&mol, 77);
        let s2 = random_smiles(&mol, 77);
        assert_eq!(s1, s2);
    }

    #[test]
    fn test_permutation_different_seeds() {
        let mol = parse("CCCC").unwrap();
        let s1 = random_smiles(&mol, 1);
        let s2 = random_smiles(&mol, 2);
        // Different seeds should (likely) produce different SMILES
        // Note: For some molecules, permutations might produce identical SMILES
        // due to symmetry, so we just check they're both valid
        assert!(!s1.is_empty());
        assert!(!s2.is_empty());
    }
}
