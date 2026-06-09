use crate::sssr::RingSet;
use chematic_core::{AtomIdx, Molecule};

/// Classification of a ring system's topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RingSystemKind {
    /// Single ring with no fusion to other rings.
    Simple,
    /// Two or more rings sharing one or more bonds (adjacent atoms).
    Fused,
    /// Two or more rings sharing exactly one atom (spiro junction).
    Spiro,
    /// Rings connected via bridgehead atoms (share 2+ non-adjacent atoms).
    Bridged,
}

/// A connected group of SSSR rings forming a ring system.
#[derive(Debug, Clone)]
pub struct RingFamily {
    /// Indices into the SSSR's ring list.
    pub ring_indices: Vec<usize>,
    /// All distinct atoms in this ring system.
    pub atoms: Vec<AtomIdx>,
    /// Classification of the ring system's topology.
    pub kind: RingSystemKind,
}

/// Find all ring systems (families of connected rings) and classify their topology.
///
/// Two SSSR rings belong to the same family if they share at least one atom.
/// The family is then classified as Simple, Fused, Spiro, or Bridged based on
/// how the rings are connected.
///
/// # Example
/// ```ignore
/// let mol = parse("C1CC2CCC1CC2").unwrap(); // norbornane (bridged bicyclic)
/// let sssr = find_sssr(&mol);
/// let families = find_ring_families(&sssr);
/// assert_eq!(families.len(), 1);
/// assert_eq!(families[0].kind, RingSystemKind::Bridged);
/// ```
pub fn find_ring_families(mol: &Molecule, sssr: &RingSet) -> Vec<RingFamily> {
    if sssr.ring_count() == 0 {
        return vec![];
    }

    let rings = sssr.rings();

    // Union-Find: group rings that share at least one atom
    let mut parent: Vec<usize> = (0..rings.len()).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], x: usize, y: usize) {
        let px = find(parent, x);
        let py = find(parent, y);
        if px != py {
            parent[px] = py;
        }
    }

    // Union rings that share atoms
    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            if rings[i].iter().any(|a| rings[j].contains(a)) {
                union(&mut parent, i, j);
            }
        }
    }

    // Group rings by their root
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for i in 0..rings.len() {
        let root = find(&mut parent, i);
        groups.entry(root).or_insert_with(Vec::new).push(i);
    }

    // Build ring families
    let mut families = Vec::new();
    for ring_indices in groups.into_values() {
        // Collect all distinct atoms in this family
        let mut atoms = std::collections::HashSet::new();
        for &idx in &ring_indices {
            for &atom in &rings[idx] {
                atoms.insert(atom);
            }
        }
        let mut atoms: Vec<AtomIdx> = atoms.into_iter().collect();
        atoms.sort_by_key(|a| a.0);

        // Classify the ring system
        let kind = classify_ring_system(mol, &ring_indices, rings);

        families.push(RingFamily {
            ring_indices,
            atoms,
            kind,
        });
    }

    // Sort by first ring index for deterministic output
    families.sort_by_key(|f| f.ring_indices.iter().min().copied().unwrap_or(0));
    families
}

/// Classify the topology of a ring system given its constituent SSSR rings.
fn classify_ring_system(mol: &Molecule, ring_indices: &[usize], rings: &[Vec<AtomIdx>]) -> RingSystemKind {
    // Single ring is always Simple
    if ring_indices.len() == 1 {
        return RingSystemKind::Simple;
    }

    let mut worst_kind = RingSystemKind::Spiro;

    // Check all pairs of rings
    for i in 0..ring_indices.len() {
        for j in (i + 1)..ring_indices.len() {
            let ring_i = &rings[ring_indices[i]];
            let ring_j = &rings[ring_indices[j]];

            // Count shared atoms
            let shared: Vec<AtomIdx> = ring_i
                .iter()
                .filter(|a| ring_j.contains(a))
                .copied()
                .collect();

            if shared.is_empty() {
                continue; // should not happen if union-find worked correctly
            }

            let kind = if shared.len() == 1 {
                RingSystemKind::Spiro
            } else if is_fused_ring(mol, &shared) {
                RingSystemKind::Fused
            } else {
                RingSystemKind::Bridged
            };

            // Keep the most complex classification
            worst_kind = worst_classification(worst_kind, kind);
        }
    }

    worst_kind
}

/// Check if shared atoms form a fused ring vs. bridged.
/// Fused: shared atoms form a closed cycle (ring). The last shared atom bonds back to the first.
/// Bridged: shared atoms form an open path. The last atom does NOT bond back to the first.
fn is_fused_ring(mol: &Molecule, atoms: &[AtomIdx]) -> bool {
    if atoms.len() <= 1 {
        return true; // Spiro: single shared atom
    }

    if atoms.len() == 2 {
        // Two atoms: bond between them (fused edge) vs. no bond (bridged)
        return mol.bond_between(atoms[0], atoms[1]).is_some();
    }

    // For 3+ atoms: check if they form a cycle (fused) or path (bridged)
    // If the last atom is bonded back to the first, it's a cycle (fused ring junction)
    // Otherwise it's a path (bridged ring junction)
    mol.bond_between(atoms[atoms.len() - 1], atoms[0]).is_some()
}

/// Return the "worst" (most complex) ring system classification.
fn worst_classification(a: RingSystemKind, b: RingSystemKind) -> RingSystemKind {
    use RingSystemKind::*;
    match (a, b) {
        (Bridged, _) | (_, Bridged) => Bridged,
        (Fused, _) | (_, Fused) => Fused,
        _ => Spiro,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_ring_families_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Simple);
        assert_eq!(families[0].ring_indices, vec![0]);
    }

    #[test]
    fn test_ring_families_naphthalene() {
        let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Fused);
        assert_eq!(families[0].ring_indices.len(), 2);
    }

    #[test]
    fn test_ring_families_biphenyl() {
        let mol = parse("c1ccccc1c1ccccc1").unwrap(); // biphenyl (two benzene rings, no connection)
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 2);
        assert!(families.iter().all(|f| f.kind == RingSystemKind::Simple));
    }

    #[test]
    fn test_ring_families_spiro() {
        let mol = parse("C1CCC2(C1)CCC2").unwrap(); // spiro[4.4]nonane
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Spiro);
    }

    #[test]
    fn test_ring_families_norbornane() {
        let mol = parse("C1CC2CCC1CC2").unwrap(); // norbornane (bridged)
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Bridged);
        assert_eq!(families[0].ring_indices.len(), 2);
    }

    #[test]
    fn test_ring_families_adamantane() {
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Bridged);
    }

    #[test]
    fn test_ring_families_caffeine() {
        // caffeine has fused rings
        let mol = parse("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 1);
        assert_eq!(families[0].kind, RingSystemKind::Fused);
    }

    #[test]
    fn test_ring_families_acyclic() {
        let mol = parse("CCCC").unwrap();
        let sssr = crate::sssr::find_sssr(&mol);
        let families = find_ring_families(&mol, &sssr);
        assert_eq!(families.len(), 0);
    }
}
