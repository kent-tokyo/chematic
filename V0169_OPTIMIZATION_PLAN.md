# v0.1.69: SSSR Caching Optimization (Design Document)

## Objective
Eliminate redundant `find_sssr(mol)` calls (11 locations) by implementing batch ring descriptor computation.

## Current State (v0.1.68)
- `find_sssr()` called independently in 11 functions
- Each call recomputes the same SSSR graph structure
- `molecule_report()` in workflow.rs invokes 3 ring descriptors sequentially

## Proposed Architecture

### 1. Batch Computation Helper (descriptors.rs)
```rust
/// Cache computed ring metrics to avoid redundant SSSR calls.
pub(crate) struct RingMetrics {
    pub ring_count: usize,
    pub aromatic_ring_count: usize,
    pub spiro_atoms: usize,
    pub bridgehead_atoms: usize,
}

pub(crate) fn compute_ring_metrics(mol: &Molecule) -> RingMetrics {
    let sssr = find_sssr(mol);
    // Compute all ring-dependent metrics once
    let ring_count = sssr.rings().len();
    let aromatic_ring_count = /* ... */;
    let spiro_atoms = /* ... */;
    let bridgehead_atoms = /* ... */;
    
    RingMetrics { ring_count, aromatic_ring_count, spiro_atoms, bridgehead_atoms }
}
```

### 2. Integration Points (workflow.rs)
Replace:
```rust
ring_count: ring_count(mol),
aromatic_ring_count: aromatic_ring_count(mol),
num_spiro_atoms: num_spiro_atoms(mol),
num_bridgehead_atoms: num_bridgehead_atoms(mol),
```

With:
```rust
let metrics = compute_ring_metrics(mol);
// Use metrics.ring_count, metrics.aromatic_ring_count, etc.
```

## Performance Impact
- **Current**: 4 × find_sssr() calls per `molecule_report()`
- **Optimized**: 1 × find_sssr() call per `molecule_report()`
- **Savings**: ~75% reduction in ring algorithm overhead for batch operations

## Implementation Notes
- Public ring descriptor functions (`ring_count()`, etc.) remain unchanged for backward compatibility
- `compute_ring_metrics()` marked as `pub(crate)` for internal workflow use
- `RingMetrics` struct is private; no API surface change
- No breaking changes; optimization is internal

## Future Enhancements
- Cache `RingMetrics` in molecule objects if molecule lifespan permits
- Extend batch computation to include `num_aliphatic_rings()`, `num_saturated_rings()`, etc.
- Measure actual performance gain on large batch screening tasks

## Scope
Low risk, isolated to descriptors.rs and workflow.rs. No public API changes.
