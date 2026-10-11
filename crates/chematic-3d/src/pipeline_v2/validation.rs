fn pipeline_config_is_invalid(config: &PipelineV2Config) -> bool {
    config.embed.materialize_implicit_h_for_chirality
        || (config.expand_implicit_h_through_pipeline && !config.embed.enforce_chirality)
}

/// True only when a declared tetrahedral center uses an implicit H. E/Z
/// declarations do not need hydrogen coordinates, so keeping E/Z-only
/// molecules on the heavy-atom path avoids the all-H distance-geometry/MMFF94
/// cost without weakening tetrahedral geometry. A narrower ring-only rule was
/// measured and rejected because four ordinary implicit-H centers developed
/// gross clashes on the 265-molecule gate.
fn mol_needs_expanded_stereo_geometry(mol: &Molecule) -> bool {
    (0..mol.atom_count()).any(|i| {
        let center = AtomIdx(i as u32);
        let atom = mol.atom(center);
        atom.chirality.is_tetrahedral() && mol.implicit_hydrogen_count(center) > 0
    })
}

/// Multi-ring junction stereocenters retain the full explicit-H treatment.
/// Their coupled ring geometry is the class for which materializing only the
/// stereocenter H can leave heavy-only final verification inconsistent with
/// the minimized internal graph. Ordinary implicit-H centers use the much
/// smaller selective expansion.
fn mol_needs_full_hydrogen_expansion(mol: &Molecule) -> bool {
    let rings = chematic_perception::find_sssr(mol);
    (0..mol.atom_count()).any(|i| {
        let center = AtomIdx(i as u32);
        let atom = mol.atom(center);
        atom.chirality.is_tetrahedral()
            && mol.implicit_hydrogen_count(center) > 0
            && rings
                .rings()
                .iter()
                .filter(|ring| ring.contains(&center))
                .take(2)
                .count()
                >= 2
    })
}

/// Whether a torsion-knowledge potential can be applied to the molecule's
/// actual geometry. Keep this predicate shared by the evidence report and
/// the fail-closed policy gate so they cannot drift apart.
fn ring_torsion_is_applicable(mol: &Molecule, potential: &TorsionPotential) -> bool {
    let (b, c) = potential.central_bond;
    mol.bond_between(b, c).is_some() && is_bridge_bond(mol, b, c)
}

fn torsion_knowledge_config(config: &PipelineV2Config) -> TorsionKnowledgeConfig {
    TorsionKnowledgeConfig {
        use_exp_torsions: config.embed.use_exp_torsions,
        use_small_ring_torsions: config.embed.use_small_ring_torsions,
        use_macrocycle_torsions: config.embed.use_macrocycle_torsions,
        use_macrocycle_14_bounds: config.embed.use_macrocycle_14_bounds,
        include_legacy_heuristic: config.include_legacy_torsion_heuristic,
    }
}

/// Verify the geometry using the coordinate view exposed by the pipeline contract.
///
/// Embedding may temporarily operate on a molecule with materialized implicit
/// hydrogens, while returned coordinates are always sized for the original molecule.
/// Keeping that rule here prevents final-stereo checks from validating a different
/// coordinate view than the one returned to the caller.
fn verify_authoritative_final_stereo(
    original_mol: &Molecule,
    internal_mol: &Molecule,
    coords: &Coords3D,
    use_expanded_geometry: bool,
    original_atom_count: usize,
) -> StereoVerification {
    if use_expanded_geometry {
        let truncated = truncate_coords(coords, original_atom_count);
        verify_stereo(original_mol, &truncated)
    } else {
        verify_stereo(internal_mol, coords)
    }
}

fn worst_bond_length(mol: &Molecule, coords: &Coords3D) -> f64 {
    mol.bonds()
        .map(|(_, bond)| coords.get(bond.atom1).distance(&coords.get(bond.atom2)))
        .fold(0.0_f64, f64::max)
}

/// Fraction of real bonds whose length deviates from `ideal_bond_length` by more than
/// `threshold` (e.g. `0.15` for the 15% rate, `0.50` for the 50% rate).
fn bond_violation_rate(mol: &Molecule, coords: &Coords3D, threshold: f64) -> f64 {
    let mut total = 0usize;
    let mut violations = 0usize;
    for (_, bond) in mol.bonds() {
        let ideal = ideal_bond_length(mol, bond.atom1, bond.atom2);
        if ideal <= 0.0 || !ideal.is_finite() {
            continue;
        }
        let actual = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
        total += 1;
        if ((actual - ideal).abs() / ideal) > threshold {
            violations += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        violations as f64 / total as f64
    }
}

fn gross_clash_count(mol: &Molecule, coords: &Coords3D) -> usize {
    let n = coords.atom_count();
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);
            if mol.bond_between(a, b).is_some() {
                continue;
            }
            if coords.get(a).distance(&coords.get(b)) < NONBONDED_CLASH_THRESHOLD_ANGSTROM {
                count += 1;
            }
        }
    }
    count
}

fn compute_final_validation(
    mol: &Molecule,
    coords: &Coords3D,
    final_stereo: &StereoVerification,
    torsion_energy_after: f64,
    ring_closure_delta: f64,
) -> FinalGeometryValidation {
    let all_finite = coords.is_finite();
    let atom_count_unchanged = coords.atom_count() == mol.atom_count();
    let worst_bond_length_v = worst_bond_length(mol, coords);
    let bounds_conformance_v = bounds_conformance(mol, coords);

    let sound = all_finite && atom_count_unchanged && worst_bond_length_v <= MAX_SANE_BOND_LENGTH;

    FinalGeometryValidation {
        all_finite,
        atom_count_unchanged,
        worst_bond_length: worst_bond_length_v,
        bond_violation_rate_15pct: bond_violation_rate(mol, coords, 0.15),
        bond_violation_rate_50pct: bond_violation_rate(mol, coords, 0.50),
        gross_clash_count: gross_clash_count(mol, coords),
        bounds_conformance: bounds_conformance_v,
        stereo_ok: final_stereo.is_fully_satisfied(),
        torsion_energy_after,
        ring_closure_delta,
        sound,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

