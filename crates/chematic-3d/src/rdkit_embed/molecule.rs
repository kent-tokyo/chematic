//! `DGeomHelpers::EmbedMultipleConfs` for one conformer of an
//! explicit-hydrogen molecule (`AllChem.EmbedMolecule(AddHs(m), ...)`).

use chematic_core::Molecule;
use chematic_smiles::RdkitMolView;

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
use rayon::prelude::*;

use super::bounds::{
    BoundsInvariant, ForcedTorsion, collect_bonds_and_angles, init_bounds_mat, set_topol_bounds,
    set_topol_bounds_with_forced_torsions,
};
use super::dg::{BoundsMatrix, ChiralSet, triangle_smooth_bounds};
use super::embed::{EmbedArgs, EmbedError, EmbedParams, embed_initial_points, embed_points};
use super::etk::{
    CrystalFfDetails, get_experimental_torsions_with_forced, minimize_with_exp_torsions,
};

/// The `EmbedParameters` fields this port honours. [`Default`] is what
/// `AllChem.EmbedMolecule(mol, randomSeed=42)` uses (ETKDGv3 settings).
#[derive(Clone, Debug)]
pub struct RdkitEmbedOptions {
    pub random_seed: i32,
    pub max_iterations: u32,
    pub use_exp_torsion_angle_prefs: bool,
    pub use_basic_knowledge: bool,
    pub et_version: u32,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_macrocycle14config: bool,
    pub enforce_chirality: bool,
    pub force_trans_amides: bool,
    pub ignore_smoothing_failures: bool,
}

impl Default for RdkitEmbedOptions {
    fn default() -> Self {
        RdkitEmbedOptions {
            random_seed: 42,
            max_iterations: 0,
            use_exp_torsion_angle_prefs: true,
            use_basic_knowledge: true,
            et_version: 2,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: true,
            use_macrocycle14config: true,
            enforce_chirality: true,
            force_trans_amides: true,
            ignore_smoothing_failures: false,
        }
    }
}

/// Why no coordinates were produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitEmbedError {
    /// The molecule is outside what the port models.
    Unsupported(String),
    /// RDKit raises an exception for this molecule.
    RdkitException(String),
    /// RDKit returns -1 (no conformer): smoothing or embedding failed.
    Failed,
}

impl core::fmt::Display for RdkitEmbedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Unsupported(w) => write!(f, "RDKit embedding: unsupported input: {w}"),
            Self::RdkitException(w) => write!(f, "RDKit embedding: RDKit raises: {w}"),
            Self::Failed => write!(f, "RDKit embedding failed (EmbedMolecule returns -1)"),
        }
    }
}

impl std::error::Error for RdkitEmbedError {}

/// Failure from the combined RDKit-compatible ETKDG + MMFF94 ensemble.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitEmbedMmffError {
    Embed(RdkitEmbedError),
    ForceField(String),
}

impl core::fmt::Display for RdkitEmbedMmffError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Embed(error) => error.fmt(f),
            Self::ForceField(reason) => write!(f, "RDKit-compatible MMFF94: {reason}"),
        }
    }
}

impl std::error::Error for RdkitEmbedMmffError {}

impl From<RdkitEmbedError> for RdkitEmbedMmffError {
    fn from(error: RdkitEmbedError) -> Self {
        Self::Embed(error)
    }
}

/// One ETKDG conformer after RDKit-compatible MMFF94 optimization.
#[derive(Clone, Debug)]
pub struct RdkitMmffConformer {
    pub coords: Vec<[f64; 3]>,
    /// RDKit status: 0 converged, 1 iteration limit, negative on failure.
    pub status: i32,
    pub energy: f64,
}

/// One optimized conformer plus diagnostics already computed by BFGS.
#[derive(Clone, Debug)]
pub struct RdkitMmffConformerDiagnostics {
    pub conformer: RdkitMmffConformer,
    pub iterations: u32,
    pub final_gradient_test: Option<f64>,
}

impl From<BoundsInvariant> for RdkitEmbedError {
    fn from(e: BoundsInvariant) -> Self {
        RdkitEmbedError::RdkitException(e.0.to_string())
    }
}

impl From<EmbedError> for RdkitEmbedError {
    fn from(e: EmbedError) -> Self {
        RdkitEmbedError::RdkitException(format!("{e:?}"))
    }
}

/// `EmbeddingOps::findChiralSets` (no coordMap, no atropisomers).
fn find_chiral_sets(v: &RdkitMolView) -> (Vec<ChiralSet>, Vec<ChiralSet>) {
    let mut chiral = Vec::new();
    let mut tet = Vec::new();
    for a in 0..v.num_atoms() {
        let at = &v.atoms[a];
        if at.atomic_num == 1 {
            continue;
        }
        let tagged = at.chiral_tag == 1 || at.chiral_tag == 2;
        if !(tagged || ((at.atomic_num == 6 || at.atomic_num == 7) && v.degree(a) == 4)) {
            continue;
        }
        let mut nbrs: Vec<usize> = v.neighbors(a).collect();
        let mut vol_lower = 5.0;
        let vol_upper = 100.0;
        if nbrs.len() < 4 {
            vol_lower = 2.0;
            nbrs.push(a);
        }
        let n_small = v.atom_ring_sizes(a).iter().filter(|&&s| s < 5).count();
        let fused = n_small > 1;
        let idx = [a, nbrs[0], nbrs[1], nbrs[2], nbrs[3]];
        if at.chiral_tag == 2 {
            chiral.push(ChiralSet {
                idx,
                vol_lower,
                vol_upper,
                in_fused_small_rings: fused,
            });
        } else if at.chiral_tag == 1 {
            chiral.push(ChiralSet {
                idx,
                vol_lower: -vol_upper,
                vol_upper: -vol_lower,
                in_fused_small_rings: fused,
            });
        } else if v.num_atom_rings(a) < 2 || v.is_atom_in_ring_of_size(a, 3) {
            // only ring atoms in two or more rings, not 3-rings
        } else {
            tet.push(ChiralSet {
                idx,
                vol_lower: 0.0,
                vol_upper: 0.0,
                in_fused_small_rings: fused,
            });
        }
    }
    (chiral, tet)
}

type DoubleBonds = (Vec<(usize, usize, usize)>, Vec<([usize; 4], i32)>);

/// `EmbeddingOps::findDoubleBonds` (no coordMap).
fn find_double_bonds(v: &RdkitMolView) -> DoubleBonds {
    let mut ends = Vec::new();
    let mut stereo = Vec::new();
    for (bi, b) in v.bonds.iter().enumerate() {
        if b.bond_type != 2 {
            continue;
        }
        for atm in [b.begin, b.end] {
            if v.degree(atm) < 2 {
                continue;
            }
            let oatm = v.other_atom(bi, atm);
            for nbr in v.neighbors(atm) {
                if nbr == oatm {
                    continue;
                }
                let obnd = v.bond_between(atm, nbr).expect("bond");
                if v.bonds[obnd].bond_type != 1 && v.degree(atm) == 2 {
                    continue;
                }
                ends.push((nbr, atm, oatm));
            }
        }
        if b.stereo > 1 {
            // Z (2) is cis-like; E (3) trans-like.
            let sign = if b.stereo == 2 { -1 } else { 1 };
            stereo.push(([b.stereo_atoms[0], b.begin, b.end, b.stereo_atoms[1]], sign));
        }
    }
    (ends, stereo)
}

/// `EmbeddingOps::setupInitialBoundsMatrix` (no coordMap).
fn setup_initial_bounds(
    v: &RdkitMolView,
    labels: &[String],
    o: &RdkitEmbedOptions,
) -> Result<Option<(BoundsMatrix, Vec<ForcedTorsion>)>, RdkitEmbedError> {
    let n = v.num_atoms();
    let mut m = BoundsMatrix::new(n);
    init_bounds_mat(&mut m);
    let forced_torsions = set_topol_bounds_with_forced_torsions(
        v,
        labels,
        &mut m,
        true,
        o.use_macrocycle14config,
        o.force_trans_amides,
    )?;
    if !triangle_smooth_bounds(&mut m, 0.0) {
        init_bounds_mat(&mut m);
        set_topol_bounds(
            v,
            labels,
            &mut m,
            false,
            o.use_macrocycle14config,
            o.force_trans_amides,
        )?;
        if !triangle_smooth_bounds(&mut m, 0.0) {
            if o.ignore_smoothing_failures {
                init_bounds_mat(&mut m);
                set_topol_bounds(
                    v,
                    labels,
                    &mut m,
                    false,
                    o.use_macrocycle14config,
                    o.force_trans_amides,
                )?;
            } else {
                return Ok(None);
            }
        }
    }
    Ok(Some((m, forced_torsions)))
}

/// `AllChem.EmbedMolecule(mh, randomSeed=..., ...)` on an explicit-hydrogen
/// molecule from `add_hydrogens` (RDKit's `AddHs` atom order): all-atom
/// coordinates, or the reason there are none.
pub fn rdkit_embed_molecule(
    mh: &Molecule,
    o: &RdkitEmbedOptions,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
    validate_options(mh, o)?;
    embed_prepared_fragments(&prepare_fragments(mh, o)?, o, o.random_seed)
}

/// RDKit 2026.09.1's `onlyInitialEmbedding=True` coordinates: the seeded
/// metric-matrix embedding before distance-geometry or ETKDG minimization.
pub fn rdkit_embed_initial_molecule(
    mh: &Molecule,
    o: &RdkitEmbedOptions,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
    validate_options(mh, o)?;
    let prepared = prepare_fragments(mh, o)?;
    embed_prepared_fragments_with(&prepared, o, o.random_seed, true)
}

/// `AllChem.EmbedMultipleConfs(mh, numConfs=..., params=...)` on an
/// explicit-hydrogen molecule from `add_hydrogens` (RDKit's `AddHs` atom
/// order).
///
/// Bounds, ring information, and experimental torsion preferences are built
/// once and reused across conformers. Failed conformer attempts are omitted,
/// matching RDKit's returned conformer-id list; input and setup failures remain
/// typed errors.
pub fn rdkit_embed_multiple_confs(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
) -> Result<Vec<Vec<[f64; 3]>>, RdkitEmbedError> {
    validate_options(mh, o)?;
    if num_confs == 0 {
        return Ok(Vec::new());
    }
    let prepared = match prepare_fragments(mh, o) {
        Ok(prepared) => prepared,
        Err(RdkitEmbedError::Failed) => return Ok(Vec::new()),
        Err(error) => return Err(error),
    };
    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    let attempts: Vec<_> = (0..num_confs)
        .into_par_iter()
        .map(|index| embed_conformer_attempt(&prepared, o, index))
        .collect();
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    let attempts: Vec<_> = (0..num_confs)
        .map(|index| embed_conformer_attempt(&prepared, o, index))
        .collect();

    let mut conformers = Vec::with_capacity(num_confs);
    for attempt in attempts {
        if let Some(coords) = attempt? {
            conformers.push(coords);
        }
    }
    Ok(conformers)
}

/// Generate an ETKDG ensemble and optimize every conformer with the
/// RDKit-compatible MMFF94 implementation.
///
/// The immutable MMFF parameterization is built once and conformers are
/// optimized in parallel on native targets. Candidate order and coordinates
/// remain deterministic.
pub fn rdkit_embed_multiple_confs_mmff(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    max_iters: u32,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    let conformers = rdkit_embed_multiple_confs(mh, num_confs, o)?;
    let Some(first) = conformers.first() else {
        return Ok(Vec::new());
    };
    let field = chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, first, 100.0, true)
        .map_err(|error| RdkitEmbedMmffError::ForceField(error.to_string()))?;

    Ok(optimize_mmff_conformers(&field, conformers, max_iters))
}

/// [`rdkit_embed_multiple_confs_mmff`] with per-conformer convergence
/// diagnostics. The coordinates, status, and energy are identical to the
/// compact API; no additional energy or gradient evaluation is performed.
pub fn rdkit_embed_multiple_confs_mmff_diagnostics(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    max_iters: u32,
) -> Result<Vec<RdkitMmffConformerDiagnostics>, RdkitEmbedMmffError> {
    let conformers = rdkit_embed_multiple_confs(mh, num_confs, o)?;
    let Some(first) = conformers.first() else {
        return Ok(Vec::new());
    };
    let field = chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, first, 100.0, true)
        .map_err(|error| RdkitEmbedMmffError::ForceField(error.to_string()))?;

    Ok(optimize_mmff_conformers_detailed(
        &field, conformers, max_iters,
    ))
}

/// Generate and optimize an ETKDG ensemble, then spend an additional MMFF94
/// budget only on the lowest-energy candidates that did not converge.
///
/// The initial pass is identical to [`rdkit_embed_multiple_confs_mmff`]. At
/// most `retry_candidates` finite, non-converged conformers are ordered by
/// `(energy, candidate index)` and continued from their current coordinates
/// for `retry_max_iters` more iterations. A retry is retained only when its
/// energy is finite and no higher than the initial result, so refinement
/// cannot make the returned ensemble worse. Candidate order remains stable.
pub fn rdkit_embed_multiple_confs_mmff_adaptive(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    initial_max_iters: u32,
    retry_max_iters: u32,
    retry_candidates: usize,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    rdkit_embed_multiple_confs_mmff_retry_policy(
        mh,
        num_confs,
        o,
        initial_max_iters,
        retry_max_iters,
        retry_candidates,
        false,
    )
}

/// Generate and optimize an ETKDG ensemble, retrying low-energy candidates
/// only when the initial pass produced no converged conformer.
///
/// This quality-oriented rescue avoids paying the additional MMFF94 budget
/// for molecules that already have at least one converged candidate. The
/// retry ordering and replacement rules are otherwise identical to
/// [`rdkit_embed_multiple_confs_mmff_adaptive`].
pub fn rdkit_embed_multiple_confs_mmff_rescue(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    initial_max_iters: u32,
    retry_max_iters: u32,
    retry_candidates: usize,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    rdkit_embed_multiple_confs_mmff_retry_policy(
        mh,
        num_confs,
        o,
        initial_max_iters,
        retry_max_iters,
        retry_candidates,
        true,
    )
}

/// Generate and optimize an ETKDG ensemble, then retry only high-gradient
/// candidates when no initial candidate converged.
///
/// Unlike the energy-only rescue, this profile uses RDKit's normalized final
/// BFGS gradient test to avoid spending extra work on candidates that are
/// already close to a stationary point. Candidate ranking remains
/// `(energy, candidate index)`, and a retry is retained only if its energy is
/// finite and no higher than the initial result.
#[allow(clippy::too_many_arguments)]
pub fn rdkit_embed_multiple_confs_mmff_gradient_rescue(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    initial_max_iters: u32,
    retry_max_iters: u32,
    retry_candidates: usize,
    minimum_gradient_test: f64,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    if !minimum_gradient_test.is_finite() || minimum_gradient_test < 0.0 {
        return Err(RdkitEmbedMmffError::ForceField(
            "minimum gradient test must be finite and non-negative".into(),
        ));
    }
    let conformers = rdkit_embed_multiple_confs(mh, num_confs, o)?;
    let Some(first) = conformers.first() else {
        return Ok(Vec::new());
    };
    let field = chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, first, 100.0, true)
        .map_err(|error| RdkitEmbedMmffError::ForceField(error.to_string()))?;
    let mut detailed = optimize_mmff_conformers_detailed(&field, conformers, initial_max_iters);

    if retry_candidates == 0
        || retry_max_iters == 0
        || detailed.iter().any(|result| result.conformer.status == 0)
    {
        return Ok(detailed
            .into_iter()
            .map(|result| result.conformer)
            .collect());
    }

    let mut retry_indices: Vec<usize> = detailed
        .iter()
        .enumerate()
        .filter_map(|(index, result)| {
            let conformer = &result.conformer;
            (conformer.status != 0
                && conformer.energy.is_finite()
                && result
                    .final_gradient_test
                    .is_some_and(|test| test >= minimum_gradient_test))
            .then_some(index)
        })
        .collect();
    retry_indices.sort_by(|&left, &right| {
        detailed[left]
            .conformer
            .energy
            .total_cmp(&detailed[right].conformer.energy)
            .then_with(|| left.cmp(&right))
    });
    retry_indices.truncate(retry_candidates);

    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    let retries: Vec<_> = retry_indices
        .par_iter()
        .map(|&index| {
            (
                index,
                optimize_mmff_conformer(
                    &field,
                    detailed[index].conformer.coords.clone(),
                    retry_max_iters,
                ),
            )
        })
        .collect();
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    let retries: Vec<_> = retry_indices
        .iter()
        .map(|&index| {
            (
                index,
                optimize_mmff_conformer(
                    &field,
                    detailed[index].conformer.coords.clone(),
                    retry_max_iters,
                ),
            )
        })
        .collect();

    for (index, retry) in retries {
        if retry.energy.is_finite() && retry.energy <= detailed[index].conformer.energy {
            detailed[index].conformer = retry;
        }
    }
    Ok(detailed
        .into_iter()
        .map(|result| result.conformer)
        .collect())
}

#[allow(clippy::too_many_arguments)]
fn rdkit_embed_multiple_confs_mmff_retry_policy(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    initial_max_iters: u32,
    retry_max_iters: u32,
    retry_candidates: usize,
    retry_only_when_none_converged: bool,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    let conformers = rdkit_embed_multiple_confs(mh, num_confs, o)?;
    let Some(first) = conformers.first() else {
        return Ok(Vec::new());
    };
    let field = chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, first, 100.0, true)
        .map_err(|error| RdkitEmbedMmffError::ForceField(error.to_string()))?;
    let mut optimized = optimize_mmff_conformers(&field, conformers, initial_max_iters);

    if retry_candidates == 0
        || retry_max_iters == 0
        || (retry_only_when_none_converged
            && optimized.iter().any(|conformer| conformer.status == 0))
    {
        return Ok(optimized);
    }
    let mut retry_indices: Vec<usize> = optimized
        .iter()
        .enumerate()
        .filter_map(|(index, conformer)| {
            (conformer.status != 0 && conformer.energy.is_finite()).then_some(index)
        })
        .collect();
    retry_indices.sort_by(|&left, &right| {
        optimized[left]
            .energy
            .total_cmp(&optimized[right].energy)
            .then_with(|| left.cmp(&right))
    });
    retry_indices.truncate(retry_candidates);

    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    let retries: Vec<_> = retry_indices
        .par_iter()
        .map(|&index| {
            (
                index,
                optimize_mmff_conformer(&field, optimized[index].coords.clone(), retry_max_iters),
            )
        })
        .collect();
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    let retries: Vec<_> = retry_indices
        .iter()
        .map(|&index| {
            (
                index,
                optimize_mmff_conformer(&field, optimized[index].coords.clone(), retry_max_iters),
            )
        })
        .collect();

    for (index, retry) in retries {
        if retry.energy.is_finite() && retry.energy <= optimized[index].energy {
            optimized[index] = retry;
        }
    }
    Ok(optimized)
}

fn optimize_mmff_conformers(
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
    conformers: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
) -> Vec<RdkitMmffConformer> {
    optimize_mmff_conformers_detailed(field, conformers, max_iters)
        .into_iter()
        .map(|result| result.conformer)
        .collect()
}

fn optimize_mmff_conformers_detailed(
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
    conformers: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
) -> Vec<RdkitMmffConformerDiagnostics> {
    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    let optimized = conformers
        .into_par_iter()
        .map(|coords| optimize_mmff_conformer_detailed(field, coords, max_iters))
        .collect();
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    let optimized = conformers
        .into_iter()
        .map(|coords| optimize_mmff_conformer_detailed(field, coords, max_iters))
        .collect();
    optimized
}

fn optimize_mmff_conformer(
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
    coords: Vec<[f64; 3]>,
    max_iters: u32,
) -> RdkitMmffConformer {
    optimize_mmff_conformer_detailed(field, coords, max_iters).conformer
}

fn optimize_mmff_conformer_detailed(
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
    mut coords: Vec<[f64; 3]>,
    max_iters: u32,
) -> RdkitMmffConformerDiagnostics {
    let result = field.optimize_detailed(coords.as_flattened_mut(), max_iters);
    RdkitMmffConformerDiagnostics {
        conformer: RdkitMmffConformer {
            coords,
            status: result.status,
            energy: result.energy,
        },
        iterations: result.iterations,
        final_gradient_test: result.final_gradient_test,
    }
}

fn embed_conformer_attempt(
    prepared: &[PreparedFragment],
    o: &RdkitEmbedOptions,
    index: usize,
) -> Result<Option<Vec<[f64; 3]>>, RdkitEmbedError> {
    let multiplier = i32::try_from(index + 1).map_err(|_| {
        RdkitEmbedError::Unsupported("numConfs exceeds the reproducible seed range".into())
    })?;
    let seed = o.random_seed.checked_mul(multiplier).ok_or_else(|| {
        RdkitEmbedError::Unsupported("conformer seed overflows a signed 32-bit integer".into())
    })?;
    match embed_prepared_fragments(prepared, o, seed) {
        Ok(coords) => Ok(Some(coords)),
        Err(RdkitEmbedError::Failed) => Ok(None),
        Err(error) => Err(error),
    }
}

fn validate_options(mh: &Molecule, o: &RdkitEmbedOptions) -> Result<(), RdkitEmbedError> {
    if o.random_seed < 0 {
        return Err(RdkitEmbedError::Unsupported(
            "randomSeed must be >= 0 (RDKit's -1 is not reproducible)".into(),
        ));
    }
    if mh.atom_count() == 0 {
        return Err(RdkitEmbedError::RdkitException(
            "molecule has no atoms".into(),
        ));
    }
    if o.use_exp_torsion_angle_prefs && !o.use_basic_knowledge {
        return Err(RdkitEmbedError::Unsupported(
            "plain ETDG (no basic knowledge) is not ported".into(),
        ));
    }
    Ok(())
}

struct PreparedFragment {
    source_atoms: Vec<usize>,
    view: RdkitMolView,
    bounds: BoundsMatrix,
    details: Option<CrystalFfDetails>,
    chiral: Vec<ChiralSet>,
    tetrahedral: Vec<ChiralSet>,
    double_bond_ends: Vec<(usize, usize, usize)>,
    stereo_double_bonds: Vec<([usize; 4], i32)>,
}

fn prepare_fragments(
    mh: &Molecule,
    o: &RdkitEmbedOptions,
) -> Result<Vec<PreparedFragment>, RdkitEmbedError> {
    let view = chematic_smiles::rdkit_mol_view(mh)
        .map_err(|e| RdkitEmbedError::Unsupported(e.to_string()))?;
    let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(mh);
    let frags = view.fragments();
    if frags.len() == 1 {
        return Ok(vec![prepare_fragment(
            (0..view.num_atoms()).collect(),
            view,
            labels,
            o,
        )?]);
    }
    frags
        .into_iter()
        .map(|atoms| {
            let piece = view.fragment(&atoms).ok_or_else(|| {
                RdkitEmbedError::Unsupported(
                    "ring perception falls back to RDKit's approximate ring finder".into(),
                )
            })?;
            let piece_labels: Vec<String> = atoms.iter().map(|&a| labels[a].clone()).collect();
            prepare_fragment(atoms, piece, piece_labels, o)
        })
        .collect()
}

fn prepare_fragment(
    source_atoms: Vec<usize>,
    view: RdkitMolView,
    labels: Vec<String>,
    o: &RdkitEmbedOptions,
) -> Result<PreparedFragment, RdkitEmbedError> {
    let etk = o.use_exp_torsion_angle_prefs || o.use_basic_knowledge;
    let Some((bounds, forced_torsions)) = setup_initial_bounds(&view, &labels, o)? else {
        return Err(RdkitEmbedError::Failed);
    };
    let details = if etk {
        let mut d = get_experimental_torsions_with_forced(
            &view,
            o.use_exp_torsion_angle_prefs,
            o.use_small_ring_torsions,
            o.use_macrocycle_torsions,
            o.use_basic_knowledge,
            o.et_version,
            &forced_torsions,
        );
        d.bonds_angles = collect_bonds_and_angles(&view);
        Some(d)
    } else {
        None
    };
    let (chiral, tetrahedral) = find_chiral_sets(&view);
    let (double_bond_ends, stereo_double_bonds) = find_double_bonds(&view);
    Ok(PreparedFragment {
        source_atoms,
        view,
        bounds,
        details,
        chiral,
        tetrahedral,
        double_bond_ends,
        stereo_double_bonds,
    })
}

fn embed_prepared_fragments(
    prepared: &[PreparedFragment],
    o: &RdkitEmbedOptions,
    seed: i32,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
    embed_prepared_fragments_with(prepared, o, seed, false)
}

fn embed_prepared_fragments_with(
    prepared: &[PreparedFragment],
    o: &RdkitEmbedOptions,
    seed: i32,
    only_initial_embedding: bool,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
    let atom_count = prepared
        .iter()
        .flat_map(|fragment| fragment.source_atoms.iter())
        .copied()
        .max()
        .map_or(0, |atom| atom + 1);
    let mut out = vec![[0.0; 3]; atom_count];
    for fragment in prepared {
        let xyz = embed_prepared_fragment(fragment, o, seed, only_initial_embedding)?;
        for (local, &source) in fragment.source_atoms.iter().enumerate() {
            out[source] = xyz[local];
        }
    }
    Ok(out)
}

/// The embedding of one prepared fragment (`EmbedMultipleConfs`' per-fragment
/// loop body and `embedHelper_`).
fn embed_prepared_fragment(
    prepared: &PreparedFragment,
    o: &RdkitEmbedOptions,
    seed: i32,
    only_initial_embedding: bool,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
    let etk = prepared.details.is_some();
    let force_tol = EmbedParams::default().optimizer_force_tol;
    let etk_stage = |p3: &mut [[f64; 3]]| -> Result<bool, EmbedError> {
        let d = prepared.details.as_ref().expect("ETKDG details");
        let mut flat: Vec<f64> = p3.iter().flatten().copied().collect();
        let planar = minimize_with_exp_torsions(
            &mut flat,
            &prepared.bounds,
            d,
            force_tol,
            o.use_basic_knowledge,
        )
        .ok_or(EmbedError::BadDirection)?;
        for (i, p) in p3.iter_mut().enumerate() {
            p.copy_from_slice(&flat[3 * i..3 * i + 3]);
        }
        Ok(planar)
    };
    let args = EmbedArgs {
        mmat: &prepared.bounds,
        chiral_centers: &prepared.chiral,
        tetrahedral_centers: &prepared.tetrahedral,
        double_bond_ends: &prepared.double_bond_ends,
        stereo_double_bonds: &prepared.stereo_double_bonds,
        exp_torsions: if etk { Some(&etk_stage) } else { None },
    };
    let params = EmbedParams {
        max_iterations: o.max_iterations,
        random_seed: seed,
        enforce_chirality: o.enforce_chirality,
        ..EmbedParams::default()
    };
    let embedded = if only_initial_embedding {
        embed_initial_points(prepared.view.num_atoms(), &args, &params, seed)?
    } else {
        embed_points(prepared.view.num_atoms(), &args, &params, seed)?
    };
    match embedded {
        Some(pos) => Ok(pos.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()),
        None => Err(RdkitEmbedError::Failed),
    }
}

/// `rdDistGeom.GetMoleculeBoundsMatrix(mh, set15bounds, scaleVDW=False,
/// doTriangleSmoothing, useMacrocycle14config)` (RDKit 2026.09.1) for an
/// explicit-hydrogen molecule from `add_hydrogens`: the N x N matrix with
/// upper bounds above the diagonal and lower bounds below it. As in RDKit,
/// a failed triangle smoothing still returns the (partially smoothed)
/// matrix.
pub fn rdkit_bounds_matrix(
    mh: &Molecule,
    set15bounds: bool,
    do_triangle_smoothing: bool,
    use_macrocycle14config: bool,
) -> Result<Vec<Vec<f64>>, RdkitEmbedError> {
    let view = chematic_smiles::rdkit_mol_view(mh)
        .map_err(|e| RdkitEmbedError::Unsupported(e.to_string()))?;
    let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(mh);
    let n = view.num_atoms();
    let mut m = BoundsMatrix::new(n);
    init_bounds_mat(&mut m);
    set_topol_bounds(
        &view,
        &labels,
        &mut m,
        set15bounds,
        use_macrocycle14config,
        true,
    )?;
    if do_triangle_smoothing {
        triangle_smooth_bounds(&mut m, 0.0);
    }
    let raw = m.raw();
    Ok((0..n).map(|i| raw[i * n..(i + 1) * n].to_vec()).collect())
}

#[cfg(test)]
mod input_diagnostic_boundary_tests {
    use super::*;

    #[test]
    fn embedding_rejects_empty_graphs_and_unsupported_modes_with_context() {
        let empty = chematic_core::MoleculeBuilder::new().build();
        let error = rdkit_embed_molecule(&empty, &RdkitEmbedOptions::default()).unwrap_err();
        assert_eq!(
            error,
            RdkitEmbedError::RdkitException("molecule has no atoms".into())
        );
        assert!(error.to_string().contains("molecule has no atoms"));
        let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
        let options = RdkitEmbedOptions {
            use_basic_knowledge: false,
            ..RdkitEmbedOptions::default()
        };
        let error = rdkit_embed_molecule(&molecule, &options).unwrap_err();
        assert!(matches!(error, RdkitEmbedError::Unsupported(_)));
        assert!(
            error
                .to_string()
                .contains("plain ETDG (no basic knowledge)")
        );
        assert_eq!(molecule.atom_count(), 9);
    }
}
