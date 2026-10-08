//! `DGeomHelpers::EmbedMultipleConfs` for one conformer of an
//! explicit-hydrogen molecule (`AllChem.EmbedMolecule(AddHs(m), ...)`).

use chematic_core::Molecule;
use chematic_smiles::RdkitMolView;

use super::bounds::{BoundsInvariant, collect_bonds_and_angles, init_bounds_mat, set_topol_bounds};
use super::dg::{BoundsMatrix, ChiralSet, triangle_smooth_bounds};
use super::embed::{EmbedArgs, EmbedError, EmbedParams, embed_points};
use super::etk::{get_experimental_torsions, minimize_with_exp_torsions};

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
) -> Result<Option<BoundsMatrix>, RdkitEmbedError> {
    let n = v.num_atoms();
    let mut m = BoundsMatrix::new(n);
    init_bounds_mat(&mut m);
    set_topol_bounds(
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
    Ok(Some(m))
}

/// `AllChem.EmbedMolecule(mh, randomSeed=..., ...)` on an explicit-hydrogen
/// molecule from `add_hydrogens` (RDKit's `AddHs` atom order): all-atom
/// coordinates, or the reason there are none.
pub fn rdkit_embed_molecule(
    mh: &Molecule,
    o: &RdkitEmbedOptions,
) -> Result<Vec<[f64; 3]>, RdkitEmbedError> {
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
    let view = chematic_smiles::rdkit_mol_view(mh)
        .map_err(|e| RdkitEmbedError::Unsupported(e.to_string()))?;
    if fragment_count(&view) > 1 {
        return Err(RdkitEmbedError::Unsupported(
            "multi-fragment molecules".into(),
        ));
    }
    if o.use_exp_torsion_angle_prefs && !o.use_basic_knowledge {
        return Err(RdkitEmbedError::Unsupported(
            "plain ETDG (no basic knowledge) is not ported".into(),
        ));
    }
    let etk = o.use_exp_torsion_angle_prefs || o.use_basic_knowledge;
    let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(mh);
    let Some(mmat) = setup_initial_bounds(&view, &labels, o)? else {
        return Err(RdkitEmbedError::Failed);
    };
    // initETKDG (getExperimentalTorsions) and setTopolBounds' bonds/angles
    let details = if etk {
        let mut d = get_experimental_torsions(
            &view,
            o.use_exp_torsion_angle_prefs,
            o.use_small_ring_torsions,
            o.use_macrocycle_torsions,
            o.use_basic_knowledge,
            o.et_version,
        );
        d.bonds_angles = collect_bonds_and_angles(&view);
        Some(d)
    } else {
        None
    };
    let (chiral, tet) = find_chiral_sets(&view);
    let (dbe, sdb) = find_double_bonds(&view);
    let force_tol = EmbedParams::default().optimizer_force_tol;
    let etk_stage = |p3: &mut [[f64; 3]]| -> Result<bool, EmbedError> {
        let d = details.as_ref().expect("ETKDG details");
        let mut flat: Vec<f64> = p3.iter().flatten().copied().collect();
        let planar =
            minimize_with_exp_torsions(&mut flat, &mmat, d, force_tol, o.use_basic_knowledge)
                .ok_or(EmbedError::BadDirection)?;
        for (i, p) in p3.iter_mut().enumerate() {
            p.copy_from_slice(&flat[3 * i..3 * i + 3]);
        }
        Ok(planar)
    };
    let args = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral,
        tetrahedral_centers: &tet,
        double_bond_ends: &dbe,
        stereo_double_bonds: &sdb,
        exp_torsions: if etk { Some(&etk_stage) } else { None },
    };
    let params = EmbedParams {
        max_iterations: o.max_iterations,
        random_seed: o.random_seed,
        enforce_chirality: o.enforce_chirality,
        ..EmbedParams::default()
    };
    // conformer 0: new_seed = (0 + 1) * randomSeed
    match embed_points(view.num_atoms(), &args, &params, o.random_seed)? {
        Some(pos) => Ok(pos.chunks(3).map(|c| [c[0], c[1], c[2]]).collect()),
        None => Err(RdkitEmbedError::Failed),
    }
}

fn fragment_count(v: &RdkitMolView) -> usize {
    let n = v.num_atoms();
    let mut seen = vec![false; n];
    let mut count = 0;
    for s in 0..n {
        if seen[s] {
            continue;
        }
        count += 1;
        seen[s] = true;
        let mut stack = vec![s];
        while let Some(a) = stack.pop() {
            for b in v.neighbors(a) {
                if !seen[b] {
                    seen[b] = true;
                    stack.push(b);
                }
            }
        }
    }
    count
}
