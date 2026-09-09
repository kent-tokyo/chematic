//! Torsion Motif Angle Distribution (TMAD) — diagnostic-only.
//!
//! Extracts rotatable-bond torsions from a molecule, builds a circular
//! (wrap-aware) histogram of their observed dihedral angles across a
//! [`ConformerEnsemble`], fits a von Mises mixture to that histogram, and
//! compares two fitted profiles with Jensen-Shannon divergence. Environment
//! labels ([`TorsionEnvironment`]) are caller-supplied metadata only — this
//! module has no notion of what "crystal" or "water" physically means, it
//! just partitions/tags observations by whatever label the caller attaches
//! to an ensemble.
//!
//! Inspired by environment-conditioned torsion-distribution analysis in
//! recent cheminformatics literature (circular-statistics histogram fitting
//! and automated distribution comparison) — this is an independent
//! from-scratch implementation of the general circular-statistics idea (von
//! Mises mixtures, Jensen-Shannon divergence), not a port of any paper's
//! code or numbers.
//!
//! **Diagnostic only**: nothing here is wired into ETKDG or any embedding/
//! sampling path. It is a measurement tool over already-generated conformer
//! ensembles.
//!
//! **Known scope limits of [`extract_torsion_motifs`]**, both inherited from
//! `chematic_chem::rotatable_bond_atom_pairs` rather than introduced here
//! (see that function's own doc comment for the full exclusion list):
//! - Amide C–N bonds (the omega torsion) are never returned — they are
//!   excluded from "rotatable" by that shared definition. A caller wanting
//!   amide torsions needs a different extraction path (e.g.
//!   `etkdg_knowledge`'s own `candidate_central_bonds`/`classify_bond`,
//!   which keeps amide bonds and flags them via `amide_like`, but that
//!   module is currently private to this crate).
//! - A biphenyl-like inter-ring single bond is only found if the input
//!   SMILES spells it as an explicit single bond (`c1ccc(-c2ccccc2)cc1`).
//!   Without the explicit `-`, both endpoints being lowercase makes this
//!   crate's SMILES parser assign `BondOrder::Aromatic` to that bond (its
//!   `implicit_bond` rule is "both atoms aromatic -> aromatic bond",
//!   independent of ring membership) even though the bond itself is not in
//!   any ring — and an aromatic-order bond is not "single", so it is
//!   invisible to the rotatable-bond filter. This is spelling-dependent,
//!   pre-existing behavior of the shared parser/descriptor, not something
//!   this diagnostic tool works around.

use std::f64::consts::PI;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, MatchConfig, QueryMolecule,
    find_matches_with_config,
};

use crate::conformer::ConformerEnsemble;
use crate::coords::Coords3D;
use crate::mol_transforms::get_dihedral_deg;

// ─── Motif identification ──────────────────────────────────────────────────

/// A candidate rotatable-bond torsion, identified by its four atoms plus
/// lightweight element/aromaticity context (element symbols and aromatic
/// flags for atoms 0..3). Atom identity is only stable within one molecule's
/// own numbering — comparing motifs across different molecules should use
/// [`TorsionMotif::signature`], not `atoms`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TorsionMotif {
    /// A-B-C-D, where B-C is the rotatable central bond.
    pub atoms: [AtomIdx; 4],
    pub elements: [&'static str; 4],
    pub aromatic: [bool; 4],
}

impl TorsionMotif {
    /// Human-readable element-pattern label, e.g. `"C-C-C-C"` or `"C-N-C-C"`.
    pub fn signature(&self) -> String {
        self.elements.join("-")
    }
}

/// Find candidate rotatable-bond torsions in `mol`.
///
/// Reuses `chematic_chem::rotatable_bond_atom_pairs` (RDKit-compatible
/// strict definition: excludes ring/amide/allene/alkyne-adjacent bonds) so
/// this module's notion of "rotatable" never drifts from the rest of the
/// codebase's. For each rotatable bond B-C, the outer atoms A and D are
/// picked deterministically as the lowest-`AtomIdx` neighbor on each side
/// (ties only affect which chemically-symmetric substituent is reported,
/// never whether a motif is found).
///
/// `coords` is used only to skip bonds with degenerate geometry (e.g.
/// collinear substituents, where no dihedral is defined); the returned
/// motifs carry no angle values of their own.
pub fn extract_torsion_motifs(mol: &Molecule, coords: &Coords3D) -> Vec<TorsionMotif> {
    let mut motifs = Vec::new();
    for (b, c) in chematic_chem::rotatable_bond_atom_pairs(mol) {
        let Some(a) = mol
            .neighbors(b)
            .filter(|&(n, _)| n != c)
            .map(|(n, _)| n)
            .min()
        else {
            continue;
        };
        let Some(d) = mol
            .neighbors(c)
            .filter(|&(n, _)| n != b)
            .map(|(n, _)| n)
            .min()
        else {
            continue;
        };
        if get_dihedral_deg(coords, a, b, c, d).is_none() {
            continue;
        }
        motifs.push(TorsionMotif {
            atoms: [a, b, c, d],
            elements: [
                mol.atom(a).element.symbol(),
                mol.atom(b).element.symbol(),
                mol.atom(c).element.symbol(),
                mol.atom(d).element.symbol(),
            ],
            aromatic: [
                mol.atom(a).aromatic,
                mol.atom(b).aromatic,
                mol.atom(c).aromatic,
                mol.atom(d).aromatic,
            ],
        });
    }
    motifs
}

/// Observed dihedral angles (degrees, `(-180, 180]`) for `motif` across
/// every conformer in `ensemble`. Conformers whose geometry makes the
/// dihedral undefined (degenerate/collinear) are skipped rather than
/// producing a bogus value.
pub fn motif_angles_deg(ensemble: &ConformerEnsemble, motif: &TorsionMotif) -> Vec<f64> {
    let [a, b, c, d] = motif.atoms;
    (0..ensemble.conformer_count())
        .filter_map(|i| ensemble.get_conformer(i))
        .filter_map(|coords| get_dihedral_deg(coords, a, b, c, d))
        .collect()
}

/// Return a symmetry-aware normalized torsion-angle distance for two
/// geometries of the same molecule.
///
/// The result is the minimum mean circular angular difference over all
/// topology-preserving self-mappings, divided by 180, so it is in `[0, 1]`.
/// This is a deterministic, bounded TFD-like diagnostic for local comparisons;
/// it is not claimed to be numerically identical to RDKit's torsion fingerprint
/// implementation.  Degenerate torsions are omitted from the mean.  If no
/// rotatable torsions are defined, the distance is zero.
pub fn torsion_distance_symmetric(mol: &Molecule, coords_a: &Coords3D, coords_b: &Coords3D) -> f64 {
    let motifs = extract_torsion_motifs(mol, coords_a);
    if motifs.is_empty() {
        return 0.0;
    }
    let query = molecule_self_query(mol);
    let matches = find_matches_with_config(
        &query,
        mol,
        &MatchConfig {
            uniquify: false,
            ..MatchConfig::default()
        },
    );
    let mut best = f64::INFINITY;
    for mapping in matches {
        let mut total = 0.0;
        let mut count = 0usize;
        for motif in &motifs {
            let [a, b, c, d] = motif.atoms;
            let mapped = [
                mapping[&(a.0 as usize)],
                mapping[&(b.0 as usize)],
                mapping[&(c.0 as usize)],
                mapping[&(d.0 as usize)],
            ];
            let Some(angle_a) = get_dihedral_deg(coords_a, a, b, c, d) else {
                continue;
            };
            let Some(angle_b) =
                get_dihedral_deg(coords_b, mapped[0], mapped[1], mapped[2], mapped[3])
            else {
                continue;
            };
            total += circular_angle_difference_deg(angle_a, angle_b);
            count += 1;
        }
        if count > 0 {
            best = best.min(total / count as f64 / 180.0);
        }
    }
    if best.is_finite() { best.min(1.0) } else { 0.0 }
}

fn circular_angle_difference_deg(a: f64, b: f64) -> f64 {
    let mut d = (a - b).abs() % 360.0;
    if d > 180.0 {
        d = 360.0 - d;
    }
    d
}

fn molecule_self_query(mol: &Molecule) -> QueryMolecule {
    let mut query = QueryMolecule::new();
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        query.add_atom(AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(
                atom.element.atomic_number(),
            ))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge))),
        ));
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        let primitive = match bond.order {
            BondOrder::Single => BondPrimitive::Single,
            BondOrder::Double => BondPrimitive::Double,
            BondOrder::Triple => BondPrimitive::Triple,
            BondOrder::Aromatic => BondPrimitive::Aromatic,
            _ => BondPrimitive::Any,
        };
        query.add_bond(
            bond.atom1.0 as usize,
            bond.atom2.0 as usize,
            BondQuery::Primitive(primitive),
        );
    }
    query
}

// ─── Environment label ─────────────────────────────────────────────────────

/// Caller-supplied context label for a conformer/ensemble. This crate does
/// not model the physics of any of these — it is purely a partition key for
/// grouping angle observations before comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TorsionEnvironment {
    Crystal,
    Vacuum,
    Water,
    Hexane,
}

impl TorsionEnvironment {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Crystal => "crystal",
            Self::Vacuum => "vacuum",
            Self::Water => "water",
            Self::Hexane => "hexane",
        }
    }
}

// ─── Circular histogram ────────────────────────────────────────────────────

/// A circular (wraparound-aware) histogram of dihedral-angle observations
/// over `[-180, 180)`, split into equal-width bins. Bin 0 covers
/// `[-180, -180+w)`; the last bin covers `[180-w, 180)`. Bins are adjacent
/// on the circle: bin 0 and the last bin are physical neighbors even though
/// they sit at opposite ends of the array.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TorsionHistogram {
    pub bin_counts: Vec<u32>,
}

/// Wrap an arbitrary angle (degrees) into the canonical `[-180, 180)` range,
/// so `+180` and `-180` (and any multiple-of-360 alias) always collapse to
/// the same representative instead of splitting a boundary cluster across
/// two floating-point spellings of the same physical angle.
fn wrap_deg(a: f64) -> f64 {
    let mut w = a % 360.0;
    if w >= 180.0 {
        w -= 360.0;
    }
    if w < -180.0 {
        w += 360.0;
    }
    w
}

impl TorsionHistogram {
    /// Build a histogram from raw angle observations (degrees, any range —
    /// each is wrapped into `[-180, 180)` before binning).
    pub fn from_angles_deg(angles_deg: &[f64], n_bins: usize) -> Self {
        assert!(n_bins > 0, "n_bins must be positive");
        let mut bin_counts = vec![0u32; n_bins];
        let width = 360.0 / n_bins as f64;
        for &a in angles_deg {
            let wrapped = wrap_deg(a);
            let idx = (((wrapped + 180.0) / width) as usize).min(n_bins - 1);
            bin_counts[idx] += 1;
        }
        Self { bin_counts }
    }

    pub fn n_bins(&self) -> usize {
        self.bin_counts.len()
    }

    pub fn bin_width_deg(&self) -> f64 {
        360.0 / self.n_bins() as f64
    }

    pub fn bin_center_deg(&self, i: usize) -> f64 {
        -180.0 + (i as f64 + 0.5) * self.bin_width_deg()
    }

    pub fn n_observations(&self) -> u32 {
        self.bin_counts.iter().sum()
    }

    /// Circular mean angle (degrees), or `None` if the histogram is empty.
    /// Computed via the sin/cos resultant vector, so a cluster straddling
    /// the ±180° wrap boundary averages correctly (unlike a naive
    /// arithmetic mean of e.g. 179° and -179°, which would wrongly give 0°).
    pub fn circular_mean_deg(&self) -> Option<f64> {
        if self.n_observations() == 0 {
            return None;
        }
        let (s, c) = self.resultant_vector();
        Some(s.atan2(c).to_degrees())
    }

    /// Mean resultant length `R` in `[0, 1]` — 1 means all mass at one
    /// angle, 0 means uniformly spread around the circle.
    pub fn mean_resultant_length(&self) -> f64 {
        let n = self.n_observations();
        if n == 0 {
            return 0.0;
        }
        let (s, c) = self.resultant_vector();
        (s * s + c * c).sqrt() / n as f64
    }

    fn resultant_vector(&self) -> (f64, f64) {
        let (mut s, mut c) = (0.0, 0.0);
        for (i, &count) in self.bin_counts.iter().enumerate() {
            let rad = self.bin_center_deg(i).to_radians();
            s += count as f64 * rad.sin();
            c += count as f64 * rad.cos();
        }
        (s, c)
    }

    /// Per-bin probabilities (sums to 1), or all-zero if empty.
    pub fn probabilities(&self) -> Vec<f64> {
        let n = self.n_observations();
        if n == 0 {
            return vec![0.0; self.n_bins()];
        }
        self.bin_counts
            .iter()
            .map(|&c| c as f64 / n as f64)
            .collect()
    }
}

// ─── Von Mises mixture fit ─────────────────────────────────────────────────

/// One von Mises component: `weight * vonmises(mu_deg, kappa)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VonMisesComponent {
    pub weight: f64,
    pub mu_deg: f64,
    pub kappa: f64,
}

/// A fitted mixture of von Mises distributions over a [`TorsionHistogram`].
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TorsionProfileFit {
    pub components: Vec<VonMisesComponent>,
}

// Guards `bessel_i0(kappa)` / `exp(kappa)` against f64 overflow (exp(709) is
// the last finite value); ponytail: fixed ceiling rather than a log-space
// density, upgrade to log-density evaluation if a caller ever needs a
// genuinely near-delta (kappa > ~700) fitted component.
const MAX_KAPPA: f64 = 700.0;

/// Modified Bessel function of the first kind, order 0 — Abramowitz &
/// Stegun 9.8.1/9.8.2 polynomial approximation (the standard "Numerical
/// Recipes bessi0" rational-polynomial form), independently implemented
/// here (not from any paper's code) since von Mises's normalizing constant
/// `1 / (2*pi*I0(kappa))` needs it and no Bessel-function crate is a
/// workspace dependency.
fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let t = (x / 3.75).powi(2);
        1.0 + t
            * (3.515_622_9
                + t * (3.089_942_4
                    + t * (1.206_749_2 + t * (0.265_973_2 + t * (0.036_076_8 + t * 0.004_581_3)))))
    } else {
        let t = 3.75 / ax;
        (ax.exp() / ax.sqrt())
            * (0.398_942_28
                + t * (0.013_285_92
                    + t * (0.002_253_19
                        + t * (-0.001_575_65
                            + t * (0.009_162_81
                                + t * (-0.020_577_06
                                    + t * (0.026_355_37
                                        + t * (-0.016_476_33 + t * 0.003_923_77))))))))
    }
}

fn von_mises_density_rad(theta: f64, mu: f64, kappa: f64) -> f64 {
    if kappa <= 0.0 {
        return 1.0 / (2.0 * PI);
    }
    (kappa * (theta - mu).cos()).exp() / (2.0 * PI * bessel_i0(kappa))
}

/// Fisher's (1993, *Statistical Analysis of Circular Data*) closed-form
/// approximate inverse of `A(kappa) = I1(kappa) / I0(kappa) = R`: a standard
/// circular-statistics method-of-moments estimator, independently
/// implemented here (three-branch polynomial/rational approximation).
fn kappa_from_mean_resultant_length(r: f64) -> f64 {
    let r = r.clamp(0.0, 0.999_999);
    let kappa = if r < 1e-12 {
        0.0
    } else if r < 0.53 {
        2.0 * r + r.powi(3) + 5.0 * r.powi(5) / 6.0
    } else if r < 0.85 {
        -0.4 + 1.39 * r + 0.43 / (1.0 - r)
    } else {
        1.0 / (r.powi(3) - 4.0 * r.powi(2) + 3.0 * r)
    };
    kappa.min(MAX_KAPPA)
}

impl TorsionProfileFit {
    /// Mixture density at `theta_deg`.
    pub fn density_deg(&self, theta_deg: f64) -> f64 {
        let theta = theta_deg.to_radians();
        self.components
            .iter()
            .map(|comp| {
                comp.weight * von_mises_density_rad(theta, comp.mu_deg.to_radians(), comp.kappa)
            })
            .sum()
    }

    /// Discretize the fitted density over `n_bins` equal-width bins spanning
    /// the circle, normalized to sum to 1 (for use as a discrete
    /// probability distribution, e.g. in [`torsion_profile_distance`]).
    pub fn probabilities(&self, n_bins: usize) -> Vec<f64> {
        if n_bins == 0 {
            return Vec::new();
        }
        let width = 360.0 / n_bins as f64;
        let raw: Vec<f64> = (0..n_bins)
            .map(|i| self.density_deg(-180.0 + (i as f64 + 0.5) * width))
            .collect();
        let sum: f64 = raw.iter().sum();
        if sum <= 0.0 {
            vec![0.0; n_bins]
        } else {
            raw.iter().map(|v| v / sum).collect()
        }
    }
}

/// Fit a mixture of `n_components` von Mises distributions to `histogram`
/// via binned EM (`n_iters` fixed iterations — no convergence-threshold
/// early exit, so the result is a pure function of the inputs and therefore
/// deterministic). The M-step uses method-of-moments (circular mean +
/// Fisher's kappa inversion above), not full MLE — reasonable-effort per
/// this tool's diagnostic scope, not a publication-grade optimizer.
///
/// Returns an empty fit (no components) if `histogram` has zero
/// observations, `n_components == 0`, or `histogram.n_bins() == 0`.
pub fn fit_von_mises_mixture(
    histogram: &TorsionHistogram,
    n_components: usize,
    n_iters: usize,
) -> TorsionProfileFit {
    let n_bins = histogram.n_bins();
    let total = histogram.n_observations() as f64;
    if n_components == 0 || n_bins == 0 || total == 0.0 {
        return TorsionProfileFit::default();
    }

    let centers_rad: Vec<f64> = (0..n_bins)
        .map(|i| histogram.bin_center_deg(i).to_radians())
        .collect();
    let weights: Vec<f64> = histogram.bin_counts.iter().map(|&c| c as f64).collect();

    let mean_deg = histogram.circular_mean_deg().unwrap_or(0.0);
    let mut mus: Vec<f64> = (0..n_components)
        .map(|k| (mean_deg + k as f64 * 360.0 / n_components as f64).to_radians())
        .collect();
    let mut kappas: Vec<f64> = vec![2.0; n_components];
    let mut pis: Vec<f64> = vec![1.0 / n_components as f64; n_components];

    for _ in 0..n_iters {
        // E-step: responsibilities gamma[i][k].
        let mut resp: Vec<Vec<f64>> = vec![vec![0.0; n_components]; n_bins];
        for i in 0..n_bins {
            let mut row_sum = 0.0;
            for k in 0..n_components {
                let d = pis[k] * von_mises_density_rad(centers_rad[i], mus[k], kappas[k]);
                resp[i][k] = d;
                row_sum += d;
            }
            if row_sum > 0.0 {
                for k in 0..n_components {
                    resp[i][k] /= row_sum;
                }
            } else {
                for r in &mut resp[i] {
                    *r = 1.0 / n_components as f64;
                }
            }
        }

        // M-step: method-of-moments per component.
        for k in 0..n_components {
            let mut nk = 0.0;
            let mut s = 0.0;
            let mut c = 0.0;
            for i in 0..n_bins {
                let wk = weights[i] * resp[i][k];
                nk += wk;
                s += wk * centers_rad[i].sin();
                c += wk * centers_rad[i].cos();
            }
            if nk > 0.0 {
                pis[k] = nk / total;
                mus[k] = s.atan2(c);
                let r = (s * s + c * c).sqrt() / nk;
                kappas[k] = kappa_from_mean_resultant_length(r);
            }
        }
    }

    TorsionProfileFit {
        components: (0..n_components)
            .map(|k| VonMisesComponent {
                weight: pis[k],
                mu_deg: mus[k].to_degrees(),
                kappa: kappas[k],
            })
            .collect(),
    }
}

// ─── Jensen-Shannon distance ────────────────────────────────────────────────

/// Number of grid points [`torsion_profile_distance`] discretizes each fit
/// over before comparing — 1° resolution, fixed so the metric is a pure
/// function of the two fits (no caller-tunable grid to accidentally make
/// two distance calls incomparable).
const JS_GRID_BINS: usize = 360;

fn kl_divergence(p: &[f64], m: &[f64]) -> f64 {
    p.iter()
        .zip(m.iter())
        .map(|(&pi, &mi)| {
            if pi <= 0.0 || mi <= 0.0 {
                0.0
            } else {
                pi * (pi / mi).ln()
            }
        })
        .sum()
}

/// Jensen-Shannon divergence (natural log, bounded by `ln(2)`) between two
/// equal-length discrete probability distributions.
fn jensen_shannon_divergence(p: &[f64], q: &[f64]) -> f64 {
    debug_assert_eq!(p.len(), q.len());
    let m: Vec<f64> = p
        .iter()
        .zip(q.iter())
        .map(|(&pi, &qi)| 0.5 * (pi + qi))
        .collect();
    0.5 * kl_divergence(p, &m) + 0.5 * kl_divergence(q, &m)
}

/// Jensen-Shannon divergence between two fitted torsion profiles: symmetric,
/// bounded in `[0, ln(2)]`, zero iff the two fits discretize to the same
/// distribution over the shared 1°-resolution grid.
pub fn torsion_profile_distance(a: &TorsionProfileFit, b: &TorsionProfileFit) -> f64 {
    let p = a.probabilities(JS_GRID_BINS);
    let q = b.probabilities(JS_GRID_BINS);
    jensen_shannon_divergence(&p, &q)
}

// ─── Deterministic JSON export ─────────────────────────────────────────────

/// Serialize one computed torsion profile (motif identity, environment
/// label, angle histogram, fitted mixture parameters) to a deterministic,
/// pretty-printed JSON string — same inputs always produce byte-identical
/// output (no `HashMap` iteration order involved: every field below comes
/// from a `Vec`/array in fixed order, and `serde_json::Value`'s default
/// `Map` is `BTreeMap`-backed so key order is alphabetical regardless of
/// insertion order).
pub fn torsion_profile_to_json(
    motif: &TorsionMotif,
    environment: TorsionEnvironment,
    histogram: &TorsionHistogram,
    fit: &TorsionProfileFit,
) -> String {
    use serde_json::json;

    let components: Vec<serde_json::Value> = fit
        .components
        .iter()
        .map(|comp| {
            json!({
                "weight": comp.weight,
                "mu_deg": comp.mu_deg,
                "kappa": comp.kappa,
            })
        })
        .collect();

    let root = json!({
        "motif": {
            "atoms": motif.atoms.iter().map(|a| a.0).collect::<Vec<u32>>(),
            "elements": motif.elements,
            "aromatic": motif.aromatic,
            "signature": motif.signature(),
        },
        "environment": environment.as_str(),
        "histogram": {
            "n_bins": histogram.n_bins(),
            "bin_width_deg": histogram.bin_width_deg(),
            "bin_counts": histogram.bin_counts,
            "n_observations": histogram.n_observations(),
        },
        "fit": {
            "components": components,
        },
    });

    serde_json::to_string_pretty(&root).unwrap_or_default()
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dg::generate_coords;
    use crate::prng::Prng;
    use chematic_smiles::parse;

    fn circular_diff_deg(a: f64, b: f64) -> f64 {
        let mut d = (a - b) % 360.0;
        if d > 180.0 {
            d -= 360.0;
        }
        if d < -180.0 {
            d += 360.0;
        }
        d.abs()
    }

    // ── extract_torsion_motifs / motif_angles_deg ──────────────────────────

    #[test]
    fn butane_has_one_rotatable_torsion_motif() {
        let mol = parse("CCCC").unwrap();
        let coords = generate_coords(&mol);
        let motifs = extract_torsion_motifs(&mol, &coords);
        assert_eq!(motifs.len(), 1, "butane has exactly one rotatable bond");
        assert_eq!(
            motifs[0].atoms,
            [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)]
        );
        assert_eq!(motifs[0].signature(), "C-C-C-C");
    }

    #[test]
    fn benzene_has_no_rotatable_torsion_motifs() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let motifs = extract_torsion_motifs(&mol, &coords);
        assert!(motifs.is_empty(), "benzene's ring bonds are not rotatable");
    }

    #[test]
    fn amide_c_n_bond_is_not_a_torsion_motif() {
        // Pins a documented scope limit (see module doc): the shared
        // rotatable-bond definition excludes amide C-N bonds, so the omega
        // torsion is never returned by this extractor.
        let mol = parse("CC(=O)NC").unwrap();
        let coords = generate_coords(&mol);
        let motifs = extract_torsion_motifs(&mol, &coords);
        assert!(
            motifs.is_empty(),
            "amide C-N bond must not be a rotatable motif"
        );
    }

    #[test]
    fn biphenyl_inter_ring_bond_is_a_motif_only_when_spelled_as_explicit_single() {
        // Pins another documented scope limit: an implicit (no "-") bond
        // between two aromatic atoms parses as BondOrder::Aromatic in this
        // crate's SMILES parser regardless of ring membership, so it is
        // invisible to the "is_single" rotatable-bond filter.
        let explicit = parse("c1ccc(-c2ccccc2)cc1").unwrap();
        let explicit_coords = generate_coords(&explicit);
        let explicit_motifs = extract_torsion_motifs(&explicit, &explicit_coords);
        assert_eq!(
            explicit_motifs.len(),
            1,
            "explicit '-' spelling must find the inter-ring torsion"
        );

        let implicit = parse("c1ccc(c2ccccc2)cc1").unwrap();
        let implicit_coords = generate_coords(&implicit);
        let implicit_motifs = extract_torsion_motifs(&implicit, &implicit_coords);
        assert!(
            implicit_motifs.is_empty(),
            "documenting current behavior: implicit spelling's Aromatic-order \
             inter-ring bond is not found -- see module doc's known scope limits"
        );
    }

    #[test]
    fn motif_angles_deg_collects_one_value_per_conformer() {
        let mol = parse("CCCC").unwrap();
        let coords = generate_coords(&mol);
        let motif = extract_torsion_motifs(&mol, &coords).remove(0);
        let mut ensemble = ConformerEnsemble::with_conformer(mol, coords.clone()).unwrap();
        ensemble.add_conformer(coords).unwrap();
        let angles = motif_angles_deg(&ensemble, &motif);
        assert_eq!(angles.len(), 2);
        assert!((angles[0] - angles[1]).abs() < 1e-9);
    }

    #[test]
    fn symmetric_torsion_distance_is_zero_for_terminal_swap() {
        let mol = parse("CCC").unwrap();
        let a = generate_coords(&mol);
        let mut b = a.clone();
        let first = a.get(AtomIdx(0));
        let last = a.get(AtomIdx(2));
        b.set(AtomIdx(0), last);
        b.set(AtomIdx(2), first);
        assert!(
            torsion_distance_symmetric(&mol, &a, &b) < 1e-12,
            "terminal-atom automorphism must make the torsion distance zero"
        );
    }

    #[test]
    fn symmetric_torsion_distance_is_bounded_and_detects_change() {
        let mol = parse("CCCC").unwrap();
        let a = generate_coords(&mol);
        let motif = extract_torsion_motifs(&mol, &a).remove(0);
        let angle = get_dihedral_deg(
            &a,
            motif.atoms[0],
            motif.atoms[1],
            motif.atoms[2],
            motif.atoms[3],
        )
        .unwrap();
        let b = crate::mol_transforms::set_dihedral(
            &a,
            &mol,
            motif.atoms[0],
            motif.atoms[1],
            motif.atoms[2],
            motif.atoms[3],
            (angle + 60.0).to_radians(),
        );
        let distance = torsion_distance_symmetric(&mol, &a, &b);
        assert!(
            distance > 0.1,
            "changed torsion should be detected: {distance}"
        );
        assert!(
            distance <= 1.0,
            "normalized distance must be bounded: {distance}"
        );
    }

    #[test]
    fn symmetric_torsion_distance_corpus_is_bounded_and_change_sensitive() {
        // Keep this corpus deliberately small and dependency-free: it is a
        // local invariant gate for the diagnostic, not an RDKit parity claim.
        let smiles = [
            "CCCCC",
            "CCOC(C)C",
            "CCNCC",
            "CCCOC",
            "CC(C)CC(C)C",
            "CCOC(=O)CC",
            "CCCCCC",
            "CC(C)COC",
            "CCOC(C)CC",
            "CCN(CC)CC",
        ];
        let mut measured = 0usize;
        for smiles in smiles {
            let mol = parse(smiles).unwrap_or_else(|err| panic!("{smiles}: {err}"));
            let coords = generate_coords(&mol);
            let Some(motif) = extract_torsion_motifs(&mol, &coords).into_iter().next() else {
                continue;
            };
            let angle = get_dihedral_deg(
                &coords,
                motif.atoms[0],
                motif.atoms[1],
                motif.atoms[2],
                motif.atoms[3],
            )
            .expect("selected corpus motif must have a defined dihedral");
            let changed = crate::mol_transforms::set_dihedral(
                &coords,
                &mol,
                motif.atoms[0],
                motif.atoms[1],
                motif.atoms[2],
                motif.atoms[3],
                (angle + 45.0).to_radians(),
            );
            let distance = torsion_distance_symmetric(&mol, &coords, &changed);
            assert!(distance.is_finite(), "{smiles}: distance must be finite");
            assert!((0.0..=1.0).contains(&distance), "{smiles}: {distance}");
            assert!(
                distance > 0.05,
                "{smiles}: changed torsion was not detected: {distance}"
            );
            measured += 1;
        }
        assert!(
            measured >= 4,
            "corpus must exercise at least four rotatable molecules"
        );
    }

    // ── TorsionHistogram: circular correctness ─────────────────────────────

    #[test]
    fn histogram_preserves_observation_count() {
        let angles = [0.0, 45.0, -170.0, 179.9, -179.9];
        let hist = TorsionHistogram::from_angles_deg(&angles, 36);
        assert_eq!(hist.n_observations(), angles.len() as u32);
        assert_eq!(hist.bin_counts.iter().sum::<u32>(), angles.len() as u32);
    }

    #[test]
    fn histogram_boundary_angles_do_not_panic_or_go_out_of_range() {
        let angles = [180.0, -180.0, 360.0, -360.0, 540.0];
        let hist = TorsionHistogram::from_angles_deg(&angles, 36);
        assert_eq!(hist.n_observations(), angles.len() as u32);
    }

    #[test]
    fn circular_mean_handles_wraparound_cluster_correctly() {
        // A tight cluster straddling the +-180 deg wrap boundary. A naive
        // arithmetic mean of e.g. -179 and 179 would give ~0 (exactly
        // wrong); the circular mean must land near +-180.
        let angles = [178.0, 179.0, -179.0, -178.0, 180.0];
        let hist = TorsionHistogram::from_angles_deg(&angles, 36);
        let mean = hist.circular_mean_deg().expect("non-empty histogram");
        let dist_to_180 = circular_diff_deg(mean, 180.0);
        assert!(
            dist_to_180 < 5.0,
            "circular mean {mean} should be within 5 deg of the wraparound cluster at 180, \
             not pulled toward 0 by a naive arithmetic average"
        );
        // Tight cluster -> high concentration, not spread across the circle.
        assert!(
            hist.mean_resultant_length() > 0.9,
            "R={} should be close to 1 for a tight cluster",
            hist.mean_resultant_length()
        );
    }

    #[test]
    fn circular_mean_of_empty_histogram_is_none() {
        let hist = TorsionHistogram::from_angles_deg(&[], 36);
        assert_eq!(hist.circular_mean_deg(), None);
        assert_eq!(hist.mean_resultant_length(), 0.0);
        assert_eq!(hist.probabilities(), vec![0.0; 36]);
    }

    #[test]
    fn probabilities_sum_to_one() {
        let angles = [-90.0, 0.0, 45.0, 170.0, -170.0];
        let hist = TorsionHistogram::from_angles_deg(&angles, 36);
        let sum: f64 = hist.probabilities().iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }

    // ── Von Mises mixture fit: recover known synthetic parameters ──────────

    /// Rejection sampler for VonMises(mu, kappa): standard sampling by
    /// rejection against a uniform proposal on the circle (accept-prob
    /// exp(kappa*(cos(theta-mu)-1)) <= 1) -- independently implemented here
    /// to synthesize ground-truth test data, not taken from any paper.
    fn sample_von_mises(rng: &mut Prng, mu: f64, kappa: f64) -> f64 {
        if kappa < 1e-8 {
            return rng.f64() * 2.0 * PI - PI;
        }
        loop {
            let theta = rng.f64() * 2.0 * PI - PI;
            let u = rng.f64();
            let accept = (kappa * ((theta - mu).cos() - 1.0)).exp();
            if u <= accept {
                return theta;
            }
        }
    }

    #[test]
    fn single_component_fit_recovers_synthetic_von_mises_parameters() {
        let true_mu_deg: f64 = 45.0;
        let true_kappa = 4.0;
        let mut rng = Prng::from_seed(20260809);
        let angles_deg: Vec<f64> = (0..5000)
            .map(|_| sample_von_mises(&mut rng, true_mu_deg.to_radians(), true_kappa).to_degrees())
            .collect();

        let hist = TorsionHistogram::from_angles_deg(&angles_deg, 72);
        let fit = fit_von_mises_mixture(&hist, 1, 25);
        assert_eq!(fit.components.len(), 1);
        let comp = fit.components[0];

        assert!(
            circular_diff_deg(comp.mu_deg, true_mu_deg) < 8.0,
            "recovered mu {} too far from true mu {}",
            comp.mu_deg,
            true_mu_deg
        );
        assert!(
            (comp.kappa - true_kappa).abs() / true_kappa < 0.4,
            "recovered kappa {} too far from true kappa {}",
            comp.kappa,
            true_kappa
        );
        assert!(
            (comp.weight - 1.0).abs() < 1e-9,
            "single component carries all the weight"
        );
    }

    #[test]
    fn two_component_fit_recovers_synthetic_von_mises_mixture_parameters() {
        // Two well-separated modes (120 deg apart, kappa=8 each -> angular
        // SD ~20 deg, comfortably resolvable) at unequal weights, so the
        // E-step actually has to do real separation work rather than the
        // K=2 default init (mean +/- 180 deg) accidentally landing on the
        // true modes.
        let (mu1_deg, kappa1, n1): (f64, f64, usize) = (-60.0, 8.0, 3000);
        let (mu2_deg, kappa2, n2): (f64, f64, usize) = (60.0, 8.0, 2000);
        let mut rng = Prng::from_seed(20260810);
        let mut angles_deg: Vec<f64> = (0..n1)
            .map(|_| sample_von_mises(&mut rng, mu1_deg.to_radians(), kappa1).to_degrees())
            .collect();
        angles_deg.extend(
            (0..n2).map(|_| sample_von_mises(&mut rng, mu2_deg.to_radians(), kappa2).to_degrees()),
        );

        let hist = TorsionHistogram::from_angles_deg(&angles_deg, 72);
        let fit = fit_von_mises_mixture(&hist, 2, 40);
        assert_eq!(fit.components.len(), 2);

        // EM has no canonical component ordering -- sort by mu before
        // comparing against the (mu1 < mu2) ground truth.
        let mut comps = fit.components.clone();
        comps.sort_by(|a, b| a.mu_deg.partial_cmp(&b.mu_deg).unwrap());

        let expected_weight1 = n1 as f64 / (n1 + n2) as f64;
        let expected_weight2 = n2 as f64 / (n1 + n2) as f64;

        assert!(
            circular_diff_deg(comps[0].mu_deg, mu1_deg) < 10.0,
            "component 0: recovered mu {} too far from true mu {}",
            comps[0].mu_deg,
            mu1_deg
        );
        assert!(
            circular_diff_deg(comps[1].mu_deg, mu2_deg) < 10.0,
            "component 1: recovered mu {} too far from true mu {}",
            comps[1].mu_deg,
            mu2_deg
        );
        assert!(
            (comps[0].weight - expected_weight1).abs() < 0.08,
            "component 0: recovered weight {} too far from expected {}",
            comps[0].weight,
            expected_weight1
        );
        assert!(
            (comps[1].weight - expected_weight2).abs() < 0.08,
            "component 1: recovered weight {} too far from expected {}",
            comps[1].weight,
            expected_weight2
        );
    }

    #[test]
    fn fit_of_empty_histogram_is_empty() {
        let hist = TorsionHistogram::from_angles_deg(&[], 36);
        let fit = fit_von_mises_mixture(&hist, 2, 10);
        assert!(fit.components.is_empty());
    }

    // ── Jensen-Shannon distance properties ──────────────────────────────────

    fn single_component_fit(mu_deg: f64, kappa: f64) -> TorsionProfileFit {
        TorsionProfileFit {
            components: vec![VonMisesComponent {
                weight: 1.0,
                mu_deg,
                kappa,
            }],
        }
    }

    #[test]
    fn distance_to_self_is_zero() {
        let fit = single_component_fit(30.0, 3.0);
        let d = torsion_profile_distance(&fit, &fit);
        assert!(d.abs() < 1e-9, "distance(a,a) = {d}, expected ~0");
    }

    #[test]
    fn distance_is_symmetric() {
        let a = single_component_fit(0.0, 2.0);
        let b = single_component_fit(90.0, 5.0);
        let d_ab = torsion_profile_distance(&a, &b);
        let d_ba = torsion_profile_distance(&b, &a);
        assert!(
            (d_ab - d_ba).abs() < 1e-9,
            "JS divergence must be symmetric: {d_ab} vs {d_ba}"
        );
    }

    #[test]
    fn distance_is_bounded_by_ln2() {
        // Two well-separated, tightly concentrated distributions -> near the
        // upper bound of ln(2).
        let a = single_component_fit(0.0, 50.0);
        let b = single_component_fit(180.0, 50.0);
        let d = torsion_profile_distance(&a, &b);
        assert!(d >= 0.0, "JS divergence must be non-negative, got {d}");
        assert!(
            d <= std::f64::consts::LN_2 + 1e-9,
            "JS divergence must be bounded by ln(2), got {d}"
        );
        assert!(
            d > 0.5,
            "well-separated concentrated distributions should be near the ln(2) bound, got {d}"
        );
    }

    // ── Deterministic JSON export ───────────────────────────────────────────

    #[test]
    fn json_export_is_byte_identical_for_same_input() {
        let motif = TorsionMotif {
            atoms: [AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)],
            elements: ["C", "C", "C", "C"],
            aromatic: [false, false, false, false],
        };
        let hist = TorsionHistogram::from_angles_deg(&[10.0, -170.0, 45.0], 36);
        let fit = fit_von_mises_mixture(&hist, 2, 15);

        let json1 = torsion_profile_to_json(&motif, TorsionEnvironment::Water, &hist, &fit);
        let json2 = torsion_profile_to_json(&motif, TorsionEnvironment::Water, &hist, &fit);
        assert_eq!(json1, json2, "same inputs must produce byte-identical JSON");

        // Sanity: it's actually valid, non-trivial JSON, not an empty fallback.
        let parsed: serde_json::Value = serde_json::from_str(&json1).unwrap();
        assert_eq!(parsed["environment"], "water");
        assert_eq!(parsed["motif"]["signature"], "C-C-C-C");
    }

    #[test]
    fn environment_as_str_labels() {
        assert_eq!(TorsionEnvironment::Crystal.as_str(), "crystal");
        assert_eq!(TorsionEnvironment::Vacuum.as_str(), "vacuum");
        assert_eq!(TorsionEnvironment::Water.as_str(), "water");
        assert_eq!(TorsionEnvironment::Hexane.as_str(), "hexane");
    }
}
