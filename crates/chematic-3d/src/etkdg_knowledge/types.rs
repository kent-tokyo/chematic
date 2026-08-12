//! Data model for the ETKDG torsion-knowledge v2 layer (3D Breakthrough
//! Program, Wave 2, Agent E, `feat/3d-torsion-knowledge-v2`).
//!
//! This is a **new, separate** type family, not an extension of the legacy
//! `TorsionPreference` (single angle + linear penalty). See
//! `docs/rfcs/3d_torsion_knowledge_audit.md` §3.1/§3.3 for concrete evidence that
//! a single-angle model cannot represent real torsion potentials (biphenyl's
//! true multi-minima shape, the gauche effect's two symmetric minima).
//!
//! Not wired into `distance_geometry_v2.rs` or `etkdg.rs` by this PR --
//! Coordinator integration is a later, separate PR (see this crate's
//! `etkdg_knowledge.rs` module docs and `docs/rfcs/3d_breakthrough_master_plan.md`).

use chematic_core::AtomIdx;

// ---------------------------------------------------------------------------
// Provenance
// ---------------------------------------------------------------------------

/// Where a [`TorsionPotential`]'s parameters came from. Never conflate
/// `LegacyHeuristic` with the others -- see
/// `docs/rfcs/3d_torsion_knowledge_audit.md` for why the pre-existing module's
/// "experimental torsion angle preferences from CSD" header comment was an
/// unverified (and, on the evidence gathered, false) label this enum exists
/// to stop from recurring.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TorsionKnowledgeSource {
    /// Translated/adapted from RDKit's `torsionPreferences_v1.in` /
    /// `torsionPreferences_v2.in` (citing Riniker & Landrum, J. Chem. Inf.
    /// Model. 2015, 55, 2562-2574, doi:10.1021/acs.jcim.5b00654, per RDKit's
    /// own in-source header comment). See
    /// `validation/manifests/etkdg_torsion_knowledge_sources.json` for the
    /// exact fetched file, commit, and per-rule source line.
    StandardExperimental,
    /// Translated/adapted from RDKit's `torsionPreferences_smallrings.in`
    /// (Wang, Witek, Landrum, Riniker, J. Chem. Inf. Model. 2020, 60,
    /// 2044-2058, doi:10.1021/acs.jcim.0c00025).
    SmallRingExperimental,
    /// Translated/adapted from RDKit's `torsionPreferences_macrocycles.in`
    /// (same 2020 paper as above).
    MacrocycleAdaptation,
    /// Structural/textbook chemical knowledge not tied to a specific fitted
    /// data source (e.g. sp centers are linear by VSEPR, flat all-sp2 rings
    /// are planar). The flat-ring rule is translated from RDKit's
    /// `useBasicKnowledge` block in `TorsionPreferences.cpp`; other rules in
    /// this category are this crate's own, clearly marked per-rule.
    BasicChemicalKnowledge,
    /// This crate's pre-existing, hand-authored, uncited heuristic layer
    /// (`get_torsion_preference` / `SMARTS_TORSION_RULES`), kept behaviorally
    /// unchanged by this PR -- see `docs/rfcs/3d_torsion_knowledge_audit.md` §2
    /// for the honest per-rule reclassification. Lowest priority tier,
    /// opt-in only.
    LegacyHeuristic,
}

// ---------------------------------------------------------------------------
// Periodic potential
// ---------------------------------------------------------------------------

/// One term of a periodic torsion potential:
/// `E(phi) = amplitude * (1 + cos(periodicity * phi - phase))`.
///
/// This is RDKit's own `V*(1 + s*cos(n*x))` form (see
/// `Code/GraphMol/ForceFieldHelpers/CrystalFF/TorsionAngleM6.cpp`, fetched
/// and hashed in the sources manifest), re-expressed with `phase_deg` instead
/// of a `±1` sign: `s=+1` maps to `phase_deg=0.0`, `s=-1` maps to
/// `phase_deg=180.0` (algebraically identical: `cos(n*x - 180deg) ==
/// -cos(n*x)` for any integer `n`, so the sign flip and the 180-degree phase
/// shift are the same transformation regardless of periodicity).
///
/// `amplitude` may be **negative** -- RDKit's own fitted `V` coefficients
/// are sometimes negative (e.g. the unsubstituted-biphenyl term translated
/// into `rules_standard.rs` has `V1=-0.7`, `V2=-8.0`), which is preserved
/// verbatim here rather than normalized away, to stay faithful to the
/// translated data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FourierTorsionTerm {
    /// `n` in `cos(n*phi)`. RDKit's own data uses `1..=6`; not enforced as a
    /// hard range here (Rust `u8` already bounds it generously) but every
    /// rule this PR ships uses a value in `1..=6`.
    pub periodicity: u8,
    /// Phase offset in degrees, conventionally `0.0` or `180.0` (see the
    /// sign-to-phase mapping above). Not restricted to those two values --
    /// a future rule could use an intermediate phase -- but none of this
    /// PR's translated rules do.
    pub phase_deg: f64,
    /// Force-constant-like coefficient (may be negative -- see above).
    pub amplitude: f64,
}

impl FourierTorsionTerm {
    pub fn new(periodicity: u8, phase_deg: f64, amplitude: f64) -> Self {
        Self {
            periodicity,
            phase_deg,
            amplitude,
        }
    }

    /// Build from RDKit's own `(periodicity, sign, V)` triple -- `sign` must
    /// be `1` or `-1` (RDKit's data files only ever use these two values;
    /// anything else is a translation bug in the caller, not a valid input,
    /// so this asserts rather than silently coercing).
    pub fn from_rdkit(periodicity: u8, sign: i8, amplitude: f64) -> Self {
        assert!(
            sign == 1 || sign == -1,
            "RDKit torsion sign must be +1 or -1, got {sign}"
        );
        let phase_deg = if sign == 1 { 0.0 } else { 180.0 };
        Self::new(periodicity, phase_deg, amplitude)
    }

    /// `E(phi) = amplitude * (1 + cos(periodicity*phi - phase))`, `phi` in
    /// degrees. Always finite for finite inputs (no division, no sqrt).
    pub fn energy(&self, phi_deg: f64) -> f64 {
        let n = f64::from(self.periodicity);
        let phi_rad = phi_deg.to_radians();
        let phase_rad = self.phase_deg.to_radians();
        self.amplitude * (1.0 + (n * phi_rad - phase_rad).cos())
    }

    /// Analytic `dE/dphi` (degrees^-1, since `phi` is in degrees):
    /// `d/dphi [A(1+cos(n*phi_rad - phase))] = -A*n*sin(n*phi_rad-phase) * (pi/180)`.
    /// Always finite for finite inputs. Verified against central-difference
    /// finite differences in this module's tests (spec §7's "derivative
    /// that is finite and verifiable" requirement).
    pub fn d_energy_d_phi_deg(&self, phi_deg: f64) -> f64 {
        let n = f64::from(self.periodicity);
        let phi_rad = phi_deg.to_radians();
        let phase_rad = self.phase_deg.to_radians();
        -self.amplitude * n * (n * phi_rad - phase_rad).sin() * (std::f64::consts::PI / 180.0)
    }
}

// ---------------------------------------------------------------------------
// Torsion potential (one central bond, possibly multiple Fourier terms)
// ---------------------------------------------------------------------------

/// A torsion potential for one A-B-C-D dihedral, where B-C is the central
/// (rotatable) bond. `terms` may hold more than one [`FourierTorsionTerm`] --
/// this is the type-level fix for the legacy `TorsionPreference`'s inability
/// to represent multi-modal potentials (spec §3, audit doc §3.1).
#[derive(Clone, Debug)]
pub struct TorsionPotential {
    pub atoms: [AtomIdx; 4],
    pub central_bond: (AtomIdx, AtomIdx),
    pub source: TorsionKnowledgeSource,
    /// Stable identifier for this specific rule (e.g.
    /// `"standard:v2.in:142:secondary_amide"`). Used for provenance tracking,
    /// conflict diagnostics, and the gap-check example's per-rule reporting.
    pub rule_id: String,
    pub terms: Vec<FourierTorsionTerm>,
    /// `Some(n)` when this potential's rule is specific to a ring of size
    /// `n` (small-ring or macrocycle rules); `None` for acyclic/basic/legacy
    /// rules.
    pub ring_size: Option<usize>,
}

impl TorsionPotential {
    /// Total energy at dihedral angle `phi_deg`, summed across all terms.
    pub fn energy(&self, phi_deg: f64) -> f64 {
        self.terms.iter().map(|t| t.energy(phi_deg)).sum()
    }

    /// Total `dE/dphi` (degrees^-1) at `phi_deg`, summed across all terms.
    pub fn d_energy_d_phi_deg(&self, phi_deg: f64) -> f64 {
        self.terms
            .iter()
            .map(|t| t.d_energy_d_phi_deg(phi_deg))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Config (mirrors EmbedParameters' 4 reserved fields exactly)
// ---------------------------------------------------------------------------

/// Which torsion-knowledge rule families are enabled.
///
/// `use_exp_torsions`, `use_small_ring_torsions`, `use_macrocycle_torsions`,
/// and `use_macrocycle_14_bounds` mirror `distance_geometry_v2::EmbedParameters`'
/// 4 reserved fields **exactly** (same names, same semantics) so Coordinator's
/// future wiring is a direct field-for-field mapping, not a translation --
/// per the Wave 2 spec's explicit instruction. `include_legacy_heuristic` is
/// a 5th, crate-local addition: `EmbedParameters` reserves no field for the
/// legacy heuristic layer (it predates this PR and was always
/// unconditionally available through its own API), so this flag has no
/// `EmbedParameters` counterpart to mirror.
///
/// Field-to-tier mapping (see `matcher.rs` for the full 6-tier priority
/// order):
/// - `use_exp_torsions` gates tier 1 (specific validated SMARTS, currently
///   empty -- see PR body), tier 4 (standard experimental), and tier 5
///   (basic chemical knowledge) -- `EmbedParameters` reserves exactly these
///   4 fields and no separate `use_basic_knowledge`, so this implementation
///   folds basic-knowledge under the same flag RDKit's own `useExpTorsions`
///   is most commonly enabled alongside (see sources manifest).
/// - `use_small_ring_torsions` gates tier 2 (small-ring).
/// - `use_macrocycle_torsions` gates tier 3 (macrocycle-specific torsion
///   potentials).
/// - `use_macrocycle_14_bounds` gates `macrocycle_14_bound_adjustments()`
///   only -- a separate API, not a torsion-potential tier at all (mirrors
///   RDKit's own `useMacrocycle14config`, which is documented in
///   `Embedder.h` as feeding `BoundsMatrixBuilder`'s 1-4 distance bounds,
///   not the ExpTorsion Fourier library).
/// - `include_legacy_heuristic` gates tier 6 (legacy heuristic opt-in).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TorsionKnowledgeConfig {
    pub use_exp_torsions: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_macrocycle_14_bounds: bool,
    pub include_legacy_heuristic: bool,
}

impl Default for TorsionKnowledgeConfig {
    /// All flags false -- matches `EmbedParameters::default()`'s 4 mirrored
    /// fields (also all false) and is the all-flags-false no-op state spec
    /// §11/§15 require (see `matcher.rs::build_torsion_knowledge`'s early
    /// return and `tests/torsion_knowledge_negative_controls.rs`).
    fn default() -> Self {
        Self {
            use_exp_torsions: false,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle_14_bounds: false,
            include_legacy_heuristic: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

/// Why a bond's torsion-knowledge lookup produced a diagnostic instead of
/// (or in addition to) a clean single potential.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TorsionKnowledgeDiagnosticKind {
    /// Two or more rules at the *same* priority tier matched the same
    /// central bond with genuinely conflicting (non-equivalent,
    /// non-composable) potentials. Spec §4: "Never arbitrarily pick one
    /// side of a genuine ambiguity."
    AmbiguousSameTierConflict,
    /// A rule's SMARTS failed to parse. The legacy
    /// `build_smarts_torsion_map` silently `continue`s on this (audit doc
    /// §3.7) -- this implementation instead always records it here.
    SmartsParseFailure,
    /// The bond spans more than one ring whose sizes fall in different
    /// small-ring/macrocycle buckets, or belongs to more than one ring of
    /// any size (fused/bridged/spiro). Spec §5: "do not naively use only
    /// the first SSSR ring's size."
    FusedOrBridgedRingBoundary,
    /// The bond was classified as non-torsional (terminal, double/triple,
    /// or otherwise not a genuine rotatable single bond) and was
    /// deliberately skipped rather than scored.
    NonRotatableBondSkipped,
}

/// One diagnostic record: which bond, what kind of issue, human-readable
/// detail, and (for conflicts) which rule ids were involved.
#[derive(Clone, Debug)]
pub struct TorsionKnowledgeDiagnostic {
    pub central_bond: (AtomIdx, AtomIdx),
    pub kind: TorsionKnowledgeDiagnosticKind,
    pub message: String,
    pub candidate_rule_ids: Vec<String>,
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

/// Full result of matching torsion knowledge against a molecule: every
/// resolved potential, which rules fired, which rotatable bonds matched
/// nothing, and every diagnosed ambiguity/skip -- never silently dropped.
#[derive(Clone, Debug, Default)]
pub struct TorsionKnowledgeReport {
    pub potentials: Vec<TorsionPotential>,
    pub matched_rule_ids: Vec<String>,
    pub unmatched_rotatable_bonds: Vec<(AtomIdx, AtomIdx)>,
    pub ambiguous_matches: Vec<TorsionKnowledgeDiagnostic>,
    pub skipped_bonds: Vec<TorsionKnowledgeDiagnostic>,
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Typed failure modes for the fallible APIs in this module family. Never a
/// silent `None`/default -- every failure mode spec §7/§13 calls out by name
/// has a variant here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TorsionKnowledgeError {
    /// An atom index in a requested torsion/pair does not exist in `mol`,
    /// or the central bond does not actually exist.
    InvalidTopology,
    /// `coords` has a different atom count than `mol`.
    CoordsAtomCountMismatch,
    /// `optimize_torsions` did not converge within
    /// `TorsionOptimizationConfig::max_iterations`. Spec §7: "non-convergence
    /// is a typed failure, not silent."
    NonConvergence,
    /// The optimizer's own internal ring-closure check detected that a
    /// rotation opened a ring (a bonded pair that should stay within
    /// bond-length range moved outside it, or a ring-closure pair's
    /// distance changed by more than the configured tolerance). Spec §7:
    /// "must never break a ring open."
    RingIntegrityViolated,
    /// The optimizer detected non-finite (`NaN`/`Inf`) energy or gradient
    /// during optimization. Spec §13: "must never accept NaN/Inf energy."
    NonFiniteEnergy,
}
