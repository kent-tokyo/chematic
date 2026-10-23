/// Force-field policy for [`minimize_with_policy`]/[`minimize_with_policy_gated`]
/// — the opt-in bridge to chematic-ff's complete implementations. Does not
/// change any existing public function's behavior or default in this crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForceFieldPolicy {
    /// Full MMFF94 (bond + angle + stretch-bend + torsion + out-of-plane +
    /// vdW + electrostatic, `chematic_ff::minimize_mmff94_lbfgs`). Refuses
    /// (typed `Err`) if even one BOND or ANGLE lacks MMFF94 parameters,
    /// rather than silently contributing zero energy/gradient for that
    /// internal coordinate — this is the exact scope named in this variant
    /// (`BondAngle`, not "every MMFF94 term"): torsion/out-of-plane coverage
    /// is always measured and reported (`coverage`/
    /// `missing_parameter_classes` on the result) but does NOT by itself
    /// trigger a refusal here. Renamed from a plain `Mmff94Strict` (found in
    /// independent review to be a scope/naming mismatch: the old name
    /// implied gating on every required MMFF94 term, but the gate only ever
    /// checked bond+angle). Use `minimize_with_policy_gated(..., true, false)` to
    /// also gate on torsion/out-of-plane, or `(..., false, true)` to also gate
    /// on stretch-bend (Priority 2 / Stage 1B).
    ///
    /// Naming/scope note for Coordinator: post-chematic-ff-#183 (bond/angle/
    /// torsion classification fixes), the 58-molecule corpus example
    /// (`examples/mmff94_bridge_coverage_report.rs`) now measures 32/58
    /// molecules pass this bond+angle-only gate vs. 31/58 under the widened
    /// (+torsion+oop) gate — nearly identical, unlike the pre-#183 numbers
    /// that originally justified keeping torsion/oop out of the gate ("fails
    /// a large fraction of ordinary organic molecules"). That justification
    /// is now stale. Widening the default gate is a small, well-scoped
    /// follow-up worth considering, but doing so also changes
    /// `Mmff94WithUffFallback`'s fallback trigger population and needs its
    /// own re-verification pass — deliberately not done in this PR to avoid
    /// conflating a naming fix with a behavior change (see PR body).
    Mmff94BondAngleStrict,
    /// Try `Mmff94BondAngleStrict` first; on any typed failure (unsupported atom
    /// element, missing bond/angle parameters, or an unsound MMFF94 result —
    /// see [`ForceFieldBridgeError::MinimizationFailed`]), fall back to
    /// `UffOnly`. When a fallback happens, it is always reported via
    /// [`PolicyMinimizeResult::fallback_reason`]/
    /// [`PolicyMinimizeResult::actual_force_field_used`]/
    /// [`PolicyMinimizeResult::missing_parameter_classes`], never silent.
    ///
    /// NOT infallible: if the `UffOnly` fallback attempt is itself unsound
    /// (measured post-chematic-ff-#183, restricted to the molecules that
    /// actually *reach* this fallback path — i.e. the ones that first fail
    /// `Mmff94BondAngleStrict`: 8 fused/conjugated polycyclic aromatics
    /// (naphthalene, quinoline, pyrene, ibuprofen, ibuprofen_S, naproxen_S,
    /// diphenhydramine, atorvastatin_fragment) plus caffeine remaining
    /// blown up even before the soundness gate — see PR #169 body), this
    /// returns `Err(MinimizationFailed)` too. This policy's contract is
    /// "never silently report success on an unsound geometry," not "always
    /// succeeds" — those are different guarantees, and only the first one
    /// is made here. This 9-molecule figure is scoped to the fallback
    /// *trigger population* specifically, not `UffOnly`'s behavior across
    /// the full corpus — see the correction on [`ForceFieldPolicy::UffOnly`]
    /// below, which a prior version of this doc comment incorrectly copied
    /// verbatim onto both variants.
    Mmff94WithUffFallback,
    /// chematic-ff's real UFF module (`chematic_ff::uff`) — generic,
    /// all-element coverage (bond lengths/angles are formula-derived from
    /// per-type constants, not a lookup table, so there is no missing-entry
    /// case to gate on). NOT infallible: can return
    /// `Err(MinimizationFailed)` if the resulting geometry is unsound (see
    /// `check_minimization_soundness`).
    ///
    /// **Correction (found during the Wave-1 C+F integration smoke test,
    /// PR #186, independently re-verified):** an earlier version of this
    /// doc comment claimed "8/58 corpus molecules, all fused/conjugated
    /// polycyclic aromatics" for `UffOnly` specifically. That figure was
    /// actually the [`ForceFieldPolicy::Mmff94WithUffFallback`] fallback
    /// *trigger population* (9 molecules, not 8, once caffeine is counted;
    /// see that variant's doc above), copied here without re-deriving it —
    /// `UffOnly` was never actually run over the full 58-molecule corpus in
    /// PR #169 itself. The real, first-ever full-corpus measurement (from
    /// PR #186, using legacy `dg::generate_coords` starting geometry, the
    /// same convention this PR's own gate-check example uses) is **17/58**
    /// blow up under `UffOnly`: the same 9 fused/conjugated aromatics above,
    /// plus 8 more with no ring fusion at all (hexane, decane,
    /// triethylene_glycol, hexanediol, hexadecane, penicillin_core,
    /// testosterone, cholesterol) that are structurally invisible to any
    /// `Mmff94WithUffFallback`-based measurement, since all 8 pass
    /// `Mmff94BondAngleStrict` cleanly and so never invoke UFF under that
    /// policy at all. Tracked as a `chematic-ff` UFF-minimizer robustness
    /// gap; see issue #185 (currently scoped only to the narrower
    /// naphthalene-vs-anthracene puzzle, not this full 17-molecule class).
    UffOnly,
    /// Existing DREIDING path (unchanged physics) — chematic-ff has no
    /// DREIDING minimizer, so this stays chematic-3d's own implementation,
    /// just routed through the same reporting result type as the other
    /// policies.
    Dreiding,
    /// No minimization — return the embedded geometry unchanged. Useful for
    /// testing/composability with a raw embedder.
    None,
}

/// Which MMFF94 internal-coordinate class a missing-parameter finding is
/// about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mmff94TermKind {
    Bond,
    Angle,
    Torsion,
    Oop,
    /// Stretch-bend cross term (Halgren MMFF.V eq. 4), keyed on the same
    /// (outer, center, outer) angle triple as [`Mmff94TermKind::Angle`].
    /// Historically never a member of this enum at all -- see
    /// `Mmff94CoverageReport::has_gate_failure`'s `include_stretch_bend`
    /// parameter for why it is opt-in, not always-gated like bond/angle.
    StretchBend,
}

/// One internal coordinate that chematic-ff's MMFF94 tables have no entry
/// for. `atoms` is the exact atom tuple (2 for bond, 3 for angle, 4 for
/// torsion/oop) so a caller can cite *which* atoms/elements were missing —
/// never just an aggregate count.
#[derive(Debug, Clone)]
pub struct Mmff94MissingTerm {
    pub kind: Mmff94TermKind,
    pub atoms: Vec<AtomIdx>,
    /// Human-readable citation, including each atom's specific numeric
    /// MMFF94 type (not just its element symbol — many distinct MMFF94
    /// types share an element, e.g. sp3 vs. aromatic vs. carbonyl carbon),
    /// e.g. `"Angle F(11)-C(1)-Cl(12) (atom indices [1, 0, 2])"`. Atom order
    /// is canonicalized (see `canonicalize_term_atoms`) so physically
    /// equivalent citations (e.g. a bond read C→N vs N→C) always produce the
    /// same string, so aggregating by this string doesn't double-count.
    pub description: String,
}

/// Full parameter-coverage measurement for one molecule under chematic-ff's
/// MMFF94 tables. Every internal coordinate chematic-ff's own energy
/// functions would evaluate is enumerated here (the enumeration mirrors
/// `chematic_ff::mmff94_minimizer`'s private energy loops exactly, so the
/// totals match what chematic-ff will actually compute over, not an
/// idealized recount) — nothing is silently excluded from the count.
#[derive(Debug, Clone, Default)]
pub struct Mmff94CoverageReport {
    pub bonds_total: usize,
    pub bonds_missing: Vec<Mmff94MissingTerm>,
    pub angles_total: usize,
    pub angles_missing: Vec<Mmff94MissingTerm>,
    pub torsions_total: usize,
    pub torsions_missing: Vec<Mmff94MissingTerm>,
    /// Torsions with NO real MMFF94 term at all, by design, not a coverage
    /// gap -- both central (j-k) atoms fully typed, but one has MMFF's
    /// `lin` flag (see `chematic_ff::torsion_no_term_by_design`'s doc, issue
    /// #227 Phase 1). Included in `torsions_total`, deliberately excluded
    /// from `torsions_missing` and from `has_gate_failure`'s torsion check.
    pub torsions_no_term_by_design: usize,
    /// Number of eligible trigonal centers, not Wilson-angle contributions.
    /// Each resolved center supplies three angles sharing one parameter lookup.
    pub oop_total: usize,
    pub oop_missing: Vec<Mmff94MissingTerm>,
    /// Stretch-bend cross-term coverage (Priority 2 / Stage 1B, issue #227).
    /// Never gated by default -- see `has_gate_failure`'s `include_stretch_bend`
    /// parameter -- but always measured/reported here regardless, same as
    /// torsion/oop.
    pub stretch_bend_total: usize,
    pub stretch_bend_missing: Vec<Mmff94MissingTerm>,
}

impl Mmff94CoverageReport {
    /// True iff no bond or angle is missing parameters — the exact scope of
    /// the RFC's mechanism-3 bug (silent zero-gradient, no restoring force
    /// while vdW repulsion still pushes atoms apart).
    pub fn bond_angle_fully_covered(&self) -> bool {
        self.bonds_missing.is_empty() && self.angles_missing.is_empty()
    }

    /// Gate check used by `Mmff94BondAngleStrict`/`Mmff94WithUffFallback`.
    /// Bond/angle always participate; torsion/out-of-plane only do when
    /// `include_torsion_oop` is set, and stretch-bend only when
    /// `include_stretch_bend` is set (see `minimize_with_policy_gated`) --
    /// each is its own independent opt-in, not bundled together, so a
    /// caller can measure the stretch-bend-only effect in isolation (Priority
    /// 2 / Stage 1B's "legacy vs. complete-term" comparison).
    pub fn has_gate_failure(&self, include_torsion_oop: bool, include_stretch_bend: bool) -> bool {
        !self.bond_angle_fully_covered()
            || (include_torsion_oop
                && (!self.torsions_missing.is_empty() || !self.oop_missing.is_empty()))
            || (include_stretch_bend && !self.stretch_bend_missing.is_empty())
    }

    /// Total count of missing internal coordinates across all 5 classes —
    /// always meaningful regardless of gate scope, since nothing here is
    /// filtered by what the gate happens to check.
    pub fn total_missing(&self) -> usize {
        self.bonds_missing.len()
            + self.angles_missing.len()
            + self.torsions_missing.len()
            + self.oop_missing.len()
            + self.stretch_bend_missing.len()
    }

    /// Flattened list of every missing internal coordinate across all 5
    /// classes — the exact "which specific element/atom-type pairs lacked
    /// coverage" citation surfaced at the top level as
    /// [`PolicyMinimizeResult::missing_parameter_classes`].
    pub fn all_missing(&self) -> Vec<Mmff94MissingTerm> {
        let mut v = Vec::with_capacity(self.total_missing());
        v.extend(self.bonds_missing.iter().cloned());
        v.extend(self.angles_missing.iter().cloned());
        v.extend(self.torsions_missing.iter().cloned());
        v.extend(self.oop_missing.iter().cloned());
        v.extend(self.stretch_bend_missing.iter().cloned());
        v
    }
}

/// Typed failure for the [`ForceFieldPolicy::Mmff94BondAngleStrict`] bridge. Never
/// silently absorbed into a zero-energy/zero-gradient term.
///
/// Naming note for Coordinator (Wave 2 reconciliation): this is deliberately
/// named for this crate's own use, not coupled to Agent C's parallel
/// `EmbedFailureCause` (distance-geometry embedding failures — a different
/// pipeline stage). The variant `MissingParameters(Mmff94CoverageReport)`
/// here is the natural candidate for an `EmbedFailureCause::
/// MissingForceFieldParameters { .. }` case if/when Coordinator unifies the
/// two vocabularies; this PR does not attempt that unification itself.
#[derive(Debug, Clone)]
pub enum ForceFieldBridgeError {
    /// `chematic_ff::assign_mmff94_numeric_types` failed outright for one or
    /// more atoms (element unsupported by the typer at all — coarser than a
    /// missing bond/angle/torsion/oop entry for an otherwise-typed atom).
    UnsupportedAtomType(String),
    /// One or more internal coordinates lack MMFF94 parameters under the
    /// active coverage gate. Boxed: `clippy::result_large_err` (the report
    /// carries per-missing-term `Vec`s / `String`s so the inline variant is
    /// large relative to the common `Ok` path).
    MissingParameters(Box<Mmff94CoverageReport>),
    /// The minimizer ran (coverage/typing succeeded) but produced a
    /// geometry that is not acceptable — NaN/Inf coordinates, a
    /// catastrophic bond blow-up, an excessive residual force, or a declared
    /// stereo contradiction. Before
    /// this variant existed, every one of these cases was returned as
    /// `Ok(PolicyMinimizeResult { converged: false, .. })` — a result a
    /// careless caller could mistake for success (independent review found
    /// exactly this on the 58-molecule corpus: UFF-fallback runs that blew
    /// worst-bond-length past 800–19,000+ Å still came back `Ok`). See
    /// [`check_minimization_soundness`] for what specifically triggers this
    /// and why plain non-convergence (`converged == false`) alone does
    /// NOT — that would be a false-failure trigger, not a safety gate (see
    /// that function's doc for the measured evidence).
    MinimizationFailed(Box<MinimizationFailureDetail>),
}

/// Which soundness check tripped inside [`check_minimization_soundness`].
#[derive(Debug, Clone, Copy)]
pub enum MinimizationFailureReason {
    /// One or more final coordinates are NaN or infinite.
    NonFiniteCoordinates,
    /// The worst bond length in the final geometry exceeds
    /// [`MAX_SANE_BOND_LENGTH`] Å — no legitimate covalent bond in the
    /// molecules this bridge handles gets remotely close to this.
    CatastrophicBondBlowup,
    /// The finite-difference max |gradient component| at the final geometry
    /// exceeds [`MAX_SANE_RESIDUAL_FORCE`] kcal/mol/Å — a backstop for a
    /// geometry that isn't bond-length-blown-up but is still nowhere near a
    /// force-balanced minimum.
    ExcessiveResidualForce,
    /// The geometry is finite and bond-length-sound, but contradicts at least
    /// one tetrahedral or E/Z declaration on the input molecule. This is used
    /// by the stereo-aware UFF bridge so a legacy starting geometry can never
    /// be returned as a successful minimization after silently inverting a
    /// declared stereo element.
    StereoConstraintViolation,
}

/// Evidence attached to [`ForceFieldBridgeError::MinimizationFailed`].
/// `converged`/`iterations` are carried here for diagnostics but are not by
/// themselves what triggered the failure; `reason` identifies the geometry
/// or stereo gate that did.
#[derive(Debug, Clone)]
pub struct MinimizationFailureDetail {
    pub policy: ForceFieldPolicy,
    pub reason: MinimizationFailureReason,
    /// The underlying minimizer's own convergence flag. Often `false` even
    /// on perfectly sound geometries that simply didn't reach a tight
    /// gradient tolerance within the default iteration budget — see
    /// `check_minimization_soundness`'s doc; this field is diagnostic only.
    pub converged: bool,
    pub iterations: usize,
    /// Optimizer stop condition for an MMFF94 attempt; absent for other
    /// force fields. The separate `reason` is the quality gate that failed.
    pub mmff94_termination: Option<Mmff94TerminationReason>,
    /// Always populated (not just when `reason` is
    /// `ExcessiveResidualForce`) so a failed molecule still reports how far
    /// from equilibrium it is.
    pub max_residual_force: f64,
    /// Always populated (not just when `reason` is `CatastrophicBondBlowup`)
    /// so a failed molecule still reports its geometry, rather than a blank
    /// — this is what stops the "measurement trap" of a fixed bug simply
    /// vanishing from a blow-up count instead of being counted as a typed
    /// failure.
    pub worst_bond_length: f64,
    /// This `MinimizationFailureDetail` always describes the ORIGINAL attempt's
    /// failure (the caller's own coordinates) — never the retry's. This flag only
    /// records whether a retry was even attempted before that original failure was
    /// returned: `true` means the original attempt failed with catastrophic/non-
    /// finite geometry or a declared-stereo violation,
    /// `embed_distance_geometry_v2` was tried as a rescue, and that rescue was **not
    /// accepted** (its own result was unsound, stereo-violating, or embedding itself
    /// errored). `false` means either no rescue was attempted at all (wrong policy,
    /// or an `ExcessiveResidualForce`-only failure) or a rescue was accepted, in
    /// which case this function returns `Ok` and this field is moot. Lets a caller
    /// distinguish "no rescue was even attempted" from "the rescue was attempted
    /// and did not help either" — never "the retry's own failure reason," which
    /// this field does not carry.
    pub distance_geometry_v2_retry_attempted: bool,
}

/// Which starting geometry actually seeded a successful UFF minimization —
/// see [`run_uff_bridge`]'s doc for the full rationale. Exists so a caller
/// can tell "my exact input coordinates were used" apart from "the input was
/// detected as catastrophically unsound after minimizing it, and
/// `embed_distance_geometry_v2` was used to re-seed instead" — never a
/// silent substitution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UffStartingGeometry {
    /// The coordinates the caller passed in were used to seed minimization
    /// and that attempt was already sound and consistent with every declared
    /// stereo element — no re-embedding was attempted.
    AsProvided,
    /// The caller-provided coordinates produced an unsound or declared-stereo-
    /// violating result; `embed_distance_geometry_v2` was used to re-seed
    /// instead, and that retry was sound and stereo-satisfied.
    ReplacedWithDistanceGeometryV2,
}

impl std::fmt::Display for ForceFieldBridgeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ForceFieldBridgeError::UnsupportedAtomType(e) => {
                write!(f, "MMFF94 atom-type assignment failed: {e}")
            }
            ForceFieldBridgeError::MissingParameters(r) => write!(
                f,
                "MMFF94 parameters missing for {} internal coordinate(s) \
                 ({} bond, {} angle, {} torsion, {} oop, {} stretch-bend)",
                r.total_missing(),
                r.bonds_missing.len(),
                r.angles_missing.len(),
                r.torsions_missing.len(),
                r.oop_missing.len(),
                r.stretch_bend_missing.len(),
            ),
            ForceFieldBridgeError::MinimizationFailed(d) => write!(
                f,
                "{:?} minimization produced an unacceptable geometry: {:?} \
                 (converged={}, iterations={}, worst_bond_length={:.2} Å, \
                 max_residual_force={:.2} kcal/mol/Å)",
                d.policy,
                d.reason,
                d.converged,
                d.iterations,
                d.worst_bond_length,
                d.max_residual_force,
            ),
        }
    }
}

impl std::error::Error for ForceFieldBridgeError {}

impl From<NumericTypeError> for ForceFieldBridgeError {
    fn from(e: NumericTypeError) -> Self {
        ForceFieldBridgeError::UnsupportedAtomType(e.to_string())
    }
}

impl From<MinimizerError> for ForceFieldBridgeError {
    fn from(e: MinimizerError) -> Self {
        // Preserve the detailed setup/calculation failure rather than
        // silently substituting an approximation in the MMFF94 layer.
        ForceFieldBridgeError::UnsupportedAtomType(e.to_string())
    }
}

/// Per-policy energy accounting. `Mmff94` reuses chematic-ff's own
/// [`EnergyBreakdown`] shape verbatim — the same struct
/// `scripts/diagnostics/etkdg_vs_rdkit_gap.py`'s "fair energy delta" metric already
/// consumes via `Mol.mmff94_energy_breakdown()` — so downstream consumers
/// never need a second energy-breakdown format.
#[derive(Debug, Clone)]
pub enum EnergyReport {
    Mmff94(EnergyBreakdown),
    Uff { total: f64 },
    Dreiding { total: f64 },
    None,
}

impl EnergyReport {
    pub fn total(&self) -> f64 {
        match self {
            EnergyReport::Mmff94(b) => b.total,
            EnergyReport::Uff { total } | EnergyReport::Dreiding { total } => *total,
            EnergyReport::None => 0.0,
        }
    }
}

/// Full result of [`minimize_with_policy`]/[`minimize_with_policy_gated`].
///
/// The `requested_force_field`/`actual_force_field_used`/`fallback_reason`/
/// `missing_parameter_classes` quartet exists specifically so that
/// `Mmff94WithUffFallback` can never silently report "MMFF94 minimized" when
/// it actually fell back to UFF — that would be the exact mechanism-3 bug
/// (silent zero/substitution) reintroduced one layer up, at the
/// policy-reporting level instead of the energy-term level.
#[derive(Debug, Clone)]
pub struct PolicyMinimizeResult {
    pub coords: Coords3D,
    /// What the caller asked for.
    pub requested_force_field: ForceFieldPolicy,
    /// What actually ran. For `Mmff94BondAngleStrict` this is always
    /// `Mmff94BondAngleStrict`; for `Mmff94WithUffFallback` this is
    /// `Mmff94BondAngleStrict` when the MMFF94 attempt succeeded (the common case —
    /// `requested_force_field == Mmff94WithUffFallback` but
    /// `actual_force_field_used == Mmff94BondAngleStrict` here is NOT a fallback, it's
    /// just the more informative "which physics actually ran" value) and
    /// `UffOnly` only when it fell back. **Check `fallback_reason.is_some()`
    /// to ask "did a fallback occur," not `actual_force_field_used !=
    /// requested_force_field`** — the latter is true on every successful
    /// `Mmff94WithUffFallback` call and does not mean a fallback happened.
    pub actual_force_field_used: ForceFieldPolicy,
    /// `Some(reason)` iff a fallback actually occurred (only possible under
    /// `Mmff94WithUffFallback`, when the MMFF94 attempt failed) — explains
    /// why. Carries the same typed reason `Mmff94BondAngleStrict` would have returned
    /// as an `Err`. This, not `actual_force_field_used != requested_force_field`,
    /// is the correct fallback-occurred check (see `actual_force_field_used`'s
    /// doc above).
    pub fallback_reason: Option<ForceFieldBridgeError>,
    /// Every specific internal coordinate (bond/angle/torsion/oop, cited by
    /// atom indices and element symbols — never just an aggregate count)
    /// that lacked MMFF94 parameter coverage. Populated whenever MMFF94
    /// typing/coverage was computed at all (`Mmff94BondAngleStrict`, or the MMFF94
    /// attempt inside `Mmff94WithUffFallback` regardless of outcome); empty
    /// (not merely absent) when MMFF94 was never attempted (`UffOnly`/
    /// `Dreiding`/`None`) or when coverage was full.
    pub missing_parameter_classes: Vec<Mmff94MissingTerm>,
    /// Full structured coverage report (all 4 term-class totals + missing
    /// lists), whenever MMFF94 typing/coverage was computed — see
    /// `missing_parameter_classes` for the flattened, always-populated view
    /// of the same data.
    pub coverage: Option<Mmff94CoverageReport>,
    pub energy_before: EnergyReport,
    pub energy_after: EnergyReport,
    pub converged: bool,
    pub iterations: usize,
    /// Terminal condition when the returned coordinates came from MMFF94.
    /// `None` for another returned force field; if an MMFF94 attempt failed
    /// before UFF fallback, its stop reason is in `fallback_reason`.
    /// A non-converged result may be a constraint-rejected stop rather than
    /// an iteration-cap outcome.
    pub mmff94_termination: Option<Mmff94TerminationReason>,
    /// Max |gradient component| (kcal/mol/Å) at the final geometry, computed
    /// by this bridge's own finite-difference pass over the exact energy
    /// function that produced `energy_after` — not read off any upstream
    /// struct (chematic-ff's `MinimizeResult`/`UffMinimizeResult` don't carry
    /// this field).
    pub max_residual_force: f64,
    /// `Some(_)` whenever UFF actually ran (`UffOnly`, or the UFF-fallback
    /// path of `Mmff94WithUffFallback`) — see [`UffStartingGeometry`] for
    /// what the two variants mean and why raw starting energy alone can't
    /// predict which one occurs. `None` for `Dreiding`/`None`, and for
    /// `Mmff94BondAngleStrict`'s own (non-UFF) success path.
    pub starting_geometry: Option<UffStartingGeometry>,
}

fn trivial_result(coords: Coords3D, policy: ForceFieldPolicy) -> PolicyMinimizeResult {
    PolicyMinimizeResult {
        coords,
        requested_force_field: policy,
        actual_force_field_used: policy,
        fallback_reason: None,
        missing_parameter_classes: Vec::new(),
        coverage: None,
        energy_before: EnergyReport::None,
        energy_after: EnergyReport::None,
        converged: true,
        iterations: 0,
        mmff94_termination: None,
        max_residual_force: 0.0,
        starting_geometry: None,
    }
}

// --- Coords3D <-> Vec<[f64; 3]> bridge helpers ------------------------------

fn coords_to_vec(coords: &Coords3D, n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect()
}

fn vec_to_coords(v: &[[f64; 3]]) -> Coords3D {
    let mut c = Coords3D::new_zeroed(v.len());
    for (i, p) in v.iter().enumerate() {
        c.set(AtomIdx(i as u32), Point3::new(p[0], p[1], p[2]));
    }
    c
}

/// Central-difference max |gradient component| over an arbitrary black-box
/// energy function. Production MMFF94 uses the matching prepared analytic
/// gradient; this remains as an independent reference for tests and for force
/// fields that do not yet expose a validated analytic residual gradient.
fn fd_max_gradient<F: Fn(&[[f64; 3]]) -> f64>(
    coords: &[[f64; 3]],
    energy_fn: F,
    delta: f64,
) -> f64 {
    let mut work = coords.to_vec();
    let mut max_g = 0.0_f64;
    for i in 0..work.len() {
        for axis in 0..3 {
            work[i][axis] += delta;
            let ep = energy_fn(&work);
            work[i][axis] -= 2.0 * delta;
            let em = energy_fn(&work);
            work[i][axis] += delta;
            let g = (ep - em) / (2.0 * delta);
            if g.abs() > max_g {
                max_g = g.abs();
            }
        }
    }
    max_g
}

fn max_gradient_component(gradient: &[[f64; 3]]) -> f64 {
    gradient
        .iter()
        .flat_map(|v| v.iter())
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max)
}

// --- Geometric/energetic soundness gate -------------------------------------
//
// Found in independent review: every policy branch below used to build its
// `PolicyMinimizeResult` directly as `Ok(..)`, including `converged: false`
// results with `max_residual_force` in the hundreds of thousands to millions
// (kcal/mol/Å) and worst bond lengths past 800–19,000+ Å (measured on the
// 58-molecule corpus's UFF-fallback cases) — a caller checking only
// `result.is_ok()` would call that "minimized successfully." This section
// converts a genuinely unsound result into a typed
// `Err(ForceFieldBridgeError::MinimizationFailed)` instead.

/// Worst bond length ceiling (Å) above which a geometry is treated as a
/// catastrophic blow-up, not merely "not yet converged." Reuses this
/// project's own pre-existing convention (`worst_bond > 3.0` is the "blown
/// up" bar used throughout `examples/mmff94_bridge_coverage_report.rs` and
/// the tests in this file) rather than inventing a new number. No
/// light-atom-only organic bond covered by MMFF94/UFF (including C-I, S-S)
/// gets remotely close to 3 Å, so this has ample margin on the "sound" side.
///
/// `pub(crate)`: Wave 2/3 Coordinator integration (`pipeline_v2.rs`'s own,
/// independent final-geometry soundness check, per that program's spec §11) reuses
/// this exact constant rather than hand-copying `3.0` a second time -- particularly
/// important because `ForceFieldPolicy::None` skips `check_minimization_soundness`
/// entirely (see `minimize_with_policy_gated`'s `None` arm), so the pipeline-level
/// check is the *only* gate against a blown-up bond length in that case.
pub(crate) const MAX_SANE_BOND_LENGTH: f64 = 3.0;

/// Residual-force ceiling (kcal/mol/Å) above which a geometry is treated as
/// unsound even if no single bond individually crossed
/// [`MAX_SANE_BOND_LENGTH`] (e.g. a purely angular/torsional distortion).
/// Measured on the 58-molecule MMFF94/UFF corpus (post chematic-ff #183):
/// every molecule that is geometrically fine but simply didn't converge
/// within the default iteration budget has `max_residual_force` ≤ 8.93
/// (cholesterol); every molecule with a real blow-up has
/// `max_residual_force` ≥ 337.99 (quinoline) — but every one of those
/// already trips the bond-length check above, so within that corpus this
/// constant has no independently-discriminating case.
///
/// An initial value of 50.0 (picked inside that 8.93–337.99 gap with no
/// corpus case to validate it) turned out to be a live tripwire, not a
/// backstop: `dreiding_policy_matches_existing_behavior_and_reports_convergence`
/// (acetic acid via [`ForceFieldPolicy::Dreiding`]) has a perfectly normal
/// worst bond length (1.496 Å) but `max_residual_force` = 55.12 after 200
/// gradient-descent iterations — DREIDING's simpler harmonic springs
/// (`BOND_SPRING_CONSTANT`/`ANGLE_SPRING_CONSTANT`) apparently settle to a
/// higher steady-state FD residual than MMFF94/UFF's more carefully
/// parameterized terms do, even on a sound geometry. Raised to 200.0 to keep
/// this measured-sound case comfortably (>3.6×) below the ceiling while
/// staying well (>40%) under the smallest measured genuine blow-up
/// (337.99) — still a backstop with only one policy's data validating the
/// "sound" side, not a value with its own discriminating corpus case.
const MAX_SANE_RESIDUAL_FORCE: f64 = 200.0;

fn worst_bond_length_vec(mol: &Molecule, coords: &[[f64; 3]]) -> f64 {
    let dist = |a: [f64; 3], b: [f64; 3]| {
        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    mol.bonds()
        .map(|(_, b)| dist(coords[b.atom1.0 as usize], coords[b.atom2.0 as usize]))
        .fold(0.0_f64, f64::max)
}

struct OptimizerStatus {
    converged: bool,
    iterations: usize,
    mmff94_termination: Option<Mmff94TerminationReason>,
}

/// Decide whether a minimization result is geometrically/energetically
/// sound enough to return as `Ok`, or must instead become a typed
/// `Err(MinimizationFailed)`.
///
/// Deliberately does NOT treat `converged == false` alone as a failure
/// trigger. Measured on the 58-molecule corpus: molecules like decane,
/// hexadecane, triethylene glycol, testosterone, cholesterol, the
/// gly-ala-gly tripeptide, and the penicillin core all report
/// `converged == false` under the default `MinimizeConfig::max_steps` — not
/// because their geometry is unsound (every one has a finite, sub-3-Å worst
/// bond and `max_residual_force` under 10 kcal/mol/Å), but because L-BFGS's
/// gradient-norm convergence tolerance is tight relative to that default
/// iteration budget. Treating that as a hard failure would make
/// `Mmff94BondAngleStrict`/`UffOnly`/`Dreiding` spuriously refuse on a large
/// fraction of ordinary, perfectly fine molecules — a false-failure
/// generator, not a safety gate. `converged`/`iterations` are still carried
/// into `MinimizationFailureDetail` for diagnostics on an actual failure,
/// and every `Ok` result with `converged == false` remains visible on
/// `PolicyMinimizeResult` for a caller who wants to distinguish
/// "definitely converged" from "sound but still iterating" — see
/// `examples/mmff94_bridge_coverage_report.rs`'s own separate count of that.
fn check_minimization_soundness(
    mol: &Molecule,
    coords: &[[f64; 3]],
    policy: ForceFieldPolicy,
    status: OptimizerStatus,
    max_residual_force: f64,
    distance_geometry_v2_retry_attempted: bool,
) -> Result<(), ForceFieldBridgeError> {
    let worst_bond_length = worst_bond_length_vec(mol, coords);
    let reason = if coords.iter().any(|p| p.iter().any(|x| !x.is_finite())) {
        Some(MinimizationFailureReason::NonFiniteCoordinates)
    } else if worst_bond_length > MAX_SANE_BOND_LENGTH {
        Some(MinimizationFailureReason::CatastrophicBondBlowup)
    } else if max_residual_force > MAX_SANE_RESIDUAL_FORCE {
        Some(MinimizationFailureReason::ExcessiveResidualForce)
    } else {
        None
    };
    match reason {
        None => Ok(()),
        Some(reason) => Err(ForceFieldBridgeError::MinimizationFailed(Box::new(
            MinimizationFailureDetail {
                policy,
                reason,
                converged: status.converged,
                iterations: status.iterations,
                mmff94_termination: status.mmff94_termination,
                max_residual_force,
                worst_bond_length,
                distance_geometry_v2_retry_attempted,
            },
        ))),
    }
}

// --- MMFF94 parameter-coverage checking -------------------------------------
//
// Post-#183 (chematic-ff MMFF94 bond/angle/torsion classification fixes),
// `bond_type_for`/`angle_type_for`/`torsion_type_for`/`MLTB_TYPES`/
// `OOP_SP2_TYPES` are now `pub` in chematic-ff (a narrow, additive visibility
// change made alongside this fix) and are called directly here — no more
// hand-copied duplicate classification. An earlier version of this bridge
// duplicated these (chematic-ff's versions were private pre-#183), which was
// deliberately kept in lockstep with chematic-ff's *pre-#183* bugs (bond
// typing via OR-logic with no bond-order check; angle type hardcoded to 0).
// That duplicate is now confirmed to have caused a real false-negative post-
// #183: a carbonyl C=O bond (types 3,7) is correctly classified bond_type=0
// by the real, order-aware `bond_type_for` (finds the `(0,3,7,...)` row —
// covered), but the old duplicate's order-blind OR-logic computed
// bond_type=1 for it (type 3 is SP2) and found no `(1,3,7,...)` row —
// falsely reporting "missing coverage" on one of the most common bonds in
// organic chemistry (every ketone/aldehyde/amide/carboxylic acid). Calling
// chematic-ff's real functions directly makes this bridge's coverage report
// track chematic-ff's actual behavior by construction, not by manually kept
// parity — see `carbonyl_c_o_double_bond_is_not_falsely_reported_missing`
// (test, below) for the pinned regression.

/// SSSR rings for `mol`, computed once per coverage pass and threaded through
/// to `angle_type_for`/`torsion_type_for` exactly as `chematic_ff`'s own
/// `mmff94_total_energy`/`angle_energy`/`torsion_energy` do internally (same
/// `find_sssr` call), so ring-size-dependent angle/torsion typing can't
/// diverge from what chematic-ff will actually compute over.
fn mmff94_rings(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    find_sssr(mol).rings().to_vec()
}

/// Reorders `atoms` into a canonical form so that physically-equivalent
/// citations (a bond read C→N vs N→C, an angle's two outer substituents in
/// either order, a torsion read forward vs backward, an out-of-plane
/// center's 3 interchangeable substituents in any order) produce the exact
/// same `atoms`/`description` — otherwise aggregating "distinct missing
/// patterns" across a corpus double-counts the same underlying gap under
/// two different string keys purely because of atom-enumeration order
/// (found in independent review: "Bond C-N" and "Bond N-C" were counted
/// separately). Ordering by ascending MMFF94 numeric type matches
/// chematic-ff's own symmetric lookups, which already sort
/// `(type_i, type_j)` before searching (see `mmff94_bond_params`/
/// `mmff94_bond_energy`), so this is the same canonical order chematic-ff
/// itself uses internally, not an arbitrary new one.
fn canonicalize_term_atoms(kind: Mmff94TermKind, atoms: &[AtomIdx], types: &[u8]) -> Vec<AtomIdx> {
    let ty = |a: AtomIdx| types[a.0 as usize];
    let mut v = atoms.to_vec();
    match kind {
        Mmff94TermKind::Bond => {
            if ty(v[0]) > ty(v[1]) {
                v.swap(0, 1);
            }
        }
        // [outer_a, center, outer_c]: center fixed, sort the outer pair.
        Mmff94TermKind::Angle => {
            if ty(v[0]) > ty(v[2]) {
                v.swap(0, 2);
            }
        }
        // [i, j, k, l]: reading backward (l, k, j, i) is the same torsion.
        Mmff94TermKind::Torsion => {
            let fwd: Vec<u8> = v.iter().map(|&a| ty(a)).collect();
            let rev: Vec<u8> = fwd.iter().rev().copied().collect();
            if rev < fwd {
                v.reverse();
            }
        }
        // [center, s1, s2, s3]: center fixed, the 3 substituents are
        // mutually interchangeable.
        Mmff94TermKind::Oop => {
            v[1..].sort_by_key(|&a| ty(a));
        }
        // Same [outer_a, center, outer_c] layout as Angle -- stretch-bend is
        // keyed on the identical angle triple.
        Mmff94TermKind::StretchBend => {
            if ty(v[0]) > ty(v[2]) {
                v.swap(0, 2);
            }
        }
    }
    v
}

fn missing_term(
    mol: &Molecule,
    types: &[u8],
    kind: Mmff94TermKind,
    atoms: &[AtomIdx],
) -> Mmff94MissingTerm {
    let atoms = canonicalize_term_atoms(kind, atoms, types);
    let labels: Vec<String> = atoms
        .iter()
        .map(|&a| format!("{}({})", mol.atom(a).element.symbol(), types[a.0 as usize]))
        .collect();
    let indices: Vec<u32> = atoms.iter().map(|a| a.0).collect();
    Mmff94MissingTerm {
        kind,
        atoms,
        description: format!("{kind:?} {} (atom indices {indices:?})", labels.join("-")),
    }
}

/// Independently measure MMFF94 parameter coverage for every bond, angle,
/// stretch-bend, torsion, and out-of-plane center that chematic-ff's own
/// energy functions would evaluate for `mol`. Mirrors
/// `chematic_ff::mmff94_minimizer`'s private energy-term enumeration loops
/// exactly, including the four-distinct-atoms constraint on torsions, so
/// the totals match the interactions evaluated by chematic-ff.
fn compute_mmff94_coverage(mol: &Molecule, types: &[u8]) -> Mmff94CoverageReport {
    let mut report = Mmff94CoverageReport::default();
    let rings = mmff94_rings(mol);

    collect_bond_coverage(mol, types, &mut report);
    collect_angle_coverage(mol, types, &rings, &mut report);
    collect_torsion_coverage(mol, types, &mut report);
    collect_oop_coverage(mol, types, &mut report);
    report
}

fn collect_bond_coverage(
    mol: &Molecule,
    types: &[u8],
    report: &mut Mmff94CoverageReport,
) {
    for (_, bond) in mol.bonds() {
        report.bonds_total += 1;
        let (a1, a2) = (bond.atom1, bond.atom2);
        let (t1, t2) = (types[a1.0 as usize], types[a2.0 as usize]);
        let bt = bond_type_for(t1, t2, bond.order);
        if mmff94_bond_energy_resolved(bt, t1, t2).is_none() {
            report
                .bonds_missing
                .push(missing_term(mol, types, Mmff94TermKind::Bond, &[a1, a2]));
        }
    }
}

fn collect_angle_coverage(
    mol: &Molecule,
    types: &[u8],
    rings: &[Vec<AtomIdx>],
    report: &mut Mmff94CoverageReport,
) {
    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                report.angles_total += 1;
                report.stretch_bend_total += 1;
                let (a, c) = (neighbors[i], neighbors[j]);
                let (ta, tc) = (types[a.0 as usize], types[c.0 as usize]);
                let at = angle_type_for(mol, rings, a.0 as usize, b_idx, c.0 as usize, types);
                // Shared by the angle-bend check below AND the stretch-bend
                // check further down -- stretch_bend_energy (chematic-ff's
                // mmff94_minimizer) iterates the identical neighbor-pair loop
                // and needs the same flanking-bond types/order (issue #227
                // Priority 2C: the STBN table is keyed by stretch-bend type,
                // not angle type -- see chematic_ff::stretch_bend_type_for's
                // doc).
                let order_ab = mol.bond_between(a, b).expect("a-b angle bond").1.order;
                let order_cb = mol.bond_between(c, b).expect("c-b angle bond").1.order;
                let bt_ab = bond_type_for(ta, types[b_idx], order_ab);
                let bt_cb = bond_type_for(tc, types[b_idx], order_cb);

                // Issue #227 Stage C: the angle-bend term itself needs both
                // flanking bonds' r0 (for the empirical rule, if the direct
                // table/eqLevel ladder don't resolve) -- matches RDKit's real
                // `getMMFFAngleBendParams`, which requires both
                // `getMMFFBondStretchParams` calls to succeed before even
                // attempting the empirical path (mirrors
                // `chematic_ff::mmff94_minimizer`'s own `angle_energy`).
                let bond_ab = mmff94_bond_energy_resolved(bt_ab, ta, types[b_idx]);
                let bond_cb = mmff94_bond_energy_resolved(bt_cb, tc, types[b_idx]);
                let angle_resolved = match (bond_ab, bond_cb) {
                    (Some((bab, _)), Some((bcb, _))) => {
                        let ring_size =
                            is_angle_in_ring_of_size_3_or_4(mol, a.0 as usize, b_idx, c.0 as usize);
                        mmff94_angle_energy_resolved(
                            at,
                            ta,
                            types[b_idx],
                            tc,
                            bab.r0,
                            bcb.r0,
                            ring_size,
                        )
                    }
                    _ => None,
                };
                if angle_resolved.is_none() {
                    report.angles_missing.push(missing_term(
                        mol,
                        types,
                        Mmff94TermKind::Angle,
                        &[a, b, c],
                    ));
                }
                let sbt = stretch_bend_type_for(at, ta, tc, bt_ab, bt_cb);
                if mmff94_stbn(
                    sbt,
                    ta,
                    types[b_idx],
                    tc,
                    mol.atom(a).element.atomic_number(),
                    mol.atom(b).element.atomic_number(),
                    mol.atom(c).element.atomic_number(),
                )
                .is_none()
                {
                    report.stretch_bend_missing.push(missing_term(
                        mol,
                        types,
                        Mmff94TermKind::StretchBend,
                        &[a, b, c],
                    ));
                }
            }
        }
    }
}

fn collect_torsion_coverage(
    mol: &Molecule,
    types: &[u8],
    report: &mut Mmff94CoverageReport,
) {
    for (_, bond) in mol.bonds() {
        let (j, k) = (bond.atom1, bond.atom2);
        let nbrs_j: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
        let nbrs_k: Vec<AtomIdx> = mol.neighbors(k).map(|(nb, _)| nb).collect();
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                if l == j || l == i {
                    continue;
                }
                report.torsions_total += 1;
                // Table row or, failing that, Halgren's empirical rule, as
                // RDKit and the energy model take it: a torsion is never
                // uncovered; all-zero barriers add no term (linear centres
                // among them).
                if mmff94_torsion_term_params(
                    mol,
                    i.0 as usize,
                    j.0 as usize,
                    k.0 as usize,
                    l.0 as usize,
                    types,
                )
                .is_none()
                {
                    report.torsions_no_term_by_design += 1;
                }
            }
        }
    }
}

fn collect_oop_coverage(
    mol: &Molecule,
    types: &[u8],
    report: &mut Mmff94CoverageReport,
) {
    for j_idx in 0..mol.atom_count() {
        let tj = types[j_idx];
        if OOP_SP2_TYPES.binary_search(&tj).is_err() {
            continue;
        }
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
        if neighbors.len() != 3 {
            continue;
        }
        report.oop_total += 1;
        let [i, k, l] = [neighbors[0], neighbors[1], neighbors[2]];
        if mmff94_oop(
            tj,
            types[i.0 as usize],
            types[k.0 as usize],
            types[l.0 as usize],
        )
        .is_none()
        {
            report
                .oop_missing
                .push(missing_term(mol, types, Mmff94TermKind::Oop, &[j, i, k, l]));
        }
    }
}
