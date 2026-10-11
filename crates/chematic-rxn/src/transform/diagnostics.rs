use chematic_core::{AtomIdx, Element};

use crate::reaction::RxnError;

/// Error type for SMIRKS transformation.
#[derive(Debug)]
pub enum TransformError {
    /// A configured match-enumeration limit was exceeded.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    SmirksParse(RxnError),
    ReactantCountMismatch {
        expected: usize,
        got: usize,
    },
}

impl core::fmt::Display for TransformError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "reaction {resource} exceeds limit {limit} (got {actual})"
            ),
            Self::SmirksParse(e) => write!(f, "SMIRKS parse error: {e}"),
            Self::ReactantCountMismatch { expected, got } => {
                write!(f, "reactant count mismatch: expected {expected}, got {got}")
            }
        }
    }
}

impl std::error::Error for TransformError {}

impl From<RxnError> for TransformError {
    fn from(error: RxnError) -> Self {
        Self::SmirksParse(error)
    }
}

/// Resource limits for reaction-template matching and product generation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionTransformLimits {
    /// Maximum number of accepted match combinations.
    pub max_matches: usize,
}

impl Default for ReactionTransformLimits {
    fn default() -> Self {
        Self {
            max_matches: 100_000,
        }
    }
}

/// Bounded accounting for one reaction-template application.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionTransformDiagnostics {
    /// Number of matches accepted by the matcher and stereo post-checks.
    pub accepted_matches: usize,
    /// Number of accepted matches that produced a valence-valid product set.
    pub applied_products: usize,
    /// Number of accepted matches rejected by product valence validation.
    pub valence_rejected_matches: usize,
    /// Always false for the current fail-on-limit policy; reserved for a future
    /// truncating mode so consumers can distinguish partial reports.
    pub truncated_matches: bool,
}

/// Why RDKit-compatible product sanitization rejected a generated product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductRejectionReason {
    /// An aromatic atom or bond is inconsistent with ring membership.
    Aromaticity,
    /// The aromatic bond system cannot be assigned a valid Kekule form.
    Kekulization,
    /// One or more atoms exceed RDKit's allowed explicit valence.
    Valence,
    /// The product was rejected but no reliable narrower cause was found.
    Unknown,
}

impl ProductRejectionReason {
    /// Stable machine-readable code shared by Rust, Python and WASM.
    pub const fn reason_code(self) -> &'static str {
        match self {
            Self::Aromaticity => "aromaticity",
            Self::Kekulization => "kekulization",
            Self::Valence => "valence",
            Self::Unknown => "unknown",
        }
    }
}

/// Reliable atom-level context for one rejected product.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductAtomRejectionDiagnostic {
    /// Atom index in the generated product before any serialization reorder.
    pub atom: AtomIdx,
    /// Product-template map label, when this atom came from a mapped slot.
    pub template_map: Option<u16>,
    /// Element of the generated product atom.
    pub element: Element,
    /// Explicit hydrogen count recorded on the generated atom.
    pub explicit_hydrogen_count: Option<u8>,
    /// Kekule explicit valence after RDKit cleanup rules, when applicable.
    pub observed_valence: Option<i16>,
    /// Largest valence RDKit permits after cleanup, when bounded.
    pub max_allowed_valence: Option<i16>,
}

/// Sanitization verdict for one generated product in a rejected product set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductSanitizationDiagnostic {
    /// Zero-based product position within the generated product set.
    pub product_index: usize,
    /// Whether this product itself is sanitizable. A set is rejected when at
    /// least one product is not, so accepted siblings remain visible here.
    pub accepted: bool,
    /// Typed rejection reason (`None` for an accepted sibling product).
    pub reason: Option<ProductRejectionReason>,
    /// At most 32 implicated atoms, in ascending product atom order.
    pub atoms: Vec<ProductAtomRejectionDiagnostic>,
}

/// Why the pinned RDKit 2026.03.6 compatibility profile cannot safely claim
/// a reaction result. Native CheMatic reaction semantics are unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionCompatibilityUnsupported {
    /// A stereo template uses syntax the RDKit reading here does not follow
    /// (component grouping, `@?`, non-tetrahedral chirality classes).
    ChiralReactantTemplateSemantics,
    /// A product's tetrahedral stereo depends on a reactant bond order the
    /// input molecule does not record (RDKit orders an atom's ring-closure
    /// bonds by ring number).
    AmbiguousStereoBondOrder,
    /// RDKit ignores reactant-template `/` `\\` when matching; this match
    /// was rejected by the native E/Z check.
    EzReactantTemplateSemantics,
    /// A reactant template asks for a ring count or ring size (`R<n>`,
    /// `r<n>`) and RDKit's symmetrized ring set for this reactant could not
    /// be settled (the ring model is ambiguous or over its budget).
    RingModelUnresolved,
}

impl ReactionCompatibilityUnsupported {
    /// Stable machine-readable refusal code shared by bindings and audit gates.
    pub const fn reason_code(self) -> &'static str {
        match self {
            Self::ChiralReactantTemplateSemantics => "chiral_reactant_template_semantics",
            Self::AmbiguousStereoBondOrder => "ambiguous_stereo_bond_order",
            Self::EzReactantTemplateSemantics => "ez_reactant_template_semantics",
            Self::RingModelUnresolved => "ring_model_unresolved",
        }
    }
}

/// Diagnostics for one normalized atomic-number variant.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionVariantDiagnostics {
    /// Stable zero-based order in the deterministic expansion.
    pub variant_index: usize,
    /// SMILES-compatible reaction template used for this variant.
    pub normalized_smirks: String,
    pub diagnostics: ReactionTransformDiagnostics,
}
