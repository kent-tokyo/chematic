#![forbid(unsafe_code)]
//! `chematic-rxn` — reaction SMILES parser and writer for chematic.
//!
//! Provides:
//! - [`Reaction`]: reactants, agents, products as `Vec<Molecule>`.
//! - [`parse_reaction`]: parse a reaction SMILES string `"R>>P"` or `"R>A>P"`.
//! - [`write_reaction`]: serialize back to reaction SMILES.
//! - [`RxnError`]: parse error type.
//! - [`run_reactants`]: apply a SMIRKS template to reactant molecules.
//! - [`PreparedReaction`]: parse/compile a SMIRKS template once for repeated
//!   application, including optional caller-provided ring perception.
//!   Variant-level diagnostics are available for both ordinary and
//!   caller-provided-ring application paths.
//! - [`ReactionRequirements`]: conservative lower bounds for fail-open
//!   candidate prefilters derived from a prepared template.
//! - [`find_reaction_matches`]/[`apply_reaction_match`]: enumerate matches and
//!   apply one of them independently, for callers that need to accept/reject
//!   individual matches rather than an entire `run_reactants` call.
//! - [`TransformError`]: error type for SMIRKS transformation.
//! - [`enumerate_library`]: combinatorial library enumeration from SMIRKS + fragment sets.

pub mod balance;
pub mod document;
pub mod enumerate;
pub mod green;
pub mod perf_counters;
pub mod query;
pub mod reaction;
pub mod requirements;
pub mod retro;
pub mod stoichiometry;
pub mod transform;

pub use balance::{BalanceResult, balance_check};
pub use document::{
    ComponentRole, ContentOrigin, ProvenanceRecord, ReactionAtomMap, ReactionComponent,
    ReactionCondition, ReactionDocument, ReactionDocumentEdit, ReactionDocumentError, ReactionLoss,
    ReactionStep,
};
pub use enumerate::{
    LibraryConfig, LibraryError, enumerate_library, enumerate_library_2way, enumerate_library_3way,
};
pub use green::{atom_economy, e_factor, pmi_rxn, reaction_mass_efficiency};
pub use perf_counters::PerfCounters;
pub use query::{
    AgentMatches, BatchQueryLimits, BatchQueryResults, ReactionPatternLibrary, ReactionQuery,
    ReactionQueryError, ReactionSmartsMatch, batch_query_reactions,
    batch_query_reactions_with_limits, batch_query_with_library,
    batch_query_with_library_with_limits, has_reaction_substructure_match, parse_reaction_query,
    query_reaction,
};
pub use reaction::{
    Reaction, ReactionCenter, ReactionParseLimits, RxnError, expand_atomic_number_primitives,
    find_reaction_center, parse_reaction, parse_reaction_with_limits, write_reaction,
};
pub use requirements::{ReactionBondKind, ReactionBondLowerBound, ReactionRequirements};
pub use retro::{DEFAULT_TEMPLATES, RetroClass, RetroResult, RetroTemplate, retro_disconnect};
pub use stoichiometry::{
    AtomInventory, ChemicalCompleteness, ComponentEvidence, DiagnosticSeverity,
    StepStoichiometryReport, StoichiometryComponent, StoichiometryDiagnostic, StoichiometryError,
    StoichiometryEvidenceScope, StoichiometryIssueCode, StoichiometryReport, StoichiometryStatus,
    StoichiometryStep, analyze_components, analyze_reaction_document, analyze_reaction_step,
};
pub use transform::{
    PreparedReaction, ReactionMatch, ReactionTransformDiagnostics, ReactionTransformLimits,
    ReactionTransformReport, ReactionVariantDiagnostics, TransformError, apply_reaction_match,
    find_reaction_matches, find_reaction_matches_with_limits, run_reactants, run_reactants_strict,
    run_reactants_strict_with_limits, run_reactants_with_diagnostics, run_reactants_with_limits,
};
