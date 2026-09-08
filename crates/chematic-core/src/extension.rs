//! Stable, host-agnostic extension points for molecule consumers.
//!
//! Extensions are deliberately small: they receive an immutable [`Molecule`]
//! and return one of the serialisable result shapes understood by the higher
//! level bindings.  Keeping this contract in `chematic-core` makes extensions
//! usable by native Rust and WASM callers without pulling in a descriptor or
//! I/O crate.

use crate::Molecule;

/// Result values produced by a [`MoleculeExtension`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ExtensionValue {
    /// One numeric result, such as a score or count.
    Scalar(f64),
    /// A fixed-order numeric vector, such as a fingerprint summary.
    Vector(Vec<f64>),
    /// A small textual result, such as a canonical identifier.
    Text(String),
}

/// Failure returned by an extension.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExtensionError {
    /// The molecule cannot be processed under the extension's contract.
    InvalidInput(String),
    /// The extension failed while calculating its result.
    Failed(String),
}

impl core::fmt::Display for ExtensionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(message) => write!(f, "invalid extension input: {message}"),
            Self::Failed(message) => write!(f, "extension failed: {message}"),
        }
    }
}

impl std::error::Error for ExtensionError {}

/// A named, read-only operation supplied by an ecosystem consumer.
pub trait MoleculeExtension: Send + Sync {
    /// Stable machine-readable identifier, for example `vendor.property.v1`.
    fn id(&self) -> &str;

    /// Extension contract version. Increment it when the result semantics change.
    fn version(&self) -> u32;

    /// Evaluate this extension for one molecule.
    fn run(&self, molecule: &Molecule) -> Result<ExtensionValue, ExtensionError>;
}

/// Deterministic in-process registry for user-supplied extensions.
#[derive(Default)]
pub struct ExtensionRegistry {
    extensions: Vec<Box<dyn MoleculeExtension>>,
}

impl ExtensionRegistry {
    /// Create an empty registry.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register an extension, rejecting a duplicate identifier.
    pub fn register<E>(&mut self, extension: E) -> Result<(), ExtensionError>
    where
        E: MoleculeExtension + 'static,
    {
        if self
            .extensions
            .iter()
            .any(|existing| existing.id() == extension.id())
        {
            return Err(ExtensionError::Failed(format!(
                "duplicate extension id: {}",
                extension.id()
            )));
        }
        self.extensions.push(Box::new(extension));
        Ok(())
    }

    /// Return registered extensions in registration order.
    pub fn extensions(&self) -> impl Iterator<Item = &dyn MoleculeExtension> {
        self.extensions.iter().map(Box::as_ref)
    }

    /// Run one extension by ID.
    pub fn run(&self, id: &str, molecule: &Molecule) -> Result<ExtensionValue, ExtensionError> {
        self.extensions
            .iter()
            .find(|extension| extension.id() == id)
            .ok_or_else(|| ExtensionError::Failed(format!("unknown extension id: {id}")))?
            .run(molecule)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Element, MoleculeBuilder};

    struct AtomCount;

    impl MoleculeExtension for AtomCount {
        fn id(&self) -> &str {
            "test.atom-count.v1"
        }

        fn version(&self) -> u32 {
            1
        }

        fn run(&self, molecule: &Molecule) -> Result<ExtensionValue, ExtensionError> {
            Ok(ExtensionValue::Scalar(molecule.atom_count() as f64))
        }
    }

    fn ethanol() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(Atom::new(Element::C));
        let c2 = builder.add_atom(Atom::new(Element::C));
        builder.add_bond(c1, c2, crate::BondOrder::Single).unwrap();
        builder.build()
    }

    #[test]
    fn registry_runs_extension_and_preserves_order() {
        let mut registry = ExtensionRegistry::new();
        registry.register(AtomCount).unwrap();
        assert_eq!(
            registry
                .extensions()
                .map(|extension| extension.id())
                .collect::<Vec<_>>(),
            ["test.atom-count.v1"]
        );
        assert_eq!(
            registry.run("test.atom-count.v1", &ethanol()),
            Ok(ExtensionValue::Scalar(2.0))
        );
    }

    #[test]
    fn registry_rejects_duplicate_and_unknown_ids() {
        let mut registry = ExtensionRegistry::new();
        registry.register(AtomCount).unwrap();
        assert!(matches!(
            registry.register(AtomCount),
            Err(ExtensionError::Failed(message)) if message.contains("duplicate")
        ));
        assert!(matches!(
            registry.run("missing.v1", &ethanol()),
            Err(ExtensionError::Failed(message)) if message.contains("unknown")
        ));
    }
}
