//! SMARTS compilation cache — avoids re-parsing repeated SMARTS strings.
//!
//! Re-parsing the same SMARTS pattern on every call is the dominant overhead
//! in high-throughput substructure screening. This LRU cache compiles each
//! pattern once and reuses the `QueryMolecule` for subsequent matches.
//!
//! ## Usage
//! ```rust,ignore
//! use chematic_smarts::SmartsCache;
//! let mut cache = SmartsCache::new(512);
//! let matches = cache.find_matches("[OH]", &mol)?;
//! ```
//!
//! ## Performance
//! Repeated SMARTS matching with a warm cache is 5–20× faster than calling
//! `parse_smarts()` + `find_matches()` each time.

use rustc_hash::FxHashMap;
use std::collections::{HashMap, VecDeque};

use crate::{
    match_vf2::{MatchConfig, find_matches_with_config},
    parser::{SmartsError, parse_smarts},
    query::QueryMolecule,
};
use chematic_core::{AtomIdx, Molecule};

/// LRU cache for compiled SMARTS patterns.
///
/// Patterns are compiled once on first use and stored as `QueryMolecule`.
/// When the cache reaches `capacity`, the least-recently-used entry is evicted.
pub struct SmartsCache {
    capacity: usize,
    // std HashMap (SipHash) intentionally — keys are user-supplied SMARTS strings;
    // FxHash has no random seed and is vulnerable to HashDoS with adversarial input.
    store: HashMap<String, QueryMolecule>,
    order: VecDeque<String>, // LRU order (front = oldest, back = newest)
}

impl SmartsCache {
    /// Create a new cache with the given maximum number of compiled patterns.
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity: capacity.max(1),
            store: HashMap::new(),
            order: VecDeque::new(),
        }
    }

    /// Compile a SMARTS pattern (or retrieve from cache) and return a reference.
    pub fn compile(&mut self, smarts: &str) -> Result<&QueryMolecule, SmartsError> {
        if !self.store.contains_key(smarts) {
            let qmol = parse_smarts(smarts)?;
            if self.store.len() >= self.capacity {
                // Evict LRU: pop_front is O(1) with VecDeque
                if let Some(oldest) = self.order.pop_front() {
                    self.store.remove(&oldest);
                }
            }
            self.store.insert(smarts.to_string(), qmol);
            self.order.push_back(smarts.to_string());
        } else {
            // Move to back (most recently used)
            if let Some(pos) = self.order.iter().position(|s| s == smarts) {
                let key = self.order.remove(pos).unwrap();
                self.order.push_back(key);
            }
        }
        Ok(self.store.get(smarts).unwrap())
    }

    /// Find all substructure matches of `smarts` in `mol` using the cache.
    pub fn find_matches(
        &mut self,
        smarts: &str,
        mol: &Molecule,
    ) -> Result<Vec<FxHashMap<usize, AtomIdx>>, SmartsError> {
        let qmol = self.compile(smarts)?.clone();
        Ok(find_matches_with_config(
            &qmol,
            mol,
            &MatchConfig::default(),
        ))
    }

    /// Find all substructure matches with custom `MatchConfig`.
    pub fn find_matches_with_config(
        &mut self,
        smarts: &str,
        mol: &Molecule,
        config: &MatchConfig,
    ) -> Result<Vec<FxHashMap<usize, AtomIdx>>, SmartsError> {
        let qmol = self.compile(smarts)?.clone();
        Ok(find_matches_with_config(&qmol, mol, config))
    }

    /// Check whether `smarts` matches at least once in `mol`.
    ///
    /// Stops the VF2 search at the first embedding (`max_matches: Some(1)`,
    /// `uniquify: false`) instead of enumerating every match and checking
    /// emptiness — existence doesn't need the full match set or dedup pass.
    pub fn has_match(&mut self, smarts: &str, mol: &Molecule) -> Result<bool, SmartsError> {
        let config = MatchConfig {
            max_matches: Some(1),
            uniquify: false,
            ..MatchConfig::default()
        };
        let matches = self.find_matches_with_config(smarts, mol, &config)?;
        Ok(!matches.is_empty())
    }

    /// Number of patterns currently in the cache.
    pub fn len(&self) -> usize {
        self.store.len()
    }

    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.store.is_empty()
    }

    /// Clear all cached patterns.
    pub fn clear(&mut self) {
        self.store.clear();
        self.order.clear();
    }
}

/// Common named SMARTS patterns for organic chemistry.
///
/// Returns the SMARTS string for a named pattern. These are preloaded into
/// `SmartsCache` via `SmartsCache::compile(named_pattern("donor"))`.
pub fn named_pattern(name: &str) -> Option<&'static str> {
    match name {
        // H-bond donors (NH, OH, SH)
        "donor" => Some("[N;!H0;v3,v4&+1]"),
        "donor_strict" => Some("[O,N;!H0]"),
        // H-bond acceptors
        "acceptor" => Some("[N;H0;v3,v4&+1]"),
        "acceptor_strict" => Some("[n;H0;+0,o,s;+0]"),
        // Aromatic
        "aromatic" => Some("[a]"),
        "aromatic_ring" => Some("[a]1:[a]:[a]:[a]:[a]:[a]:1"),
        // Hydrophobic
        "hydrophobic" => Some("[c,s,S&H0&v2,F,Cl,Br,I,#6&!$([#6]~[#7,#8,#15,#16,F,Cl,Br,I])]"),
        // Positively charged
        "positive" => Some("[#7&+,NH2&+0&$(N-[#6]),NH1&+0&$(N(-[#6])-[#6])]"),
        // Negatively charged
        "negative" => Some("[C,S](=[O,S,P])-[OH,O-]"),
        // Functional groups
        "carboxylic_acid" => Some("[C;X3;$(C(-O)-O)](=O)"),
        "aldehyde" => Some("[CX3H1](=O)[#6]"),
        "ketone" => Some("[#6][CX3](=O)[#6]"),
        "alcohol" => Some("[OX2H][CX4;!$(C([OX2H])[O,S,#7,#15])]"),
        "phenol" => Some("[OX2H][cX3]:[c]"),
        "amine_primary" => Some("[NH2;$(N-[#6])]"),
        "amine_secondary" => Some("[NH1;$(N(-[#6])-[#6])]"),
        "amine_tertiary" => Some("[NH0;$(N(-[#6])(-[#6])-[#6])]"),
        "amide" => Some("[NX3][CX3](=[OX1])[#6]"),
        "ester" => Some("[#6][CX3](=O)[OX2H0][#6]"),
        "ether" => Some("[OD2]([#6])[#6]"),
        "halide" => Some("[F,Cl,Br,I]"),
        "aromatic_n" => Some("[n]"),
        "sulfonamide" => {
            Some("[$([#16X4](=[OX1])(=[OX1])([#6])[NX3]),$([#16X4+2]([OX1-])([OX1-])([#6])[NX3])]")
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn cache_compiles_once() {
        let mut cache = SmartsCache::new(10);
        let mol = parse("CCO").expect("parse");
        let m1 = cache.find_matches("[OH]", &mol).expect("m1");
        let m2 = cache.find_matches("[OH]", &mol).expect("m2");
        assert_eq!(m1.len(), m2.len());
        assert_eq!(cache.len(), 1);
    }

    #[test]
    fn cache_evicts_lru() {
        let mut cache = SmartsCache::new(2);
        let mol = parse("C").expect("parse");
        cache.find_matches("[C]", &mol).unwrap();
        cache.find_matches("[#6]", &mol).unwrap();
        assert_eq!(cache.len(), 2);
        // Adding a third evicts the LRU ("[C]")
        cache.find_matches("[H]", &mol).unwrap();
        assert_eq!(cache.len(), 2);
        assert!(!cache.store.contains_key("[C]"));
    }

    #[test]
    fn named_pattern_donor_matches_alcohol() {
        let mol = parse("CCO").expect("parse ethanol");
        let smarts = named_pattern("donor_strict").expect("pattern");
        let mut cache = SmartsCache::new(50);
        let matches = cache.find_matches(smarts, &mol).expect("matches");
        assert!(!matches.is_empty(), "ethanol O-H should match donor");
    }

    #[test]
    fn named_patterns_all_parse() {
        let mut cache = SmartsCache::new(50);
        for name in &[
            "donor",
            "acceptor",
            "aromatic",
            "hydrophobic",
            "positive",
            "negative",
            "carboxylic_acid",
            "aldehyde",
            "ketone",
            "alcohol",
            "phenol",
            "amine_primary",
            "amine_secondary",
            "amide",
            "ester",
            "halide",
        ] {
            if let Some(pat) = named_pattern(name) {
                cache
                    .compile(pat)
                    .unwrap_or_else(|e| panic!("pattern '{}' failed: {}", name, e));
            }
        }
    }

    #[test]
    fn has_match_works() {
        let mut cache = SmartsCache::new(10);
        let mol = parse("c1ccccc1").expect("benzene");
        assert!(cache.has_match("[a]", &mol).unwrap());
        assert!(!cache.has_match("[OH]", &mol).unwrap());
    }
}
