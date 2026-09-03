//! Descriptor caching for improved performance.
//!
//! Simple LRU-like cache for descriptors keyed by molecule canonical SMILES.
//! Improves performance when computing same molecules repeatedly.

use std::collections::{HashMap, VecDeque};
use std::sync::{Arc, Mutex};

/// Descriptor cache entry: stores computed descriptor values.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct DescriptorEntry {
    /// Molecular weight (Ø = computed on demand)
    pub mw: Option<f64>,
    /// LogP (Ø)
    pub logp: Option<f64>,
    /// TPSA (Ø)
    pub tpsa: Option<f64>,
    /// HBA count (Ø)
    pub hba: Option<usize>,
    /// HBD count (Ø)
    pub hbd: Option<usize>,
    /// Rotatable bond count (Ø)
    pub rotb: Option<usize>,
}

/// Thread-safe descriptor cache with max_size limit.
#[derive(Clone, Debug)]
pub struct DescriptorCache {
    state: Arc<Mutex<DescriptorCacheState>>,
    max_size: usize,
}

#[derive(Debug, Default)]
struct DescriptorCacheState {
    entries: HashMap<String, DescriptorEntry>,
    order: VecDeque<String>,
}

impl DescriptorCache {
    /// Create new cache with given max size.
    pub fn new(max_size: usize) -> Self {
        Self {
            state: Arc::new(Mutex::new(DescriptorCacheState::default())),
            max_size,
        }
    }

    /// Get cached entry for molecule (keyed by canonical SMILES).
    pub fn get(&self, smiles: &str) -> Option<DescriptorEntry> {
        let mut state = self.state.lock().ok()?;
        let entry = state.entries.get(smiles).cloned();
        if entry.is_some() {
            if let Some(pos) = state.order.iter().position(|key| key == smiles) {
                state.order.remove(pos);
            }
            state.order.push_back(smiles.to_owned());
        }
        entry
    }

    /// Store/update cached entry.
    pub fn put(&self, smiles: String, entry: DescriptorEntry) {
        if self.max_size == 0 {
            return;
        }
        if let Ok(mut state) = self.state.lock() {
            let is_new = !state.entries.contains_key(&smiles);
            if !state.entries.contains_key(&smiles)
                && state.entries.len() >= self.max_size
                && let Some(oldest) = state.order.pop_front()
            {
                state.entries.remove(&oldest);
            }
            state.entries.insert(smiles.clone(), entry);
            if !is_new && let Some(pos) = state.order.iter().position(|key| key == &smiles) {
                state.order.remove(pos);
            }
            state.order.push_back(smiles);
        }
    }

    /// Clear all cached entries.
    pub fn clear(&self) {
        if let Ok(mut state) = self.state.lock() {
            state.entries.clear();
            state.order.clear();
        }
    }

    /// Get cache size.
    pub fn len(&self) -> usize {
        self.state
            .lock()
            .map(|state| state.entries.len())
            .unwrap_or(0)
    }

    /// Check if cache is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cache_put_get() {
        let cache = DescriptorCache::new(100);
        let entry = DescriptorEntry {
            mw: Some(46.0),
            ..DescriptorEntry::default()
        };

        cache.put("CC".to_string(), entry.clone());
        let retrieved = cache.get("CC").unwrap();
        assert_eq!(retrieved.mw, Some(46.0));
    }

    #[test]
    fn cache_miss() {
        let cache = DescriptorCache::new(100);
        assert_eq!(cache.get("CC"), None);
    }

    #[test]
    fn cache_eviction() {
        let cache = DescriptorCache::new(2);
        let entry = DescriptorEntry::default();

        cache.put("C1".to_string(), entry.clone());
        cache.put("C2".to_string(), entry.clone());
        cache.put("C3".to_string(), entry.clone());

        assert_eq!(cache.len(), 2);
        assert!(cache.get("C1").is_none());
    }

    #[test]
    fn zero_capacity_disables_storage() {
        let cache = DescriptorCache::new(0);
        cache.put("CC".to_string(), DescriptorEntry::default());
        assert_eq!(cache.len(), 0);
        assert!(cache.is_empty());
    }

    #[test]
    fn cache_eviction_is_lru() {
        let cache = DescriptorCache::new(2);
        let entry = DescriptorEntry::default();
        cache.put("C1".to_string(), entry.clone());
        cache.put("C2".to_string(), entry.clone());
        assert!(cache.get("C1").is_some());
        cache.put("C3".to_string(), entry);
        assert!(cache.get("C1").is_some());
        assert!(cache.get("C2").is_none());
    }

    #[test]
    fn cache_clear() {
        let cache = DescriptorCache::new(100);
        let entry = DescriptorEntry::default();

        cache.put("CC".to_string(), entry);
        assert!(!cache.is_empty());

        cache.clear();
        assert!(cache.is_empty());
    }
}
