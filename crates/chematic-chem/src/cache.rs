//! Descriptor caching for improved performance.
//!
//! Simple LRU-like cache for descriptors keyed by molecule canonical SMILES.
//! Improves performance when computing same molecules repeatedly.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
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
    access_generation: HashMap<String, u64>,
    order: BinaryHeap<Reverse<(u64, String)>>,
    next_generation: u64,
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
            touch(&mut state, smiles);
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
                && let Some(oldest) = pop_oldest(&mut state)
            {
                state.entries.remove(&oldest);
                state.access_generation.remove(&oldest);
            }
            state.entries.insert(smiles.clone(), entry);
            if !is_new {
                state.access_generation.remove(&smiles);
            }
            touch(&mut state, &smiles);
        }
    }

    /// Clear all cached entries.
    pub fn clear(&self) {
        if let Ok(mut state) = self.state.lock() {
            state.entries.clear();
            state.access_generation.clear();
            state.order.clear();
            state.next_generation = 0;
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

fn touch(state: &mut DescriptorCacheState, key: &str) {
    state.next_generation = state.next_generation.wrapping_add(1);
    let generation = state.next_generation;
    state.access_generation.insert(key.to_owned(), generation);
    state.order.push(Reverse((generation, key.to_owned())));
    let rebuild_at = state.entries.len().saturating_mul(4).max(16);
    if state.order.len() > rebuild_at {
        state.order = state
            .access_generation
            .iter()
            .map(|(key, &generation)| Reverse((generation, key.clone())))
            .collect();
    }
}

fn pop_oldest(state: &mut DescriptorCacheState) -> Option<String> {
    while let Some(Reverse((generation, key))) = state.order.pop() {
        if state.access_generation.get(&key) == Some(&generation) {
            return Some(key);
        }
    }
    None
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
    fn repeated_hits_bound_recency_heap_growth() {
        let cache = DescriptorCache::new(1);
        cache.put("CC".to_string(), DescriptorEntry::default());
        for _ in 0..1000 {
            assert!(cache.get("CC").is_some());
        }
        let state = cache.state.lock().expect("cache state");
        assert!(state.order.len() <= 16);
        assert_eq!(state.access_generation.len(), 1);
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
