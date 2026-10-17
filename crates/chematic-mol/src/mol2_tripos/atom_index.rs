use chematic_core::AtomIdx;

use super::Mol2Error;

// Tiny records avoid sorting. Larger records pay one in-place sort instead
// of scanning every atom twice for every bond. AtomIdx follows input order,
// so the secondary key retains the old first-match rule for duplicate IDs.
const LINEAR_LIMIT: usize = 32;

pub(super) fn prepare_atom_ids(ids: &mut [(u32, AtomIdx)]) {
    if ids.len() > LINEAR_LIMIT {
        ids.sort_unstable_by_key(|&(id, index)| (id, index.0));
    }
}

pub(super) fn mol2_atom_index(
    ids: &[(u32, AtomIdx)],
    id: u32,
    line: usize,
) -> Result<AtomIdx, Mol2Error> {
    let found = if ids.len() <= LINEAR_LIMIT {
        ids.iter().find(|&&(candidate, _)| candidate == id)
    } else {
        let first = ids.partition_point(|&(candidate, _)| candidate < id);
        ids.get(first).filter(|&&(candidate, _)| candidate == id)
    };
    found
        .map(|&(_, index)| index)
        .ok_or_else(|| Mol2Error::InvalidBondLine {
            line,
            detail: format!("atom_id {id} not found"),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prepared_lookup_matches_frozen_first_match_for_every_id() {
        for len in [0, 1, 3, 31, 32, 33, 64, 257, 4096] {
            for pattern in 0..6 {
                let original: Vec<_> = (0..len)
                    .map(|i| {
                        let id = match pattern {
                            0 => i as u32 + 1,
                            1 => (len - i) as u32,
                            2 => (i % 13) as u32,
                            3 => u32::MAX - (i as u32 * 17),
                            4 => (i as u32 * 7919) % 104729,
                            _ => 0,
                        };
                        (id, AtomIdx(i as u32))
                    })
                    .collect();
                let mut prepared = original.clone();
                prepare_atom_ids(&mut prepared);
                // Existing IDs, absent IDs and the unsigned boundaries.
                for id in original
                    .iter()
                    .map(|&(id, _)| id)
                    .chain([0, 1, 999_999, u32::MAX])
                {
                    let expected = original
                        .iter()
                        .find(|&&(candidate, _)| candidate == id)
                        .map(|&(_, idx)| idx)
                        .ok_or_else(|| Mol2Error::InvalidBondLine {
                            line: 123,
                            detail: format!("atom_id {id} not found"),
                        });
                    assert_eq!(
                        mol2_atom_index(&prepared, id, 123),
                        expected,
                        "len={len} pattern={pattern} id={id}"
                    );
                }
            }
        }
    }
}
