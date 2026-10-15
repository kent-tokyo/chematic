//! Preserve lexical neighbour order while resolving ring closures.
use super::*;

impl Parser<'_> {
    // -- ring closures -------------------------------------------------------

    /// Handle a ring closure: either close an existing open ring or open a new one.
    pub(super) fn handle_ring_closure(
        &mut self,
        mol: &mut QueryMolecule,
        current: usize,
        ring_num: u8,
        ring_bond: Option<BondQuery>,
        open_rings: &mut HashMap<u8, (usize, Option<BondQuery>, usize)>,
    ) -> Result<(), SmartsError> {
        if let Some((open_atom, open_bond, open_position)) = open_rings.remove(&ring_num) {
            // Resolve the bond type.
            let bond = match (open_bond, ring_bond) {
                (Some(a), Some(b)) if a == b => a,
                (Some(_), Some(_)) => {
                    // Conflicting ring bond specifications — use Any as a fallback.
                    BondQuery::Any
                }
                (Some(b), None) | (None, Some(b)) => b,
                (None, None) => BondQuery::Any,
            };
            mol.add_bond(open_atom, current, bond);
            // `add_bond` appends at both endpoints. At the closing endpoint
            // that is the lexical position of the ring digit. At the opening
            // endpoint the digit was seen earlier, so move only that adjacency
            // entry back to the recorded position. Graph semantics do not
            // change; this preserves the order needed for stereo parity.
            let open_adj = &mut mol.adj[open_atom];
            // safe: add_bond appended the new adjacency entry.
            let entry = open_adj
                .pop()
                .expect("the just-added ring bond has an adjacency entry");
            let pending_before = open_rings
                .values()
                .filter(|&&(atom, _, position)| atom == open_atom && position < open_position)
                .count();
            let position = open_position.saturating_sub(pending_before);
            open_adj.insert(position.min(open_adj.len()), entry);
        } else {
            let pending = open_rings
                .values()
                .filter(|&&(atom, _, _)| atom == current)
                .count();
            open_rings.insert(
                ring_num,
                (current, ring_bond, mol.adj[current].len() + pending),
            );
        }
        Ok(())
    }
}
