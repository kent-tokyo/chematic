use super::MesoStereoGroupError as Error;
use super::graph::Graph;
use chematic_core::{AtomIdx, STEREO_H_SENTINEL};

pub(super) struct Budget {
    pub visits: u64,
    pub limit: u64,
}
impl Budget {
    fn probe(&mut self) -> Result<(), Error> {
        if self.visits == self.limit {
            return Err(Error::BudgetExceeded {
                visits: self.visits,
                limit: self.limit,
            });
        }
        self.visits += 1;
        Ok(())
    }
}
fn mapped(index: u32, map: &[usize]) -> Option<u32> {
    if index == STEREO_H_SENTINEL {
        Some(index)
    } else {
        map.get(index as usize)
            .copied()
            .filter(|&v| v != usize::MAX)
            .map(|v| v as u32)
    }
}
fn even(from: &[u32; 4], to: &[u32; 4]) -> bool {
    let mut positions = [0; 4];
    for i in 0..4 {
        positions[i] = to
            .iter()
            .position(|&v| v == from[i])
            .expect("mapped full frame has the same ligands");
    }
    let inversions = (0..4)
        .map(|i| (i + 1..4).filter(|&j| positions[i] > positions[j]).count())
        .sum::<usize>();
    inversions.is_multiple_of(2)
}
fn tetrahedral(a: &Graph<'_>, b: &Graph<'_>, map: &[usize]) -> bool {
    for (source, frame) in a.frames.iter().enumerate() {
        let Some(frame) = frame else {
            continue;
        };
        let target = map[source];
        if target == usize::MAX {
            continue;
        }
        let Some(values) = frame
            .iter()
            .map(|&i| mapped(i, map))
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        let from: [u32; 4] = values.try_into().expect("four slots");
        let to = b.frames[target].expect("tetrahedral presence is an atom color");
        // Every frame slot is a validated neighbor or the unique H sentinel.
        // Adjacency and H colors make this a full permutation, not raw @/@@.
        let same = a.mol.atom(AtomIdx(source as u32)).chirality
            == b.mol.atom(AtomIdx(target as u32)).chirality;
        if same != even(&from, &to) {
            return false;
        }
    }
    true
}
fn geometric(a: &Graph<'_>, b: &Graph<'_>, map: &[usize]) -> bool {
    for fact in &a.facts {
        let mut mapped_ends = Vec::new();
        for end in &fact.ends {
            let outer = a
                .mol
                .bond(end.reference)
                .other(end.atom)
                .expect("reference is incident");
            let (Some(center), Some(outer)) = (mapped(end.atom.0, map), mapped(outer.0, map))
            else {
                break;
            };
            mapped_ends.push((center as usize, outer as usize));
        }
        if mapped_ends.len() != 2 {
            continue;
        }
        let edge = b.adj[mapped_ends[0].0][mapped_ends[1].0].expect("mapped double bond");
        let target = &b.facts
            [b.fact_by_bond[edge.index.0 as usize].expect("geometric presence is an edge color")];
        let mut parity = target.ends[0].up ^ target.ends[1].up;
        for (center, outer) in mapped_ends {
            let end = target
                .ends
                .iter()
                .find(|e| e.atom.0 as usize == center)
                .expect("matching double endpoint");
            let reference = b
                .mol
                .bond(end.reference)
                .other(end.atom)
                .expect("reference is incident");
            // Choosing the other explicit substituent reverses the reference
            // side; reversing both cancels. Endpoint order itself is irrelevant.
            parity ^= outer != reference.0 as usize;
        }
        if parity != (fact.ends[0].up ^ fact.ends[1].up) {
            return false;
        }
    }
    true
}
fn compatible(a: &Graph<'_>, b: &Graph<'_>, map: &[usize], source: usize, target: usize) -> bool {
    if a.colors[source] != b.colors[target] {
        return false;
    }
    for (s, &t) in map.iter().enumerate() {
        if t == usize::MAX {
            continue;
        }
        if a.adj[source][s].map(|e| e.signature()) != b.adj[target][t].map(|e| e.signature()) {
            return false;
        }
    }
    true
}
fn search(
    a: &Graph<'_>,
    b: &Graph<'_>,
    map: &mut [usize],
    used: &mut [bool],
    budget: &mut Budget,
) -> Result<bool, Error> {
    let mut choice = None;
    let mut candidates = Vec::new();
    for source in 0..map.len() {
        if map[source] != usize::MAX {
            continue;
        }
        let mut possible = Vec::new();
        for (target, &taken) in used.iter().enumerate() {
            if taken {
                continue;
            }
            budget.probe()?;
            if compatible(a, b, map, source, target) {
                possible.push(target);
            }
        }
        if possible.is_empty() {
            return Ok(false);
        }
        if choice.is_none() || possible.len() < candidates.len() {
            choice = Some(source);
            candidates = possible;
        }
    }
    let Some(source) = choice else {
        return Ok(true);
    };
    for target in candidates {
        map[source] = target;
        used[target] = true;
        if tetrahedral(a, b, map) && geometric(a, b, map) && search(a, b, map, used, budget)? {
            return Ok(true);
        }
        map[source] = usize::MAX;
        used[target] = false;
    }
    Ok(false)
}
pub(super) fn find_mapping(
    a: &Graph<'_>,
    b: &Graph<'_>,
    budget: &mut Budget,
) -> Result<Option<Vec<usize>>, Error> {
    let mut map = vec![usize::MAX; a.colors.len()];
    let mut used = vec![false; map.len()];
    if search(a, b, &mut map, &mut used, budget)? {
        Ok(Some(map))
    } else {
        Ok(None)
    }
}

pub(super) fn find_mapping_anchored(
    a: &Graph<'_>,
    b: &Graph<'_>,
    budget: &mut Budget,
    anchor: usize,
) -> Result<Option<Vec<usize>>, Error> {
    let mut map = vec![usize::MAX; a.colors.len()];
    let mut used = vec![false; map.len()];
    if !compatible(a, b, &map, anchor, anchor) {
        return Ok(None);
    }
    map[anchor] = anchor;
    used[anchor] = true;
    if tetrahedral(a, b, &map)
        && geometric(a, b, &map)
        && search(a, b, &mut map, &mut used, budget)?
    {
        Ok(Some(map))
    } else {
        Ok(None)
    }
}
