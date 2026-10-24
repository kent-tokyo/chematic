use super::*;
use std::collections::{BTreeMap, HashSet};

fn number(text: &str) -> Result<u32, ReactionCxError> {
    if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(ReactionCxError::Invalid(format!(
            "invalid unsigned index/number '{text}'"
        )));
    }
    text.parse()
        .map_err(|_| ReactionCxError::Invalid("index/number exceeds u32".into()))
}

fn head(text: &str) -> Result<StereoGroupKind, ReactionCxError> {
    if text == "a" {
        return Ok(StereoGroupKind::Absolute);
    }
    let kind = if let Some(n) = text.strip_prefix('&') {
        StereoGroupKind::And(number(n)?)
    } else if let Some(n) = text.strip_prefix('o') {
        StereoGroupKind::Or(number(n)?)
    } else {
        return Err(ReactionCxError::Unsupported(format!("CX field '{text}'")));
    };
    group_key(&kind)?;
    Ok(kind)
}

fn fields(extension: &str, total: usize) -> Result<Vec<StereoGroup>, ReactionCxError> {
    if extension.is_empty() {
        return Ok(Vec::new());
    }
    let mut groups: Vec<StereoGroup> = Vec::new();
    let mut keys = HashSet::new();
    let mut members = HashSet::new();
    for token in extension.split(',') {
        let index = if let Some((label, index)) = token.split_once(':') {
            let kind = head(label)?;
            if !keys.insert(group_key(&kind)?) {
                return Err(ReactionCxError::Invalid(
                    "repeated group kind/number".into(),
                ));
            }
            groups.push(StereoGroup {
                kind,
                atom_indices: Vec::new(),
            });
            index
        } else {
            token
        };
        let index = number(index)?;
        if index as usize >= total || !members.insert(index) {
            return Err(ReactionCxError::Invalid(format!(
                "out-of-range or repeated atom index {index}"
            )));
        }
        let group = groups
            .last_mut()
            .ok_or_else(|| ReactionCxError::Invalid("atom index has no group".into()))?;
        group.atom_indices.push(AtomIdx(index));
    }
    Ok(groups)
}

fn split(input: &str) -> Result<(&str, Option<&str>), ReactionCxError> {
    let input = input.trim();
    let (base, extension) = if let Some((base, tail)) = input.split_once('|') {
        let extension = tail.strip_suffix('|').ok_or_else(|| {
            ReactionCxError::Invalid("CX block must end with | and have no title".into())
        })?;
        if extension.contains('|') {
            return Err(ReactionCxError::Invalid("multiple CX blocks".into()));
        }
        (base.trim_end(), Some(extension))
    } else {
        (input, None)
    };
    if base.chars().any(char::is_whitespace) {
        return Err(ReactionCxError::Unsupported(
            "reaction titles/embedded whitespace".into(),
        ));
    }
    Ok((base, extension))
}

pub(super) fn parse(
    input: &str,
    limits: &ReactionParseLimits,
) -> Result<Reaction, ReactionCxError> {
    if input.len() > limits.max_input_bytes {
        return Err(RxnError::ResourceLimit {
            resource: "input bytes including CX",
            actual: input.len(),
            limit: limits.max_input_bytes,
        }
        .into());
    }
    let (base, extension) = split(input)?;
    let mut reaction = crate::parse_reaction_with_limits(base, limits)?;
    let Some(extension) = extension else {
        return Ok(reaction);
    };
    let mut molecules: Vec<_> = reaction
        .reactants
        .iter_mut()
        .chain(&mut reaction.agents)
        .chain(&mut reaction.products)
        .collect();
    let mut total = 0usize;
    let starts: Vec<_> = molecules
        .iter()
        .map(|mol| {
            let start = total;
            total += mol.atom_count();
            start
        })
        .collect();
    let groups = fields(extension, total)?;
    let mut local = vec![Vec::new(); molecules.len()];
    for group in groups {
        let mut by_component = BTreeMap::<usize, Vec<AtomIdx>>::new();
        for atom in group.atom_indices {
            let component = starts.partition_point(|&start| start <= atom.0 as usize) - 1;
            let index = AtomIdx(atom.0 - starts[component] as u32);
            if !molecules[component].atom(index).chirality.is_tetrahedral() {
                return Err(ReactionCxError::Invalid(format!(
                    "group atom {} lacks tetrahedral stereo",
                    atom.0
                )));
            }
            by_component.entry(component).or_default().push(index);
        }
        for (component, indices) in by_component {
            local[component].push(StereoGroup::new(group.kind.clone(), indices));
        }
    }
    for (mol, groups) in molecules.iter_mut().zip(local) {
        mol.set_stereo_groups(groups);
    }
    Ok(reaction)
}
