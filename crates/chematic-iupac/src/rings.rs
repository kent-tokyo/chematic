//! Aromatic, saturated-heterocycle, and carbocyclic (cycloalkane/spiro/bicyclo) naming.

use crate::helpers::{alkane_base, alkane_suffix, best_benzene_locants, find_bridge_sizes};
use crate::{IupacError, Namer};
use chematic_core::{AtomIdx, BondOrder, implicit_hcount};
use chematic_perception::{RingSystemKind, find_ring_families, find_sssr};
use std::collections::{HashSet, VecDeque};

impl<'a> Namer<'a> {
    // -----------------------------------------------------------------------
    // Aromatic ring naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_aromatic_ring(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // All ring atoms must be aromatic.
        if !ring_atoms.iter().all(|&i| mol.atom(i).aromatic) {
            return Err(IupacError::NotSupported);
        }

        let n_n = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 7)
            .count();
        let n_o = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 8)
            .count();
        let n_s = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 16)
            .count();
        let sz = ring_atoms.len();

        // Case 1: Pure aromatic ring (no substituents).
        if ring_atoms.len() == mol.atom_count() {
            return match (sz, n_n, n_o, n_s) {
                (6, 0, 0, 0) => Ok("benzene".into()),
                (6, 1, 0, 0) => Ok("pyridine".into()),
                (6, 2, 0, 0) => Ok("pyrimidine".into()),
                (5, 0, 1, 0) => Ok("furan".into()),
                (5, 0, 0, 1) => Ok("thiophene".into()),
                (5, 1, 0, 0) => Ok("pyrrole".into()),
                (5, 2, 0, 0) => Ok("imidazole".into()),
                (10, 0, 0, 0) => Ok("naphthalene".into()),
                _ => Err(IupacError::NotSupported),
            };
        }

        // Case 2: Monosubstituted benzene (phenol, toluene, aniline, etc.)
        // Only support pure benzene ring (6 C, no N/O/S in ring).
        if sz == 6 && n_n == 0 && n_o == 0 && n_s == 0 {
            let sub_atoms: Vec<AtomIdx> = mol
                .atoms()
                .filter(|(i, _)| !ring_atoms.contains(i))
                .map(|(i, _)| i)
                .collect();
            return self.name_monosubstituted_benzene(ring_atoms, &sub_atoms);
        }

        Err(IupacError::NotSupported)
    }

    // -----------------------------------------------------------------------
    // Monosubstituted benzene naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_monosubstituted_benzene(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        sub_atoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Count how many ring C have substituents.
        let attach_count = ring_atoms
            .iter()
            .filter(|&&r| mol.neighbors(r).any(|(nb, _)| !ring_atoms.contains(&nb)))
            .count();
        if attach_count == 2 {
            return self.name_disubstituted_benzene(ring_atoms, sub_atoms);
        }
        if attach_count == 3 {
            return self.name_trisubstituted_benzene(ring_atoms);
        }
        if attach_count != 1 {
            return Err(IupacError::NotSupported);
        }

        // Classify substituent by element counts + bond types.
        let mut n_c = 0usize;
        let mut n_n = 0usize;
        let mut n_o = 0usize;
        let mut n_hal = 0usize;
        let mut halogen_an = 0u8;
        for &a in sub_atoms {
            match mol.atom(a).element.atomic_number() {
                6 => n_c += 1,
                7 => n_n += 1,
                8 => n_o += 1,
                1 => {}
                an @ (9 | 17 | 35 | 53) => {
                    n_hal += 1;
                    halogen_an = an;
                }
                _ => return Err(IupacError::NotSupported),
            }
        }

        let sub_set: HashSet<AtomIdx> = sub_atoms.iter().copied().collect();
        let has_triple = mol.bonds().any(|(_, b)| {
            b.order == BondOrder::Triple
                && (sub_set.contains(&b.atom1) || sub_set.contains(&b.atom2))
        });
        let has_double = mol.bonds().any(|(_, b)| {
            b.order == BondOrder::Double
                && (sub_set.contains(&b.atom1) || sub_set.contains(&b.atom2))
        });

        match (n_c, n_n, n_o, n_hal, has_double, has_triple) {
            // Phenol: c1ccccc1O
            (0, 0, 1, 0, false, false) => Ok("phenol".into()),
            // Aniline: c1ccccc1N
            (0, 1, 0, 0, false, false) => Ok("aniline".into()),
            // Halo-benzenes
            (0, 0, 0, 1, false, false) => {
                let prefix = match halogen_an {
                    9 => "fluoro",
                    17 => "chloro",
                    35 => "bromo",
                    53 => "iodo",
                    _ => return Err(IupacError::NotSupported),
                };
                Ok(format!("{prefix}benzene"))
            }
            // Toluene: c1ccccc1C (one CH3)
            (1, 0, 0, 0, false, false) => Ok("toluene".into()),
            // Benzaldehyde: c1ccccc1C=O (n_c=1, n_o=1, has_double)
            (1, 0, 1, 0, true, false) => Ok("benzaldehyde".into()),
            // Benzoic acid: c1ccccc1C(=O)O (n_c=1, n_o=2, has_double)
            (1, 0, 2, 0, true, false) => Ok("benzoic acid".into()),
            // Benzonitrile: c1ccccc1C#N (n_c=1, n_n=1, has_triple)
            (1, 1, 0, 0, false, true) => Ok("benzonitrile".into()),
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Disubstituted benzene naming (e.g., "4-chlorophenol")
    // -----------------------------------------------------------------------

    pub(crate) fn name_disubstituted_benzene(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        _sub_atoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;

        // Identify the two ring C attachment points and their substituent sets.
        let attach_points: Vec<AtomIdx> = ring_atoms
            .iter()
            .filter(|&&r| mol.neighbors(r).any(|(nb, _)| !ring_atoms.contains(&nb)))
            .copied()
            .collect();
        if attach_points.len() != 2 {
            return Err(IupacError::NotSupported);
        }

        // Compute ring distance (shortest path within ring) between the two attachment points.
        let ring_dist = {
            let ring_vec: Vec<AtomIdx> = ring_atoms.iter().copied().collect();
            let mut dist = usize::MAX;
            // BFS within the ring
            let mut queue = VecDeque::new();
            let mut visited: HashSet<AtomIdx> = HashSet::new();
            queue.push_back((attach_points[0], 0usize));
            visited.insert(attach_points[0]);
            while let Some((cur, d)) = queue.pop_front() {
                if cur == attach_points[1] {
                    dist = d;
                    break;
                }
                for (nb, _) in mol.neighbors(cur) {
                    if ring_atoms.contains(&nb) && visited.insert(nb) {
                        queue.push_back((nb, d + 1));
                    }
                }
            }
            // Take minimum of this and the longer path
            dist.min(ring_vec.len() - dist)
        };

        // Classify each substituent group.
        let classify_sub = |attach: AtomIdx| -> Option<(&str, bool)> {
            // Returns (substituent_name, is_principal)
            // is_principal: true if this substituent determines the compound root name
            // Collect sub atoms for this attachment (unused for now, just for documentation)
            // Simple: just look at atoms directly bonded to attach that are not in ring
            let direct: Vec<AtomIdx> = mol
                .neighbors(attach)
                .filter(|(nb, _)| !ring_atoms.contains(nb))
                .map(|(nb, _)| nb)
                .collect();
            if direct.is_empty() {
                return None;
            }
            let first = direct[0];
            let an = mol.atom(first).element.atomic_number();
            match an {
                8 if !mol
                    .neighbors(first)
                    .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double) =>
                {
                    Some(("hydroxy", true)) // -OH → phenol as principal
                }
                7 if implicit_hcount(mol, first) > 0 => Some(("amino", true)), // -NH2 → aniline
                6 => Some(("methyl", false)), // -CH3 → toluene substituent
                17 => Some(("chloro", false)),
                35 => Some(("bromo", false)),
                9 => Some(("fluoro", false)),
                53 => Some(("iodo", false)),
                _ => None,
            }
        };

        let sub_a = classify_sub(attach_points[0]);
        let sub_b = classify_sub(attach_points[1]);

        let (sub_a, sub_b) = match (sub_a, sub_b) {
            (Some(a), Some(b)) => (a, b),
            _ => return Err(IupacError::NotSupported),
        };

        // Determine locant prefix (1,2= ortho, 1,3= meta, 1,4= para for 6-ring).
        let pos2 = ring_dist + 1; // position of the second substituent from first

        // Build name: principal group determines root, non-principal is prefix.
        let (prefix_sub, root_name) = if sub_a.1 {
            // sub_a is principal (phenol/aniline): prefix comes from sub_b
            let root = match sub_a.0 {
                "hydroxy" => "phenol",
                "amino" => "aniline",
                _ => return Err(IupacError::NotSupported),
            };
            (sub_b.0, root)
        } else if sub_b.1 {
            let root = match sub_b.0 {
                "hydroxy" => "phenol",
                "amino" => "aniline",
                _ => return Err(IupacError::NotSupported),
            };
            (sub_a.0, root)
        } else {
            // Neither is principal — both are substituents on benzene.
            // Alphabetically first substituent gets locant 1.
            let (s1, s2) = if sub_a.0 <= sub_b.0 {
                (sub_a.0, sub_b.0)
            } else {
                (sub_b.0, sub_a.0)
            };
            return if s1 == s2 {
                Ok(format!("1,{}-di{}benzene", pos2, s1))
            } else {
                Ok(format!("1-{}-{}-{}benzene", s1, pos2, s2))
            };
        };

        Ok(format!("{}-{}{}", pos2, prefix_sub, root_name))
    }

    // -----------------------------------------------------------------------
    // Trisubstituted benzene naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_trisubstituted_benzene(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let attach_points: Vec<AtomIdx> = ring_atoms
            .iter()
            .filter(|&&r| mol.neighbors(r).any(|(nb, _)| !ring_atoms.contains(&nb)))
            .copied()
            .collect();
        if attach_points.len() != 3 {
            return Err(IupacError::NotSupported);
        }

        let locant_map = best_benzene_locants(mol, ring_atoms, &attach_points);

        // Classify each substituent.
        let mut sub_list: Vec<(usize, String)> = Vec::new();
        for &(locant, attach) in &locant_map {
            let sub = self
                .classify_benzene_sub_simple(attach, ring_atoms)
                .ok_or(IupacError::NotSupported)?;
            sub_list.push((locant, sub));
        }

        // Sort alphabetically by substituent name, then numerically by locant.
        sub_list.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

        // Group identical substituents for di/tri multiplier.
        let mut groups: Vec<(String, Vec<usize>)> = Vec::new();
        for (locant, name) in sub_list {
            if let Some(last) = groups.last_mut()
                && last.0 == name
            {
                last.1.push(locant);
                continue;
            }
            groups.push((name, vec![locant]));
        }

        let mut parts: Vec<String> = Vec::new();
        for (name, mut locs) in groups {
            locs.sort_unstable();
            let locant_str = locs
                .iter()
                .map(|l| l.to_string())
                .collect::<Vec<_>>()
                .join(",");
            let mult = match locs.len() {
                1 => String::new(),
                2 => "di".to_string(),
                3 => "tri".to_string(),
                _ => return Err(IupacError::NotSupported),
            };
            parts.push(format!("{}-{}{}", locant_str, mult, name));
        }

        Ok(format!("{}benzene", parts.join("-")))
    }

    /// Classify a single benzene substituent by element type (for trisubstituted naming).
    pub(crate) fn classify_benzene_sub_simple(
        &self,
        attach: AtomIdx,
        ring_atoms: &HashSet<AtomIdx>,
    ) -> Option<String> {
        let mol = self.mol;
        let direct: Vec<AtomIdx> = mol
            .neighbors(attach)
            .filter(|(nb, _)| !ring_atoms.contains(nb))
            .map(|(nb, _)| nb)
            .collect();
        if direct.is_empty() {
            return None;
        }
        let first = direct[0];
        match mol.atom(first).element.atomic_number() {
            6 => Some("methyl".to_string()),
            7 => Some("amino".to_string()),
            8 => Some("hydroxy".to_string()),
            9 => Some("fluoro".to_string()),
            17 => Some("chloro".to_string()),
            35 => Some("bromo".to_string()),
            53 => Some("iodo".to_string()),
            _ => None,
        }
    }

    // -----------------------------------------------------------------------
    // Cycloalkane naming
    // -----------------------------------------------------------------------

    pub(crate) fn name_cycloalkane(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        if carbons.iter().any(|&c| mol.atom(c).aromatic) {
            return Err(IupacError::NotSupported);
        }
        // All carbons in ring: unsubstituted cycloalkane.
        if ring_atoms.len() == carbons.len() {
            return Ok(format!("cyclo{}", alkane_suffix(ring_atoms.len())));
        }
        let outside: Vec<AtomIdx> = carbons
            .iter()
            .filter(|&&c| !ring_atoms.contains(&c))
            .copied()
            .collect();

        let is_terminal_methyl = |sub: AtomIdx| -> bool {
            mol.neighbors(sub)
                .filter(|(nb, _)| {
                    mol.atom(*nb).element.atomic_number() == 6 && !ring_atoms.contains(nb)
                })
                .count()
                == 0
        };

        if outside.len() == 1 && is_terminal_methyl(outside[0]) {
            return Ok(format!("methylcyclo{}", alkane_suffix(ring_atoms.len())));
        }

        if outside.len() == 2 && is_terminal_methyl(outside[0]) && is_terminal_methyl(outside[1]) {
            let att_a = mol
                .neighbors(outside[0])
                .find(|(nb, _)| ring_atoms.contains(nb))
                .map(|(nb, _)| nb)
                .ok_or(IupacError::NotSupported)?;
            let att_b = mol
                .neighbors(outside[1])
                .find(|(nb, _)| ring_atoms.contains(nb))
                .map(|(nb, _)| nb)
                .ok_or(IupacError::NotSupported)?;
            // BFS shortest path within ring.
            let raw_dist = {
                let mut dist = 0usize;
                let mut queue = VecDeque::new();
                let mut visited: HashSet<AtomIdx> = HashSet::new();
                queue.push_back((att_a, 0usize));
                visited.insert(att_a);
                'bfs: while let Some((cur, d)) = queue.pop_front() {
                    if cur == att_b {
                        dist = d;
                        break 'bfs;
                    }
                    for (nb, _) in mol.neighbors(cur) {
                        if ring_atoms.contains(&nb) && visited.insert(nb) {
                            queue.push_back((nb, d + 1));
                        }
                    }
                }
                dist
            };
            let ring_dist = raw_dist.min(ring_atoms.len() - raw_dist);
            return Ok(format!(
                "1,{}-dimethylcyclo{}",
                ring_dist + 1,
                alkane_suffix(ring_atoms.len())
            ));
        }

        Err(IupacError::NotSupported)
    }

    // -----------------------------------------------------------------------
    // Polycyclic naming: spiro and bicyclo (bridged)
    // -----------------------------------------------------------------------

    pub(crate) fn name_polycyclic(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let sssr = find_sssr(mol);
        let families = find_ring_families(mol, &sssr);

        // Must be a single ring family.
        if families.len() != 1 {
            return Err(IupacError::NotSupported);
        }
        let family = &families[0];

        match family.kind {
            RingSystemKind::Spiro => self.name_spiro(ring_atoms, carbons, &sssr),
            RingSystemKind::Bridged => self.name_bicyclo(ring_atoms, carbons, &sssr),
            _ => Err(IupacError::NotSupported),
        }
    }

    /// Name a simple spiro compound: spiro[a.b]alkane.
    pub(crate) fn name_spiro(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
        sssr: &chematic_perception::RingSet,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let rings = sssr.rings();

        // Must have exactly 2 rings, no substituents outside the ring system.
        if rings.len() != 2 {
            return Err(IupacError::NotSupported);
        }
        if ring_atoms.len() != carbons.len() {
            return Err(IupacError::NotSupported);
        }
        // All ring atoms must be sp3 carbons.
        if ring_atoms.iter().any(|&a| {
            mol.atom(a).aromatic
                || mol.bonds().any(|(_, b)| {
                    (b.atom1 == a || b.atom2 == a)
                        && !matches!(b.order, chematic_core::BondOrder::Single)
                })
        }) {
            return Err(IupacError::NotSupported);
        }

        // Find spiro atom (shared by both rings).
        let shared: Vec<AtomIdx> = rings[0]
            .iter()
            .filter(|a| rings[1].contains(a))
            .copied()
            .collect();
        if shared.len() != 1 {
            return Err(IupacError::NotSupported);
        }
        let _spiro_atom = shared[0];

        // Bridge sizes = ring_size - 1 (excluding the spiro atom).
        let mut bridges = [rings[0].len() - 1, rings[1].len() - 1];
        bridges.sort_unstable();
        let total_c = ring_atoms.len();

        Ok(format!(
            "spiro[{}.{}]{}",
            bridges[0],
            bridges[1],
            alkane_suffix(total_c)
        ))
    }

    /// Name a bicyclic bridged compound: bicyclo[x.y.z]alkane.
    pub(crate) fn name_bicyclo(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
        _sssr: &chematic_perception::RingSet,
    ) -> Result<String, IupacError> {
        let mol = self.mol;

        // Must have no substituents outside the ring system.
        if ring_atoms.len() != carbons.len() {
            return Err(IupacError::NotSupported);
        }
        // All ring atoms must be sp3 (no double bonds, no aromatic).
        if ring_atoms.iter().any(|&a| {
            mol.atom(a).aromatic
                || mol.bonds().any(|(_, b)| {
                    (b.atom1 == a || b.atom2 == a)
                        && !matches!(b.order, chematic_core::BondOrder::Single)
                })
        }) {
            return Err(IupacError::NotSupported);
        }

        // Bridgehead atoms: ring atoms with degree ≥ 3 in the ring subgraph.
        let bridgeheads: Vec<AtomIdx> = ring_atoms
            .iter()
            .copied()
            .filter(|&a| {
                mol.neighbors(a)
                    .filter(|(nb, _)| ring_atoms.contains(nb))
                    .count()
                    >= 3
            })
            .collect();
        if bridgeheads.len() != 2 {
            return Err(IupacError::NotSupported);
        }
        let bh0 = bridgeheads[0];
        let bh1 = bridgeheads[1];

        // Find bridge sizes: simple paths between bridgeheads through non-bridgehead ring atoms.
        let bridge_sizes = find_bridge_sizes(mol, bh0, bh1, ring_atoms);
        if bridge_sizes.len() != 3 {
            return Err(IupacError::NotSupported);
        }

        let mut bridges = bridge_sizes;
        bridges.sort_unstable_by(|a, b| b.cmp(a)); // descending
        let total_c = ring_atoms.len();

        Ok(format!(
            "bicyclo[{}.{}.{}]{}",
            bridges[0],
            bridges[1],
            bridges[2],
            alkane_suffix(total_c)
        ))
    }

    // -----------------------------------------------------------------------
    // Cycloalkanol naming (cyclopentanol, cyclohexanol, ...)
    // -----------------------------------------------------------------------

    pub(crate) fn name_cycloalkanol(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
        carbons: &[AtomIdx],
        o_atoms: &[AtomIdx],
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        // Only one OH substituent.
        if o_atoms.len() != 1 {
            return Err(IupacError::NotSupported);
        }
        let o_idx = o_atoms[0];
        // O must be single-bond –OH (not carbonyl).
        if mol
            .neighbors(o_idx)
            .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double)
        {
            return Err(IupacError::NotSupported);
        }
        // O must have implicit H.
        if implicit_hcount(mol, o_idx) == 0 {
            return Err(IupacError::NotSupported);
        }
        // No exocyclic carbons (unsubstituted ring + OH only).
        let exo_c = carbons
            .iter()
            .filter(|&&c| !ring_atoms.contains(&c))
            .count();
        if exo_c > 0 {
            return Err(IupacError::NotSupported);
        }
        Ok(format!("cyclo{}ol", alkane_base(ring_atoms.len())))
    }

    // -----------------------------------------------------------------------
    // Saturated N-heterocycles (piperidine, pyrrolidine, azetidine)
    // -----------------------------------------------------------------------

    pub(crate) fn name_aza_ring(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let sz = ring_atoms.len();
        let n_n = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 7)
            .count();
        // Only support pure monocyclic rings (no substituents).
        if mol.atom_count() != sz {
            return Err(IupacError::NotSupported);
        }
        match (sz, n_n) {
            (4, 1) => Ok("azetidine".into()),
            (5, 1) => Ok("pyrrolidine".into()),
            (6, 1) => Ok("piperidine".into()),
            (6, 2) => Ok("piperazine".into()),
            (7, 1) => Ok("azepane".into()),
            _ => Err(IupacError::NotSupported),
        }
    }

    // -----------------------------------------------------------------------
    // Saturated N+O heterocycles (morpholine) and N+N (piperazine)
    // -----------------------------------------------------------------------

    pub(crate) fn name_oxaaza_ring(
        &self,
        ring_atoms: &HashSet<AtomIdx>,
    ) -> Result<String, IupacError> {
        let mol = self.mol;
        let sz = ring_atoms.len();
        if mol.atom_count() != sz {
            return Err(IupacError::NotSupported);
        }
        let n_n = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 7)
            .count();
        let n_o = ring_atoms
            .iter()
            .filter(|&&i| mol.atom(i).element.atomic_number() == 8)
            .count();
        // Morpholine: 6-membered ring, 4C + 1N + 1O
        if sz == 6 && n_n == 1 && n_o == 1 {
            return Ok("morpholine".into());
        }
        // Piperazine: 6-membered ring, 4C + 2N
        if sz == 6 && n_n == 2 && n_o == 0 {
            return Ok("piperazine".into());
        }
        Err(IupacError::NotSupported)
    }
}
