use super::*;

impl EFrag {
    pub(super) fn merge_embedded_atoms(&mut self, embedded: &EFrag, common: &[usize]) {
        for (&atom_id, other_atom) in &embedded.ea {
            if !common.contains(&atom_id) {
                assign_eatom(self.at(atom_id), other_atom);
                if !other_atom.neighs.is_empty() && !self.attach.contains(&atom_id) {
                    self.attach.push(atom_id);
                }
                continue;
            }
            let atom = self.at(atom_id);
            if other_atom.cis_trans_nbr >= 0 {
                atom.cis_trans_nbr = other_atom.cis_trans_nbr;
                atom.normal = other_atom.normal;
                atom.ccw = other_atom.ccw;
            }
            if other_atom.angle > 0.0 {
                atom.angle = other_atom.angle;
                atom.nbr1 = other_atom.nbr1;
                atom.nbr2 = other_atom.nbr2;
            }
        }
    }

    pub(super) fn apply_open_angle_rotation(
        &mut self,
        atom1: usize,
        atom2: usize,
        pivot1: usize,
        pivot2: usize,
        kind: u8,
        cross: f64,
    ) {
        match kind {
            1 => {
                let angle = if cross < 0.0 { -ANGLE_OPEN } else { ANGLE_OPEN };
                let first = T2::rotation_about(self.at(pivot1).loc, angle);
                let second = T2::rotation_about(self.at(pivot2).loc, -angle);
                first.apply(&mut self.at(atom1).loc);
                second.apply(&mut self.at(atom2).loc);
            }
            2 => {
                let angle = if cross < 0.0 {
                    -2.0 * ANGLE_OPEN
                } else {
                    2.0 * ANGLE_OPEN
                };
                T2::rotation_about(self.at(pivot1).loc, angle).apply(&mut self.at(atom1).loc);
            }
            _ => {
                let angle = if cross < 0.0 {
                    2.0 * ANGLE_OPEN
                } else {
                    -2.0 * ANGLE_OPEN
                };
                T2::rotation_about(self.at(pivot2).loc, angle).apply(&mut self.at(atom2).loc);
            }
        }
    }
}
