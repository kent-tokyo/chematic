use super::*;
use fields::{cell, data, integer, layout, number, rows};
pub(super) struct Frame {
    pub positions: Vec<[f64; 3]>,
    pub identities: Option<Vec<Gromos96Identity>>,
    pub cell: Option<[[f64; 3]; 3]>,
    velocity: bool,
}
#[derive(Default)]
struct Frames {
    ready: Vec<Frame>,
    current: Option<Frame>,
    pending_time: bool,
    atoms: usize,
}
impl Frames {
    fn finish(&mut self) {
        if let Some(frame) = self.current.take() {
            self.ready.push(frame);
        }
    }
    fn accept(
        &mut self,
        block: &Gromos96Block,
        limits: &Gromos96ParseLimits,
    ) -> Result<(), Gromos96Error> {
        match block.kind.as_str() {
            "TIMESTEP" => {
                if self.pending_time {
                    return Err(error(0, "TIMESTEP without positions"));
                }
                self.finish();
                let lines = data(block);
                if lines.len() != 1 || !lines[0].is_ascii() || lines[0].len() != 30 {
                    return Err(error(0, "expected one 30-column TIMESTEP row"));
                }
                integer(&lines[0][..15], true)?;
                number(&lines[0][15..])?;
                self.pending_time = true;
            }
            "POSITION" | "POSITIONRED" => {
                self.finish();
                if self.ready.len() >= limits.max_frames {
                    return Err(error(0, "frame limit exceeded"));
                }
                let (positions, identities) = rows(block, limits)?;
                self.atoms = self
                    .atoms
                    .checked_add(positions.len())
                    .ok_or_else(|| error(0, "atom count overflow"))?;
                if self.atoms > limits.max_atoms {
                    return Err(error(0, "total position-row limit exceeded"));
                }
                self.current = Some(Frame {
                    positions,
                    identities,
                    cell: None,
                    velocity: false,
                });
                self.pending_time = false;
            }
            "VELOCITY" | "VELOCITYRED" => self.velocity(block, limits)?,
            "BOX" => {
                let frame = self
                    .current
                    .as_mut()
                    .ok_or_else(|| error(0, "BOX before positions"))?;
                if frame.cell.is_some() {
                    return Err(error(0, "duplicate BOX"));
                }
                frame.cell = Some(cell(block)?);
            }
            _ => return Err(error(0, "unsupported or repeated block")),
        }
        Ok(())
    }
    fn velocity(
        &mut self,
        block: &Gromos96Block,
        limits: &Gromos96ParseLimits,
    ) -> Result<(), Gromos96Error> {
        let frame = self
            .current
            .as_mut()
            .ok_or_else(|| error(0, "velocity before positions"))?;
        if frame.velocity || frame.cell.is_some() {
            return Err(error(0, "duplicate or out-of-order velocity block"));
        }
        let (values, identities) = rows(block, limits)?;
        if values.len() != frame.positions.len() || identities != frame.identities {
            return Err(error(
                0,
                "velocity row count/identity does not match positions",
            ));
        }
        frame.velocity = true;
        Ok(())
    }
}
pub(super) fn collect(
    document: &Gromos96Document,
    limits: &Gromos96ParseLimits,
) -> Result<Vec<Frame>, Gromos96Error> {
    if document.blocks.first().map(|b| b.kind.as_str()) != Some("TITLE")
        || document.blocks.len() > limits.max_frames.saturating_mul(4).saturating_add(1)
    {
        return Err(error(0, "missing TITLE or block limit exceeded"));
    }
    let mut frames = Frames::default();
    for (i, block) in document.blocks.iter().enumerate() {
        for line in &block.leading {
            check_line(line, 0, limits)?;
            if !layout(line) {
                return Err(error(0, "non-layout leading record"));
            }
        }
        for line in &block.lines {
            check_line(line, 0, limits)?;
            if line.trim() == "END" {
                return Err(error(0, "END cannot occur in edited block payload"));
            }
        }
        if i != 0 {
            frames.accept(block, limits)?;
        }
    }
    for line in &document.trailing {
        check_line(line, 0, limits)?;
        if !layout(line) {
            return Err(error(0, "non-layout trailing record"));
        }
    }
    if frames.pending_time {
        return Err(error(0, "trailing TIMESTEP without positions"));
    }
    frames.finish();
    if frames.ready.is_empty() {
        return Err(error(0, "no coordinate frames"));
    }
    Ok(frames.ready)
}
