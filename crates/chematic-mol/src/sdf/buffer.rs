//! Borrow complete, bounded SDF records already exposed by a BufRead source.
//! Incomplete records and resource-limit cases retain the original line reader.

use super::{MolParseError, SdfFileReader, SdfParseLimits, SdfRecord};
use std::io::{BufRead, ErrorKind};

fn record_span(
    input: &[u8],
    limits: SdfParseLimits,
    mut bytes_read: usize,
    records_read: usize,
) -> Option<(usize, usize)> {
    let mut start = 0;
    for newline in memchr::memchr_iter(b'\n', input) {
        let end = newline + 1;
        let line = &input[start..end];
        if line.len() > limits.max_line_bytes {
            return None;
        }
        bytes_read = bytes_read.saturating_add(line.len());
        if bytes_read > limits.max_input_bytes {
            return None;
        }
        let content = &line[..line.len() - 1];
        if content.strip_suffix(b"\r").unwrap_or(content) == b"$$$$" {
            if !input[..start].iter().all(u8::is_ascii_whitespace)
                && (start > limits.max_record_bytes || records_read >= limits.max_records)
            {
                return None;
            }
            return Some((start, end));
        }
        start = end;
    }
    None
}

impl<R: BufRead> SdfFileReader<R> {
    // None asks the caller to use its original line reader. Empty complete
    // records are consumed iteratively without allocating a result wrapper.
    pub(super) fn next_buffered(&mut self) -> Option<Result<SdfRecord, MolParseError>> {
        loop {
            let available = match self.reader.fill_buf() {
                Ok(available) => available,
                Err(error) if error.kind() == ErrorKind::Interrupted => continue,
                Err(error) => {
                    self.done = true;
                    return Some(Err(MolParseError::Io(error.to_string())));
                }
            };
            let (end, consumed) =
                record_span(available, self.limits, self.bytes_read, self.records_read)?;
            self.bytes_read = self.bytes_read.saturating_add(consumed);
            if available[..end].iter().all(u8::is_ascii_whitespace) {
                self.reader.consume(consumed);
                continue;
            }
            self.records_read += 1;
            let result = Self::parse_block(&available[..end], self.diagnostics);
            // The borrowed record is no longer needed after the complete parse.
            self.reader.consume(consumed);
            return Some(result);
        }
    }
}

#[cfg(test)]
#[path = "buffer_tests.rs"]
mod tests;
