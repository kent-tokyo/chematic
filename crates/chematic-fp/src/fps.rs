//! FPS ("Fingerprint file format") read/write — the plain-text fingerprint
//! interchange format popularized by chemfp/OpenBabel.
//!
//! A `.fps` file is a `#`-prefixed metadata header followed by one
//! tab-separated `<hex-fingerprint>\t<identifier>` record per line:
//!
//! ```text
//! #FPS1
//! #num_bits=8
//! #type=Test/1
//! #software=chematic-fp/0.15.0
//! 01\tfirst
//! 0a\tsecond
//! c0\tthird
//! ```
//!
//! # Hex bit-ordering convention
//!
//! Verified against the real chemfp FPS format specification
//! (<https://chemfp.com/fps_format/>, fetched 2026-08-14):
//!
//! - **Byte order**: the left-most 2 hex characters encode fingerprint byte 0
//!   (bits 0-7), the next 2 encode byte 1 (bits 8-15), and so on.
//! - **Bit order within a byte**: "big-endian" in chemfp's own terminology,
//!   meaning the byte's hex value equals the standard binary value of that
//!   byte with bit 0 as the *lowest* (rightmost) bit — hex `"01"` is byte
//!   value 1 (bit 0 set), `"0a"` is byte value 10 = `0b1010` (bits 1 and 3
//!   set), `"c0"` is byte value 192 = `0b11000000` (bits 6 and 7 set). These
//!   three exact examples from the spec are used as positive-control test
//!   fixtures below.
//! - **Padding**: fingerprints whose bit count isn't a multiple of 8 are
//!   padded with zero bits at the *most significant* end of the last byte
//!   (i.e. at the high end of the highest-index byte).
//!
//! This lines up exactly with [`BitVecN`]/[`BitVec2048`]'s own bit
//! numbering (`bit n` lives at `(word n/64, offset n%64)`, tested via each
//! word's `1u64 << (bit % 64)`): global fingerprint bit `n` maps to byte
//! `n/8`, bit position `n%8` within that byte, LSB-first — so encoding walks
//! `get(byte*8..byte*8+8)` into a `u8` with `bit_in_byte` as the shift
//! amount, which is precisely the chemfp convention above.

use std::io::{BufRead, Write};

use crate::bitvec::{BitVec2048, BitVecN};

/// Header metadata for an FPS file.
///
/// `num_bits`, `type`, `software`, `source`, and `comment` are modeled
/// explicitly (each is a standard/conventional FPS header field). Every
/// other `#`-prefixed header line — a custom field, a stray `#FPS1` version
/// line, or anything this reader doesn't specifically recognize — is
/// preserved verbatim in `extra` (without its leading `#`) so round-tripping
/// a file produced by another tool never silently drops header content.
///
/// Field order is *not* preserved relative to `extra` on write-back:
/// recognized fields are always written first, in a fixed canonical order,
/// followed by `extra` lines in their original relative order — with one
/// exception. A version line (`#FPS<digits>`, e.g. `#FPS1`) is required by
/// the spec to be the very first line of the file, so any such line found
/// in `extra` is hoisted ahead of `#num_bits=` on write, regardless of
/// where it appeared in `extra`'s own order; every other `extra` line stays
/// after the recognized fields. This does not affect fingerprint round-trip
/// fidelity (only the informational header), and keeps this struct a plain
/// field bag instead of a single ordered enum-of-header-line-kinds list.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct FpsHeader {
    /// Bits per fingerprint. Required — reading a file with no `#num_bits=`
    /// line fails with [`FpsError::MissingNumBits`].
    pub num_bits: usize,
    /// Fingerprint algorithm/parameters label, e.g. `"ECFP4"`.
    pub fp_type: Option<String>,
    /// Generating software, e.g. `"chematic-fp/0.15.0"`.
    pub software: Option<String>,
    /// Input file/dataset the fingerprints were computed from. Repeatable
    /// per the FPS spec, so every `#source=` line is kept.
    pub source: Vec<String>,
    /// Free-text `#comment=` lines, in file order.
    pub comments: Vec<String>,
    /// Every other `#`-prefixed header line, verbatim minus the leading
    /// `#` (covers unrecognized custom keys and non-`key=value` lines like
    /// `#FPS1`).
    pub extra: Vec<String>,
}

impl FpsHeader {
    /// A bare header with just the required bit count set.
    pub fn new(num_bits: usize) -> Self {
        Self {
            num_bits,
            ..Default::default()
        }
    }

    /// A header for fingerprints computed by this crate: stamps
    /// `software=chematic-fp/<crate version>` and `type=<fp_type>` so
    /// provenance isn't left blank when chematic-fp itself is the producer.
    pub fn for_chematic(num_bits: usize, fp_type: impl Into<String>) -> Self {
        Self {
            num_bits,
            fp_type: Some(fp_type.into()),
            software: Some(format!("chematic-fp/{}", env!("CARGO_PKG_VERSION"))),
            ..Default::default()
        }
    }

    /// Set (overwrite) the fingerprint type label.
    pub fn with_type(mut self, fp_type: impl Into<String>) -> Self {
        self.fp_type = Some(fp_type.into());
        self
    }

    /// Set (overwrite) the software label.
    pub fn with_software(mut self, software: impl Into<String>) -> Self {
        self.software = Some(software.into());
        self
    }

    /// Append a `#source=` line.
    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.source.push(source.into());
        self
    }

    /// Append a `#comment=` line.
    pub fn with_comment(mut self, comment: impl Into<String>) -> Self {
        self.comments.push(comment.into());
        self
    }

    fn absorb_line(&mut self, rest: &str) -> Result<(), FpsError> {
        match split_header_kv(rest) {
            Some(("num_bits", value)) => {
                self.num_bits = value
                    .parse()
                    .map_err(|_| FpsError::InvalidNumBits(value.to_string()))?;
            }
            Some(("type", value)) => self.fp_type = Some(value.to_string()),
            Some(("software", value)) => self.software = Some(value.to_string()),
            Some(("source", value)) => self.source.push(value.to_string()),
            Some(("comment", value)) => self.comments.push(value.to_string()),
            _ => self.extra.push(rest.to_string()),
        }
        Ok(())
    }

    fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        // The spec requires the version line (`#FPS1`) to be the file's
        // very first line, so hoist any such line ahead of #num_bits= even
        // though it's carried in `extra` -- everything else in `extra`
        // keeps its position after the recognized fields.
        for extra in &self.extra {
            if is_version_line(extra) {
                writeln!(w, "#{extra}")?;
            }
        }
        writeln!(w, "#num_bits={}", self.num_bits)?;
        if let Some(t) = &self.fp_type {
            writeln!(w, "#type={t}")?;
        }
        if let Some(s) = &self.software {
            writeln!(w, "#software={s}")?;
        }
        for src in &self.source {
            writeln!(w, "#source={src}")?;
        }
        for c in &self.comments {
            writeln!(w, "#comment={c}")?;
        }
        for extra in &self.extra {
            if !is_version_line(extra) {
                writeln!(w, "#{extra}")?;
            }
        }
        Ok(())
    }
}

/// Is `s` (a header line's content after the leading `#`) a version line
/// like `FPS1`? The spec's version line has no `=` and must be the file's
/// first line.
fn is_version_line(s: &str) -> bool {
    s.strip_prefix("FPS")
        .is_some_and(|rest| !rest.is_empty() && rest.bytes().all(|b| b.is_ascii_digit()))
}

/// Split a header line's post-`#` content into `(key, value)` if it matches
/// the FPS spec's `key=value` grammar (`key` matches `[A-Za-z_][A-Za-z0-9_]*`).
/// Lines that don't match (no `=`, or an invalid key — including the bare
/// `#FPS1` version line) return `None` and are preserved verbatim instead.
fn split_header_kv(rest: &str) -> Option<(&str, &str)> {
    let eq = rest.find('=')?;
    let (key, value) = (&rest[..eq], &rest[eq + 1..]);
    let mut chars = key.chars();
    let first_ok = chars
        .next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == '_');
    let rest_ok = chars.all(|c| c.is_ascii_alphanumeric() || c == '_');
    (first_ok && rest_ok).then_some((key, value))
}

/// One parsed FPS record: an identifier and its fingerprint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FpsRecord {
    pub id: String,
    pub fingerprint: BitVecN,
}

/// Why reading or parsing an FPS stream failed.
#[derive(Debug)]
pub enum FpsError {
    /// Underlying I/O failure.
    Io(std::io::Error),
    /// No `#num_bits=` header line was seen before data started (or before
    /// EOF).
    MissingNumBits,
    /// `#num_bits=` value was not a valid positive integer.
    InvalidNumBits(String),
    /// A data line has no tab separator between the hex fingerprint and the
    /// identifier.
    MalformedRecord(String),
    /// The hex fingerprint's character count doesn't match what
    /// `#num_bits=` requires (`ceil(num_bits / 8) * 2` hex chars).
    HexLengthMismatch {
        expected_hex_chars: usize,
        got_hex_chars: usize,
        num_bits: usize,
    },
    /// A non-hex-digit character appeared in the fingerprint field.
    InvalidHexDigit(char),
}

impl std::fmt::Display for FpsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FpsError::Io(e) => write!(f, "fps: io error: {e}"),
            FpsError::MissingNumBits => {
                write!(f, "fps: missing required #num_bits= header line")
            }
            FpsError::InvalidNumBits(v) => write!(f, "fps: invalid #num_bits= value: {v:?}"),
            FpsError::MalformedRecord(line) => {
                write!(
                    f,
                    "fps: record line has no <hex>\\t<id> separator: {line:?}"
                )
            }
            FpsError::HexLengthMismatch {
                expected_hex_chars,
                got_hex_chars,
                num_bits,
            } => write!(
                f,
                "fps: hex fingerprint has {got_hex_chars} chars, expected {expected_hex_chars} for num_bits={num_bits}"
            ),
            FpsError::InvalidHexDigit(c) => write!(f, "fps: invalid hex digit {c:?}"),
        }
    }
}

impl std::error::Error for FpsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            FpsError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for FpsError {
    fn from(e: std::io::Error) -> Self {
        FpsError::Io(e)
    }
}

fn strip_newline(s: &str) -> &str {
    s.trim_end_matches(['\n', '\r'])
}

/// Streaming FPS reader over any [`BufRead`] source.
///
/// Reads and parses the `#`-prefixed header block on construction, then
/// yields one [`FpsRecord`] per data line via [`Iterator`] — records are
/// parsed one at a time as the caller pulls them, never buffering the whole
/// file. Iteration stops for good after the first error (a later line in a
/// corrupt file can't be trusted either).
#[derive(Debug)]
pub struct FpsReader<R: BufRead> {
    inner: R,
    header: FpsHeader,
    pending: Option<String>,
    done: bool,
}

impl<R: BufRead> FpsReader<R> {
    /// Read and parse the header block, leaving the stream positioned at
    /// the first data line.
    pub fn new(mut inner: R) -> Result<Self, FpsError> {
        let mut header = FpsHeader::default();
        let mut pending = None;
        loop {
            let mut buf = String::new();
            let n = inner.read_line(&mut buf)?;
            if n == 0 {
                break;
            }
            let line = strip_newline(&buf);
            match line.strip_prefix('#') {
                Some(rest) => header.absorb_line(rest)?,
                None => {
                    pending = Some(line.to_string());
                    break;
                }
            }
        }
        if header.num_bits == 0 {
            return Err(FpsError::MissingNumBits);
        }
        Ok(Self {
            inner,
            header,
            pending,
            done: false,
        })
    }

    /// The parsed header block.
    pub fn header(&self) -> &FpsHeader {
        &self.header
    }
}

impl<R: BufRead> Iterator for FpsReader<R> {
    type Item = Result<FpsRecord, FpsError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.done {
                return None;
            }
            let line = match self.pending.take() {
                Some(line) => line,
                None => {
                    let mut buf = String::new();
                    match self.inner.read_line(&mut buf) {
                        Ok(0) => {
                            self.done = true;
                            return None;
                        }
                        Ok(_) => strip_newline(&buf).to_string(),
                        Err(e) => {
                            self.done = true;
                            return Some(Err(FpsError::Io(e)));
                        }
                    }
                }
            };
            // Tolerate stray/trailing blank lines rather than erroring —
            // some writers leave one at EOF.
            if line.trim().is_empty() {
                continue;
            }
            return Some(match parse_record(&line, self.header.num_bits) {
                Ok(record) => Ok(record),
                Err(e) => {
                    self.done = true;
                    Err(e)
                }
            });
        }
    }
}

fn parse_record(line: &str, num_bits: usize) -> Result<FpsRecord, FpsError> {
    let (hex, id) = line
        .split_once('\t')
        .ok_or_else(|| FpsError::MalformedRecord(line.to_string()))?;
    let fingerprint = hex_to_bitvecn(hex, num_bits)?;
    Ok(FpsRecord {
        id: id.to_string(),
        fingerprint,
    })
}

fn hex_nibble(c: u8) -> Result<u8, FpsError> {
    match c {
        b'0'..=b'9' => Ok(c - b'0'),
        b'a'..=b'f' => Ok(c - b'a' + 10),
        b'A'..=b'F' => Ok(c - b'A' + 10),
        _ => Err(FpsError::InvalidHexDigit(c as char)),
    }
}

/// Decode a hex fingerprint string into a [`BitVecN`] of `num_bits` bits,
/// per the FPS convention documented at the top of this module.
fn hex_to_bitvecn(hex: &str, num_bits: usize) -> Result<BitVecN, FpsError> {
    let expected_hex_chars = num_bits.div_ceil(8) * 2;
    let hex_bytes = hex.as_bytes();
    if hex_bytes.len() != expected_hex_chars {
        return Err(FpsError::HexLengthMismatch {
            expected_hex_chars,
            got_hex_chars: hex_bytes.len(),
            num_bits,
        });
    }
    let mut bv = BitVecN::new(num_bits);
    for (byte_idx, pair) in hex_bytes.as_chunks::<2>().0.iter().enumerate() {
        let byte = (hex_nibble(pair[0])? << 4) | hex_nibble(pair[1])?;
        for bit_in_byte in 0..8usize {
            let bit_idx = byte_idx * 8 + bit_in_byte;
            if bit_idx < num_bits && (byte >> bit_in_byte) & 1 == 1 {
                bv.set(bit_idx);
            }
            // Set padding bits (bit_idx >= num_bits) are tolerated and
            // silently ignored rather than rejected — the spec says
            // writers must zero them, but a lenient reader shouldn't choke
            // on a nonconformant producer over bits that don't exist.
        }
    }
    Ok(bv)
}

/// Encode a [`BitVecN`] as an FPS hex fingerprint string.
fn bitvecn_to_hex(bv: &BitVecN) -> String {
    let bits = bv.bit_width();
    let num_bytes = bits.div_ceil(8);
    let mut hex = String::with_capacity(num_bytes * 2);
    for byte_idx in 0..num_bytes {
        let mut byte = 0u8;
        for bit_in_byte in 0..8usize {
            let bit_idx = byte_idx * 8 + bit_in_byte;
            if bit_idx < bits && bv.get(bit_idx) {
                byte |= 1 << bit_in_byte;
            }
        }
        use std::fmt::Write as _;
        write!(hex, "{byte:02x}").expect("String write is infallible");
    }
    hex
}

/// Streaming FPS writer over any [`Write`] sink.
///
/// Writes the header block immediately on construction, then one data line
/// per [`write_record`](Self::write_record) call — no buffering of the
/// whole output.
pub struct FpsWriter<W: Write> {
    inner: W,
}

impl<W: Write> FpsWriter<W> {
    /// Write the header block and return a writer ready for
    /// [`write_record`](Self::write_record) calls.
    pub fn new(mut inner: W, header: &FpsHeader) -> std::io::Result<Self> {
        header.write_to(&mut inner)?;
        Ok(Self { inner })
    }

    /// Write one `<hex>\t<id>` data line.
    ///
    /// # Errors
    /// Returns an [`std::io::ErrorKind::InvalidInput`] error if `id`
    /// contains a tab or newline (the FPS identifier field forbids both).
    pub fn write_record(&mut self, id: &str, fingerprint: &BitVecN) -> std::io::Result<()> {
        if id.contains(['\t', '\n', '\r']) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "FPS identifier must not contain a tab or newline",
            ));
        }
        let hex = bitvecn_to_hex(fingerprint);
        writeln!(self.inner, "{hex}\t{id}")
    }

    /// Convenience wrapper for the common case of writing a fixed 2048-bit
    /// fingerprint (ECFP4/ECFP6/MACCS/etc. in this crate all return
    /// [`BitVec2048`]).
    pub fn write_record_2048(&mut self, id: &str, fingerprint: &BitVec2048) -> std::io::Result<()> {
        self.write_record(id, &fingerprint.to_bitvecn())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;
    use std::io::Cursor;

    // ===== Header parsing =====

    #[test]
    fn header_recognizes_standard_fields_and_keeps_extras() {
        let text = "#FPS1\n#num_bits=8\n#type=Test/1\n#software=chematic-fp/0.15.0\n#source=in.smi\n#source=in2.smi\n#comment=hello world\n#custom_key=custom_value\n01\tfirst\n";
        let reader = FpsReader::new(Cursor::new(text)).unwrap();
        let header = reader.header();
        assert_eq!(header.num_bits, 8);
        assert_eq!(header.fp_type.as_deref(), Some("Test/1"));
        assert_eq!(header.software.as_deref(), Some("chematic-fp/0.15.0"));
        assert_eq!(header.source, vec!["in.smi", "in2.smi"]);
        assert_eq!(header.comments, vec!["hello world"]);
        assert_eq!(header.extra, vec!["FPS1", "custom_key=custom_value"]);
    }

    #[test]
    fn header_write_back_keeps_version_line_first_and_preserves_extras() {
        // Read a file with a #FPS1 version line plus an unrecognized custom
        // key, write the *parsed* header straight back out, and confirm:
        // (a) #FPS1 is still line 1 (the spec requires the version line to
        // be first -- it must not just be preserved *somewhere*), and
        // (b) the custom extra line survives too.
        let text = "#FPS1\n#num_bits=8\n#custom_key=custom_value\n01\tfirst\n";
        let reader = FpsReader::new(Cursor::new(text)).unwrap();
        let header = reader.header().clone();

        let mut out = Vec::new();
        header.write_to(&mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = out.lines().collect();

        assert_eq!(lines[0], "#FPS1", "version line must be the first line");
        assert!(
            lines.contains(&"#custom_key=custom_value"),
            "unrecognized extra header line must survive write-back: {lines:?}"
        );
        assert!(lines.contains(&"#num_bits=8"));

        // Re-parsing the written-back header must reproduce the same
        // extras (fully round-trip-safe, not just "line present").
        let reparsed = FpsHeader::default().tap_absorb_all(&lines);
        assert_eq!(reparsed.extra, header.extra);
    }

    // Small test-only helper: replay a set of raw "#..." lines through
    // `absorb_line` the same way `FpsReader::new` does, to check write_to's
    // output re-parses to the same header fields.
    impl FpsHeader {
        fn tap_absorb_all(mut self, lines: &[&str]) -> Self {
            for line in lines {
                if let Some(rest) = line.strip_prefix('#') {
                    self.absorb_line(rest).unwrap();
                }
            }
            self
        }
    }

    #[test]
    fn missing_num_bits_errors() {
        let text = "#type=Test\n01\tfirst\n";
        let err = FpsReader::new(Cursor::new(text)).unwrap_err();
        assert!(matches!(err, FpsError::MissingNumBits));
    }

    // ===== Spec bit-order positive controls =====
    // These three hex values and their bit patterns are taken directly from
    // the chemfp FPS format spec (https://chemfp.com/fps_format/): "the hex
    // value "01" encodes a byte value of 1, the hex value "0a" encodes the
    // byte value 10, and the hex value "c0" encodes the byte value 192."

    #[test]
    fn spec_example_bit_patterns_round_trip() {
        let text = "#num_bits=8\n01\tone\n0a\tten\nc0\tonehundredninetytwo\n";
        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();

        let r1 = reader.next().unwrap().unwrap();
        assert_eq!(r1.id, "one");
        assert_bits(&r1.fingerprint, &[0]);

        let r2 = reader.next().unwrap().unwrap();
        assert_eq!(r2.id, "ten");
        assert_bits(&r2.fingerprint, &[1, 3]);

        let r3 = reader.next().unwrap().unwrap();
        assert_eq!(r3.id, "onehundredninetytwo");
        assert_bits(&r3.fingerprint, &[6, 7]);

        assert!(reader.next().is_none());
    }

    fn assert_bits(bv: &BitVecN, expected_set: &[usize]) {
        for i in 0..bv.bit_width() {
            let want = expected_set.contains(&i);
            assert_eq!(bv.get(i), want, "bit {i} mismatch");
        }
    }

    #[test]
    fn write_matches_spec_example_hex() {
        let mut bv = BitVecN::new(8);
        bv.set(6);
        bv.set(7);
        assert_eq!(bitvecn_to_hex(&bv), "c0");

        let mut bv = BitVecN::new(8);
        bv.set(1);
        bv.set(3);
        assert_eq!(bitvecn_to_hex(&bv), "0a");
    }

    // ===== Non-byte-aligned width (166-bit MACCS-shaped external file) =====

    #[test]
    fn reads_non_byte_aligned_166_bit_external_fixture() {
        // 166 bits -> ceil(166/8) = 21 bytes -> 42 hex chars, last byte
        // padded with 2 zero bits at the high end (bits 166, 167 unused).
        // Set bit 0 (first byte, LSB) and bit 165 (last valid bit, i.e. bit
        // 5 of byte 20) by hand.
        let mut bytes = [0u8; 21];
        bytes[0] = 0b0000_0001; // bit 0
        bytes[20] = 0b0010_0000; // byte 20 covers bits 160..168; bit 165 = bit_in_byte 5
        let hex: String = bytes.iter().map(|b| format!("{b:02x}")).collect();
        let text = format!("#num_bits=166\n#type=MACCS-166\n{hex}\texternal_maccs\n");

        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();
        assert_eq!(reader.header().num_bits, 166);
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.id, "external_maccs");
        assert_eq!(record.fingerprint.bit_width(), 166);
        assert!(record.fingerprint.get(0));
        assert!(record.fingerprint.get(165));
        assert_eq!(record.fingerprint.popcount(), 2);
    }

    #[test]
    fn hex_length_mismatch_is_rejected() {
        let text = "#num_bits=16\nff\tshort\n"; // 16 bits needs 4 hex chars, only 2 given
        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();
        let err = reader.next().unwrap().unwrap_err();
        assert!(matches!(err, FpsError::HexLengthMismatch { .. }));
    }

    #[test]
    fn invalid_hex_digit_is_rejected() {
        let text = "#num_bits=8\nzz\tbad\n";
        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();
        let err = reader.next().unwrap().unwrap_err();
        assert!(matches!(err, FpsError::InvalidHexDigit('z')));
    }

    #[test]
    fn malformed_record_missing_tab_is_rejected() {
        let text = "#num_bits=8\nc0nolabel\n";
        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();
        let err = reader.next().unwrap().unwrap_err();
        assert!(matches!(err, FpsError::MalformedRecord(_)));
    }

    #[test]
    fn write_record_rejects_tab_in_identifier() {
        let header = FpsHeader::new(8);
        let mut out = Vec::new();
        let mut writer = FpsWriter::new(&mut out, &header).unwrap();
        let bv = BitVecN::new(8);
        let err = writer.write_record("bad\tid", &bv).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    }

    #[test]
    fn trailing_blank_line_is_tolerated() {
        let text = "#num_bits=8\n01\tfirst\n\n";
        let mut reader = FpsReader::new(Cursor::new(text)).unwrap();
        assert!(reader.next().unwrap().is_ok());
        assert!(reader.next().is_none());
    }

    // ===== Round trips against this crate's own fingerprints =====

    fn benzene() -> chematic_core::Molecule {
        parse("c1ccccc1").unwrap()
    }

    #[test]
    fn ecfp4_round_trip_read_write_read_byte_identical_hex() {
        let fp = crate::ecfp4(&benzene());
        let header = FpsHeader::for_chematic(2048, "ECFP4");

        let mut buf1 = Vec::new();
        {
            let mut w = FpsWriter::new(&mut buf1, &header).unwrap();
            w.write_record_2048("benzene", &fp).unwrap();
        }

        // Read pass 1.
        let mut r1 = FpsReader::new(Cursor::new(&buf1)).unwrap();
        assert_eq!(r1.header().num_bits, 2048);
        assert_eq!(r1.header().fp_type.as_deref(), Some("ECFP4"));
        let rec1 = r1.next().unwrap().unwrap();
        assert_eq!(rec1.fingerprint, fp.to_bitvecn());

        // Write pass 2 from the just-parsed header + record (not the
        // original in-memory `header` -- this exercises the header's own
        // write-back path, not just the fingerprint's).
        let mut buf2 = Vec::new();
        {
            let mut w = FpsWriter::new(&mut buf2, r1.header()).unwrap();
            w.write_record(&rec1.id, &rec1.fingerprint).unwrap();
        }

        // The hex fingerprint field must be byte-identical across both
        // write passes.
        let hex1 = std::str::from_utf8(&buf1)
            .unwrap()
            .lines()
            .last()
            .unwrap()
            .split_once('\t')
            .unwrap()
            .0;
        let hex2 = std::str::from_utf8(&buf2)
            .unwrap()
            .lines()
            .last()
            .unwrap()
            .split_once('\t')
            .unwrap()
            .0;
        assert_eq!(hex1, hex2);

        // Read pass 2 must reproduce the exact same BitVecN.
        let mut r2 = FpsReader::new(Cursor::new(&buf2)).unwrap();
        let rec2 = r2.next().unwrap().unwrap();
        assert_eq!(rec2.fingerprint, fp.to_bitvecn());
    }

    #[test]
    fn maccs_round_trip_read_write_read_byte_identical_hex() {
        // chematic's maccs() returns a full BitVec2048 (2048 bits, only
        // bits 0..165 ever set), not the conventional 166-bit-only
        // representation some other toolkits use -- verified by reading
        // maccs.rs (`pub fn maccs(mol: &Molecule) -> BitVec2048`). The FPS
        // header's num_bits must match what's actually serialized (2048),
        // and the type label spells out the mismatch explicitly
        // ("MACCS166/2048") rather than writing `type=MACCS166` next to
        // `num_bits=2048`, which would read as self-contradictory to an
        // external FPS consumer expecting a real 166-bit vector.
        let fp = crate::maccs(&benzene());
        let header = FpsHeader::for_chematic(2048, "MACCS166/2048");

        let mut buf1 = Vec::new();
        {
            let mut w = FpsWriter::new(&mut buf1, &header).unwrap();
            w.write_record_2048("benzene", &fp).unwrap();
        }
        let mut r1 = FpsReader::new(Cursor::new(&buf1)).unwrap();
        assert_eq!(r1.header().num_bits, 2048);
        let rec1 = r1.next().unwrap().unwrap();
        assert_eq!(rec1.fingerprint, fp.to_bitvecn());

        let mut buf2 = Vec::new();
        {
            let mut w = FpsWriter::new(&mut buf2, r1.header()).unwrap();
            w.write_record(&rec1.id, &rec1.fingerprint).unwrap();
        }
        assert_eq!(
            buf1, buf2,
            "MACCS FPS output must be byte-identical across read/write/read"
        );

        let mut r2 = FpsReader::new(Cursor::new(&buf2)).unwrap();
        let rec2 = r2.next().unwrap().unwrap();
        assert_eq!(rec2.fingerprint, fp.to_bitvecn());
    }

    #[test]
    fn multi_record_stream_iterates_in_order() {
        let header = FpsHeader::for_chematic(2048, "ECFP4");
        let fps = [
            crate::ecfp4(&benzene()),
            crate::ecfp4(&parse("CC").unwrap()),
            crate::ecfp4(&parse("C").unwrap()),
        ];
        let mut buf = Vec::new();
        {
            let mut w = FpsWriter::new(&mut buf, &header).unwrap();
            for (i, fp) in fps.iter().enumerate() {
                w.write_record_2048(&format!("mol{i}"), fp).unwrap();
            }
        }
        let reader = FpsReader::new(Cursor::new(&buf)).unwrap();
        let records: Vec<_> = reader.collect::<Result<_, _>>().unwrap();
        assert_eq!(records.len(), 3);
        for (i, (record, fp)) in records.iter().zip(fps.iter()).enumerate() {
            assert_eq!(record.id, format!("mol{i}"));
            assert_eq!(record.fingerprint, fp.to_bitvecn());
        }
    }

    #[test]
    fn for_chematic_stamps_software_and_type() {
        let header = FpsHeader::for_chematic(2048, "ECFP4");
        assert_eq!(header.fp_type.as_deref(), Some("ECFP4"));
        assert_eq!(
            header.software.as_deref(),
            Some(concat!("chematic-fp/", env!("CARGO_PKG_VERSION")))
        );
    }

    #[test]
    fn builder_methods_set_arbitrary_caller_supplied_fields() {
        // General case: caller supplies its own type/software (e.g.
        // round-tripping a file produced by a different tool).
        let header = FpsHeader::new(1024)
            .with_type("RDKit-Fingerprint/2")
            .with_software("RDKit/2025.09.4")
            .with_source("compounds.smi")
            .with_comment("example");
        assert_eq!(header.num_bits, 1024);
        assert_eq!(header.fp_type.as_deref(), Some("RDKit-Fingerprint/2"));
        assert_eq!(header.software.as_deref(), Some("RDKit/2025.09.4"));
        assert_eq!(header.source, vec!["compounds.smi"]);
        assert_eq!(header.comments, vec!["example"]);
    }
}
