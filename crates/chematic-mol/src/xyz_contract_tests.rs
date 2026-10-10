use crate::xyz::*;
use std::io::{self, BufRead, Cursor, Read};

#[test]
fn malformed_extxyz_fields_report_the_rejected_field() {
    for (header, atom, diagnostic) in [
        ("=bad", "C 0 0 0", "info field"),
        ("name=\"unfinished", "C 0 0 0", "unterminated"),
        ("name=a name=b", "C 0 0 0", "duplicate"),
        ("Lattice=\"1 2\"", "C 0 0 0", "Lattice"),
        ("Lattice=\"1 0 0 0 x 0 0 0 1\"", "C 0 0 0", "Lattice"),
        ("Lattice=\"1 0 0 0 NaN 0 0 0 1\"", "C 0 0 0", "Lattice"),
        ("Properties=species:S", "C 0 0 0", "Properties"),
        ("Properties=species:Q:1:pos:R:3", "C 0 0 0", "Properties"),
        ("Properties=species:S:x:pos:R:3", "C 0 0 0", "Properties"),
        ("Properties=species:S:0:pos:R:3", "C 0 0 0", "Properties"),
        ("Properties=pos:R:3:species:S:1", "C 0 0 0", "Properties"),
        ("Properties=species:S:1:pos:R:3", "C 0 0", "columns"),
        ("", "Qq 0 0 0", "element"),
        ("", "C invalid 0 0", "coordinate"),
        ("", "C NaN 0 0", "non-finite"),
        (
            "Properties=species:S:1:pos:R:3:force:R:1",
            "C 0 0 0 bad",
            "force",
        ),
        (
            "Properties=species:S:1:pos:R:3:force:R:1",
            "C 0 0 0 inf",
            "force",
        ),
        (
            "Properties=species:S:1:pos:R:3:tag:I:1",
            "C 0 0 0 1.5",
            "tag",
        ),
        (
            "Properties=species:S:1:pos:R:3:fixed:L:1",
            "C 0 0 0 maybe",
            "fixed",
        ),
    ] {
        let input = format!("1\n{header}\n{atom}\n");
        let error = parse_extxyz(&input).unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{input}: {error}");
        let mut stream = ExtxyzFileReader::new(Cursor::new(input));
        assert_eq!(
            stream.next().unwrap().unwrap_err().to_string(),
            error.to_string()
        );
        assert!(stream.next().is_none(), "must stop after a bad frame");
    }
}

#[test]
fn xyz_streams_fail_closed_on_truncation_and_limits() {
    let valid = "1\nframe\nC 0 0 0\n";
    let defaults = XyzParseLimits::default();
    for (input, limits, diagnostic) in [
        ("bad\n", defaults, "atom-count"),
        ("1\n", defaults, "found 0"),
        ("2\nframe\nC 0 0 0\n", defaults, "found 1"),
        ("1\nframe\nQq 0 0 0\n", defaults, "element"),
        (
            valid,
            XyzParseLimits {
                max_atoms_per_frame: 0,
                ..defaults
            },
            "atoms per frame",
        ),
        (
            valid,
            XyzParseLimits {
                max_input_bytes: 0,
                ..defaults
            },
            "input bytes",
        ),
        (
            valid,
            XyzParseLimits {
                max_input_bytes: 3,
                ..defaults
            },
            "input bytes",
        ),
        (
            valid,
            XyzParseLimits {
                max_line_bytes: 2,
                ..defaults
            },
            "line bytes",
        ),
        (
            valid,
            XyzParseLimits {
                max_frames: 0,
                ..defaults
            },
            "frames",
        ),
    ] {
        let mut plain = XyzFileReader::with_limits(Cursor::new(input), limits);
        let error = plain.next().unwrap().unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{error}");
        assert!(plain.next().is_none());
        let mut ext = ExtxyzFileReader::with_limits(Cursor::new(input), limits);
        let error = ext.next().unwrap().unwrap_err();
        assert!(error.to_string().contains(diagnostic), "{error}");
        assert!(ext.next().is_none());
    }
    for input in ["", "\n\n", "\r\n"] {
        assert!(XyzFileReader::new(Cursor::new(input)).next().is_none());
        assert!(ExtxyzFileReader::new(Cursor::new(input)).next().is_none());
    }
}

/// Fail at any frame boundary, including in the count, comment and atom lines.
struct FailingReader {
    good: Cursor<Vec<u8>>,
}
impl Read for FailingReader {
    fn read(&mut self, _out: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::other("injected read failure"))
    }
}
impl BufRead for FailingReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.good.position() as usize == self.good.get_ref().len() {
            Err(io::Error::other("injected read failure"))
        } else {
            self.good.fill_buf()
        }
    }
    fn consume(&mut self, amount: usize) {
        self.good.consume(amount);
    }
}

#[test]
fn io_failures_are_distinguished_from_end_of_stream() {
    for prefix in ["", "1\n", "1\nframe\n"] {
        let reader = || FailingReader {
            good: Cursor::new(prefix.as_bytes().to_vec()),
        };
        let mut plain = XyzFileReader::new(reader());
        assert!(
            plain
                .next()
                .unwrap()
                .unwrap_err()
                .to_string()
                .contains("I/O error")
        );
        assert!(plain.next().is_none());
        let mut ext = ExtxyzFileReader::new(reader());
        assert!(
            ext.next()
                .unwrap()
                .unwrap_err()
                .to_string()
                .contains("I/O error")
        );
        assert!(ext.next().is_none());
    }
}
