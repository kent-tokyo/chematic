//! Per-input compression selection for native streaming SDF readers.

use flate2::read::MultiGzDecoder;
use std::io::{self, BufRead, BufReader, Chain, Cursor, Read};

/// Compression policy for one input, independent of all previous readers.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SdfCompression {
    /// Read bytes unchanged (the compatibility default).
    #[default]
    Plain,
    /// Require gzip; corrupt headers, payloads and trailers are IO errors.
    Gzip,
    /// Inspect the first two bytes for gzip magic, without using the filename.
    Auto,
}

type Replay<R> = Chain<Cursor<Vec<u8>>, R>;

enum Input<R: Read> {
    Plain(BufReader<Replay<R>>),
    Gzip(BufReader<MultiGzDecoder<Replay<R>>>),
}

/// Buffered plain/gzip input suitable for [`super::SdfFileReader`].
///
/// Auto detection reads at most two bytes on construction, replaying them into
/// the selected reader. Short/non-seekable reads are supported. Concatenated
/// gzip members are decoded as one byte stream, including a record split
/// across members. A fresh reader selects its own compression mode.
///
/// Integrity failures may be detected after records have already been yielded;
/// consume the SDF iterator through EOF to check the entire compressed input.
pub struct SdfInput<R: Read> {
    inner: Input<R>,
}

impl<R: Read> SdfInput<R> {
    /// Construct one independent input. Header-probe IO failures return here;
    /// decompression failures are returned by subsequent reads.
    pub fn new(mut source: R, compression: SdfCompression) -> io::Result<Self> {
        let mut prefix = Vec::new();
        let selected = if compression == SdfCompression::Auto {
            (&mut source).take(2).read_to_end(&mut prefix)?;
            if prefix == [0x1f, 0x8b] {
                SdfCompression::Gzip
            } else {
                SdfCompression::Plain
            }
        } else {
            compression
        };
        let replay = Cursor::new(prefix).chain(source);
        let inner = match selected {
            SdfCompression::Gzip => Input::Gzip(BufReader::new(MultiGzDecoder::new(replay))),
            _ => Input::Plain(BufReader::new(replay)),
        };
        Ok(Self { inner })
    }

    /// Resolved mode (`Plain` or `Gzip`), never `Auto`.
    pub fn compression(&self) -> SdfCompression {
        match self.inner {
            Input::Plain(_) => SdfCompression::Plain,
            Input::Gzip(_) => SdfCompression::Gzip,
        }
    }
}

impl<R: Read> Read for SdfInput<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        match &mut self.inner {
            Input::Plain(input) => input.read(output),
            Input::Gzip(input) => input.read(output),
        }
    }
}

impl<R: Read> BufRead for SdfInput<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match &mut self.inner {
            Input::Plain(input) => input.fill_buf(),
            Input::Gzip(input) => input.fill_buf(),
        }
    }

    fn consume(&mut self, amount: usize) {
        match &mut self.inner {
            Input::Plain(input) => input.consume(amount),
            Input::Gzip(input) => input.consume(amount),
        }
    }
}
