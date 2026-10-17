"""Bounded, source-preserving file-backed SDF partitions (Open Babel #2147)."""

from __future__ import annotations

from array import array
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import operator
import os
import re
import stat
from typing import TYPE_CHECKING, BinaryIO, Iterator, SupportsIndex, Tuple, Union, cast

if TYPE_CHECKING:
    _Offsets = array[int]
else:
    _Offsets = array

SdfSplitPath = Union[str, bytes, os.PathLike[str], os.PathLike[bytes]]
_Identity = Tuple[int, int, int, int, int]
_MAX_BYTES = 1 << 30
_MAX_RECORDS = 1_000_000
_MAX_PARTS = 10_000
_MAX_LINE = 1 << 20
_MAX_RECORD = 64 << 20
_MAX_CHUNK = 1 << 20


def _integer(value: SupportsIndex, name: str, lower: int, upper: int) -> int:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer, excluding bool")
    n = operator.index(value)
    if not lower <= n <= upper:
        raise ValueError(f"{name} must be in {lower}..{upper}")
    return n


def _identity(file: BinaryIO) -> _Identity:
    s = os.fstat(file.fileno())
    if not stat.S_ISREG(s.st_mode):
        raise ValueError("SDF split source must be a regular file")
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


@contextmanager
def _open(path: Union[str, bytes]) -> Iterator[BinaryIO]:
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NONBLOCK", 0)
    fd = os.open(path, flags)
    try:
        file = cast(BinaryIO, os.fdopen(fd, "rb"))
    except BaseException:
        os.close(fd)
        raise
    with file:
        _identity(file)
        yield file


@dataclass(frozen=True)
class SdfSplitPart:
    """A contiguous byte interval and zero-based source record interval."""

    filename: str
    byte_offset: int
    byte_len: int
    first_record: int
    record_count: int
    sha256: str


@dataclass(frozen=True)
class SdfSplitPlan:
    """Immutable split metadata; iterators read/verify bounded chunks on demand."""

    source_path: Union[str, bytes]
    source_bytes: int
    source_sha256: str
    record_count: int
    unterminated_final_record: bool
    parts: Tuple[SdfSplitPart, ...]
    _source_identity: _Identity

    def iter_part(
        self, index: SupportsIndex, *, chunk_bytes: SupportsIndex = _MAX_CHUNK
    ) -> Iterator[bytes]:
        """Read a part in at most 1 MiB chunks, verifying its hash at exhaustion.

        An index/buffer error is immediate. Opening/reading is lazy. Consume to
        exhaustion to verify the complete part; a read/mutation error can occur
        after chunks were yielded. The API writes no output files.
        """
        if isinstance(index, bool):
            raise TypeError("part index must be an integer, excluding bool")
        i = operator.index(index)
        if not 0 <= i < len(self.parts):
            raise IndexError("SDF split part index outside plan")
        size = _integer(chunk_bytes, "chunk_bytes", 1, _MAX_CHUNK)
        part = self.parts[i]
        if (
            part.byte_len <= 0
            or part.byte_offset < 0
            or part.byte_offset + part.byte_len > self._source_identity[2]
        ):
            raise ValueError("invalid byte range in SDF split plan")
        return _chunks(self.source_path, self._source_identity, part, size)


def _chunks(
    path: Union[str, bytes], identity: _Identity, part: SdfSplitPart, size: int
) -> Iterator[bytes]:
    with _open(path) as file:
        if _identity(file) != identity:
            raise ValueError("SDF split source changed since planning")
        file.seek(part.byte_offset)
        remaining = part.byte_len
        digest = hashlib.sha256()
        while remaining:
            data = file.read(min(size, remaining))
            if not data:
                raise ValueError("SDF split source truncated while reading")
            remaining -= len(data)
            digest.update(data)
            yield data
        if digest.hexdigest() != part.sha256 or _identity(file) != identity:
            raise ValueError(
                "SDF split source changed while reading; part hash not verified"
            )


def _scan(
    file: BinaryIO, max_bytes: int, max_records: int, max_line: int, max_record: int
) -> tuple[_Offsets, str, bool]:
    if file.read(2) == b"\x1f\x8b":
        raise ValueError("compressed SDF is outside this plain-file split profile")
    file.seek(0)
    ends = array("Q")
    digest = hashlib.sha256()
    offset = start = 0
    content = False
    while True:
        line = file.readline(max_line + 1)
        if not line:
            break
        offset += len(line)
        if len(line) > max_line or offset > max_bytes or offset - start > max_record:
            raise ValueError("SDF split input/line/record byte limit exceeded")
        digest.update(line)
        body = line.removesuffix(b"\n").removesuffix(b"\r")
        if body == b"$$$$":
            if not content:
                raise ValueError("empty SDF record is unsupported")
            if len(ends) >= max_records:
                raise ValueError("SDF split record count limit exceeded")
            ends.append(offset)
            start = offset
            content = False
        else:
            content = content or bool(body.strip(b" \t\r\n\v\f"))
    unterminated = content
    if content:
        if len(ends) >= max_records:
            raise ValueError("SDF split record count limit exceeded")
        ends.append(offset)
    elif ends:
        # Trailing whitespace belongs to the last part, not a phantom record.
        last_start = ends[-2] if len(ends) > 1 else 0
        if offset - last_start > max_record:
            raise ValueError(
                "SDF split record byte limit exceeded including trailing whitespace"
            )
        ends[-1] = offset
    if not ends:
        raise ValueError("SDF split requires a nonempty record")
    return ends, digest.hexdigest(), unterminated


def _parts(
    file: BinaryIO, ends: _Offsets, count: int, prefix: str
) -> tuple[Tuple[SdfSplitPart, ...], str]:
    quotient, remainder = divmod(len(ends), count)
    width = len(str(count))
    first = start = 0
    result = []
    whole = hashlib.sha256()
    for i in range(count):
        records = quotient + int(i < remainder)
        stop = ends[first + records - 1]
        digest = hashlib.sha256()
        file.seek(start)
        remaining = stop - start
        while remaining:
            data = file.read(min(_MAX_CHUNK, remaining))
            if not data:
                raise ValueError("SDF split source truncated during planning")
            digest.update(data)
            whole.update(data)
            remaining -= len(data)
        result.append(
            SdfSplitPart(
                f"{prefix}{i+1:0{width}d}.sdf",
                start,
                stop - start,
                first,
                records,
                digest.hexdigest(),
            )
        )
        first += records
        start = stop
    return tuple(result), whole.hexdigest()


def plan_sdf_split(
    path: SdfSplitPath,
    parts: SupportsIndex,
    *,
    prefix: str = "split",
    max_input_bytes: SupportsIndex = _MAX_BYTES,
    max_records: SupportsIndex = _MAX_RECORDS,
    max_line_bytes: SupportsIndex = _MAX_LINE,
    max_record_bytes: SupportsIndex = _MAX_RECORD,
) -> SdfSplitPlan:
    """Plan exactly `parts` balanced nonempty groups, preserving source bytes.

    This is raw SDF record framing, not molecular validation: invalid MOL blocks
    and unknown/duplicate SD properties are retained. A delimiter is an exact
    $$$$ line with LF, CRLF or EOF. Empty records are refused; an unterminated
    nonempty last record is retained. Larger groups occur first and differ in
    record count by at most one. Planning keeps 8-byte offsets per record and
    hashes each part in a second pass; the input file is never loaded in full.
    """
    count = _integer(parts, "parts", 1, _MAX_PARTS)
    byte_limit = _integer(max_input_bytes, "max_input_bytes", 1, _MAX_BYTES)
    record_limit = _integer(max_records, "max_records", 1, _MAX_RECORDS)
    line_limit = _integer(max_line_bytes, "max_line_bytes", 1, _MAX_LINE)
    record_bytes = _integer(max_record_bytes, "max_record_bytes", 1, _MAX_RECORD)
    if (
        not isinstance(prefix, str)
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}", prefix) is None
    ):
        raise ValueError(
            "prefix must be a 1..64 character ASCII file stem without paths"
        )
    source = os.path.abspath(os.fspath(path))
    with _open(source) as file:
        identity = _identity(file)
        if identity[2] > byte_limit:
            raise ValueError("SDF split input byte limit exceeded")
        ends, digest, unterminated = _scan(
            file, byte_limit, record_limit, line_limit, record_bytes
        )
        if count > len(ends):
            raise ValueError("parts must not exceed the source record count")
        planned, rechecked_digest = _parts(file, ends, count, prefix)
        if rechecked_digest != digest or _identity(file) != identity:
            raise ValueError("SDF split source changed during planning")
    return SdfSplitPlan(
        source, identity[2], digest, len(ends), unterminated, planned, identity
    )
