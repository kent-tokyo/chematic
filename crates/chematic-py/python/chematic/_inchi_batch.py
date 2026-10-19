"""Position-preserving, identity-checked Standard InChIKey batches."""

from __future__ import annotations

from dataclasses import dataclass
import operator
import re
from typing import TYPE_CHECKING, Iterable, Literal, Optional, SupportsIndex

import numpy as np

if TYPE_CHECKING:
    from . import ChematicInputError, from_inchi
else:
    from .chematic import ChematicInputError, from_inchi

_MAX_BYTES = 16 << 20
_MAX_RECORD_BYTES = 8192
_MAX_FORMULA_ATOMS = 4096
_KEY = re.compile(r"[A-Z]{14}-[A-Z]{10}-[A-Z]\Z")
_FORMULA = re.compile(r"(?:[A-Z][a-z]?[0-9]*)+\Z")
_ELEMENT_COUNT = re.compile(r"[A-Z][a-z]?([0-9]*)")
InchiKeyStatus = Literal[
    "ok", "invalid", "unsupported", "resource_limit", "conversion_error"
]


@dataclass(frozen=True)
class InchiKeyRecord:
    """One input position, retaining failures instead of dropping rows.

    key is populated only for status='ok'; code/message describe other states.
    source_text is exact. A failed row does not establish that its input is
    chemically invalid: 'unsupported' and 'conversion_error' are distinct.
    """

    index: int
    source_text: str
    status: InchiKeyStatus
    key: Optional[str]
    code: Optional[str]
    message: Optional[str]


def _stage(values: Iterable[str], limit: int) -> list[tuple[str, Optional[int]]]:
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError("values must be an iterable of individual InChI strings")
    result: list[tuple[str, Optional[int]]] = []
    total = 0
    for value in values:
        if len(result) >= limit:
            raise ValueError("InChI batch exceeds max_records")
        if type(value) is not str:
            raise TypeError("each InChI input must be a plain str")
        if len(value) > _MAX_BYTES:
            raise ValueError("InChI batch exceeds 16 MiB")
        try:
            size: Optional[int] = len(value.encode("utf-8"))
        except UnicodeError:
            size = None
        total += (
            size if size is not None else len(value.encode("utf-8", "surrogatepass"))
        )
        if total > _MAX_BYTES:
            raise ValueError("InChI batch exceeds 16 MiB")
        result.append((value, size))
    return result


def _failure(
    index: int, text: str, status: InchiKeyStatus, code: str, message: str
) -> InchiKeyRecord:
    return InchiKeyRecord(index, text, status, None, code, message[:1024])


def _preflight(index: int, text: str, size: Optional[int]) -> Optional[InchiKeyRecord]:
    if size is None:
        return _failure(
            index, text, "invalid", "invalid_utf8", "Input is not valid UTF-8"
        )
    if size > _MAX_RECORD_BYTES:
        return _failure(
            index,
            text,
            "resource_limit",
            "input_byte_limit",
            "InChI exceeds 8192 bytes",
        )
    if not text.startswith("InChI="):
        return _failure(
            index, text, "invalid", "inchi_prefix", "Input must start with InChI="
        )
    if not text.startswith("InChI=1S/"):
        return _failure(
            index,
            text,
            "unsupported",
            "standard_inchi_required",
            "A Standard InChI=1S/ identifier is required",
        )
    if not re.fullmatch(r"[\x21-\x7e]+", text):
        return _failure(
            index,
            text,
            "invalid",
            "invalid_characters",
            "Whitespace/control/non-ASCII characters are not supported",
        )
    # Bound numeric ranges before either the native or fallback parser runs.
    for token in re.findall(r"[0-9]+", text):
        significant = token.lstrip("0")
        if len(significant) > 4 or int(significant or "0") > _MAX_FORMULA_ATOMS:
            return _failure(
                index,
                text,
                "resource_limit",
                "numeric_limit",
                "Numeric fields must not exceed 4096",
            )
    formula = text.split("/", 2)[1]
    atoms = 0
    for component in formula.split("."):
        match = re.fullmatch(r"([0-9]*)(.*)", component)
        assert match is not None
        multiplier, body = match.groups()
        if not _FORMULA.fullmatch(body):
            return _failure(
                index,
                text,
                "invalid",
                "invalid_formula",
                "Unsupported or malformed formula layer",
            )
        counts = [
            int(value.lstrip("0") or "0") if value else 1
            for value in _ELEMENT_COUNT.findall(body)
        ]
        copies = int(multiplier.lstrip("0") or "0") if multiplier else 1
        if copies == 0 or 0 in counts:
            return _failure(
                index,
                text,
                "invalid",
                "invalid_formula",
                "Formula counts and multipliers must be positive",
            )
        atoms += copies * sum(counts)
        if atoms > _MAX_FORMULA_ATOMS:
            return _failure(
                index,
                text,
                "resource_limit",
                "formula_atom_limit",
                "Formula exceeds 4096 atoms including H and component multipliers",
            )
    return None


def inchi_keys_checked(
    values: Iterable[str], *, max_records: SupportsIndex = 1000
) -> tuple[InchiKeyRecord, ...]:
    """Return one indexed success/failure record per Standard InChI input.

    Args:
        values: Plain str identifiers, in source order. Generators and duplicate
            values are supported. No whitespace stripping or normalization.
        max_records: Integer in 1..10000, excluding bool. Input text totals at
            most 16 MiB. All source items are staged before native conversion.

    Returns:
        Frozen records, including rejected rows. A standard IUPAC key is returned
        only when native regeneration exactly matches the original identifier.
        No approximate Mol.inchikey or unverified readback key is used. The native
        InChI feature is needed; its RuntimeError becomes conversion_error.

    Raises:
        TypeError: Invalid batch/row/limit type.
        ValueError: Invalid batch limit or aggregate size. Exceeding max_records
            consumes one extra item. No partial tuple is returned.

    Individual oversize/numeric/formula errors are resource_limit rows. Other
    parser categories remain distinct; changed identity is unsupported, not
    proof of invalid input. MemoryError and unexpected exceptions propagate.
    Iterator exceptions propagate and native work has not started at that point.
    Caller iterator side effects are not rolled back. No files are written.
    """
    if isinstance(max_records, (bool, np.bool_)):
        raise TypeError("max_records must be an integer, excluding bool")
    limit = operator.index(max_records)
    if not 1 <= limit <= 10000:
        raise ValueError("max_records must be in 1..10000")
    staged = _stage(values, limit)
    output = []
    for index, (text, size) in enumerate(staged):
        rejected = _preflight(index, text, size)
        if rejected is not None:
            output.append(rejected)
            continue
        try:
            molecule = from_inchi(text)
            if molecule.standard_inchi != text:
                output.append(
                    _failure(
                        index,
                        text,
                        "unsupported",
                        "identity_changed",
                        "Regenerated Standard InChI differs from input; no key returned",
                    )
                )
                continue
            key = molecule.standard_inchikey
        except ChematicInputError as error:
            category = getattr(error, "category", "malformed")
            status: InchiKeyStatus = "invalid"
            if category == "resource_limit":
                status = "resource_limit"
            elif category in ("unsupported", "ambiguous"):
                status = "unsupported"
            code = getattr(error, "code", "inchi_parse")
            output.append(
                _failure(
                    index,
                    text,
                    status,
                    code if isinstance(code, str) else "inchi_parse",
                    str(error),
                )
            )
            continue
        except RuntimeError as error:
            output.append(
                _failure(
                    index, text, "conversion_error", "native_inchi_error", str(error)
                )
            )
            continue
        if not _KEY.fullmatch(key):
            output.append(
                _failure(
                    index,
                    text,
                    "conversion_error",
                    "invalid_key",
                    "Native output is not an InChIKey",
                )
            )
        else:
            output.append(InchiKeyRecord(index, text, "ok", key, None, None))
    return tuple(output)
