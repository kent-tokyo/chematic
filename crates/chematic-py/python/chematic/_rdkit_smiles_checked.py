"""An opt-in fixed-point check for the legacy RDKit-compatible writer."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Optional

if TYPE_CHECKING:
    from . import Mol, from_smiles
else:
    from .chematic import Mol, from_smiles


class RdkitSmilesStabilityError(ValueError):
    """A candidate could not be read and rewritten to the identical string.

    The candidate strings are diagnostic evidence, not accepted identifiers.
    ``code`` is ``unreadable_output``, ``rewrite_failed`` or ``non_idempotent``.
    """

    def __init__(
        self,
        code: Literal["unreadable_output", "rewrite_failed", "non_idempotent"],
        first_smiles: str,
        roundtrip_smiles: Optional[str],
    ) -> None:
        self.code = code
        self.first_smiles = first_smiles
        self.roundtrip_smiles = roundtrip_smiles
        super().__init__(
            f"RDKit-compatible canonical SMILES failed stability check: {code}"
        )


def rdkit_smiles_checked(mol: Mol, /) -> str:
    """Return default RDKit-compatible SMILES only if one rewrite is identical.

    Write ``mol.rdkit_smiles``, parse that string, then write it again. Reject
    unequal strings with :class:`RdkitSmilesStabilityError`; never select a
    member of a cycle, remove stereo tags or retry until convergence. The
    input molecule is not modified. Initial writer failures propagate.

    This checks textual idempotence only. It does not prove preservation of
    the input's chemical identity or completeness of its stereochemistry.
    RDKit issue #9391 is one known legacy writer failure caught by this check.
    """
    if not isinstance(mol, Mol):
        raise TypeError("mol must be a chematic.Mol")
    first = mol.rdkit_smiles
    try:
        readback = from_smiles(first)
    except (ValueError, RuntimeError) as error:
        raise RdkitSmilesStabilityError("unreadable_output", first, None) from error
    try:
        rewritten = readback.rdkit_smiles
    except (ValueError, RuntimeError) as error:
        raise RdkitSmilesStabilityError("rewrite_failed", first, None) from error
    if first != rewritten:
        raise RdkitSmilesStabilityError("non_idempotent", first, rewritten)
    return first
