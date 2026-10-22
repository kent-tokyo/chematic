"""Typed names for the native MCS modes; legacy string values remain accepted."""
from enum import Enum
from typing import Literal, Union


class McsAtomCompare(str, Enum):
    """Native MCS atom comparison modes.

    ELEMENTS requires the same atomic number. ANY_HEAVY_ATOM ignores element
    identity but excludes hydrogen. ANY allows all atoms, including hydrogen.
    RDKIT_ANY_HEAVY_ATOM also permits explicit H-to-H matches, following
    RDKit's CompareAnyHeavyAtom. These are native comparator choices, not
    RDKit comparator subclass adapters.
    """

    ELEMENTS = "elements"
    ANY_HEAVY_ATOM = "any_heavy_atom"
    RDKIT_ANY_HEAVY_ATOM = "rdkit_any_heavy_atom"
    ANY = "any"


class McsBondCompare(str, Enum):
    """Native MCS bond comparison modes.

    ORDER_OR_AROMATIC requires equal normalized bond order (Up/Down become
    Single). ANY ignores bond order while retaining required connectivity.
    """

    ORDER_OR_AROMATIC = "order_or_aromatic"
    ANY = "any"


McsAtomComparison = Union[
    McsAtomCompare,
    Literal["elements", "any_heavy_atom", "rdkit_any_heavy_atom", "any"],
]
McsBondComparison = Union[McsBondCompare, Literal["order_or_aromatic", "any"]]
