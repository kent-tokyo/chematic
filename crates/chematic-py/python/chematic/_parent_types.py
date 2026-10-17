"""Public dictionary contracts for native Parent transformation reports."""

from typing import List, Literal, TypedDict


class ParentStage(TypedDict):
    """One intermediate structure in the native composed Parent workflow."""

    name: Literal["fragment", "charge", "isotope", "stereo", "tautomer"]
    smiles: str


class SuperParentReport(TypedDict):
    """Final structure/status and ordered intermediate structures.

    Status remains str because native status strings can contain reason payloads
    and the underlying Rust enum is non-exhaustive. This is a structural typing
    contract; native report dictionaries are ordinary mutable Python dicts.
    """

    smiles: str
    status: str
    stages: List[ParentStage]
