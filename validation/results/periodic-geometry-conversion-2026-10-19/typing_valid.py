from chematic import (
    PeriodicGeometryConversion,
    poscar_to_turbomole_checked,
    turbomole_to_poscar_checked,
)


def convert(text: str) -> tuple[str, tuple[int, ...]]:
    result: PeriodicGeometryConversion = turbomole_to_poscar_checked(text)
    reverse: PeriodicGeometryConversion = poscar_to_turbomole_checked(result.text)
    return reverse.text, result.source_indices
