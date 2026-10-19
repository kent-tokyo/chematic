from typing import Optional
from chematic import InchiKeyRecord, InchiKeyStatus, inchi_keys_checked


def convert(inputs: list[str]) -> list[Optional[str]]:
    rows: tuple[InchiKeyRecord, ...] = inchi_keys_checked(iter(inputs))
    for row in rows:
        status: InchiKeyStatus = row.status
        assert status in (
            "ok",
            "invalid",
            "unsupported",
            "resource_limit",
            "conversion_error",
        )
    return [row.key for row in rows]
