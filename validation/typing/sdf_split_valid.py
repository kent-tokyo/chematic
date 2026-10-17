from pathlib import Path
from typing import Iterator
from chematic import SdfSplitPart, SdfSplitPlan, plan_sdf_split


def split(path: Path) -> Iterator[bytes]:
    plan: SdfSplitPlan = plan_sdf_split(path, 3, prefix="batch")
    part: SdfSplitPart = plan.parts[0]
    assert part.record_count > 0
    return plan.iter_part(0, chunk_bytes=1024)
