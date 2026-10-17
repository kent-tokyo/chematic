"""Compare the checked Python alignment API with RDKit on original point clouds."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
import rdkit
from rdkit.Numerics import rdAlignment as reference
from rdkit.Geometry import Point3D
from chematic.rdkit_compat import rdAlignment as candidate
import chematic


def main():
    rng = np.random.default_rng(9601)
    cases = 0
    metamorphic = 0
    max_ssd_difference = 0.0
    max_transform_difference = 0.0
    checksum = hashlib.sha256()
    for n in [1, 2, 3, 4, 8, 16]:
        for sample in range(16):
            ref = rng.normal(size=(n, 3))
            probe = rng.normal(size=(n, 3)) + rng.uniform(-3.0, 3.0, 3)
            for weights in [None, [], rng.uniform(0.05, 4.0, n)]:
                w = np.ones(n) if weights is None or len(weights) == 0 else weights
                for reflect in [False, True]:
                    kind = sample % 3
                    if kind == 0:
                        a, b = ref.tolist(), probe.tolist()
                    elif kind == 1:
                        a, b = ref, probe
                    else:
                        a = [Point3D(*p) for p in ref]
                        b = [Point3D(*p) for p in probe]
                    args = (a, b) if weights is None else (a, b, weights)
                    ssd, matrix = candidate.GetAlignmentTransform(
                        *args, reflect=reflect
                    )
                    r_ssd, r_matrix = reference.GetAlignmentTransform(
                        *args, reflect=reflect
                    )
                    np.testing.assert_allclose(ssd, r_ssd, rtol=1e-8, atol=1e-9)
                    delta = probe @ matrix[:3, :3].T + matrix[:3, 3] - ref
                    direct = float(np.sum(w[:, None] * delta**2))
                    np.testing.assert_allclose(ssd, direct, rtol=1e-10, atol=1e-10)
                    np.testing.assert_allclose(
                        matrix[:3, :3] @ matrix[:3, :3].T, np.eye(3), atol=1e-12
                    )
                    np.testing.assert_allclose(
                        np.linalg.det(matrix[:3, :3]),
                        -1.0 if reflect else 1.0,
                        atol=1e-12,
                    )
                    max_ssd_difference = max(max_ssd_difference, abs(ssd - r_ssd))
                    if n >= 4:
                        # Eigen/SVD tolerances differ; matrix agreement is diagnostic.
                        max_transform_difference = max(
                            max_transform_difference,
                            float(np.abs(matrix - r_matrix).max()),
                        )
                    checksum.update(
                        ref.tobytes() + probe.tobytes() + w.tobytes() + bytes([reflect])
                    )
                    scaled, _ = candidate.GetAlignmentTransform(
                        ref * 3.0, probe * 3.0, weights, reflect
                    )
                    np.testing.assert_allclose(
                        scaled, ssd * 9.0, rtol=1e-10, atol=1e-10
                    )
                    multiplied, _ = candidate.GetAlignmentTransform(
                        ref, probe, w * 7.0, reflect
                    )
                    np.testing.assert_allclose(
                        multiplied, ssd * 7.0, rtol=1e-10, atol=1e-10
                    )
                    metamorphic += 2
                    cases += 1
    extension = Path(chematic.chematic.__file__)
    print(
        json.dumps(
            {
                "scope": "paired point-cloud weighted SSD; no molecule matching or conformer generation",
                "rdkit_version": rdkit.__version__,
                "numpy_version": np.__version__,
                "reference_cases": cases,
                "metamorphic_checks": metamorphic,
                "corpus_sha256": checksum.hexdigest(),
                "maximum_absolute_ssd_difference": max_ssd_difference,
                "maximum_nondegenerate_transform_difference_diagnostic": max_transform_difference,
                "native_fixture_sha256": hashlib.sha256(
                    extension.read_bytes()
                ).hexdigest(),
                "native_fixture_source_commit": "e2041226b67030108acb00309da5f0659f636c1d",
                "native_fixture_scope": "reused import fixture; alignment implementation is Python/NumPy only; no fresh Rust build",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
