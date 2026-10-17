#!/usr/bin/env python3
"""Original bit corpora; compare checked adapter outputs to live RDKit.

Only ordinary writable contiguous arrays are passed into RDKit. Native strided
outputs are checked against the independent on-bit reference and RDKit's normal
output values, not against RDKit's behavior on strided/invalid destinations.
"""
import gzip
import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

import chematic
import numpy as np
from chematic.rdkit_compat import DataStructs as adapter, ExplicitBitVect
from rdkit import DataStructs as rdkit, rdBase


class Source:
    def __init__(self, width, bits):
        self.width, self.bits = width, bits

    def GetNumBits(self):
        return self.width

    def GetOnBits(self):
        return self.bits


def vector(cls, width, bits):
    out = cls(width)
    for bit in bits:
        out.SetBit(bit)
    return out


def sources(width, bits):
    yield "chematic_explicit", vector(ExplicitBitVect, width, bits)
    yield "rdkit_explicit", vector(rdkit.ExplicitBitVect, width, bits)
    yield "rdkit_sparse", vector(rdkit.SparseBitVect, width, bits)
    yield "lazy_ids", Source(width, (i for i in bits + bits[:3]))
    view = np.zeros(len(bits) * 2, dtype=np.uint32)
    view[::2] = bits
    yield "numpy_ids", Source(np.int64(width), view[::2])


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    records, failures = [], []
    runs = oracle_runs = 0
    dtypes = [np.bool_, np.int8, np.uint16, np.int64, np.float32, np.float64]
    for width, density, seed, dtype in itertools.product(
        [1, 7, 8, 9, 31, 64, 128, 2048], [0, 0.2, 1], range(4), dtypes
    ):
        rng = random.Random(seed)
        bits = [i for i in range(width) if rng.random() < density]
        on = set(bits)
        expected = np.array([int(i in on) for i in range(width)], dtype=dtype)
        rdkit_outputs = {}
        for name, cls in [
            ("explicit", rdkit.ExplicitBitVect),
        ]:
            result = np.full(width, 7, dtype=dtype)
            returned = rdkit.ConvertToNumpyArray(vector(cls, width, bits), result)
            oracle_runs += 1
            if not np.array_equal(result, expected):
                failures.append([len(records), "rdkit", name])
            rdkit_outputs[name] = dict(
                values=result.tolist(), return_is_none=returned is None
            )
        outputs = {}
        for mode in ["contiguous", "strided", "reversed"]:
            for name, fp in sources(width, bits):
                backing = np.full(width * 2 + 1, 7, dtype=dtype)
                dest = backing[1::2]
                if mode == "contiguous":
                    dest = np.full(width, 7, dtype=dtype)
                if mode == "reversed":
                    dest = dest[::-1]
                returned = adapter.ConvertToNumpyArray(fp, dest)
                runs += 1
                if returned is not dest or not np.array_equal(dest, expected):
                    failures.append([len(records), mode, name])
                if mode != "contiguous" and not np.all(
                    backing[::2] == np.asarray(7, dtype=dtype)
                ):
                    failures.append([len(records), "padding overwritten", name])
                outputs[mode + "/" + name] = dest.tolist()
        for name, fp in sources(width, bits):
            result = adapter.ConvertToNumpyArray(fp)
            runs += 1
            if (
                result.dtype != np.int8
                or result.tolist() != expected.astype(np.int8).tolist()
            ):
                failures.append([len(records), "optional destination", name])
        records.append(
            dict(
                width=width,
                density=density,
                seed=seed,
                dtype=np.dtype(dtype).str,
                on_bits=bits,
                reference=expected.tolist(),
                rdkit=rdkit_outputs,
                adapter=outputs,
            )
        )
    try:
        rdkit.ConvertToNumpyArray(
            vector(rdkit.SparseBitVect, 8, [2]), np.zeros(8, dtype=np.int8)
        )
        sparse_diagnostic = {"accepted": True}
    except Exception as error:
        sparse_diagnostic = {
            "accepted": False,
            "error_type": type(error).__name__,
            "error": str(error),
        }
    summary = dict(
        schema_version=1,
        scope="binary bit projection; ordinary RDKit outputs plus independently checked adapter strided destinations",
        datasets=len(records),
        adapter_calls=runs,
        rdkit_calls=oracle_runs,
        failures=failures,
        failure_count=len(failures),
        rdkit_version=rdBase.rdkitVersion,
        rdkit_sparse_binary_conversion=sparse_diagnostic,
        numpy_version=np.__version__,
        python_version=sys.version,
        adapter_sha256=hashlib.sha256(
            Path(chematic.rdkit_compat.__file__).read_bytes()
        ).hexdigest(),
        helper_sha256=hashlib.sha256(
            (
                Path(chematic.rdkit_compat.__file__).parent / "_numpy_fingerprint.py"
            ).read_bytes()
        ).hexdigest(),
        reused_extension_sha256=hashlib.sha256(
            Path(chematic.chematic.__file__).read_bytes()
        ).hexdigest(),
        extension_scope="Only reused to import the existing package. Modified conversion code is pure Python/NumPy and does not call chematic Rust APIs; no main-wheel build claim.",
    )
    with gzip.open(out / "corpus-results.json.gz", "wt") as stream:
        json.dump(records, stream, separators=(",", ":"))
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
