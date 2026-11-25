"""Run the experimental timed conformer API from an importable program.

Install a chematic build that contains ``timed_conformer_batch`` first; the
published 1.2.4 wheel does not contain this example API. For example::

    python timed_conformer_batch.py CCO 'C1=CC=CC=C1' invalid
    python timed_conformer_batch.py --allow-unoptimized 'B(O)O'

Each molecule has a 0.2-second generation deadline after worker startup.
Timeout terminates the worker and preserves candidates already received;
process cleanup can extend elapsed time beyond 0.2 seconds. An optional MMFF
fallback is explicitly unoptimized: its status is -1 and energy is None.

On one macOS host the 2,029-molecule experimental candidate evaluation gave
1,280 recoveries versus RDKit 1,173 with 0.2-second deadlines. The candidate's
median elapsed time including cleanup was about 0.211 seconds. This describes
that fixed input and host, not a general speed or runtime guarantee. See
benchmarks/2026-11-24-candidate-v3-timed-summary.json for the full result.
"""

import argparse
import json

from chematic import timed_conformer_batch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="+", help="SMILES strings to process")
    parser.add_argument("--allow-unoptimized", action="store_true", help="permit explicitly marked fallback after MMFF type-assignment failure")
    args = parser.parse_args()
    results = timed_conformer_batch(args.smiles, seconds=0.2, num_confs=20,
                                    repair_mmff_iters=400,
                                    allow_unoptimized=args.allow_unoptimized)
    for smiles, result in zip(args.smiles, results):
        conformers = result["conformers"]
        print(json.dumps({
            "smiles": smiles,
            "status": result["status"],
            "elapsed_seconds": result["elapsed"],
            "conformers": len(conformers),
            "converged": sum(status == 0 for status, _, _ in conformers),
            "unoptimized_fallback": sum(status == -1 for status, _, _ in conformers),
            "fallback_reason": result["fallback_reason"],
            "error": result["error"],
        }, ensure_ascii=False))


if __name__ == "__main__":
    main()
