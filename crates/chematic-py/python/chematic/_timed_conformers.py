"""Process-isolated, wall-clock-limited conformer generation."""

import json
import math
import os
import subprocess
import sys
import time


def generate_conformers_timed(smiles, num_confs, *, timeout_seconds,
                              random_seed=42, max_iterations=0,
                              max_mmff_iters=200, repair_mmff_iters=0):
    """Generate an explicit-hydrogen mixed ensemble in an isolated process.

    The timeout includes interpreter startup and terminates an ongoing attempt.
    Timeout results contain no partial ensemble. Coordinates use the atom order
    obtained by parsing ``smiles`` and adding hydrogens in the worker.
    """
    if not isinstance(smiles, str):
        raise TypeError("smiles must be a string")
    if not isinstance(num_confs, int) or isinstance(num_confs, bool) or num_confs < 0:
        raise ValueError("num_confs must be a nonnegative integer")
    if not isinstance(timeout_seconds, (int, float)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive and finite")
    request = {"smiles": smiles, "num_confs": num_confs, "random_seed": random_seed,
               "max_iterations": max_iterations, "max_mmff_iters": max_mmff_iters,
               "repair_mmff_iters": repair_mmff_iters,
               "package_file": os.path.realpath(sys.modules[__package__].__file__)}
    started = time.monotonic()
    process = subprocess.Popen([sys.executable, "-c", "from chematic._timed_conformers import _main; _main()"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    try:
        try:
            remaining = max(0.0, timeout_seconds - (time.monotonic() - started))
            stdout, stderr = process.communicate(json.dumps(request), timeout=remaining)
        except subprocess.TimeoutExpired:
            process.kill()
            process.communicate()
            return {"status": "timeout", "elapsed_seconds": time.monotonic() - started}
        if process.returncode != 0:
            return {"status": "worker_error", "error": stderr[-2000:],
                    "elapsed_seconds": time.monotonic() - started}
        try:
            response = json.loads(stdout)
        except (ValueError, TypeError) as error:
            return {"status": "worker_error", "error": str(error),
                    "elapsed_seconds": time.monotonic() - started}
        response["elapsed_seconds"] = time.monotonic() - started
        return response
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()


def _main():
    try:
        import chematic
        request = json.load(sys.stdin)
        if os.path.realpath(chematic.__file__) != request["package_file"]:
            raise RuntimeError("worker imported a different chematic package")
        mol = chematic.from_smiles(request["smiles"]).add_hydrogens()
        rows, embedded, force_field_error, violated, unevaluable, invalid = (
            mol.etkdg_embed_multiple_confs_mixed_checked_diagnostics(
                request["num_confs"], random_seed=request["random_seed"],
                max_iterations=request["max_iterations"], max_mmff_iters=request["max_mmff_iters"],
                repair_mmff_iters=request["repair_mmff_iters"])
        )
        response = {"status": "completed", "conformers": rows, "embedded": embedded,
                    "force_field_error": force_field_error, "stereo_violated": violated,
                    "stereo_unevaluable": unevaluable, "invalid_coordinates": invalid}
    except Exception as error:
        response = {"status": "error", "error_type": type(error).__name__, "error": str(error)}
    print(json.dumps(response, allow_nan=False))


if __name__ == "__main__":
    _main()
