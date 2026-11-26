"""Process-isolated, wall-clock-limited conformer generation."""

import contextlib
import json
import math
import os
import subprocess
import sys
import tempfile
import time


def generate_conformers_timed(smiles, num_confs, *, timeout_seconds,
                              random_seed=42, max_iterations=0,
                              max_mmff_iters=200, repair_mmff_iters=0,
                              return_partial=False):
    """Generate an explicit-hydrogen mixed ensemble in an isolated process.

    The timeout includes interpreter startup and terminates an ongoing attempt.
    By default timeout results contain no partial ensemble. When return_partial
    is true, completed individual corrected-gradient attempts are retained.
    Coordinates use the atom order
    obtained by parsing ``smiles`` and adding hydrogens in the worker.
    """
    if not isinstance(smiles, str):
        raise TypeError("smiles must be a string")
    if not isinstance(num_confs, int) or isinstance(num_confs, bool) or num_confs < 0:
        raise ValueError("num_confs must be a nonnegative integer")
    if not isinstance(timeout_seconds, (int, float)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive and finite")
    if not isinstance(return_partial, bool):
        raise TypeError("return_partial must be a bool")
    request = {"smiles": smiles, "num_confs": num_confs, "random_seed": random_seed,
               "max_iterations": max_iterations, "max_mmff_iters": max_mmff_iters,
               "repair_mmff_iters": repair_mmff_iters,
               "package_file": os.path.realpath(sys.modules[__package__].__file__)}
    with tempfile.TemporaryDirectory(prefix="chematic-conformers-") if return_partial else contextlib.nullcontext() as directory:
        if directory is not None:
            request["checkpoint"] = os.path.join(directory, "completed.json")
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
                response = {"status": "timeout"}
                if directory is not None:
                    response.update(_read_checkpoint(request["checkpoint"]))
                response["elapsed_seconds"] = time.monotonic() - started
                return response
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


def _read_checkpoint(path):
    try:
        with open(path, encoding="utf-8") as stream:
            result = json.load(stream)
        if not isinstance(result, dict) or not isinstance(result.get("conformers"), list):
            raise ValueError("invalid checkpoint")
        return result
    except FileNotFoundError:
        return {"conformers": [], "attempts_completed": 0, "embedded": 0,
                "force_field_error": None, "stereo_violated": None,
                "stereo_unevaluable": None, "invalid_coordinates": None}


def _run_partial(mol, request):
    # Individual mixed-ensemble calls with one attempt use the corrected
    # gradient. This is intentionally a separate profile from a mixed batch.
    result = {"conformers": [], "attempts_completed": 0, "embedded": 0,
              "force_field_error": None, "stereo_violated": None,
              "stereo_unevaluable": None, "invalid_coordinates": None}
    counts_known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates")}
    for index in range(request["num_confs"]):
        base = request["random_seed"]
        seed = (index + 1 if base == 0 and index else base * (index + 1))
        if not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1:
            raise ValueError("conformer seed outside signed 32-bit nonnegative range")
        rows, embedded, error, violated, unevaluable, invalid = (
            mol.etkdg_embed_multiple_confs_mixed_checked_diagnostics(
                1, random_seed=seed, max_iterations=request["max_iterations"],
                max_mmff_iters=request["max_mmff_iters"],
                repair_mmff_iters=request["repair_mmff_iters"])
        )
        result["conformers"].extend(rows)
        result["embedded"] += embedded
        result["force_field_error"] = result["force_field_error"] or error
        for key, value in (("stereo_violated", violated), ("stereo_unevaluable", unevaluable),
                           ("invalid_coordinates", invalid)):
            if value is None:
                counts_known[key] = False
                result[key] = None
            elif counts_known[key]:
                result[key] = (result[key] or 0) + value
        result["attempts_completed"] += 1
        temporary = request["checkpoint"] + ".tmp"
        with open(temporary, "w", encoding="utf-8") as stream:
            json.dump(result, stream, allow_nan=False)
        os.replace(temporary, request["checkpoint"])
    return result


def _main():
    try:
        import chematic
        request = json.load(sys.stdin)
        if os.path.realpath(chematic.__file__) != request["package_file"]:
            raise RuntimeError("worker imported a different chematic package")
        mol = chematic.from_smiles(request["smiles"]).add_hydrogens()
        if "checkpoint" in request:
            response = {"status": "completed", **_run_partial(mol, request)}
        else:
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
