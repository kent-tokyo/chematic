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
                              return_partial=False, return_preview=False,
                              preview_require_regular_geometry=False):
    """Generate an explicit-hydrogen mixed ensemble in an isolated process.

    The timeout includes interpreter startup and terminates an ongoing attempt.
    By default timeout results contain no partial ensemble. When return_partial
    is true, completed individual corrected-gradient attempts are retained on
    timeout or worker failure when their checkpoint is readable. With both
    return_partial and return_preview, a stereo-checked, unoptimized embedding
    for an interrupted attempt may be returned separately as in_progress.
    The optional regular-geometry screen applies only to that preview, not to
    the completed conformers.
    Coordinates use the atom order
    obtained by parsing ``smiles`` and adding hydrogens in the worker.
    """
    if not isinstance(smiles, str):
        raise TypeError("smiles must be a string")
    if not isinstance(num_confs, int) or isinstance(num_confs, bool) or num_confs < 0:
        raise ValueError("num_confs must be a nonnegative integer")
    if not isinstance(timeout_seconds, (int, float)) or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive and finite")
    if not isinstance(return_partial, bool) or not isinstance(return_preview, bool):
        raise TypeError("return_partial and return_preview must be bools")
    if return_preview and not return_partial:
        raise ValueError("return_preview requires return_partial=True")
    if not isinstance(preview_require_regular_geometry, bool):
        raise TypeError("preview_require_regular_geometry must be a bool")
    if preview_require_regular_geometry and not return_preview:
        raise ValueError("preview_require_regular_geometry requires return_preview=True")
    request = {"smiles": smiles, "num_confs": num_confs, "random_seed": random_seed,
               "max_iterations": max_iterations, "max_mmff_iters": max_mmff_iters,
               "repair_mmff_iters": repair_mmff_iters, "return_preview": return_preview,
               "preview_require_regular_geometry": preview_require_regular_geometry,
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
                _attach_checkpoint(response, request.get("checkpoint"))
                response["elapsed_seconds"] = time.monotonic() - started
                return response
            if process.returncode != 0:
                response = {"status": "worker_error", "error": stderr[-2000:]}
                _attach_checkpoint(response, request.get("checkpoint"))
                response["elapsed_seconds"] = time.monotonic() - started
                return response
            try:
                response = json.loads(stdout)
            except (ValueError, TypeError) as error:
                response = {"status": "worker_error", "error": str(error)}
            if not isinstance(response, dict) or response.get("status") not in ("completed", "error", "worker_error"):
                response = {"status": "worker_error", "error": "invalid worker response"}
            if response["status"] != "completed":
                _attach_checkpoint(response, request.get("checkpoint"))
            response["elapsed_seconds"] = time.monotonic() - started
            return response
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()


def _attach_checkpoint(response, path):
    if path is None:
        return
    try:
        response.update(_read_checkpoint(path, recover_prefix=True))
    except (OSError, ValueError) as error:
        # Preserve the worker's original failure instead of turning an
        # unreadable checkpoint into an unrelated exception.
        response["checkpoint_error"] = f"{type(error).__name__}: {error}"


def _empty_checkpoint():
    return {"conformers": [], "attempts_completed": 0, "embedded": 0,
            "attempt_diagnostics": [],
            "force_field_error": None, "stereo_violated": None,
            "stereo_unevaluable": None, "invalid_coordinates": None,
            "unclassified_rejections": None}


def _attempt_diagnostic(entry):
    counts = {"converged": 0, "unfinished": 0, "minimizer_failed": 0,
              "force_field_unavailable": 0, "other": 0}
    for row in entry["rows"]:
        if (not isinstance(row, (list, tuple)) or len(row) != 4
                or type(row[0]) is not int
                or (row[3] is not None and not isinstance(row[3], str))):
            raise ValueError("invalid checkpoint conformer")
        status, _, _, error = row
        if status == 0:
            counts["converged"] += 1
        elif status == 1:
            counts["unfinished"] += 1
        elif status == -1 and error is not None:
            counts["force_field_unavailable"] += 1
        elif status == -1:
            counts["minimizer_failed"] += 1
        else:
            counts["other"] += 1
    rejected = entry["embedded"] - len(entry["rows"])
    rejection_counts = (entry["violated"], entry["unevaluable"], entry["invalid"])
    known_rejections = sum(value for value in rejection_counts if value is not None)
    if rejected < known_rejections:
        raise ValueError("inconsistent checkpoint candidate counts")
    remaining = rejected - known_rejections
    unclassified = (remaining if all(value is not None for value in rejection_counts)
                    or remaining == 0 else None)
    return {"attempt": entry["attempt"], "embedded": entry["embedded"],
            "returned": len(entry["rows"]), "force_field_error": entry["error"],
            "stereo_violated": entry["violated"], "stereo_unevaluable": entry["unevaluable"],
            "invalid_coordinates": entry["invalid"], "unclassified_rejections": unclassified,
            "conformer_status_counts": counts}


def _apply_checkpoint_entry(result, counts_known, entry):
    if (not isinstance(entry, dict)
            or not {"attempt", "rows", "embedded", "error", "violated", "unevaluable", "invalid"} <= entry.keys()
            or type(entry.get("attempt")) is not int
            or entry["attempt"] != result["attempts_completed"] + 1
            or not isinstance(entry.get("rows"), list)
            or type(entry.get("embedded")) is not int or entry["embedded"] < 0
            or (entry.get("error") is not None and not isinstance(entry["error"], str))):
        raise ValueError("invalid checkpoint entry")
    values = (("stereo_violated", entry.get("violated")),
              ("stereo_unevaluable", entry.get("unevaluable")),
              ("invalid_coordinates", entry.get("invalid")))
    if any(value is not None and (type(value) is not int or value < 0) for _, value in values):
        raise ValueError("invalid checkpoint rejection count")
    diagnostic = _attempt_diagnostic(entry)
    values += (("unclassified_rejections", diagnostic["unclassified_rejections"]),)
    result["conformers"].extend(entry["rows"])
    result["attempt_diagnostics"].append(diagnostic)
    result["embedded"] += entry["embedded"]
    result["force_field_error"] = result["force_field_error"] or entry["error"]
    for key, value in values:
        if value is None:
            counts_known[key] = False
            result[key] = None
        elif counts_known[key]:
            result[key] = (result[key] or 0) + value
    result["attempts_completed"] += 1
    result.pop("in_progress", None)


def _apply_preview(result, entry):
    if (type(entry.get("preview_attempt")) is not int
            or entry["preview_attempt"] != result["attempts_completed"] + 1
            or "in_progress" in result or not isinstance(entry.get("coordinates"), list)):
        raise ValueError("invalid preview checkpoint")
    if "stereo_screen" in entry and entry["stereo_screen"] != "regular_geometry":
        raise ValueError("invalid preview stereo screen")
    if "geometry" in entry:
        geometry = entry["geometry"]
        if geometry is not None:
            if (not isinstance(geometry, dict)
                    or set(geometry) != {"bond_count", "worst_bond_length", "bond_violation_rate_50pct", "gross_clash_count"}
                    or any(type(geometry[key]) is not int or geometry[key] < 0 for key in ("bond_count", "gross_clash_count"))
                    or any(value is not None and (type(value) not in (int, float) or not math.isfinite(value))
                           for value in (geometry["worst_bond_length"], geometry["bond_violation_rate_50pct"]))):
                raise ValueError("invalid preview geometry diagnostics")
        if entry.get("geometry_error") is not None and not isinstance(entry["geometry_error"], str):
            raise ValueError("invalid preview geometry error")
    result["in_progress"] = {"attempt": entry["preview_attempt"],
                             "coordinates": entry["coordinates"]}
    if "stereo_screen" in entry:
        result["in_progress"]["stereo_screen"] = entry["stereo_screen"]
    if "geometry" in entry:
        result["in_progress"]["geometry"] = entry["geometry"]
        result["in_progress"]["geometry_error"] = entry.get("geometry_error")


def _read_checkpoint(path, *, recover_prefix=False):
    result = _empty_checkpoint()
    counts_known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates",
                                         "unclassified_rejections")}
    records_read = 0
    try:
        with open(path, "rb") as stream:
            for line_number, line in enumerate(stream, 1):
                # A killed process can leave a partial final write. Only a
                # newline-terminated record is a completed checkpoint.
                if not line.endswith(b"\n"):
                    break
                try:
                    entry = json.loads(line)
                    if isinstance(entry, dict) and "preview_attempt" in entry:
                        _apply_preview(result, entry)
                    else:
                        _apply_checkpoint_entry(result, counts_known, entry)
                except ValueError as error:
                    if not recover_prefix or records_read == 0:
                        raise
                    result["checkpoint_error"] = f"line {line_number}: {type(error).__name__}: {error}"
                    break
                records_read += 1
    except OSError as error:
        if isinstance(error, FileNotFoundError) and records_read == 0:
            pass
        elif recover_prefix and records_read:
            result["checkpoint_error"] = f"{type(error).__name__}: {error}"
        else:
            raise
    return result


def _run_partial(mol, request):
    # Individual mixed-ensemble calls with one attempt use the corrected
    # gradient. This is intentionally a separate profile from a mixed batch.
    result = _empty_checkpoint()
    counts_known = {key: True for key in ("stereo_violated", "stereo_unevaluable", "invalid_coordinates",
                                         "unclassified_rejections")}
    prepared = None
    if request.get("return_preview") and request["num_confs"]:
        base = request["random_seed"]
        if not isinstance(base, int) or not 0 <= base <= 2**31 - 1:
            raise ValueError("conformer seed outside signed 32-bit nonnegative range")
        prepared = mol._timed_prepare_corrected(max_iterations=request["max_iterations"])
    for index in range(request["num_confs"]):
        base = request["random_seed"]
        seed = (index + 1 if base == 0 and index else base * (index + 1))
        if not isinstance(seed, int) or not 0 <= seed <= 2**31 - 1:
            raise ValueError("conformer seed outside signed 32-bit nonnegative range")
        if request.get("return_preview"):
            coords = prepared.embed(seed)
            if coords is None:
                rows, embedded, error, violated, unevaluable, invalid = [], 0, None, None, None, None
            else:
                regular_screen = request.get("preview_require_regular_geometry", False)
                valid = (mol._timed_stereo_valid(coords, require_regular_geometry=True)
                         if regular_screen else mol._timed_stereo_valid(coords))
                if valid:
                    preview = {"preview_attempt": index + 1, "coordinates": coords}
                    if regular_screen:
                        preview["stereo_screen"] = "regular_geometry"
                    try:
                        count, worst, violation_rate, clashes = mol._timed_preview_geometry(coords)
                        preview["geometry"] = {"bond_count": count, "worst_bond_length": worst,
                                               "bond_violation_rate_50pct": violation_rate,
                                               "gross_clash_count": clashes}
                        preview["geometry_error"] = None
                    except Exception as geometry_error:
                        preview["geometry"] = None
                        preview["geometry_error"] = f"{type(geometry_error).__name__}: {geometry_error}"
                    _apply_preview(result, preview)
                    _append_checkpoint(request["checkpoint"], preview)
                rows, embedded, error, violated, unevaluable, invalid = mol._timed_optimize_embedded(
                    coords, max_mmff_iters=request["max_mmff_iters"],
                    repair_mmff_iters=request["repair_mmff_iters"])
        else:
            rows, embedded, error, violated, unevaluable, invalid = (
                mol.etkdg_embed_multiple_confs_mixed_checked_diagnostics(
                    1, random_seed=seed, max_iterations=request["max_iterations"],
                    max_mmff_iters=request["max_mmff_iters"],
                    repair_mmff_iters=request["repair_mmff_iters"])
            )
        entry = {"attempt": index + 1, "rows": rows, "embedded": embedded, "error": error,
                 "violated": violated, "unevaluable": unevaluable, "invalid": invalid}
        _apply_checkpoint_entry(result, counts_known, entry)
        _append_checkpoint(request["checkpoint"], entry)
    return result


def _append_checkpoint(path, entry):
    payload = (json.dumps(entry, allow_nan=False, separators=(",", ":")) + "\n").encode("utf-8")
    with open(path, "ab") as stream:
        stream.write(payload)


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
