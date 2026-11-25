"""Run checked conformer generation in a replaceable worker process."""

import math
import multiprocessing
import threading
import time


def _volume(a, b, c, d):
    x = [a[i] - d[i] for i in range(3)]
    y = [b[i] - d[i] for i in range(3)]
    z = [c[i] - d[i] for i in range(3)]
    return (x[0] * (y[1] * z[2] - y[2] * z[1])
            - x[1] * (y[0] * z[2] - y[2] * z[0])
            + x[2] * (y[0] * z[1] - y[1] * z[0]))


def _stereo_checked(mol, coords):
    if not all(math.isfinite(x) for point in coords for x in point):
        return False
    report = mol.check_geometry(coords)
    stereo = report["stereo"]
    if (report["stereo_groups_not_evaluated"] or stereo["n_violations"]
            or stereo["n_unevaluable"]):
        return False
    bonds = mol.bond_table
    for center in stereo["tetrahedral"]:
        atom = center["atom"]
        neighbors = [b if a == atom else a for a, b, _, _ in bonds if a == atom or b == atom]
        if len(neighbors) != 4:
            continue
        points = [coords[i] for i in neighbors]
        total = _volume(*points)
        if not math.isfinite(total) or abs(total) <= 1e-12:
            return False
        for i in range(4):
            vertex = list(points)
            vertex[i] = coords[atom]
            if _volume(*vertex) / total <= 0:
                return False
    return True


def _worker(connection):
    from . import from_smiles

    connection.send(None)
    while True:
        request = connection.recv()
        if request is None:
            return
        smiles, candidates, seed, iterations, repair_iterations, allow_unoptimized = request
        try:
            mol = from_smiles(smiles).add_hydrogens()
            fallback_reason = None
            for attempt in range(candidates):
                current_seed = seed * (attempt + 1)
                if fallback_reason is None:
                    try:
                        result = mol.etkdg_embed_multiple_confs_mmff_corrected(
                            1, random_seed=current_seed, max_mmff_iters=iterations,
                            repair_mmff_iters=repair_iterations,
                        )
                    except ValueError as exc:
                        if not allow_unoptimized or "MMFF94 type assignment failed" not in str(exc):
                            raise
                        fallback_reason = str(exc)
                        connection.send(("fallback", fallback_reason))
                if fallback_reason is not None:
                    result = []
                    for coords in mol.rdkit_embed_multiple_confs(1, random_seed=current_seed):
                        if _stereo_checked(mol, coords):
                            result.append((-1, None, coords))
                else:
                    result = [candidate for candidate in result if _stereo_checked(mol, candidate[2])]
                connection.send(("candidate", result))
            connection.send(("ok", None))
        except Exception as exc:
            connection.send(("error", f"{type(exc).__name__}: {exc}"))


def timed_conformer_batch(smiles, *, seconds=0.2, num_confs=1, random_seed=42,
                          max_mmff_iters=200, repair_mmff_iters=0,
                          allow_unoptimized=False):
    """Generate checked conformers with a per-molecule process deadline.

    Returns one dictionary per input with status, conformers, error, and elapsed
    seconds. Each candidate uses seed ``random_seed * (attempt + 1)``. Worker
    startup is excluded from each molecule's deadline. On timeout, completed
    candidates received within the deadline are retained, and the worker is
    terminated and replaced before the next molecule. Elapsed
    time includes process cleanup, which can extend beyond the deadline.
    Call from an importable Python program under the usual __main__ guard::

        from chematic import timed_conformer_batch

        if __name__ == "__main__":
            for result in timed_conformer_batch(["CCO", "c1ccccc1"], seconds=0.2):
                print(result["status"], len(result["conformers"]))

    A candidate is ``(status, energy, coordinates)``. Status 0 indicates
    convergence; status 1 may contain a stereo-consistent unoptimized geometry.
    With ``allow_unoptimized=True``, an MMFF type-assignment failure permits
    stereo-checked embedded coordinates with status -1 and energy None.
    The result's ``fallback_reason`` records the type-assignment error.
    """
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("seconds must be finite and positive")
    if num_confs < 1 or max_mmff_iters < 0 or repair_mmff_iters < 0:
        raise ValueError("invalid conformer or iteration count")
    context = multiprocessing.get_context("spawn")
    process = None
    connection = None
    results = []
    try:
        for text in smiles:
            if not isinstance(text, str):
                raise TypeError("each SMILES input must be a string")
            if process is None:
                connection, child = context.Pipe()
                process = context.Process(target=_worker, args=(child,))
                process.start()
                child.close()
                if not connection.poll(30):
                    raise RuntimeError("conformer worker did not start")
                connection.recv()
            start = time.monotonic()
            deadline = start + seconds
            response = []
            candidates = []
            fallback_reasons = []
            transport_failed = []
            finished = threading.Event()

            def receive():
                try:
                    connection.send((text, num_confs, random_seed, max_mmff_iters,
                                     repair_mmff_iters, allow_unoptimized))
                    while True:
                        kind, value = connection.recv()
                        received = time.monotonic()
                        if kind == "candidate":
                            if received <= deadline:
                                candidates.extend(value)
                        elif kind == "fallback":
                            if received <= deadline:
                                fallback_reasons.append(value)
                        else:
                            response.append(((kind, value), received))
                            break
                except Exception as exc:
                    transport_failed.append(True)
                    response.append((("error", f"worker exited: {exc}"), time.monotonic()))
                finally:
                    finished.set()

            reader = threading.Thread(target=receive, daemon=True)
            reader.start()
            finished.wait(max(0.0, deadline - time.monotonic()))
            if finished.is_set() and response[0][1] <= deadline:
                status, value = response[0][0]
            else:
                status, value = "timeout", None
            if process is not None and (status == "timeout" or transport_failed or not process.is_alive()):
                if process is not None and process.is_alive():
                    process.terminate()
                if process is not None:
                    process.join(timeout=1)
                    if process.is_alive():
                        process.kill()
                        process.join()
                    reader.join()
                    connection.close()
                process = connection = None
            results.append({"status": status, "conformers": candidates,
                            "error": value if status == "error" else None,
                            "fallback_reason": fallback_reasons[0] if fallback_reasons else None,
                            "elapsed": time.monotonic() - start})
    finally:
        if process is not None:
            if process.is_alive():
                try:
                    connection.send(None)
                except (BrokenPipeError, EOFError, OSError):
                    pass
            process.join(timeout=1)
            if process.is_alive():
                process.terminate()
                process.join(timeout=1)
            if process.is_alive():
                process.kill()
                process.join()
            connection.close()
    return results
