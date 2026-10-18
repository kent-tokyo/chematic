"""Independent ASE reads of positive-NVAL inputs and every scalar export."""

import hashlib
from io import StringIO
import itertools
import json
from pathlib import Path

import ase
from ase.io.cube import read_cube
from ase.units import Bohr
import numpy as np

from chematic import read_cube_datasets

HERE = Path(__file__).parent
corpus = hashlib.sha256()
documents = selected_fields = selected_values = positive_fields = 0
negative_accepted = negative_refused = 0
max_value_error = max_geometry_error = 0.0

for case in range(96):
    shape = (case % 4 + 1, case % 3 + 2, case % 5 + 1)
    count = case % 4 + 1
    ids = [11, 29, 5, 83][:count]
    points = list(itertools.product(*(range(n) for n in shape)))
    fields = [
        [100 * d + 10 * x + 2 * y + 0.125 * z for x, y, z in points]
        for d in range(count)
    ]
    interleaved = [fields[d][i] for i in range(len(points)) for d in range(count)]
    data = "\n".join(
        " ".join(format(v, ".17g") for v in interleaved[i : i + 7])
        for i in range(0, len(interleaved), 7)
    )
    for negative in (False, True):
        count_line = "-2 .125 -.25 .375" if negative else f"2 .125 -.25 .375 {count}"
        header = f"synthetic case {case}\nopaque scalar fields\n{count_line}\n{shape[0]} .5 .125 0\n{shape[1]} 0 .75 .25\n{shape[2]} .375 0 1\n6 4.0 .125 .25 .375\n8 8.0 -.5 .625 .75\n"
        labels = f"{count}\n" + "\n".join(map(str, ids)) + "\n" if negative else ""
        text = header + labels + data + "\n"
        corpus.update(text.encode())
        doc = read_cube_datasets(text)
        assert doc.source_text == text
        assert doc.dataset_ids == (tuple(ids) if negative else None)
        np.testing.assert_array_equal(doc.datasets, fields)
        documents += 1
        if not negative:
            original = read_cube(StringIO(text))
            np.testing.assert_array_equal(original["datas"].reshape(count, -1), fields)
            positive_fields += count
        else:
            try:
                original = read_cube(StringIO(text))
            except ValueError:
                negative_refused += 1
            else:
                np.testing.assert_array_equal(
                    original["datas"].reshape(count, -1), fields
                )
                negative_accepted += 1
        for d in range(count):
            selection = doc.select(d)
            assert selection.dataset_id == (ids[d] if negative else None)
            assert selection.omitted_indices == tuple(i for i in range(count) if i != d)
            assert int(selection.text.splitlines()[2].split()[0]) == 2
            reference = read_cube(StringIO(selection.text))
            actual = selection.to_grid()
            values = reference["data"].reshape(-1)
            np.testing.assert_array_equal(values, fields[d])
            np.testing.assert_array_equal(actual.values, values)
            max_value_error = max(
                max_value_error, float(np.max(np.abs(actual.values - values)))
            )
            for a, b in [
                (reference["origin"], np.array(actual.origin) * Bohr),
                (reference["spacing"], actual.axes * Bohr),
                (
                    reference["atoms"].positions,
                    np.array([atom[2] for atom in actual.atoms]) * Bohr,
                ),
            ]:
                np.testing.assert_allclose(a, b, atol=1e-14, rtol=0)
                max_geometry_error = max(
                    max_geometry_error, float(np.max(np.abs(a - b)))
                )
            selected_fields += 1
            selected_values += len(points)

summary = dict(
    ase=ase.__version__,
    numpy=np.__version__,
    input_documents=documents,
    positive_nval_original_fields=positive_fields,
    selected_single_fields=selected_fields,
    selected_scalar_values=selected_values,
    max_value_error=max_value_error,
    max_geometry_error_angstrom=max_geometry_error,
    ase_negative_original_accepted=negative_accepted,
    ase_negative_original_refused=negative_refused,
    inputs_sha256=corpus.hexdigest(),
    scope="ASE positive-NVAL original reads and all scalar exports; negative multi-dataset original refusals reported separately; exact authored values independently checked",
)
(HERE / "oracle-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
