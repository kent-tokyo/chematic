"""The geometry comparison must not count crashes or usage text as safe results."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("gjf_validation", Path(__file__).parents[1] / "validation/validate_gjf_cartesian.py")
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


class ProbeObservation(unittest.TestCase):
    def test_only_the_explicit_parser_error_exit_is_a_refusal(self):
        for code, expected in [(2, "refused"), (0, "process_error"), (1, "process_error"), (-11, "process_error")]:
            self.assertEqual(validation.probe_result({"stdout": "", "exit_code": code}, None), expected)

    def test_usage_output_is_never_a_geometry(self):
        self.assertEqual(validation.probe_result({"stdout": "Usage: obabel", "exit_code": 1}, None), "process_error")
        self.assertEqual(validation.probe_result({"stdout": "Usage: obabel", "exit_code": 0}, None), "unreadable_output")

    def test_nonfinite_and_partial_json_are_not_matches(self):
        for output in ['{"coords":[[NaN,0,0]]}', '{"coords":[]}']:
            self.assertEqual(validation.probe_result({"stdout": output, "exit_code": 0}, {}), "unreadable_output")

    def test_emission_for_a_refused_profile_is_not_a_match(self):
        output = '{"coords":[[1,2,3]],"atom_count":1,"charge":0,"multiplicity":1}'
        self.assertEqual(validation.probe_result({"stdout": output, "exit_code": 0}, None), "emitted_mismatch")


if __name__ == "__main__":
    unittest.main()
