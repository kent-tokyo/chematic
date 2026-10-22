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


writer_spec = importlib.util.spec_from_file_location("gjf_write_validation", Path(__file__).parents[1] / "validation/validate_gjf_write.py")
writer = importlib.util.module_from_spec(writer_spec)
writer_spec.loader.exec_module(writer)


class WriteObservation(unittest.TestCase):
    def test_five_fields_are_observed_as_flag_and_three_coordinates(self):
        text = "# HF Units=Ang\n\noriginal\n\n-1 2\nC -1 1.25 2.5 3.75\n\n"
        result = writer.observe(text)
        self.assertEqual(result["coords"], [[1.25, 2.5, 3.75]])
        self.assertEqual(result["freeze_codes"], [-1])

    def test_corrupted_or_partial_output_cannot_be_a_geometry_match(self):
        base = "# HF Units=Ang\n\noriginal\n\n0 1\nC 1 2 3\n\n"
        for bad in [base.replace("C 1 2 3", "C 1 2"), base.replace("C 1 2 3", "C -1 1 2 3 H"),
                    base.replace("C 1 2 3", "C 1 2 NaN"), base.replace("0 1", "0 0"),
                    base + base, base.replace("Units=Ang", "Units=Unknown"),
                    base.replace("Units=Ang", "Units=Ang Units=AU")]:
            with self.assertRaises(ValueError): writer.observe(bad)

    def test_signed_zero_cannot_hide_as_a_bitwise_match(self):
        self.assertFalse(writer.same_numbers([[0.,0.,0.]], [[-0.,0.,0.]]))
        self.assertFalse(writer.same_numbers([[1.,2.,3.]], [[1.,2.,3.],[4.,5.,6.]]))


if __name__ == "__main__":
    unittest.main()
