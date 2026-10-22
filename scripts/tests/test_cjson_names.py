import json
import unittest

from scripts.validation.check_cjson_names import inspect, packet


class NamePacketTests(unittest.TestCase):
    def test_legacy_location_is_not_counted_as_name_loss(self):
        value = {"atoms": {"elements": {"number": [2]}, "coords": {"3d": [0, 0, 0]}},
                 "properties": {"name": "title"}}
        fields, name = inspect(json.dumps(value), "cjson", {"name": "title"})
        self.assertEqual(name, "title")
        self.assertTrue(all(fields.values()))
        value["name"] = "different"
        with self.assertRaises(ValueError):
            inspect(json.dumps(value), "cjson", {"name": "title"})

    def test_changed_title_graph_or_coordinate_fails(self):
        value = {"atoms": {"elements": {"number": [2]}, "coords": {"3d": [0, 0, 0]}},
                 "name": "title"}
        fields, _ = inspect(json.dumps(value), "cjson", {"name": "other"})
        self.assertFalse(fields["name"])
        value["atoms"]["coords"]["3d"][0] = 1
        self.assertFalse(inspect(json.dumps(value), "cjson", {"name": "title"})[0]["explicit_graph_and_coordinates"])
        value["atoms"]["coords"]["3d"][0] = 0
        value["atoms"]["elements"]["number"] = [1]
        self.assertFalse(inspect(json.dumps(value), "cjson", {"name": "title"})[0]["explicit_graph_and_coordinates"])

    def test_refusal_cases_are_not_positive_coverage(self):
        cases = packet()
        self.assertEqual(len(cases), 11)
        self.assertEqual(sum(c["expected"]["refuse"] for c in cases), 4)
        self.assertEqual(len({c["id"] for c in cases}), len(cases))


if __name__ == "__main__":
    unittest.main()
