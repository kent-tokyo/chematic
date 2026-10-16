"""Syntax and real archive failures must stop the gate without importing code."""
from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
import warnings
import zipfile
from pathlib import Path
from scripts import check_python_stubs as gate


class StubGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.package = self.root / 'chematic'
        self.package.mkdir()
        self.valid = b'from typing import Sequence\ndef from_smiles(smiles: str) -> object: ...\n'
        (self.package / '__init__.pyi').write_bytes(self.valid)
        (self.package / 'py.typed').write_bytes(b'')

    def wheel(self, members=None, name='chematic-test.whl'):
        path = self.root / name
        members = [('chematic/__init__.pyi', self.valid), ('chematic/py.typed', b'')] if members is None else members
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(path, 'w') as z:
                for key, data in members:
                    z.writestr(key, data)
        return path

    def test_source_and_wheel_accept_valid_bytes_and_nested_stubs(self):
        sub = self.package / 'nested'
        sub.mkdir()
        (sub / 'types.pyi').write_bytes(b'class Atom: ...\n')
        p = self.wheel([('chematic/__init__.pyi',self.valid),('chematic/nested/types.pyi',b'class Atom: ...\n'),('chematic/py.typed',b''),('chematic/native.so',b'not executed')])
        r = gate.audit(self.package,[str(p)])
        self.assertEqual(r['status'],'passed')
        self.assertEqual(len(r['source']['stubs']),2)
        self.assertEqual(len(r['wheels'][0]['stubs']),2)

    def test_all_reported_syntax_shapes_are_rejected_with_locations(self):
        samples = [b'None: int\n', b'class Example:\n    value: unsigned int\n',
                   b'def f(default=0, required): ...\n', b'def f(x, x): ...\n']
        for data in samples:
            with self.subTest(data=data):
                issue=gate.check_syntax(data,'original-negative.pyi')
                self.assertIsNotNone(issue)
                self.assertEqual(issue['kind'],'syntax')
                self.assertGreater(issue['line'],0)

    def test_parse_and_compile_failures_are_distinguished(self):
        issue=gate.check_syntax(b'None: int\n','parse.pyi')
        self.assertEqual(issue['phase'],'parse')
        for data in [b'def f(x, x): ...\n',b'import typing\nfrom __future__ import annotations\n']:
            issue=gate.check_syntax(data,'compile.pyi')
            self.assertEqual(issue['phase'],'compile')

    def test_invalid_encoding_and_null_bytes_are_rejected(self):
        for data in [b'bad: \xff\n', b'x: int\0\n']:
            with self.subTest(data=data):
                self.assertIsNotNone(gate.check_syntax(data,'encoding.pyi'))

    def test_python_39_grammar_rejects_newer_statement_syntax(self):
        self.assertIsNotNone(gate.check_syntax(b'match x:\n    case 1: ...\n','newer.pyi'))
        self.assertIsNone(gate.check_syntax(b'from typing import Optional\ndef f(x: Optional[int]) -> list[str]: ...\n','portable.pyi'))

    def test_does_not_execute_stub_imports_annotations_or_bodies(self):
        marker=self.root/'executed'
        data=f'from pathlib import Path\nx: Path({str(marker)!r}).write_text("bad")\nraise RuntimeError("never execute")\n'.encode()
        (self.package/'__init__.pyi').write_bytes(data)
        p=self.wheel([('chematic/__init__.pyi',data),('chematic/py.typed',b'')])
        self.assertEqual(gate.audit(self.package,[str(p)])['status'],'passed')
        self.assertFalse(marker.exists())

    def test_missing_root_stub_and_marker_fail(self):
        (self.package/'__init__.pyi').unlink()
        (self.package/'py.typed').unlink()
        r=gate.audit(self.package,[])
        self.assertEqual(r['status'],'failed')
        self.assertEqual(len(r['errors']),2)

    def test_missing_package_and_unmatched_glob_fail(self):
        r=gate.audit(self.root/'absent',['missing-*.whl'])
        self.assertEqual(r['status'],'failed')
        self.assertTrue(any('matched no files' in e['message'] for e in r['errors']))

    def test_missing_wheel_stub_and_marker_fail(self):
        r=gate.audit(self.package,[str(self.wheel([]))])
        self.assertEqual(r['status'],'failed')
        self.assertEqual({e['kind'] for e in r['errors']},{'missing'})

    def test_stale_but_valid_stub_and_marker_bytes_fail(self):
        p=self.wheel([('chematic/__init__.pyi',b'x: int\n'),('chematic/py.typed',b'partial\n')])
        r=gate.audit(self.package,[str(p)])
        self.assertEqual(r['status'],'failed')
        self.assertEqual(len(r['errors']),2)
        self.assertTrue(all(e['kind']=='mismatch' for e in r['errors']))

    def test_unexpected_packaged_stub_and_other_package_syntax_are_checked(self):
        p=self.wheel([('chematic/__init__.pyi',self.valid),('chematic/py.typed',b''),('chematic/extra.pyi',b'x: int\n'),('other/broken.pyi',b'None: int\n')])
        r=gate.audit(self.package,[str(p)])
        self.assertEqual({e['kind'] for e in r['errors']},{'unexpected','syntax'})

    def test_duplicate_zip_members_never_select_one_silently(self):
        p=self.wheel([('chematic/__init__.pyi',self.valid),('chematic/__init__.pyi',b'None: int\n'),('chematic/py.typed',b'')])
        r=gate.audit(self.package,[str(p)])
        self.assertEqual(r['status'],'failed')
        self.assertIn('duplicate',r['errors'][0]['message'])

    def test_corrupt_archive_and_directory_in_wheel_glob_fail(self):
        p=self.root/'broken.whl';p.write_bytes(b'not zip')
        d=self.root/'directory.whl';d.mkdir()
        for candidate in [p,d]:
            with self.subTest(candidate=candidate):
                self.assertEqual(gate.audit(self.package,[str(candidate)])['status'],'failed')

    def test_entry_and_total_byte_budgets_reject_before_reading(self):
        old_entry,old_total=gate.MAX_STUB_BYTES,gate.MAX_TOTAL_BYTES
        self.addCleanup(setattr,gate,'MAX_STUB_BYTES',old_entry)
        self.addCleanup(setattr,gate,'MAX_TOTAL_BYTES',old_total)
        p=self.wheel()
        gate.MAX_STUB_BYTES=8
        self.assertEqual(gate.audit(self.package,[str(p)])['status'],'failed')
        gate.MAX_STUB_BYTES=old_entry;gate.MAX_TOTAL_BYTES=8
        self.assertEqual(gate.audit(self.package,[str(p)])['status'],'failed')

    def test_multiple_files_errors_are_all_reported(self):
        (self.package/'__init__.pyi').write_bytes(b'None: int\n')
        (self.package/'other.pyi').write_bytes(b'def f(x, x): ...\n')
        r=gate.audit(self.package,[])
        self.assertEqual(len(r['errors']),2)

    def test_cli_json_and_exit_status_do_not_silently_accept_failures(self):
        p=self.wheel();report=self.root/'report.json'
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            good=gate.main(['--source-package',str(self.package),'--wheel',str(p),'--report',str(report)])
            self.assertEqual(good,0)
            self.assertEqual(json.loads(report.read_text())['status'],'passed')
            bad=gate.main(['--source-package',str(self.package),'--wheel',str(self.root/'missing.whl'),'--report',str(report)])
        self.assertEqual(bad,1)
        self.assertEqual(json.loads(report.read_text())['status'],'failed')

    def test_quoted_glob_checks_all_wheels_and_deduplicates_paths(self):
        one=self.wheel(name='one.whl');self.wheel(name='two.whl')
        r=gate.audit(self.package,[str(self.root/'*.whl'),str(one)])
        self.assertEqual(r['status'],'passed')
        self.assertEqual(len(r['wheels']),2)


if __name__=='__main__':
    unittest.main()
