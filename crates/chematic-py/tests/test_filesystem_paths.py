"""Original path-protocol controls for the filesystem part of RDKit #9601."""
import errno
import os
from pathlib import Path

import chematic
import pytest

READERS = ['iter_sdf', 'iter_sdf_checked', 'iter_sdf_batched', 'SdfDataset',
           'SdfDataset.open', 'SDMolSupplier', 'SmilesMolSupplier', 'TDTMolSupplier',
           'iter_xyz_batched', 'iter_extxyz_batched']
WRITERS = ['SDWriter', 'SmilesWriter', 'TDTWriter', 'Report.save', 'report.output']
OPERATIONS = READERS + WRITERS
EXPECTED_SMILES = chematic.from_smiles('CCO').smiles


class PathValue:
    def __init__(self, value):
        self.value = value
        self.calls = 0

    def __fspath__(self):
        self.calls += 1
        return self.value


class StringSubclass(str):
    pass


class BytesSubclass(bytes):
    pass


class StringOnly:
    def __str__(self):
        raise AssertionError('arbitrary __str__ must not be called')


class BrokenPath:
    def __fspath__(self):
        raise RuntimeError('path callback failed')


@pytest.fixture
def files(tmp_path):
    folder = tmp_path / '日本語 paths'
    folder.mkdir()
    mol = chematic.from_smiles('CCO')
    sdf = mol.to_mol_block() + '\n$$$$\n'
    result = {}
    for name in OPERATIONS:
        path = folder / (name + '.input')
        if name in READERS:
            if name == 'SmilesMolSupplier':
                data = 'SMILES Name\nCCO original\n'
            elif name == 'TDTMolSupplier':
                data = '$SMI<CCO>\nNAME<original>\n|\n'
            elif 'xyz' in name:
                data = '1\noriginal\nC 0 0 0\n'
            else:
                data = sdf
            path.write_text(data)
        result[name] = path
    return result


def invoke(name, path):
    if name == 'SdfDataset.open':
        return chematic.SdfDataset.open(path)
    if name == 'Report.save':
        report = chematic.report([chematic.from_smiles('CCO')])
        report.save(path)
        return report
    if name == 'report.output':
        return chematic.report([chematic.from_smiles('CCO')], output=path)
    if name in ['SDWriter', 'SmilesWriter', 'TDTWriter']:
        writer = getattr(chematic, name)(path)
        writer.write(chematic.from_smiles('CCO'))
        writer.close()
        return None
    return getattr(chematic, name)(path)


def check_result(name, value, path):
    if name.startswith('SdfDataset'):
        assert len(value) == 1
        assert value[0].mol.smiles == EXPECTED_SMILES
    elif name == 'iter_sdf_checked':
        rows = list(value)
        assert len(rows) == 1 and rows[0][0] == 0 and rows[0][2] is None
        assert rows[0][1].mol.smiles == EXPECTED_SMILES
    elif name == 'iter_sdf_batched':
        assert [[record.mol.smiles for record in batch] for batch in value] == [[EXPECTED_SMILES]]
    elif 'xyz' in name:
        rows = list(value)
        assert len(rows) == 1 and len(rows[0]) == 1
        assert rows[0][0]['coords'] == [[0.0, 0.0, 0.0]]
    elif name == 'iter_sdf':
        assert [record.mol.smiles for record in value] == [EXPECTED_SMILES]
    elif name in READERS:
        assert [mol.smiles for mol in value] == [EXPECTED_SMILES]
    elif name == 'SDWriter':
        assert [record.mol.smiles for record in chematic.iter_sdf(os.fsdecode(path))] == [EXPECTED_SMILES]
    elif name == 'SmilesWriter':
        assert [mol.smiles for mol in chematic.SmilesMolSupplier(os.fsdecode(path))] == [EXPECTED_SMILES]
    elif name == 'TDTWriter':
        assert [mol.smiles for mol in chematic.TDTMolSupplier(os.fsdecode(path))] == [EXPECTED_SMILES]
    else:
        with open(path, encoding='utf8') as stream:
            assert stream.read() == str(value)


FORMS = [str, os.fsencode, Path, lambda p: PathValue(str(p)),
         lambda p: PathValue(os.fsencode(p)), StringSubclass,
         lambda p: BytesSubclass(os.fsencode(p))]


@pytest.mark.parametrize('name', OPERATIONS)
@pytest.mark.parametrize('form', FORMS, ids=['str', 'bytes', 'Path', 'PathLike-str',
                                           'PathLike-bytes', 'str-subclass', 'bytes-subclass'])
def test_native_path_protocol_and_file_contents(files, name, form):
    path = form(files[name])
    value = invoke(name, path)
    if isinstance(path, PathValue):
        assert path.calls == 1
    check_result(name, value, files[name])


@pytest.mark.parametrize('name', OPERATIONS)
@pytest.mark.parametrize('bad', [None, 3, True, 1.5, [], {}, bytearray(b'x'),
                                StringOnly(), PathValue(3), PathValue(Path('x'))])
def test_invalid_path_types_do_not_touch_files(files, name, bad):
    # report(output=None) intentionally means no file output.
    if name == 'report.output' and bad is None:
        assert chematic.report([], output=None) is not None
        return
    files[name].write_bytes(b'sentinel')
    with pytest.raises(TypeError):
        invoke(name, bad)
    assert files[name].read_bytes() == b'sentinel'


@pytest.mark.parametrize('name', OPERATIONS)
def test_missing_path_parent_is_an_io_error(tmp_path, name):
    with pytest.raises(OSError):
        invoke(name, tmp_path / 'absent' / 'file')


@pytest.mark.parametrize('name', OPERATIONS)
def test_callback_exception_is_preserved(name):
    with pytest.raises(RuntimeError, match='path callback failed'):
        invoke(name, BrokenPath())


@pytest.mark.parametrize('name', OPERATIONS)
def test_nul_paths_are_io_errors(name):
    with pytest.raises(OSError):
        invoke(name, b'bad\x00path')


@pytest.mark.skipif(os.name != 'posix', reason='raw POSIX filename bytes')
@pytest.mark.parametrize('name', OPERATIONS)
@pytest.mark.parametrize('as_surrogate', [False, True], ids=['bytes', 'surrogate-str'])
def test_non_utf8_filename_bytes_are_lossless(files, name, as_surrogate):
    raw = os.fsencode(files[name]) + b'\xff'
    # APFS on this macOS host rejects such names even through Python's open.
    # Only skip that explicit OS limitation; other I/O failures still fail.
    try:
        with open(raw, 'wb'):
            pass
    except OSError as error:
        if error.errno != errno.EILSEQ:
            raise
        path = os.fsdecode(raw) if as_surrogate else PathValue(raw)
        with pytest.raises(OSError):
            invoke(name, path)
        pytest.skip('filesystem rejects raw non-UTF-8 filenames (EILSEQ)')
    if name in READERS:
        os.rename(files[name], raw)
    path = os.fsdecode(raw) if as_surrogate else PathValue(raw)
    value = invoke(name, path)
    if isinstance(path, PathValue):
        assert path.calls == 1
    check_result(name, value, raw)


def test_file_api_stubs_cover_all_path_parameters():
    import ast
    stub = Path(__file__).resolve().parents[1] / 'python/chematic/__init__.pyi'
    tree = ast.parse(stub.read_text(encoding='utf-8'))
    declarations = {node.name: node for node in tree.body if hasattr(node, 'name')}
    for name in OPERATIONS:
        argument = 'output' if name == 'report.output' else 'path'
        if name == 'report.output':
            function = declarations['report']
        elif '.' in name:
            owner, method = name.split('.')
            function = next(n for n in declarations[owner].body if isinstance(n, ast.FunctionDef) and n.name == method)
        elif isinstance(declarations[name], ast.ClassDef):
            function = next(n for n in declarations[name].body if isinstance(n, ast.FunctionDef) and n.name == '__init__')
        else:
            function = declarations[name]
        parameter = next(arg for arg in function.args.args if arg.arg == argument)
        assert '_FileSystemPath' in ast.unparse(parameter.annotation), name
