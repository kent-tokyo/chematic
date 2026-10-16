"""Binary stream contracts motivated by RDKit #9601 / #9615."""
import gc
import gzip
import io
import weakref

import pytest

import chematic
from chematic import rdkit_compat as Chem


@pytest.fixture(params=[chematic, Chem], ids=["native", "compat"])
def api(request):
    return request.param


def record(smiles="CCO", title="ethanol", value="7.2"):
    block = chematic.from_smiles(smiles).to_mol_block()
    block = title + "\n" + block.split("\n", 1)[1]
    return (block + f"> <Activity>\n{value}\n\n$$$$\n").encode("utf-8")


def atom_count(mol):
    return mol.atom_count if isinstance(mol, chematic.Mol) else mol.GetNumAtoms()


@pytest.mark.parametrize("source", [42, None, object(), b"not a stream"])
def test_missing_read_rejected_at_construction(api, source):
    with pytest.raises(ValueError, match="read"):
        api.ForwardSDMolSupplier(source)


@pytest.mark.parametrize("read", [None, 42, "read"])
def test_noncallable_read_rejected(api, read):
    source = type("Source", (), {"read": read})()
    with pytest.raises(ValueError, match="callable"):
        api.ForwardSDMolSupplier(source)


def test_getattr_error_preserved(api):
    error = RuntimeError("read lookup failed")

    class Source:
        @property
        def read(self):
            raise error

    with pytest.raises(RuntimeError) as caught:
        api.ForwardSDMolSupplier(Source())
    assert caught.value is error


def test_properties_unicode_and_unterminated_record(api):
    source = io.BytesIO(record(title="エタノール", value="日本語\nsecond line")
                        + record("CN", "amine").removesuffix(b"$$$$\n"))
    supplier = api.ForwardSDMolSupplier(source)
    assert iter(supplier) is supplier
    molecules = list(supplier)
    assert len(molecules) == 2
    assert [atom_count(m) for m in molecules] == [3, 2]
    assert molecules[0].GetProp("_Name") == "エタノール"
    assert molecules[0].GetProp("Activity") == "日本語\nsecond line"
    assert molecules[1].GetProp("_Name") == "amine"
    assert list(supplier) == []
    assert not source.closed


@pytest.mark.parametrize("chunk", [1, 3, 17, 8192])
def test_short_reads_crlf_and_no_eager_consumption(api, chunk):
    data = record(title="日本語").replace(b"\n", b"\r\n") * 40

    class Source:
        def __init__(self):
            self.position = 0
            self.requests = []

        def read(self, size):
            assert 0 < size <= 8192
            self.requests.append(size)
            result = data[self.position:self.position + min(size, chunk)]
            self.position += len(result)
            return result

    source = Source()
    supplier = api.ForwardSDMolSupplier(source)
    assert source.requests == []
    assert next(supplier).GetProp("_Name") == "日本語"
    assert source.position < len(data)
    assert len(list(supplier)) == 39
    requests = len(source.requests)
    assert list(supplier) == []
    assert len(source.requests) == requests


def test_gzip_stream_and_ownership(api):
    source = gzip.GzipFile(fileobj=io.BytesIO(gzip.compress(record() * 2)))
    supplier = api.ForwardSDMolSupplier(source)
    assert len(list(supplier)) == 2
    del supplier
    gc.collect()
    assert not source.closed
    source.close()


def test_source_retained_and_released(api):
    source = io.BytesIO(record())
    reference = weakref.ref(source)
    supplier = api.ForwardSDMolSupplier(source)
    del source
    gc.collect()
    assert reference() is not None
    mol = next(supplier)
    del supplier
    gc.collect()
    assert reference() is None
    assert mol.GetProp("_Name") == "ethanol"


def test_source_cycle_collected(api):
    source = io.BytesIO(record())
    reference = weakref.ref(source)
    source.supplier = api.ForwardSDMolSupplier(source)
    del source
    gc.collect()
    assert reference() is None


def test_static_read_still_retains_source(api):
    class Source:
        read = staticmethod(io.BytesIO(record()).read)

    source = Source()
    reference = weakref.ref(source)
    supplier = api.ForwardSDMolSupplier(source)
    del source
    gc.collect()
    assert reference() is not None
    assert atom_count(next(supplier)) == 3
    del supplier
    gc.collect()
    assert reference() is None


@pytest.mark.parametrize("strict", [True, False])
def test_parse_errors_keep_record_position(api, strict):
    source = io.BytesIO(record() + b"malformed\n$$$$\n" + record("CN"))
    supplier = api.ForwardSDMolSupplier(source, strictParsing=strict)
    assert atom_count(next(supplier)) == 3
    if strict:
        with pytest.raises(ValueError, match="malformed SDF"):
            next(supplier)
    else:
        assert next(supplier) is None
    assert atom_count(next(supplier)) == 2
    assert list(supplier) == []


@pytest.mark.parametrize("result", ["", "text", None, bytearray(b""), memoryview(b"")])
def test_invalid_read_result_is_terminal(api, result):
    class Source:
        def read(self, size):
            return result

    supplier = api.ForwardSDMolSupplier(Source(), strictParsing=False)
    with pytest.raises(TypeError, match="bytes"):
        next(supplier)
    assert list(supplier) == []


def test_oversized_read_rejected(api):
    class Source:
        def read(self, size):
            return b"x" * (size + 1)

    supplier = api.ForwardSDMolSupplier(Source())
    with pytest.raises(ValueError, match="more than size"):
        next(supplier)
    assert list(supplier) == []


@pytest.mark.parametrize("strict", [True, False])
def test_read_exception_preserved_after_valid_record(api, strict):
    error = OSError("custom read failure")

    class Source:
        first = True

        def read(self, size):
            if self.first:
                self.first = False
                return record()
            raise error

    supplier = api.ForwardSDMolSupplier(Source(), strictParsing=strict)
    assert atom_count(next(supplier)) == 3
    with pytest.raises(OSError) as caught:
        next(supplier)
    assert caught.value is error
    assert list(supplier) == []


def test_invalid_utf8_is_not_silently_skipped(api):
    supplier = api.ForwardSDMolSupplier(io.BytesIO(record(title="x").replace(b"x\n", b"\xff\n", 1)),
                                      strictParsing=False)
    with pytest.raises(OSError):
        next(supplier)


def test_v3000_matches_native_path_reader(api, tmp_path):
    path = tmp_path / "stereo.sdf"
    with chematic.SDWriter(path) as writer:
        writer.SetForceV3000(True)
        mol = chematic.from_smiles("F[C@H](Cl)Br")
        mol.SetProp("_Name", "stereo")
        writer.write(mol)
    expected = list(chematic.SDMolSupplier(path))
    with path.open("rb") as source:
        actual = list(api.ForwardSDMolSupplier(source))
    assert len(actual) == len(expected) == 1
    assert actual[0].GetProp("_Name") == expected[0].GetProp("_Name")
    if api is Chem:
        assert Chem.MolToSmiles(actual[0]) == expected[0].smiles
    else:
        assert actual[0].smiles == expected[0].smiles


def test_v2000_wedge_perception_is_preserved(api):
    # RDKit SDWriter-style bond wedge, with the stereocentre as bond start.
    block = """stereo
     RDKit          2D

  4  3  0  0  0  0  0  0  0  0999 V2000
    1.2990    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -1.5000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  2  1  1  1
  2  3  1  0
  2  4  1  0
M  END
"""
    expected = chematic.from_mol_block(block).smiles
    assert "@" in expected
    mol = next(api.ForwardSDMolSupplier(io.BytesIO((block + "$$$$\n").encode())))
    assert (mol.smiles if api is chematic else api.MolToSmiles(mol)) == expected


def test_empty_and_blank_stream(api):
    assert list(api.ForwardSDMolSupplier(io.BytesIO())) == []
    assert list(api.ForwardSDMolSupplier(io.BytesIO(b"\n$$$$\n\n"))) == []


def test_rdkit_record_identifier_keeps_sd_property(api, tmp_path):
    data = record().replace(b"> <Activity>", b">  <Activity>  (1) ")
    mol = next(api.ForwardSDMolSupplier(io.BytesIO(data)))
    assert mol.GetProp("Activity") == "7.2"
    path = tmp_path / "rdkit.sdf"
    path.write_bytes(data)
    assert next(chematic.SDMolSupplier(path)).GetProp("Activity") == "7.2"
    assert next(chematic.iter_sdf_str(data.decode())).get("Activity") == "7.2"


def test_resource_limit_is_not_lenient_parse_failure(api):
    supplier = api.ForwardSDMolSupplier(
        io.BytesIO(b"malformed\n$$$$\n" * 100_001), strictParsing=False,
    )
    for _ in range(100_000):
        assert next(supplier) is None
    with pytest.raises(chematic.ChematicInputError) as caught:
        next(supplier)
    assert caught.value.category == "resource_limit"
    assert caught.value.code == "sdf_resource_limit"
    assert list(supplier) == []
