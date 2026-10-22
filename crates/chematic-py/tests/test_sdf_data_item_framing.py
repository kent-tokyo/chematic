"""Original SDF field values, through native and compatibility stream readers."""
import io
import json
from pathlib import Path

import pytest
import chematic
from chematic import rdkit_compat as Chem

PACKET = json.loads((Path(__file__).parents[3]/'validation/sdf_data_item_framing/cases.json').read_text(encoding="utf-8"))


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c:c['id'])
@pytest.mark.parametrize('fmt', ['v2000','v3000'])
@pytest.mark.parametrize('ending', ['\n','\r\n'], ids=['lf','crlf'])
@pytest.mark.parametrize('api', [chematic,Chem], ids=['native','compat'])
def test_data_item_boundaries_and_exact_values(case, fmt, ending, api):
    text = (PACKET['ctabs'][fmt]+case['data']+'$$$$\n').replace('\n',ending)
    records = list(api.ForwardSDMolSupplier(io.BytesIO(text.encode())))
    assert len(records) == 1 and records[0] is not None
    mol = records[0]
    assert {k:mol.GetProp(k) for k in mol.GetPropNames() if not k.startswith('_')} == case['expected_properties']


@pytest.mark.parametrize('case', PACKET['cases'], ids=lambda c:c['id'])
@pytest.mark.parametrize('v3000', [False,True], ids=['v2000','v3000'])
def test_property_writer_readback_preserves_literal_header_lines_and_whitespace(case, v3000, tmp_path):
    source = PACKET['ctabs']['v2000']+case['data']+'$$$$\n'
    mol = next(chematic.ForwardSDMolSupplier(io.BytesIO(source.encode())))
    path = tmp_path/'result.sdf'
    writer = chematic.SDWriter(path,compute2d=False)
    writer.SetForceV3000(v3000)
    writer.write(mol)
    writer.close()
    again = next(chematic.ForwardSDMolSupplier(io.BytesIO(path.read_bytes())))
    assert {k:again.GetProp(k) for k in again.GetPropNames() if not k.startswith('_')} == case['expected_properties']
