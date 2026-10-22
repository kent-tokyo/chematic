"""Original opaque XML retention; chemical interpretation is deliberately separate."""
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).parents[3] / 'validation/cml_source_preservation/cases.json').read_text(encoding="utf-8"))['cases']


@pytest.mark.parametrize('case', PACKET, ids=lambda c: c['id'])
@pytest.mark.parametrize('policy', ['strict', 'allow_loss', 'preserve_source'])
def test_single_document_retains_exact_source_and_opaque_scope(case, policy):
    source = case['source']
    report = chematic.convert_with_report(source, ' .CML ', '.cml', policy=policy)
    assert report == json.loads(chematic.convert_with_report_json(source, 'cml', 'cml', policy))
    assert report['output'].encode('utf8') == source.encode('utf8')
    assert report['status'] == 'completed'
    assert report['verification_scope'] == 'exact_source_bytes'
    assert report['candidate_graph_verified'] is False
    assert report['candidate_atom_mapping'] == report['candidate_bond_mapping'] == []
    assert {o['field']: o['status'] for o in report['observations']} == {
        'source_text': 'preserved', 'source_semantics': 'opaque_retained',
    }


@pytest.mark.parametrize('source', [
    '<molecule/>',
    '<molecule xmlns="http://www.xml-cml.org/schema"><molecule/></molecule>',
    '<molecule xmlns="http://www.xml-cml.org/schema">&#0;</molecule>',
    '<!DOCTYPE molecule><molecule xmlns="http://www.xml-cml.org/schema"/>',
    '<?xml version="1.1"?><molecule xmlns="http://www.xml-cml.org/schema"/>',
    '<molecule xmlns="http://www.xml-cml.org/schema"><atom></molecule>',
])
def test_invalid_or_unsupported_framing_is_an_error(source):
    with pytest.raises(ValueError):
        chematic.convert_with_report(source, 'cml', 'cml')


def test_cross_format_source_retention_is_unsupported():
    with pytest.raises(ValueError):
        chematic.convert_with_report(PACKET[0]['source'], 'cml', 'smiles', policy='preserve_source')
