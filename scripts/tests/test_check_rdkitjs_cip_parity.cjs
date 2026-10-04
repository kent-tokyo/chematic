const assert = require('node:assert/strict');
const test = require('node:test');
const { parseV2000, indexCorresponds, differingLabels, chematicLabels, rdkitLabels } = require('../check_rdkitjs_cip_parity.cjs');

test('V2000 index proof ignores coordinates and bond orientation but not element order', () => {
  const text = [
    '', '  test', '', '  2  1  0  0  0  0  0  0  0  0999 V2000',
    '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0   0  0',
    '    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0   0  0',
    '  1  2  1  0  0  0  0', 'M  END',
  ].join('\n');
  const parsed = parseV2000(text);
  assert.deepEqual(parsed.atoms, ['C', 'O']);
  assert.deepEqual(parsed.bonds, [{ a: 0, b: 1, order: 1 }]);
  assert.equal(indexCorresponds(parsed, { atoms: ['C', 'O'], bonds: [{ a: 1, b: 0, order: 2 }] }), true);
  assert.equal(indexCorresponds(parsed, { atoms: ['O', 'C'], bonds: parsed.bonds }), false);
});

test('unspecified RDKit stereo tags are not assigned CIP labels', () => {
  const rdkit = rdkitLabels({ CIP_atoms: [[1, '(S)'], [2, '(?)']], CIP_bonds: [[3, 4, '(E)']] });
  assert.deepEqual(rdkit, { atoms: { 1: 'S' }, bonds: { '3-4': 'E' } });
  const chematic = chematicLabels([
    { atomIdx: 1, cipCode: 'S' }, { atomIdx: 3, cipCode: 'E' },
  ], { bonds: [{ a: 3, b: 4, order: 2 }] });
  assert.deepEqual(differingLabels(chematic, rdkit), []);
  assert.throws(() => chematicLabels([{ atomIdx: 3, cipCode: 'E' }], { bonds: [] }), /ambiguous E\/Z/);
});
