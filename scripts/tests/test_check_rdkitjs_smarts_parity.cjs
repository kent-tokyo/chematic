const assert = require('node:assert/strict');
const test = require('node:test');
const { chematicAtomSets } = require('../check_rdkitjs_smarts_parity.cjs');

test('CheMatic match order and duplicate embeddings do not affect atom sets', () => {
  assert.deepEqual(chematicAtomSets('[]'), []);
  assert.deepEqual(chematicAtomSets('[[2,1],[1,2],[4]]'), ['[1,2]', '[4]']);
  assert.throws(() => chematicAtomSets('{"error":"bad"}'), /unexpected CheMatic/);
});
