const assert = require('node:assert/strict');
const test = require('node:test');
const { atomSets, options } = require('../compare_rdkitjs_release_chemistry.cjs');

test('normalizes empty MinimalLib result and order-independent atom sets', () => {
  assert.deepEqual(atomSets('{}'), []);
  assert.deepEqual(
    atomSets('[{"atoms":[2,1],"bonds":[0]},{"atoms":[1,2],"bonds":[0]}]'),
    ['[1,2]'],
  );
  assert.throws(() => atomSets('{"error":"unexpected"}'), /unexpected substructure output/);
});

test('requires explicit old/new packages and output paths', () => {
  assert.throws(() => options([]), /missing --old-package/);
  assert.deepEqual(options([
    '--old-package', 'old', '--new-package', 'new', '--corpus', 'rows',
    '--queries', 'queries', '--summary', 'summary', '--rows-output', 'rows.gz',
  ])['new-package'], 'new');
});
