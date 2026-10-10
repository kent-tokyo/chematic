# CDK and Indigo comparison

This gate compares chematic with the official stable CDK 2.13 and Indigo
1.46.0 distributions. It is intentionally narrower than either competitor's
total feature set. A row is rankable only when chematic and the competitor
both complete equivalent work and match the pinned reference on every corpus
row.

The adapters use public APIs only. No CDK or Indigo implementation code is
copied into chematic. See `contract.json` for versions, licenses, operations,
and exit criteria.

## Build and run

Build the Java adapter with JDK 17 or newer and Maven:

```bash
cd validation/cdk_indigo_comparison/cdk-adapter
mvn -q package
```

Install `chematic`, `rdkit==2026.3.6`, and `epam.indigo==1.46.0` into the
Python used below. From the repository root:

```bash
CORPUS=validation/cdk_indigo_comparison/corpus-v1.jsonl
MANIFEST=validation/cdk_indigo_comparison/corpus-v1.manifest.json
PY=python3
JAVA=java
ADAPTER=validation/cdk_indigo_comparison/python_adapter.py
CDK=validation/cdk_indigo_comparison/cdk-adapter/target/cdk-comparison-adapter.jar

$PY $ADAPTER --engine rdkit --corpus $CORPUS > /tmp/rdkit.jsonl
$PY $ADAPTER --engine chematic --corpus $CORPUS > /tmp/chematic.jsonl
$PY $ADAPTER --engine indigo --corpus $CORPUS > /tmp/indigo.jsonl
$JAVA -jar $CDK --corpus $CORPUS > /tmp/cdk.jsonl

$PY validation/cosmolkit_comparison/scorecard.py \
  --result rdkit=/tmp/rdkit.jsonl --result chematic=/tmp/chematic.jsonl \
  --result indigo=/tmp/indigo.jsonl --result cdk=/tmp/cdk.jsonl \
  --reference rdkit --corpus $CORPUS --manifest $MANIFEST \
  --output /tmp/cdk-indigo-scorecard.json
```

The persistent-process benchmark keeps JVM and native-library startup outside
the hot lanes, rotates engine order between blocks, and records every raw
sample. CheMatic's prepared SMARTS lane uses the public
`SmartsQuery.matches_many()` batch API; the report records this in engine
metadata. Pass the same scorecard so mismatched operations are never ranked:

```bash
$PY validation/cdk_indigo_comparison/benchmark.py \
  --corpus $CORPUS --scorecard /tmp/cdk-indigo-scorecard.json \
  --adapter "rdkit=$PY $ADAPTER --engine rdkit" \
  --adapter "chematic=$PY $ADAPTER --engine chematic" \
  --adapter "indigo=$PY $ADAPTER --engine indigo" \
  --adapter "cdk=$JAVA -jar $CDK" \
  --blocks 21 --output /tmp/cdk-indigo-benchmark.json
```

The 24-row corpus checks adapter wiring and obvious chemical-class regressions
only. Accuracy or performance claims require a larger frozen stratified corpus,
exact artifact hashes, a second host, and a separate peak-RSS lane.

## Current source-candidate baseline

On the 24-row smoke corpus, the current source candidate matches the pinned
RDKit reference on all 15 measured operations. The MOL V2000 round trip now
retains `M  ISO` isotopes and coordinate-derived E/Z; those defects were found
by this comparison and are covered by Rust regressions.

The paired macOS arm64 smoke benchmark records source-candidate median and 95%
interval wins over CDK 2.13 and Indigo 1.46.0 for every rankable lane. HBD is
2.3x faster than CDK in the parse-inclusive lane and 6.6x in the prepared
lane; the direct implementation remains 10,000/10,000 exact against RDKit.
One high-repetition `[#8]` lane wins 20/21 blocks, so the contract's stricter
every-block exit remains open. These are development signals, not a broad or
published-package superiority claim. See the dated benchmark record.

`rings` is diagnostic only. RDKit/CheMatic report the pinned RDKit ring model,
whereas CDK and Indigo expose SSSR counts; cubane therefore differs without
showing that either graph is chemically wrong.
