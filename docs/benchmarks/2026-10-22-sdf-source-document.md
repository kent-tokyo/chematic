# Complete SDF retention and single-record conversion gates

The Python/WASM common bridges previously normalized `sdf` to `mol`. SDF → SDF
therefore parsed one molecule and regenerated a graph-only MOL block. A two-
record carbon/oxygen source returned only carbon, without a record delimiter or
SD fields. This loses collection membership as well as source information.

The common bridges now distinguish SDF documents from MOL blocks:

```python
retained = chematic.convert_format(source, "sdf", "sdf")
assert retained == source
assert chematic.roundtrip_sdf_source(source.encode("utf-8")) == source
```

```javascript
const retained = wasm.convert_common_format(source, ".SDF", "sdf");
if (retained !== source) throw new Error("document changed");
// wasm.roundtrip_sdf_source(source) uses the same native implementation.
```

Every record is validated before any source string is returned. One failed
record or exceeded limit returns no output, rather than a valid prefix. Original
coordinate tokens, version/dimension declarations, SD fields (including duplicate
names), comments, titles, record order and line endings survive without rewriting.
This is immutable source retention; metadata, query/polymer and stereochemical
meaning is not independently interpreted or verified. It is not an editable
multi-molecule model or general CTfile validator.

For other common-bridge destinations, `mol` and `sdf` input must contain exactly
one physical SDF/MOL record. Extra records fail instead of selecting the first.
SD fields and delimiter are intentionally omitted, and the source version selects
the existing V2000/V3000 graph reader. Python uses physical LF lines, including
CRLF, rather than treating a Unicode line separator inside a UTF-8 header as an
extra record/header line. Other graph/coordinate/stereo limits remain applicable.
Use the existing SDF iterators for batch conversion with record-level diagnostics.
Non-SDF → SDF now writes one MOL block followed by a `$$$$` delimiter.

Rust exports `roundtrip_sdf_source[_with_limits]`,
`validate_sdf_source_with_limits`, `SdfSourceLimits` and `SdfSourceSummary`
(records/total atoms). `single_sdf_mol_block[_with_limits]` extracts one unchanged
MOL block through its CTAB terminator, excluding SD data/delimiter. Python exposes
`roundtrip_sdf_source` and `single_sdf_mol_block` for str/UTF-8 bytes. Existing
low-level graph readers and SDF iterators keep their behavior.

## Bounds and compatibility

Defaults are 16 MiB input and per record, 1 MiB per physical line (including its
ending), 10,000 records, 100,000 atoms per record and 100,000 atoms in total.
WASM tightens input/record to 1,000,000 bytes and per-record/total atoms to 10,000.
Byte/line/record bounds precede native projection; atom counts are checked on
each bounded projection. Native CTAB allocation/count limits also apply. Graphs
are dropped after validation of each record; the returned source has input size.

The local profile requires a nonempty document and an exact `M  END` in every
record after its four header lines. Exact physical `$$$$` lines delimit records.
Whitespace-only records between delimiters refuse; a zero-atom MOL record with
SD metadata is valid and counted. Whitespace after the final delimiter is retained.
The final record may end at EOF without `$$$$`; this is an explicitly tolerated
reader profile, not a claim of strict SDF conformance. SD field/opaque-extension
syntax is not fully validated. Malformed native cores and nonfinite coordinates
fail, including failures in a later record.

Python SDF → SDF rejects coordinate/charge overrides, nondefault name and
nonempty comment. Single-molecule routes also tighten previous acceptance of
multi-record input and unterminated CTABs. Native callers can select resource
limits; single-block extraction always applies a one-record cap. No partial
output, generated coordinates or implicit record selection is added.

## Development comparison

`validation/sdf_source_document/cases.json` contains eight original synthetic MIT
documents × LF/CRLF/no-final-newline: **24** inputs. They cover two/three records,
mixed V2000/V3000, single V3000 with SD fields, duplicate titles, absent terminal
delimiter, blank title, header `M  END` text, duplicate SD fields and literal
`M  END`/dollar substrings in values. Declared ordered titles, record counts and
total atom counts are independent oracles. The packet was used for development;
it is not a holdout or upstream issue reproduction.

The independent observer recognizes physical record framing, version-marker
headers, ordered titles, exact source and unexpected trailing bytes. It is not
a CTAB/graph/SD-field parser or chemistry-equivalence test. Exact source is a
separate metric: regenerated text can differ while retaining the tested record
count/titles. Unrecognized text and NUL tails are recorded separately from
recognized record count.

| Route | Exact source | Recognized record count matches | Ordered titles match | Recognized framing without extra text |
| --- | ---: | ---: | ---: | ---: |
| Baseline Python common | 0/24 | 0/24 | 0/24 | 21/24 |
| Baseline WASM common | 0/24 | 0/24 | 0/24 | 21/24 |
| Candidate Python common | 24/24 | 24/24 | 24/24 | 24/24 |
| Candidate WASM common | 24/24 | 24/24 | 24/24 | 24/24 |
| Open Babel 3.2.1 `-isdf -osdf` | 0/24 | 24/24 | 24/24 | 24/24 |
| Open Babel 3.2.1 with `-x3` | 0/24 | 24/24 | 24/24 | 24/24 |
| Open Babel raw-copy, normal read, framed carrier | 16/24 | 17/24 | 17/24 | 16/24 |
| Open Babel raw-copy, title-only, framed carrier | 14/24 | 15/24 | 15/24 | 16/24 |

Open Babel's ordinary SDF conversion retains all measured record counts/titles;
this fix catches up with that existing batch behavior. Its changed bytes do not
by themselves establish chemical information loss. Baseline common bridges emit
one graph-only V2000 record in 21 cases; the three single-V3000 variants fail.
The candidate retains the original documents in all cases.

Raw-copy controls use seekable files with terminal SDF framing, a completed last
line and, when needed, a blank final SD-value separator. The carrier's bytes,
not the original EOF-tolerant document, are that route's exact-source oracle.
`-ocopy` uses normal molecular reading; `-ocopy -aT` uses title-only reading.
The latter does not establish graph validation. The eight CRLF controls produce
normalized line endings/NULs; the single-record CRLF control is still recognized
but has a NUL tail. Multi-record CRLF controls have unrecognized text/framing.
Title-only reading also omits the blank-title record on both LF carrier variants;
normal raw-copy retains it. These are local workflow observations, not proof
that Open Babel lacks raw-copy functionality. Its own documentation warns about
Windows endings on that utility format.

The first carrier implementation appended `$$$$` after an EOF SD value without
restoring the required blank separator. That control was corrected; scored
product fixture inputs did not change. The earlier complete baseline/candidate
outputs and carrier function remain in
`sdf-source-document-development-carrier-v1-2026-10-22.json.gz` rather than being
silently discarded. The table uses the corrected carrier.

Baseline Python is main `b9412a06` package source with its attested extension;
WASM is attested `c713b21a` with an unrelated CJSON-source branch. Candidate
incorporates main `b7817c64`'s unrelated force-field/conformer updates; relevant
SDF/core/common-route baseline paths are verified unchanged across intervening
main updates. Earlier I/O PRs are not dependencies. Source and binary hashes,
complete outputs and terminal logs are indexed in
`validation/results/sdf-source-document-validation-2026-10-22.json`.

This establishes complete immutable-document retention and clearer record
accounting on a small local profile. It does not establish toolkit-wide
superiority, editable/cross-format fidelity, timing, release availability or
complete CTfile semantics.

## Reproduction

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_sdf_source_document.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/sdf-source-document.json
```

`--baseline` records prior bindings without requiring candidate success.
Subprocess time/output are bounded. Native tests cover source equality/counts,
later-record failure, empty versus zero-atom records, and inclusive byte/line/
record/per-record-atom/total-atom limits. Python and real Node/WASM test both
public paths, multi-record single-molecule refusals, framing and UTF-8 header
classification; Node tests 10,000/10,001 total/per-record atoms. Observer tests
separate count, order, empty title, header/value text, EOF framing and NUL tails.
CI discovers these native/Python/Node suites and runs the observer unit tests.

References: [BIOVIA CTfile Formats](https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf),
[Open Babel raw-copy format](https://openbabel.org/docs/FileFormats/Copy_raw_text.html).
