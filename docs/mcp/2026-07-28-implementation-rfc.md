# chematic-mcp: 2026-07-28 stateless tools-only implementation RFC

Status: implemented (this PR). Scope: `chematic-mcp` protocol core only —
Remote HTTP, OAuth, Tasks, and MCP Apps are explicitly out of scope (see
"Descoped" below).

This document records the primary-source basis for every non-obvious design
decision in this PR, and — per the grounding note in the original task —
is explicit about where this implementation departs from the task brief's
assumptions, with citations, rather than silently reconciling them.

## 0. What this server is (and is not)

**MCP 2026-07-28 tools-only stateless stdio server**, with byte-compatible
legacy (`2024-11-05`-style) support on the same stdio connection. It is
*not* "MCP 2026-07-28 complete" — resources, prompts, sampling, roots,
logging, tasks, MCP Apps, subscriptions, and any HTTP transport are all
unimplemented and not advertised (see `docs/specification/draft/server/discover.mdx`
at the pinned commit: a server "MUST NOT" advertise a capability it doesn't
support).

## 1. Schema provenance

Primary source: `modelcontextprotocol/modelcontextprotocol`, tag
`2026-07-28-RC`, commit `9d700ed62dcf86cb77475c9b81930611a9182f46` (verified:
`git/refs/tags/2026-07-28-RC` resolves to exactly this SHA — the short hash
given in the task brief is the real tagged RC, not an arbitrary commit).
Full fetch manifest with URLs and SHA-256 hashes:
`validation/manifests/mcp_2026_07_28_rc.json`.

Everything below was read from that fetched schema — none of it is
reconstructed from training-data memory of 2024-11-05/2025-03-26/2025-06-18/2025-11-25.

## 2. B0 — SDK adoption decision: **declined**

Checklist from the task brief, evaluated against the *real*, currently
released `rmcp` (official Rust SDK, `modelcontextprotocol/rust-sdk`)
v3.0.0-beta.2 (crates.io, published 2026-07-24):

| Checklist item | rmcp v3.0.0-beta.2 status | Evidence |
|---|---|---|
| `server/discover` | Implemented | PR #973 "add server discovery and negotiation (SEP-2575)" |
| 2026-07-28 stateless era | Implemented | PR #999 "serve draft-version requests statelessly per SEP-2567"; `StreamableHttpServerConfig::legacy_session_mode` doc: "the `2026-07-28` draft version is always served statelessly" |
| dual-era stdio negotiation | Implemented | PR #995 "add modern client lifecycle modes (SEP-2575)" |
| per-request MCP metadata envelope | Implemented | PR #993 "align metadata models with draft schema" |
| `CacheableResult` ttlMs/cacheScope | Implemented | PR #889/#1025 "implement SEP-2549 cache hints" / "client-side TTL-honoring response cache" |
| 2026 error codes | Implemented, but **-32020/-32021/-32022**, not the brief's exact guess validated against the wrong source (see §4) | `crates/rmcp/src/model.rs` `ErrorCode` constants (fetched at rust-sdk HEAD `14298b72e0b25473ea79d5465fe186e22eb86397`) |
| full JSON Schema 2020-12 tool schemas | Implemented | PR #895/#933 "relax outputSchema/structuredContent to accept non-object JSON Schema types (SEP-2106)" |
| structuredContent | Implemented | as above |

**Every checklist item is substantively met.** So why not adopt it?

1. **It's a beta, not a stable release** (`3.0.0-beta.2`). `chematic-mcp`
   is a published, versioned crate on crates.io; depending on a pre-1.0-of-a-major
   beta for wire-protocol-critical logic is a stability risk this repo's
   dependency policy (pinned `0.7.0` path deps, no other crate in the
   workspace depends on an unstable beta) does not otherwise accept.
2. **Scope bloat.** Adopting `rmcp` pulls in Streamable HTTP, OAuth/DCR
   (device/client-credentials flows), the Tasks extension, and a
   distributed SSE event store — all explicitly deferred to a later PR by
   the task brief (§11). A hand-rolled core has zero surface for any of
   these; `rmcp` has all of it compiled in.
3. **Full rewrite risk against byte-compat fixtures.** All 20 tools would
   need to be re-authored as `rmcp` tool handlers. The existing legacy wire
   fixtures (§13) need to stay byte-identical; a full SDK migration is a
   much larger, harder-to-audit diff than a ~1200-line hand-rolled codec
   for exactly the subset of MCP this server needs (tools-only,
   stateless, stdio).

Given the task brief's explicit "don't make this a giant PR" instruction,
declining SDK adoption for policy reasons — not because the SDK is
incapable — is the correct call. Corroborating detail: rmcp's own PR #1038
("omit resultType for legacy protocol sessions") independently arrived at
the same legacy/modern `resultType` split this implementation uses (§6).

**Decision: hand-rolled, typed codec.** No new runtime dependency added to
the published crate (see §9 for the one *dev*-dependency added for tests).

## 3. Error codes: a real, sourced discrepancy — not a guess

The task brief's §9 named `-32020`/`-32021`/`-32022` for header-mismatch /
missing-capability / unsupported-version. The **pinned RC tag** (`9d700ed`)
does **not** define these — its `schema.ts` defines:

```
MISSING_REQUIRED_CLIENT_CAPABILITY = -32003
UNSUPPORTED_PROTOCOL_VERSION        = -32004
```

with **no** header-mismatch code at all (`grep -n "header" schema.ts` at
`9d700ed` returns nothing but a hyperlink to the HTTP `MCP-Protocol-Version`
header in a doc comment for a different field).

Fetching the spec repo's **untagged `main` branch** (no `2026-07-28` final
tag exists yet — confirmed via `git/refs/tags/2026-07-28`, 404) at commit
`7634684382c3d14cf7e9f14073fe40a2d8ace3fa` (2026-07-23, four days before
this PR was written) shows the codes were renumbered and a new code added:

```
HEADER_MISMATCH                     = -32020   (new)
MISSING_REQUIRED_CLIENT_CAPABILITY  = -32021   (was -32003)
UNSUPPORTED_PROTOCOL_VERSION        = -32022   (was -32004)
```

with a new documented allocation scheme: `-32000..-32019` is
implementation-defined (never assigned by the spec), `-32020..-32099` is
reserved for the spec, allocated sequentially starting at `-32020`.

Three independent sources agree on the renumbered values:

| Source | Value used |
|---|---|
| Spec `main` @ `7634684` (2026-07-23, untagged) | `-32020`/`-32021`/`-32022` |
| `rmcp` v3.0.0-beta.2 `ErrorCode` constants (rust-sdk HEAD `14298b72e0b25473ea79d5465fe186e22eb86397`) | `-32020`/`-32021`/`-32022` |
| Task brief §9 (this repo) | `-32020`/`-32021`/`-32022` |

**Decision: ship `-32020`/`-32021`/`-32022`**, not the RC tag's own
`-32003`/`-32004`. The RC tag's values are demonstrably dead code by the
time this PR ships; shipping them would mean speaking a wire contract no
real client (including `rmcp`) will ever send. See
`validation/manifests/mcp_2026_07_28_rc.json`'s
`post_rc_untagged_drift_checked_for_reconciliation` block for the full
diff and fetch provenance.

Emission status in this server:

- `-32700`/`-32600`/`-32601`/`-32602`/`-32603`: standard JSON-RPC, all
  actively emitted (parse errors, malformed requests, unknown methods,
  argument/schema violations, panics).
- `-32022` `UNSUPPORTED_PROTOCOL_VERSION`: actively emitted — the only one
  of the three MCP-specific codes this server can actually trigger, since
  it's the only one of the three tied to an input this server checks
  (`_meta.protocolVersion` against a fixed one-element supported list).
- `-32020` `HEADER_MISMATCH`: **defined, never emitted.** This is a
  stdio-only transport; there are no HTTP headers to mismatch. Reserved as
  a typed constant/test fixture for a future HTTP adapter PR, per the task
  brief's explicit request.
- `-32021` `MISSING_REQUIRED_CLIENT_CAPABILITY`: **defined, never
  emitted.** All 20 tools are plain request/response with no elicitation,
  sampling, or roots dependency; this server never requires *any* client
  capability, so there is no legitimate input that triggers this code. A
  codec-level test proves the shape can be serialized correctly (see
  `crates/chematic-mcp/src/server.rs`'s
  `header_mismatch_and_missing_capability_codes_are_serializable`); the
  conformance matrix (§14) marks the *wire scenario* `not_applicable`
  rather than wiring a fake trigger to turn a row green.

## 4. Dual-era negotiation vs. the protocol's own statelessness

`docs/specification/draft/basic/lifecycle.mdx` (pinned RC commit), lines
17–22: servers **MUST NOT** rely on prior requests over the same
connection to establish protocol version, capabilities, or client
identity — "an open connection or STDIO process is not a conversation or
session."

The task brief's §4 asks for exactly the opposite-sounding thing: pin the
era on the connection's first request, and reject a request that tries to
switch dialects mid-connection.

These are reconciled, not in tension, once the pin is understood narrowly:

- The pin decides **which dialect's grammar** a given method/params shape
  belongs to (legacy vs. modern) — it never supplies data a request didn't
  provide itself.
- Every modern request still carries, and is validated against, its own
  full `_meta` triple (`protocolVersion`/`clientInfo`/`clientCapabilities`)
  — `RequestContext` is built fresh per request from that request's own
  `params`, never from a cached previous value (see
  `transport::Connection::handle_line`).
- §10 of the task brief explicitly whitelists "connection-pinned protocol
  era" as allowed transport state, alongside "stdio process lifecycle" —
  this is precisely that carve-out, not an invented exception to
  statelessness.
- The pin rejects a **dialect switch** (e.g., a legacy `initialize` after a
  modern `server/discover` already answered on the same connection) with a
  typed `-32600 Invalid Request` error — there is no dedicated MCP error
  code for this (the spec doesn't officially model a stateful connection
  noticing a dialect switch, since it's specified to be irrelevant), so the
  standard JSON-RPC "this request is invalid in context" code is used, with
  a `data.reason` explaining the mismatch.

Notifications (no `id`) never pin — they can't be replied to or classified
against, and per `lifecycle.mdx`'s own STDIO backward-compatibility probe
description, the entire negotiation flow is keyed off id-bearing requests
(`server/discover` first, `Method not found` fallback to `initialize`).

## 5. `resultType`: a requirement the brief didn't mention

`Result.resultType` (`docs/specification/draft/basic/index.mdx` §"ResultType",
pinned RC commit) is **required** on every modern-era result: `"complete"`
for a finished request, with backward-compat text stating that a client
**MUST** treat an *absent* `resultType` as `"complete"` when talking to a
server on an earlier protocol revision. This server:

- Includes `"resultType": "complete"` on every modern `DiscoverResult`,
  `ListToolsResult`, and `CallToolResult`.
- Never includes it on any legacy-era result (legacy fixtures assert its
  absence explicitly — see `lib.rs`'s `test_initialize`/`test_tools_list`/
  `test_tools_call_parse_smiles`).

This server implements no Multi Round-Trip Requests (elicitation/sampling),
so `"input_required"` is never produced.

## 6. `server/discover` result shape: a three-way hybrid, not a pure pick

`DiscoverResult` differs between the RC tag and the untagged `main` drift
(§3's methodology, same three sources):

| Field | RC tag (`9d700ed`) | Spec `main` @ `7634684` (untagged) | `rmcp` v3.0.0-beta.2 `DiscoverResult` struct | This server |
|---|---|---|---|---|
| `serverInfo` (top-level) | Required field | **Deleted** (moved to an optional `_meta.io.modelcontextprotocol/serverInfo`, self-reported, non-security-use) | Kept as a required top-level field | **Kept**, top-level |
| `ttlMs`/`cacheScope` | Absent (`extends Result`) | Present (`extends CacheableResult`) | Present (both fields on the struct) | **Present** |

Two of three real sources (a shipping SDK and the untagged drift) agree
`ttlMs`/`cacheScope` belong on `DiscoverResult`; `rmcp` is the only source
that actually ships code, and it keeps `serverInfo` top-level even while
adopting the cache fields. This server mirrors `rmcp`'s pragmatic hybrid
rather than either doc snapshot in isolation, on the reasoning that a real
SDK's shipped behavior is stronger evidence of what interoperates than a
prose diff in an unreleased branch.

`capabilities` declares only `{"tools": {}, "extensions": {}}` — no
`resources`/`prompts`/`sampling`/`roots`/`logging`/`completions` keys are
present at all (not even as `false`/empty — `ServerCapabilities` makes
every key optional, so omission is the correct way to say "not
implemented," per task brief §5's explicit prohibition list).

## 7. `ping` and `initialize`: removed in the modern era, not "kept"

The task brief's §4 lists `ping` as a method the modern era should
implement. `docs/specification/draft/changelog.mdx` at the pinned RC
commit, "Major changes" item 5: **"Remove `ping`, `logging/setLevel`, and
`notifications/roots/list_changed`."** — not deprecated, removed; confirmed
absent from `schema.ts` at `9d700ed` (`grep -i ping` returns nothing), and
absent from the "Deprecated" registry in `deprecated.mdx` (it isn't in a
transition window; it's gone).

This server: `ping` and `initialize` remain fully supported in the legacy
era (byte-compatible), and in the modern era both simply fall through to
`-32601 Method not found` like any other unrecognized method — this is not
a dialect violation (see §4), just an accurate reflection that these
methods do not exist in the modern vocabulary. Test:
`server::tests::modern_ping_is_method_not_found`.

## 8. Tool-call error taxonomy (§9 of the brief)

`ToolCallError` (`tools.rs`) splits every tool failure into exactly two
kinds, matching the brief's own three examples:

- `InvalidArgs`: argument-shape/schema violations (missing/wrong-typed
  argument, unknown tool name, a value that fails the tool's declared
  `inputSchema` — e.g. `find_mcs`'s `minItems: 2`). Modern era: rejected
  before any chemistry runs, as `-32602 Invalid Params`. Legacy era:
  unchanged wire behavior (implementation-defined `-32000`, see below).
- `Domain`: the arguments were well-formed but the requested chemistry
  failed (unparseable SMILES/SMARTS, molecule too large/disconnected for a
  bounded algorithm, a failed PubChem lookup). Modern era: a **successful**
  `tools/call` response with `isError: true` and a machine-readable
  `structuredContent.error.{code,message,details}` object — never a
  JSON-RPC transport error, so an LLM client can see and self-correct from
  it. Legacy era: unchanged wire behavior.
- Genuine Rust panics (neither variant — caught via `catch_unwind` around
  every tool dispatch, in both eras) map to `-32603 Internal Error`.
  Verified `catch_unwind` is not a silent no-op here: the workspace root
  `Cargo.toml`'s only `[profile]` section (`[profile.release]`) does not
  set `panic = "abort"`, so the default `unwind` strategy applies.

**Legacy-era byte-compatibility note:** before this PR, *every* tool
failure (both kinds) surfaced as a JSON-RPC-level error with the
implementation-defined code `-32000` and the tool's raw message string —
arguably not fully spec-correct even under 2024-11-05 (which also defines
`isError`-based tool results), but this is the server's *actual prior wire
behavior*, and no existing fixture pins the old shape as intentional
protocol conformance rather than an accident. Per the task brief's own
explicit allowance ("legacy eraでは既存のtext JSON contractを維持して構いません",
§8), this PR leaves the legacy-era error shape untouched rather than
"fixing" it — the fix only applies to the new modern era, where it's
required by spec, not retrofitted onto a path this task never asked to
touch.

## 9. `outputSchema` in the legacy `tools/list`: one deliberate, additive change to the "byte-identical" rule

Every one of the 20 tools now carries an `outputSchema` (task brief §7,
unconditional — no era carve-out stated there), and every `inputSchema`
gained `additionalProperties: false` plus tightened `minLength`/`maxLength`/
`minItems`/`maxItems` bounds. These are emitted in **both** eras' `tools/list`
responses, because they're tool metadata (part of the `Tool` object), not a
new protocol-level envelope construct — unlike `ttlMs`/`cacheScope`/
`resultType`, which are genuinely new wire-level fields section 6 of the
brief says legacy must not carry.

This is the one place this PR does not produce byte-identical legacy
output relative to the pre-refactor server: every `Tool` object in a
legacy `tools/list` response now has more keys than before. This is
intentional (required unconditionally by §7) and low-risk (JSON's
additive-field tolerance; no legacy fixture asserted an exact-key-set
match, only tool count).

## 10. `structuredContent`/`content` design: matches the RC's own example, not the brief's illustrative one

The task brief's §8 illustrative example shows `content[0].text` as
hand-written prose ("Parsed molecule: 6 atoms, 6 bonds") distinct from
`structuredContent`. The **actual RC schema's own example fixture**
(`schema/draft/examples/CallToolResult/result-with-structured-content.json`
at `9d700ed`) instead uses the JSON-stringified payload as `content[0].text`,
identical in content to `structuredContent`:

```json
{
  "content": [{ "type": "text", "text": "{\"temperature\": 22.5, ...}" }],
  "structuredContent": { "temperature": 22.5, ... }
}
```

This server follows the real fixture: `content` is generated from the same
one `Value` the tool computed (`payload.to_string()`), and `structuredContent`
is that same `Value` — never two independently-authored representations
that could drift apart, and no bespoke per-tool prose to invent and
maintain for 20 tools. Each `tool_*` function in `tools.rs` computes its
chemistry exactly once; `server.rs`'s presentation layer (`content_only`
for legacy, `success_result`/`domain_error_result` for modern) is the only
place that decides how to wrap it per era.

## 11. Caching semantics and staleness

`tools/list`'s cache hint: `ttlMs: 86400000` (24h), `cacheScope: "public"`.
Per `docs/specification/draft/server/utilities/caching.mdx` (pinned RC
commit), TTL is a freshness *hint*, not a correctness guarantee — the
underlying data can change before it expires. This server's registry is a
compiled-in, static `json!` literal with no runtime mutation path, so
"can the registry go stale mid-TTL" reduces to "can the registry change
while the process is running," which it cannot (a change ships as a new
`chematic-mcp` release = a new process). A client that restarts the server
(new binary, new version) gets a fresh, uncached response on its first
request regardless of any TTL window a previous process advertised —
process/version boundaries are the versioning mechanism here, not
`notifications/tools/list_changed` (unimplemented — this server never
declares `listChanged: true`).

## 12. Runtime vs. test-time schema validation

Two independent validators exist, deliberately:

- **Runtime** (`src/schema.rs`): a small hand-rolled subset validator
  covering exactly the keywords the 20 tools' schemas use (`type`,
  `properties`, `required`, `additionalProperties`, `items`, `minLength`,
  `maxLength`, `minItems`, `maxItems`, `minimum`, `maximum`, `enum`,
  `const`, `oneOf`). Runs before every modern `tools/call` dispatch to
  produce `-32602` on a schema violation. No new dependency in the
  published crate. Bounded by the same `MAX_JSON_DEPTH` used for
  adversarial-input limits (§13), so a schema/data pair engineered to
  recurse cannot make this unbounded.
- **Test-time** (`crates/chematic-mcp/tests/schema_conformance.rs`): the
  real `jsonschema` crate (MIT, `default-features = false` to drop its
  optional remote-`$ref`-fetching features — this project never resolves
  external `$ref`s, so even a dev-dependency shouldn't carry that
  capability) as a **dev-dependency only**. Confirms every `inputSchema`/
  `outputSchema` is itself valid JSON Schema 2020-12
  (`jsonschema::meta::validate`), and that every tool's actual
  `structuredContent` conforms to its declared `outputSchema`
  (`jsonschema::validator_for(...).is_valid(...)`). This independently
  verifies the hand-rolled runtime validator's target schemas are correct,
  rather than trusting one implementation's self-consistency.

## 13. Adversarial-input limits (§12)

`protocol.rs`: `MAX_REQUEST_BYTES` (1 MiB, checked on the raw line before
parsing), `MAX_JSON_DEPTH` (64, checked *twice*: a byte-level bracket-depth
prescan before `serde_json::from_str` ever recurses into the input, and
again on the parsed `Value` tree as an independent second check),
`MAX_ARRAY_LEN` (10,000 elements), `MAX_STRING_LEN` (256 KiB). The prescan
respects string literals (tracks `\"` escapes) so brackets inside a JSON
string payload are never miscounted as structural nesting — verified by
`protocol::tests::prescan_accepts_brackets_inside_strings`.

## 14. Descoped (explicitly, per task brief §11)

Streamable HTTP, `Mcp-Method`/`Mcp-Name`/`Mcp-Param-*` HTTP headers,
OAuth/OIDC, MCP Apps, the Tasks extension, `subscriptions/listen`,
resources, prompts, an OpenTelemetry exporter, remote deployment. The
transport/protocol-codec/server-core/tool-registry layering exists
precisely so an HTTP adapter can be added later without touching
`tools.rs` or `server.rs`'s dispatch logic — but no HTTP code exists in
this PR.

`name_to_smiles` remains the sole tool making a live network call
(PubChem REST API) — unchanged, documented in `crates/chematic-mcp/README.md`
and excluded from the automated all-20-tools smoke test to avoid CI
flakiness under network restrictions (its schema/error-taxonomy are still
checked).

## 15. Final-vs-RC reconciliation (re-checked at PR time)

Re-verified on 2026-07-27, the day this PR is being opened (one day before
the `2026-07-28` name the spec repo's own tag references):

```
$ gh api repos/modelcontextprotocol/modelcontextprotocol/tags --jq '.[] | select(.name | test("2026-07-28")) | .name'
2026-07-28-RC

$ gh api repos/modelcontextprotocol/modelcontextprotocol/commits/main --jq '{sha: .sha, date: .commit.committer.date}'
{"date":"2026-07-23T23:49:30Z","sha":"7634684382c3d14cf7e9f14073fe40a2d8ace3fa"}
```

Both results are unchanged from the check recorded in §3 and in
`validation/manifests/mcp_2026_07_28_rc.json`'s
`post_rc_untagged_drift_checked_for_reconciliation` block: still only the
`-RC` tag exists, `main` is still sitting at the same untagged commit it was
at four days ago. **No final, non-RC `2026-07-28` tag or release exists as
of this PR.** Per the task dispatcher's explicit instruction, this is
recorded as "not yet applicable" rather than fabricated — there is no final
spec text to diff against, only the RC (used as primary source, §1) and the
untagged `main` drift already reconciled into the error-code decision (§3).

If a final tag appears after this PR is opened, the reconciliation work is:
re-fetch `schema.ts`/`schema.json` at the new tag, diff against the RC files
already hashed in the manifest, and re-run `tests/schema_conformance.rs`
(which validates against the schemas embedded in `tools.rs`, not a fetched
copy, so a real behavioral drift would need a follow-up PR regardless of
what a text diff shows).
