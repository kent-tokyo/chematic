//! chematic-mcp binary — MCP stdio server for chematic.
//!
//! Reads newline-delimited JSON-RPC 2.0 requests from stdin, writes
//! responses to stdout. Stderr is used for diagnostics only. Supports both
//! the legacy (`2024-11-05`-style `initialize` handshake) and modern
//! (`2026-07-28` stateless, per-request `_meta`) protocol eras on the same
//! connection — see `chematic_mcp::Connection` for the era-pinning logic.

#![forbid(unsafe_code)]

fn main() {
    chematic_mcp::run_stdio();
}
