//! chematic-mcp binary — MCP stdio server for chematic.
//!
//! Reads newline-delimited JSON-RPC 2.0 requests from stdin,
//! writes responses to stdout. Stderr is used for diagnostics only.

#![forbid(unsafe_code)]

use std::io::{BufRead, Write};

fn main() {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    for line in stdin.lock().lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("chematic-mcp: read error: {e}");
                break;
            }
        };

        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        if let Some(response) = chematic_mcp::handle_line(trimmed) {
            let serialized = match serde_json::to_string(&response) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("chematic-mcp: serialization error: {e}");
                    continue;
                }
            };
            if let Err(e) = writeln!(out, "{serialized}") {
                eprintln!("chematic-mcp: write error: {e}");
                break;
            }
            if let Err(e) = out.flush() {
                eprintln!("chematic-mcp: flush error: {e}");
                break;
            }
        }
    }
}
