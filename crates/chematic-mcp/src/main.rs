//! chematic-mcp binary — MCP server for chematic.
//!
//! By default (no arguments, or `--transport stdio`) it reads
//! newline-delimited JSON-RPC 2.0 requests from stdin and writes responses
//! to stdout; stderr is used for diagnostics only. It supports both the
//! legacy (`2024-11-05`-style `initialize` handshake) and modern
//! (`2026-07-28` stateless, per-request `_meta`) protocol eras on the same
//! connection — see `chematic_mcp::Connection` for the era-pinning logic.
//!
//! `--transport streamable-http` serves the modern dialect over HTTP
//! instead (#779), bound to loopback unless `--allow-non-loopback` is given;
//! see `chematic_mcp::http`.

#![forbid(unsafe_code)]

use std::net::IpAddr;
use std::process::ExitCode;

const USAGE: &str = "\
usage: chematic-mcp [--transport stdio]
       chematic-mcp --transport streamable-http [--bind ADDR] [--port PORT]
                    [--allow-non-loopback] [--allowed-origin ORIGIN]...
                    [--max-concurrency N]

  --transport stdio            newline-delimited JSON-RPC on stdin/stdout (default)
  --transport streamable-http  MCP 2026-07-28 Streamable HTTP at POST /mcp
  --bind ADDR                  address to bind (default 127.0.0.1)
  --port PORT                  port to bind, 0 for any free port (default 3000)
  --allow-non-loopback         allow a non-loopback --bind; the HTTP transport has
                               no authentication, so anyone who reaches the port can
                               call every tool
  --allowed-origin ORIGIN      a browser Origin allowed to call the endpoint
                               (repeatable; by default any Origin is refused)
  --max-concurrency N          connections served at once (default 16)";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        chematic_mcp::run_stdio();
        return ExitCode::SUCCESS;
    }
    let mut transport = "stdio".to_string();
    let mut config = chematic_mcp::HttpConfig::default();
    let mut http_only = Vec::new();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let mut value = |flag: &str| -> Result<String, String> {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{flag} needs a value"))
        };
        let result: Result<(), String> = (|| {
            match arg.as_str() {
                "-h" | "--help" => {
                    println!("{USAGE}");
                    std::process::exit(0);
                }
                "-V" | "--version" => {
                    println!("chematic-mcp {}", env!("CARGO_PKG_VERSION"));
                    std::process::exit(0);
                }
                "--transport" => transport = value("--transport")?,
                "--bind" => {
                    http_only.push("--bind");
                    config.bind = value("--bind")?
                        .parse::<IpAddr>()
                        .map_err(|e| format!("--bind: {e}"))?;
                }
                "--port" => {
                    http_only.push("--port");
                    config.port = value("--port")?
                        .parse()
                        .map_err(|e| format!("--port: {e}"))?;
                }
                "--allow-non-loopback" => {
                    http_only.push("--allow-non-loopback");
                    config.allow_non_loopback = true;
                }
                "--allowed-origin" => {
                    http_only.push("--allowed-origin");
                    config.allowed_origins.push(value("--allowed-origin")?);
                }
                "--max-concurrency" => {
                    http_only.push("--max-concurrency");
                    config.max_concurrency = value("--max-concurrency")?
                        .parse()
                        .map_err(|e| format!("--max-concurrency: {e}"))?;
                }
                other => return Err(format!("unknown argument: {other}")),
            }
            Ok(())
        })();
        if let Err(message) = result {
            eprintln!("chematic-mcp: {message}\n\n{USAGE}");
            return ExitCode::from(2);
        }
    }
    match transport.as_str() {
        "stdio" => {
            if let Some(flag) = http_only.first() {
                eprintln!("chematic-mcp: {flag} applies to --transport streamable-http only");
                return ExitCode::from(2);
            }
            chematic_mcp::run_stdio();
            ExitCode::SUCCESS
        }
        "streamable-http" => match chematic_mcp::run_http(config) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("chematic-mcp: {e}");
                ExitCode::FAILURE
            }
        },
        other => {
            eprintln!(
                "chematic-mcp: unknown transport '{other}' (stdio or streamable-http)\n\n{USAGE}"
            );
            ExitCode::from(2)
        }
    }
}
