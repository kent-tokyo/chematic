#!/usr/bin/env python3
"""Query the published MCP wire surface and record a runtime tool inventory.

This deliberately speaks JSON-RPC over stdio instead of importing server
internals. It is therefore a packaging/runtime complement to the Rust schema
conformance suite and detects a binary whose registry differs from source.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = ["cargo", "run", "--offline", "-q", "-p", "chematic-mcp", "--"]
META = {
    "io.modelcontextprotocol/protocolVersion": "2026-07-28",
    "io.modelcontextprotocol/clientInfo": {"name": "chematic-runtime-inventory", "version": "1"},
    "io.modelcontextprotocol/clientCapabilities": {},
}


def raw_request(
    process: subprocess.Popen[str], request_id: int, method: str, params: dict[str, object]
) -> dict[str, object]:
    payload = {"jsonrpc": "2.0", "id": request_id, "method": method, "params": {**params, "_meta": META}}
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps(payload) + "\n")
    process.stdin.flush()
    line = process.stdout.readline()
    if not line:
        stderr = process.stderr.read() if process.stderr else ""
        raise RuntimeError(f"MCP process ended before response: {stderr}")
    return json.loads(line)


def request(process: subprocess.Popen[str], request_id: int, method: str, params: dict[str, object]) -> dict[str, object]:
    response = raw_request(process, request_id, method, params)
    if "error" in response:
        raise RuntimeError(f"{method} failed: {response['error']}")
    return response["result"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", nargs="+", default=DEFAULT_BINARY, help="MCP server command")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()
    process = subprocess.Popen(args.binary, cwd=ROOT, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        discovery = request(process, 1, "server/discover", {})
        server_info = discovery.get("serverInfo")
        if (
            discovery.get("resultType") != "complete"
            or discovery.get("supportedVersions") != ["2026-07-28"]
            or not isinstance(server_info, dict)
            or server_info.get("name") != "chematic-mcp"
            or not isinstance(server_info.get("version"), str)
            or discovery.get("cacheScope") != "public"
            or not isinstance(discovery.get("ttlMs"), int)
        ):
            raise RuntimeError(f"server/discover did not return the modern discovery contract: {discovery!r}")
        result = request(process, 2, "tools/list", {})
        tools = result.get("tools")
        if not isinstance(tools, list):
            raise RuntimeError("tools/list result has no tools array")
        names = [tool.get("name") for tool in tools if isinstance(tool, dict)]
        if len(tools) != 20 or len(set(names)) != 20 or any(not isinstance(name, str) for name in names):
            raise RuntimeError(f"expected exactly 20 uniquely named tools, got {names}")
        missing_schemas = [name for name, tool in zip(names, tools, strict=True) if not isinstance(tool.get("inputSchema"), dict) or not isinstance(tool.get("outputSchema"), dict)]
        if missing_schemas:
            raise RuntimeError(f"tools missing input/output schemas: {missing_schemas}")
        canonical = request(process, 3, "tools/call", {"name": "canonical_smiles", "arguments": {"smiles": "OCC"}})
        structured = canonical.get("structuredContent")
        if not isinstance(structured, dict) or not isinstance(structured.get("canonical"), str) or not structured["canonical"]:
            raise RuntimeError(f"canonical_smiles runtime output lacks a canonical string: {structured!r}")
        invalid_smiles = request(process, 4, "tools/call", {"name": "canonical_smiles", "arguments": {"smiles": "C1CC"}})
        invalid_error = invalid_smiles.get("structuredContent")
        if (
            invalid_smiles.get("isError") is not True
            or not isinstance(invalid_error, dict)
            or not isinstance(invalid_error.get("error"), dict)
            or invalid_error["error"].get("code") != "INVALID_SMILES"
        ):
            raise RuntimeError(f"invalid SMILES was not reported as a structured domain error: {invalid_smiles!r}")
        oversized = raw_request(
            process,
            5,
            "tools/call",
            {"name": "canonical_smiles", "arguments": {"smiles": "C" * 100_001}},
        )
        oversized_error = oversized.get("error")
        if (
            not isinstance(oversized_error, dict)
            or oversized_error.get("code") != -32602
            or "maxLength 100000" not in str(oversized_error.get("message", ""))
        ):
            raise RuntimeError(f"oversized SMILES was not rejected as invalid parameters: {oversized!r}")
        recovery = request(process, 6, "tools/call", {"name": "canonical_smiles", "arguments": {"smiles": "CCN"}})
        recovery_structured = recovery.get("structuredContent")
        if not isinstance(recovery_structured, dict) or recovery_structured.get("canonical") != "C(C)N":
            raise RuntimeError(f"server did not recover after a bounded-input rejection: {recovery!r}")
        inventory = {
            "schema": "chematic.mcp-runtime-inventory.v3",
            "execution_boundary": "MCP JSON-RPC stdio subprocess",
            "discovery": {
                "protocol_version": "2026-07-28",
                "result_type": discovery["resultType"],
                "server_info": server_info,
                "cache_scope": discovery["cacheScope"],
                "ttl_ms": discovery["ttlMs"],
            },
            "tool_count": len(tools),
            "tools": [{"name": tool["name"], "has_input_schema": True, "has_output_schema": True} for tool in tools],
            "representative_call": {"tool": "canonical_smiles", "input": "OCC", "structured_output": structured},
            "error_contracts": {
                "invalid_smiles": {
                    "tool": "canonical_smiles",
                    "input": "C1CC",
                    "transport_error": False,
                    "structured_error_code": "INVALID_SMILES",
                },
                "oversized_smiles": {
                    "tool": "canonical_smiles",
                    "input_bytes": 100_001,
                    "transport_error_code": -32602,
                    "message_contains": "maxLength 100000",
                },
                "recovery": {"tool": "canonical_smiles", "input": "CCN", "canonical": "C(C)N"},
            },
            "status": "pass",
        }
    finally:
        if process.stdin:
            process.stdin.close()
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
    rendered = json.dumps(inventory, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        output = args.json_out if args.json_out.is_absolute() else ROOT / args.json_out
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
