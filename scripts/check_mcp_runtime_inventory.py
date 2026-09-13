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


def request(process: subprocess.Popen[str], request_id: int, method: str, params: dict[str, object]) -> dict[str, object]:
    payload = {"jsonrpc": "2.0", "id": request_id, "method": method, "params": {**params, "_meta": META}}
    assert process.stdin is not None and process.stdout is not None
    process.stdin.write(json.dumps(payload) + "\n")
    process.stdin.flush()
    line = process.stdout.readline()
    if not line:
        stderr = process.stderr.read() if process.stderr else ""
        raise RuntimeError(f"MCP process ended before response: {stderr}")
    response = json.loads(line)
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
        result = request(process, 1, "tools/list", {})
        tools = result.get("tools")
        if not isinstance(tools, list):
            raise RuntimeError("tools/list result has no tools array")
        names = [tool.get("name") for tool in tools if isinstance(tool, dict)]
        if len(tools) != 20 or len(set(names)) != 20 or any(not isinstance(name, str) for name in names):
            raise RuntimeError(f"expected exactly 20 uniquely named tools, got {names}")
        missing_schemas = [name for name, tool in zip(names, tools, strict=True) if not isinstance(tool.get("inputSchema"), dict) or not isinstance(tool.get("outputSchema"), dict)]
        if missing_schemas:
            raise RuntimeError(f"tools missing input/output schemas: {missing_schemas}")
        canonical = request(process, 2, "tools/call", {"name": "canonical_smiles", "arguments": {"smiles": "OCC"}})
        structured = canonical.get("structuredContent")
        if not isinstance(structured, dict) or not isinstance(structured.get("canonical"), str) or not structured["canonical"]:
            raise RuntimeError(f"canonical_smiles runtime output lacks a canonical string: {structured!r}")
        inventory = {
            "schema": "chematic.mcp-runtime-inventory.v1",
            "execution_boundary": "MCP JSON-RPC stdio subprocess",
            "tool_count": len(tools),
            "tools": [{"name": tool["name"], "has_input_schema": True, "has_output_schema": True} for tool in tools],
            "representative_call": {"tool": "canonical_smiles", "input": "OCC", "structured_output": structured},
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
