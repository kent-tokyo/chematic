use chematic_mcp::{
    ClientMeta, INVALID_PARAMS, JsonRpcRequest, McpServer, ProtocolEra, RequestContext,
};
use serde_json::{Value, json};

fn call(ctx: &RequestContext, method: &str, params: &Value) -> Value {
    McpServer::new().handle_request(
        JsonRpcRequest {
            id: &json!(42),
            method,
            params,
        },
        ctx,
    )
}

#[test]
fn legacy_missing_tool_name_returns_invalid_params_with_request_id() {
    let ctx = RequestContext {
        era: ProtocolEra::Legacy,
        client_meta: None,
    };
    let response = call(&ctx, "tools/call", &json!({}));
    assert_eq!(response["id"], 42);
    assert_eq!(response["error"]["code"], INVALID_PARAMS);
    assert_eq!(response["error"]["message"], "Missing tool name in params");
}

#[test]
fn modern_dispatch_rejects_missing_context_name_and_unknown_tool() {
    let mut ctx = RequestContext {
        era: ProtocolEra::Modern20260728,
        client_meta: None,
    };
    let response = call(&ctx, "tools/call", &json!({}));
    assert_eq!(response["error"]["code"], INVALID_PARAMS);
    assert!(
        response["error"]["message"]
            .as_str()
            .unwrap()
            .contains("missing required _meta")
    );
    ctx.client_meta = Some(ClientMeta {
        protocol_version: "2026-07-28".into(),
        client_info: None,
        client_capabilities: json!({}),
        log_level: None,
    });
    let missing = call(&ctx, "tools/call", &json!({}));
    assert_eq!(missing["error"]["code"], INVALID_PARAMS);
    assert_eq!(missing["error"]["message"], "Missing tool name in params");
    let unknown = call(&ctx, "tools/call", &json!({"name": "unknown-tool"}));
    assert_eq!(unknown["error"]["code"], INVALID_PARAMS);
    assert_eq!(unknown["error"]["message"], "Unknown tool: unknown-tool");
}
