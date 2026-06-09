use serde_json::json;
use std::io::{self, BufRead};

mod tools;
use tools::ChematicTools;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tools = ChematicTools::new();

    eprintln!("chematic-mcp server started");

    let stdin = io::stdin();
    let reader = stdin.lock();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        // Parse incoming MCP request
        let request: serde_json::Value = match serde_json::from_str(&line) {
            Ok(v) => v,
            Err(_) => {
                eprintln!("Failed to parse JSON: {}", line);
                continue;
            }
        };

        // Handle request (simplified)
        let response = if let Some(tool_name) = request.get("tool").and_then(|v| v.as_str()) {
            let args = request.get("args").unwrap_or(&json!({}));
            match tools.call_tool(tool_name, args).await {
                Ok(result) => json!({ "status": "ok", "result": result }),
                Err(e) => json!({ "status": "error", "error": e }),
            }
        } else {
            json!({ "status": "error", "error": "No tool specified" })
        };

        println!("{}", serde_json::to_string(&response)?);
    }

    Ok(())
}
