use crate::{http, protocol, schema, tools, transport};
use serde_json::{Value, json};
use std::io::{self, BufReader, Cursor, Read};

#[test]
fn schema_boundaries_match_independent_json_schema_validator() {
    for (schema, values) in [
        (json!({"const":7}), vec![json!(7), json!(8)]),
        (json!({"enum":["a",2,null]}), vec![json!("a"), json!(3)]),
        (
            json!({"type":"string","minLength":2,"maxLength":4}),
            vec![json!(""), json!("a"), json!("ab"), json!("abcde")],
        ),
        (
            json!({"type":"integer","minimum":-1,"maximum":2}),
            vec![json!(-2), json!(-1), json!(2), json!(3), json!(1.5)],
        ),
        (
            json!({"type":["null","number"]}),
            vec![json!(null), json!(3.5), json!("x")],
        ),
        (
            json!({"type":"object","properties":{"x":{"type":"boolean"}},"additionalProperties":false}),
            vec![json!({}), json!({"x":true}), json!({"x":1}), json!({"y":1})],
        ),
        (
            json!({"type":"array","items":{"type":"string"},"minItems":1,"maxItems":2}),
            vec![json!([]), json!(["a"]), json!([1]), json!(["a", "b", "c"])],
        ),
    ] {
        for value in values {
            assert_eq!(
                schema::validate(&schema, &value).is_ok(),
                jsonschema::is_valid(&schema, &value),
                "{schema}: {value}"
            );
        }
    }
    let mut shape = json!({"type":"null"});
    let mut value = Value::Null;
    for _ in 0..protocol::MAX_JSON_DEPTH + 2 {
        shape = json!({"type":"array","items":shape});
        value = json!([value]);
    }
    assert!(
        schema::validate(&shape, &value)
            .unwrap_err()
            .contains("depth")
    );
    assert!(
        protocol::validate_value_limits(&value)
            .unwrap_err()
            .contains("depth")
    );
}

#[test]
fn tool_resource_and_argument_errors_keep_structured_codes() {
    for (name, args, code) in [
        (
            "parse_smiles",
            json!({"smiles":"C".repeat(100001)}),
            "INVALID_ARGUMENTS",
        ),
        (
            "parse_smiles",
            json!({"smiles":"C".repeat(10001)}),
            "MOLECULE_TOO_LARGE",
        ),
        (
            "moljson_to_smiles",
            json!({"moljson":"x".repeat(100001)}),
            "INVALID_ARGUMENTS",
        ),
        (
            "smarts_match",
            json!({"smarts":"C".repeat(100001),"smiles":"C"}),
            "INVALID_ARGUMENTS",
        ),
        (
            "smarts_match",
            json!({"smarts":"[","smiles":"C"}),
            "INVALID_SMARTS",
        ),
        (
            "find_mcs",
            json!({"smiles_list":vec!["C";21]}),
            "INVALID_ARGUMENTS",
        ),
        (
            "find_mcs",
            json!({"smiles_list":["C".repeat(201),"C".to_string()]}),
            "MOLECULE_TOO_LARGE",
        ),
        (
            "name_to_smiles",
            json!({"name":"x".repeat(501)}),
            "INVALID_ARGUMENTS",
        ),
        (
            "retrosynthesis",
            json!({"smiles":"C".repeat(501)}),
            "MOLECULE_TOO_LARGE",
        ),
    ] {
        let error = tools::call_tool(name, &args).unwrap_err();
        assert_eq!(
            error.to_structured_error()["code"],
            code,
            "{name}: {error:?}"
        );
        assert!(!error.legacy_message().is_empty());
    }
    let result = tools::call_tool("find_mcs", &json!({"smiles_list":["[He]","CC"]})).unwrap();
    assert!(result["mcs"].is_null());
    assert_eq!(result["atom_count"], 0);
    for source in ["C=C", "C#C"] {
        assert!(tools::call_tool("find_mcs",&json!({"smiles_list":[source,source]})).unwrap()["bond_count"].as_u64().unwrap()>0);
    }
}

#[test]
fn http_reader_rejects_malformed_framing_before_dispatch() {
    for (request, status) in [
        ("POST /", 400),
        ("GET / HTTP/2.0\r\n\r\n", 400),
        ("GET / HTTP/1.1\r\nBroken\r\n\r\n", 400),
        ("GET / HTTP/1.1\r\nBad Name: x\r\n\r\n", 400),
        ("POST / HTTP/1.1\r\nTransfer-Encoding: chunked\r\n\r\n", 411),
        (
            "POST / HTTP/1.1\r\nContent-Length: 0\r\nContent-Length: 0\r\n\r\n",
            400,
        ),
        ("POST / HTTP/1.1\r\nContent-Length: bad\r\n\r\n", 400),
        ("POST / HTTP/1.1\r\nContent-Length: 3\r\n\r\nx", 400),
    ] {
        let result = http::read_request(&mut Cursor::new(request));
        assert_eq!(result.unwrap_err().status, status, "{request}");
    }
    let headers = (0..101)
        .map(|i| format!("X-{i}: x\r\n"))
        .collect::<String>();
    assert_eq!(
        http::read_request(&mut Cursor::new(format!("GET / HTTP/1.1\r\n{headers}\r\n")))
            .unwrap_err()
            .status,
        431
    );
    assert_eq!(
        http::read_request(&mut Cursor::new(format!(
            "GET / HTTP/1.1\r\nX: {}\r\n\r\n",
            "x".repeat(33000)
        )))
        .unwrap_err()
        .status,
        431
    );
    let response = http::HttpResponse {
        status: 400,
        headers: vec![],
        body: b"{}".to_vec(),
    };
    assert_eq!(response.json_body(), Some(json!({})));
    assert!(
        http::HttpServer::bind(http::HttpConfig {
            max_concurrency: 0,
            port: 0,
            ..Default::default()
        })
        .err()
        .unwrap()
        .to_string()
        .contains("max_concurrency")
    );
    let mut connection = transport::Connection::default();
    assert!(connection.handle_line(" \n").is_none());
    assert_eq!(
        connection.handle_line("{\"id\":1}").unwrap()["error"]["code"],
        protocol::INVALID_REQUEST
    );
}

struct FailedRead {
    prefix: Cursor<Vec<u8>>,
    kind: io::ErrorKind,
}
impl Read for FailedRead {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.prefix.position() == self.prefix.get_ref().len() as u64 {
            Err(io::Error::from(self.kind))
        } else {
            self.prefix.read(out)
        }
    }
}
#[test]
fn http_reader_distinguishes_timeout_from_short_body() {
    for (kind, status) in [
        (io::ErrorKind::TimedOut, 408),
        (io::ErrorKind::WouldBlock, 408),
        (io::ErrorKind::Other, 400),
    ] {
        let prefix = Cursor::new(b"POST / HTTP/1.1\r\nContent-Length: 4\r\n\r\nx".to_vec());
        let reader = FailedRead { prefix, kind };
        assert_eq!(
            http::read_request(&mut BufReader::with_capacity(1, reader))
                .unwrap_err()
                .status,
            status
        );
    }
}
