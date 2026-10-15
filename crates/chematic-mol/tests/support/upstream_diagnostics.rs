/// Rejection fixtures also verify the public human-readable error contract.
pub fn diagnostic<T, E: core::fmt::Display>(result: Result<T, E>) -> Result<T, E> {
    result.inspect_err(|error| {
        assert!(
            !error.to_string().trim().is_empty(),
            "rejected input must provide a diagnostic to native and binding callers"
        );
    })
}
