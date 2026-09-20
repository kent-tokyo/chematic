// Explicit accounting for a known-length Explorer input.
//
// `unprocessedCount` is deliberately not called `skippedCount`: a client cap
// or cancellation means those rows were never given a terminal operation
// outcome.  Consumers must not infer batch success from `failedCount === 0`.

export function summarizeKnownLengthBatch({ inputCount, completedCount, terminalReason = null }) {
  for (const [name, value] of Object.entries({ inputCount, completedCount })) {
    if (!Number.isInteger(value) || value < 0) {
      throw new TypeError(`${name} must be a non-negative integer`);
    }
  }
  if (completedCount > inputCount) {
    throw new RangeError("completedCount cannot exceed inputCount");
  }

  const complete = completedCount === inputCount;
  if (complete && terminalReason !== null) {
    throw new RangeError("a complete batch cannot have a terminal reason");
  }
  if (!complete && terminalReason === null) {
    throw new RangeError("an incomplete batch must declare a terminal reason");
  }

  return {
    inputCount,
    completedCount,
    unprocessedCount: inputCount - completedCount,
    complete,
    terminalReason,
  };
}

export function formatBatchOutcome(outcome, { acceptedCount, rejectedCount }) {
  if (acceptedCount + rejectedCount !== outcome.completedCount) {
    throw new RangeError("accepted and rejected counts must account for completed records");
  }
  const loaded = `${acceptedCount} loaded, ${rejectedCount} failed to parse`;
  if (outcome.complete) return `${loaded}.`;
  return `${loaded}; ${outcome.unprocessedCount} not processed (${outcome.terminalReason}; complete=false).`;
}

// Streaming input has no trustworthy total until its reader reaches EOF.  Keep
// rows already observed separate from the unread suffix: an aborted Worker
// request may have observed rows for which no terminal result exists, while a
// cancelled reader can still contain an unknown number of later rows.
//
// This deliberately does not pretend that an aborted stream has `skipped`
// records. `unreadInput` is an explicit protocol value, not a count estimate.
export function summarizeStreamingBatch({
  observedInputCount,
  completedCount,
  streamComplete,
  terminalReason = null,
}) {
  for (const [name, value] of Object.entries({ observedInputCount, completedCount })) {
    if (!Number.isInteger(value) || value < 0) {
      throw new TypeError(`${name} must be a non-negative integer`);
    }
  }
  if (completedCount > observedInputCount) {
    throw new RangeError("completedCount cannot exceed observedInputCount");
  }
  if (typeof streamComplete !== "boolean") {
    throw new TypeError("streamComplete must be a boolean");
  }

  const complete = streamComplete && completedCount === observedInputCount;
  if (complete && terminalReason !== null) {
    throw new RangeError("a complete stream cannot have a terminal reason");
  }
  if (!complete && terminalReason === null) {
    throw new RangeError("an incomplete stream must declare a terminal reason");
  }

  return {
    inputKind: "unknown_length_stream",
    observedInputCount,
    completedCount,
    unprocessedObservedCount: observedInputCount - completedCount,
    unreadInput: streamComplete ? 0 : "unknown",
    complete,
    terminalReason,
  };
}

export function formatStreamingBatchOutcome(outcome, { acceptedCount, rejectedCount }) {
  if (acceptedCount + rejectedCount !== outcome.completedCount) {
    throw new RangeError("accepted and rejected counts must account for completed records");
  }
  const loaded = `${acceptedCount} loaded, ${rejectedCount} failed to parse`;
  // Keep the longstanding success wording: downstream smoke tests and, more
  // importantly, people watching a large local import can distinguish a
  // completed stream from a partial one at a glance.
  if (outcome.complete && rejectedCount === 0) {
    return `${acceptedCount} molecule${acceptedCount === 1 ? "" : "s"} loaded.`;
  }
  if (outcome.complete) return `${loaded}.`;
  const observed = outcome.unprocessedObservedCount === 0
    ? "no observed rows pending"
    : `${outcome.unprocessedObservedCount} observed row${outcome.unprocessedObservedCount === 1 ? "" : "s"} not processed`;
  return `${loaded}; ${observed}; unread input=${outcome.unreadInput} (${outcome.terminalReason}; complete=false).`;
}
