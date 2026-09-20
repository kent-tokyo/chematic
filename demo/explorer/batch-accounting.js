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
