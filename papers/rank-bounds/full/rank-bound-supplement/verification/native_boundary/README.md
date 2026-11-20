# Native arithmetic boundary checks

From the supplement root:

```sh
python3 verification/native_boundary/check.py
```

This requires Python 3.10+, Linux signal support and g++ with C++17 and
undefined-behavior sanitizer support. The script compiles each native backend
with its production lcm ceiling and with a deliberately small ceiling. It
uses two positive near-terminal input rows from the covering census, at square
rank 39 and nonsquare rank 41. Their long fixed prefixes force a handoff before
the assigned prefix ends, exercising the rule that the entire assigned prefix
must reach the Python continuation.

Each native task and exact continuation has a five-second limit. A timeout is
a failed check. Each compiler invocation has a 30-second limit. No rank census
or full partition is generated. Every handoff is exhausted; both settings
produce the same nonempty completed row set, and each resulting degree list
fails a necessary terminal test. `certificate.json` records the fixtures,
counts, completion trailers, source hashes and equality checks. `--write`
regenerates this receipt after an intentional, reviewed source change.

The source also guards lcm multiplication before entering a fixed-width
state. At rank at most 44, a native state has `L <= 2^58`; the interval bound
limits the next denominator by `sqrt(h+1/36)L < 7L`. Nonsquare states additionally
check `(36h+1)L^2` against the unsigned 128-bit range and hand off before the
terminal scale could exceed the 64-bit range. The square backend carries the
terminal scale and emitted degrees in 128 bits. In both backends the exact
handoff preserves the full denominator prefix (and the nonsquare support).

These tests exercise specific transitions, not every arithmetic operation or
every discarded prefix. Completion accounting is tested separately in
`census/tests/test_driver_accounting.py`, including malformed/truncated output,
a nonzero native exit despite a matching trailer, and an exact generator
interrupted after emitting a row. None may publish a completed task ledger.
