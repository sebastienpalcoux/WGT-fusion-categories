# Stored odd execution records

```sh
python3 verification/odd_receipts/check.py
```

This short check regenerates the rank-209 prime profile census and all 1,461
square reservations. It compares the order of those reservations with the
stored outcomes, checks every exclusion status and node total, and validates
the nonsquare receipt's internal accounting. `certificate.json` binds all
core `odd_square/` and `odd_nonsquare/` code and result files by SHA-256.
These 25 files are byte-identical to the supplied code and output.

The checker does not execute either final odd continuation or independently
prove the recorded exclusions. The exact executions are available through
`odd_square/check.py` and `odd_nonsquare/check.py`; their supplied complete run
is recorded in `verification/replay.log`. The separate independent finite
checker verifies the initial counts and final nonsquare equation by different
arithmetic algorithms. Hashes establish source/result identity, not execution
correctness. `--write` regenerates the binding/accounting receipt only.
