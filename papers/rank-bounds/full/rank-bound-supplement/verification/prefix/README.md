# Prefix-capacity verification

The production cut uses Corollary `new:prime-free-square-trace`. At a fixed
three-prime support with squarefree kernel3, a completed denominator-p class
cannot contain an equal ternary block, for either p=2 or the third prime.
Every member therefore needs two degree-half partners (denominator2p).
A prefix can be rejected when its current p-count exceeds half the largest
possible final2p-count. The check includes every unfinished partner slot,
and stops adding capacity only after the ordered search has passed2p.

Run `python3 verification/prefix/check.py` from the supplement directory.
It extracts the reference sources to a temporary directory, compiles the producer and optimized
native programs, compares four small tasks, and verifies every removed
complete row with the current trace or prime-cut condition. Eleven boundary
fixtures and 5,756 finite-completion implications check that unfinished
partner classes remain conservative. Two arbitrary-precision continuations
compare the reference mass helper with the current helper, which also uses
the prime-cut mass reservation. Each task and continuation has a five-second
guard; a timeout fails the check. No rank census is run.

The [bounded comparison receipt](current-prefix-check.json) binds the
checker, both producer/optimized mass helpers and both terminal conditions.
It records 11 boundary cases, 5,756 partner-capacity implications, four native
task comparisons and two arbitrary-precision continuations. It does not
assert equivalence of every discarded prefix.

The prime-cut mass theorem has a separate small-completion audit in
[test_prime_cut.py](../../census/tests/test_prime_cut.py). Its 12,110
completion-prefix projections include odd dual-pair multiplicities, open
degree classes, and allowed supports containing primes that need not occur.
Run it through `python3 -m unittest discover -s census/tests` from the
supplement directory, or run the default `python3 verify.py` audit. The
12,110 prime-cut projections are not the 5,756 implications recorded by
the native partner-capacity comparison.
