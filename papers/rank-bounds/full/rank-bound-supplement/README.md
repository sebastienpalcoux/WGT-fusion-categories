# Computational supplement

The rank censuses concern non-pointed based-simple integral rotation-admissible
data. The mathematical implications are proved in the manuscript; this folder
contains census inputs, code, finite arithmetic checks and execution records.
Section 10 also uses small finite counts for perfect categories.

## General rank: the fixed covering corpus

Table 1 counts the lists emitted by the source-bound producer archived in
`census/reference-census.zip`: 12,273 completed tasks, 41,915 exact
continuations and 18,199 distinct degree multisets.

| Rank | Square | Nonsquare | Mass | Gcd | Tiling | Block | Trace |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1–38 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 39 | 1 | 0 | 1 | 1 | 1 | 0 | 0 |
| 40 | 2 | 0 | 2 | 2 | 2 | 0 | 0 |
| 41 | 17 | 4 | 17 | 15 | 15 | 0 | 0 |
| 42 | 111 | 65 | 119 | 98 | 98 | 0 | 0 |
| 43 | 914 | 723 | 955 | 801 | 801 | 30 | 0 |
| 44 | 6269 | 10093 | 7113 | 5277 | 5249 | 214 | 0 |
| Total | 7314 | 10885 | 8207 | 6194 | 6166 | 244 | 0 |

The [terminal replay](census/five_tests.py) verifies the
[table](census/rank-counts.json) and
[first-exclusion ledger](census/rank-dispositions.jsonl.gz).
The [independent verifier](verification/independent/check.py) imports no
production modules and supplies an additional
[exact witness for every list](verification/independent/terminal-witnesses.jsonl.gz)
and [source-bound receipt](verification/independent/certificate.json).
The [seven-test comparison](census/recheck.py) excludes the same corpus using
separate degree-slice, overlap and fusion conditions.

The optimized generator has stronger necessary prefix tests and can emit
fewer lists than the source-bound producer. Table 1 always uses the fixed
covering corpus. A fresh optimized run has its own task and execution records;
its emitted-list counts need not coincide with the covering-corpus counts.
See the [census instructions](census/README.txt) for fresh generation.
The [rank-45 prefix replay](verification/sector_frontier/check.py) counts
possible starts at the next rank; it asserts no exclusion at rank 45.

## Odd dimension

The nonsquare [generator](odd_nonsquare/census.py) uses ternary rank ceiling
211. On fourteen supports the multiplicity caps and positive residual mass
permit 4,189,686 tuples. For `s` unfilled positions, `M_min = 3 product(P)`
and residual reciprocal-square mass `E`, the additional floor-15 bound
`E <= 1/M_min^2 + s/225` leaves 685,040 starting profiles. Four survive the
structural tests. Three fail ordinary mass; the last nine-pair equation has
no solution. The [production checker](odd_nonsquare/check.py) regenerates its
[receipt](odd_nonsquare/receipt.json), with 72 recursion nodes. The independent
implementation uses 172 memoized states for the same residual equation.
These are different recursions, not conflicting counts. The rank congruence
gives the nonsquare lower bound 219.

The square [profile generator](odd_square/profiles.py) at rank ceiling 209
produces 16,512 prime profiles and [1,151 edgeless plus 310 single-edge
reservations](odd_square/census.json). The [production checker](odd_square/check.py)
replays both complete recursions and compares all outcomes with the
[certificate](odd_square/certificate.json): 6,721,348 and 2,753,944 nodes,
respectively, with no retained reservation. The independent checker verifies
the initial profile count; it does not repeat these recursions. The
[stored-record audit](verification/odd_receipts/check.py) regenerates all
1,461 reservations, validates the outcome/node accounting and checks the
[source/result bindings](verification/odd_receipts/certificate.json), without
replaying either final odd continuation.

The [rank-217 frontier](verification/frontier217/check.py) leaves 227 prime
profiles and 1,628 graph reservations (1,133 edgeless and 495 single-edge).
This asserts neither realizability nor nonexistence at rank 217.

## Perfect categories through rank 44

The [Section 10 arithmetic](arithmetic/perfect_counts.py) checks the prime-power
lists, seven neutral-orbit patterns, reciprocal-square obstruction, cohomology
orbit arithmetic and deleted-gauging sector degrees. The independent verifier
also constructs the small permutation groups and all 59 subgroups of A5.
The categorical reductions, cohomology inputs and classification exhaustiveness
are mathematical arguments in Section 10; finite checks do not replace them.

## Commands and scope

Use Python 3.10 or newer. The default suite uses the standard library. The
independent suite additionally requires SymPy. Run from this directory:

```sh
# Bounded audit, including independent terminal and finite checks:
python3 verify.py --quick --independent

# Also exercise native prefix comparisons and forced-width handoffs:
python3 verify.py --quick --independent --native-prefix

# Replay both complete production odd certificates as well:
python3 verify.py --jobs 3
```

The native option requires g++ with C++17 and undefined-behavior sanitizer
support. Each selected native task and exact continuation has a five-second
guard; timeouts fail the check. The [boundary receipt](verification/native_boundary/certificate.json)
records equality of completed row sets at production and forced handoff
ceilings. Accounting tests reject truncated output, nonzero exits and
interrupted exact continuations.

The [short audit log](verification/revision-audit.log) records the bounded
suite. The supplied [full replay log](verification/replay.log) records the
complete odd recursions. `--independent-odd` adds separate rank-201 and rank-209
comparison checks; `--full` adds auxiliary arithmetic replays. The optional
Gauss backends are described in [native odd instructions](ODD-NATIVE.md).

A terminal replay does not regenerate the native task partition or certify
every discarded prefix independently. The archive audit checks sources,
declared-task completion and record accounting. Full mathematical coverage
also uses the necessary structural reductions and producer implementation,
indexed in the [proof map](proof-map.md). No source hash or execution receipt
is represented as a proof-assistant certificate. Resource limits and
interrupted continuations never certify a mathematical exclusion.
