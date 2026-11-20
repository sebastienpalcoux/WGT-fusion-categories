GENERAL RANK BOUND AND TABLE 1

The reference census covers nonpointed based-simple integral, rotation-
admissible modular data in the manuscript, through rank 44. It holds all
12,273 completed tasks, 41,915 completed exact handoffs, 18,199 degree
multisets (7,314 square and 10,885 nonsquare), and the generation sources.
The five terminal stages eliminate all these multisets.

REPRODUCE THE TABLE

Run from this directory:

    python3 five_tests.py

This verifies the rank-counts.json table and rank-dispositions.jsonl.gz
ledger against the lists in reference-census.zip. The primary mass,
prime-square gcd, defect-one tiling, principal-block and ternary trace
conditions leave 8,207, 6,194, 6,166, 244 and 0 lists respectively.

SEVEN-TEST COMPARISON

    python3 recheck.py

The program checks the archived coverage and completion records, validates
the degree lists, and evaluates the seven comparison stages. The source bindings, full per-rank table and hash of the uncompressed
row ledger are in ../verification/rank_tables/table1.json. Each entry in
../verification/rank_tables/dispositions.jsonl.gz gives the rank, task index,
reference-row digest, and first exclusion. The cumulative numbers
are 393 after degree slices, 284 after principal overlap and 215 after the
fusion product test. Four-prime isolation and the prime-square tests leave
143, principal blocks leave
124, simultaneous block traces leave 12 and defect-one tilings leave zero.
Counts use the first applicable obstruction in exactly this order.

Run the independent checks from the supplement's parent code directory:

    cd ..
    python3 verify.py --quick --independent

This also compares both degree-slice decisions at every archived profile,
rederives principal-overlap and fusion-product rejection certificates,
checks finite odd starting counts and the tail, and exercises accounting
boundaries. The independent checker requires SymPy. The full production odd
recursions remain available with python3 verify.py. To run a
single independent general test from the supplement root:

    python3 verification/rank_tables/audit_slices.py census/reference-census.zip
    python3 verification/rank_tables/audit_overlap.py census/reference-census.zip
    python3 verification/rank_tables/audit_fusion.py census/reference-census.zip

A general terminal replay verifies the complete emitted records without
rerunning native generation. Independent tests consume the same reference
archive. They do not infer coverage outside its declared rank interval.
Optional resource limits on individual mathematical tests retain candidates.

FRESH ENUMERATION

From this census directory, on a system with Python 3.9+ and a C++17 compiler:

    python3 run.py --through 44 --jobs 8 --work fresh-all-ranks

For rank 44 only:

    python3 run.py --rank 44 --jobs 8 --work fresh-rank44

The fresh enumerator rebuilds the native task partitions and emits every
complete multiset reaching its current terminal classifier, including
exact Python continuations of fixed-width native cutoffs. The immutable
covering archive binds its own producer sources. The optimized driver has
stronger necessary prefix rules, so a fresh run may emit fewer multisets.
Table 1 always uses the fixed covering corpus, not the optimized output.
A fresh optimized run has its own task and execution records; its emitted-list
counts need not coincide with the covering-corpus counts. Use the same --rank
or --through and --work selection to resume or verify a fresh work directory:

    python3 run.py --through 44 --jobs 8 --work fresh-all-ranks --resume
    python3 run.py --through 44 --work fresh-all-ranks --verify --recheck

The output file is fresh-all-ranks/rank-results.zip, or the corresponding
rank-only result. A completed mathematical exclusion requires completed
tasks and exact continuations. The run and verify programs report unfinished
work or inconsistent records as errors. Run proof code without -O or -OO.

INDEPENDENT TERMINAL WITNESSES

    python3 ../verification/independent/check.py

This imports no production code and regenerates every Table 1 cell and an
exact witness for each emitted row. The receipt and witness schema are in
../verification/independent. These witnesses certify terminal arithmetic on
the fixed corpus; they do not independently certify native prefix coverage.
Task/receipt consistency, bounded native comparisons and mathematical
necessity of the pruning predicates are separate parts of verification.
