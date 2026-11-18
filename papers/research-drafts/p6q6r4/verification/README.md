# Compact exact checks and source records

This companion deliberately contains no large row ledgers, compressed certificate
shards, generated field caches, or global remaining-state catalogue. The **total
size of all delivered supporting material**, excluding the manuscript TeX and
PDF, must be below 10,000,000 bytes. No large file has been divided into parts.
The complete theory remains in the manuscript.

## Run

From this directory, with Python 3.10 or later and assertions enabled:

```sh
python3 -B verify.py --output /tmp/p6q6r4-checks
```

The default bounded checks require NumPy and SymPy for two independent scripts.
No floating-point calculation or numerical optimizer is used: NumPy integer
arrays have a proven overflow bound; the dyadic matrices use arbitrary-precision
integer coefficients with exact common denominators. The exact versions used
are recorded in `reference/validation_summary.json`.

To run only checks using the Python standard library:

```sh
python3 -B verify.py standard-library --output /tmp/p6q6r4-standard
```

The separate small structural generators and their independently implemented
auditors can be run with:

```sh
python3 -B verify.py structural --output /tmp/p6q6r4-structural
```

Generated results and runtime copies remain outside this immutable directory.
Do not use `python -O` or `PYTHONOPTIMIZE`. No Internet access is needed.

## What is checked

| Scope | Files and exact result | Limit |
| --- | --- | --- |
| Certificate-set enlargement | `checks/check_certificate_scope.py`, `checks/exact_certificate_checks.py`: linear, bounded and modular witnesses valid on the old set are rejected on an enlarged set | Synthetic regression; not a coverage record for every historical certificate |
| First odd scalar frontier | `checks/first_odd_family_certificate.py`: 41 ordered triples, 40 excluded, explicit cyclic mask at (5,7,23) | Arithmetic necessary condition; no modular datum or category is constructed |
| Canonical length ten | `checks/verify_length10.py`: 26 profiles, 208 assignments, final skew alternatives excluded | Uses the stated palette, inherited source-rank bounds and categorical determinant constraints |
| Two degree-two blocks at rank fourteen | `checks/verify_two_blocks.py`: 32 placements, five determinant survivors, all five forbidden-pattern cases excluded | Only this algebra shape; does not eliminate the other rank-fourteen shape |
| Initial input counts and one explicit obstruction | `checks/verify_p6q6r4_review.py`: initial generating functions; thirteen necessary rows, four orthogonality edges, no triangle; three printed terminal square tests | No subsequent field/involution inventory, row coverage or uniqueness in the global catalogue |
| Dyadic compression | `checks/verify_dyadic_compression.py`: eight gamma values and seven nonidentity twists; projector, adjoint, polynomial and equality identities | Finite modulus-eight computation; not the preceding general representation-theoretic reduction |
| Coefficient controls | `checks/check_ch_pair_independent.py` and four small `controls/` pools; S3/Q8 identities also checked | Predicates on supplied pools, not independent completeness of those pools or new global exclusions |
| Small structural enumerations | `checks/structural/`: orthogonal codegree-six catalogue, cubic degree-four E3 and E4 branches, and symplectic profiles on at most seven slots, each with independent auditor | Exact displayed finite systems and profile norms 1 through 15, 18, 25 and 27, under their mathematical hypotheses; not other cubic branches |

`theorem_status.json` maps these checks to the manuscript labels and identifies
which large claims remain conditional on the manuscript's verification
hypothesis. Hashes authenticate bytes; they do not certify mathematical validity.

## Correct treatment of certificate enlargement

A certificate checked on a larger necessary row set also excludes every subset.
It cannot be carried from a subset to a larger set without rechecking all newly
admitted rows. Bounded separators must recompute the complete negative budget;
modular separators must check all rows with positive admissible bounds.

The current exact checker enforces a nonzero integral target and explicit
nonnegative-budget preconditions. The regression controls include zero targets,
unbounded zero rows and invalid negative budget coordinates. These safeguards
are local implementation tests. They do not demonstrate that all admissible
fields, involutions, repeated roots, quaternion alternatives or dimension ranges
were included in any historical global run.

## Preserved sources, not a full replay package

`historical_sources.zip` preserves all 520 original Python source files, byte for
byte (about 3.09 MB when expanded). `historical_sources.json` records their
original byte counts and SHA-256 hashes. Their historical directory names are
kept to preserve imports and provenance. The archive is a source record, not an
installation or one-command reproduction of the former large computations.
Many of these programs require data and certificates deliberately omitted here;
optional numerical proposal scripts may require additional packages. They are
not executed by the compact verifier. The revised exact separator checker is a
separate, hardened copy; the archived original remains unchanged.

In particular, this package does **not** validate the large rank exclusions,
canonical length-twelve/thirteen sieves, complete field/quaternion catalogues,
600,558-state inventory, disjointness of excluded sets, or global uniqueness
claims. Initial counts, small control examples, saved aggregate numbers and
source checksums cannot substitute for those missing verifications. These
computational conclusions are conditional in the revised manuscript. No
existence or categorification assertion is made.

## Integrity

```sh
python3 -B verify.py integrity --output /tmp/p6q6r4-integrity
sha256sum -c SHA256SUMS.txt
```

The integrity command also checks and parses all archived Python sources.
`SHA256SUMS.txt` covers every delivered file except itself. The small reference
outputs are comparisons for the bounded checks, not substitutes for running them.
