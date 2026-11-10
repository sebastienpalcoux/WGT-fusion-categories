# Rank bounds and proof calculations

The [calculation manifest](../certificates/manifest.json) gives hypotheses, inputs, digests, regeneration commands and terminal counts for each rank calculation R01–R12. The [current statement index](../certificates/statement-index.tsv) maps the permanent labels to the revised manuscript. See the [release guide](../certificates/README.md) for resources and the distinction between testing candidates and regenerating them. Run these commands from this directory.

| Category | Lower bound | Calculation |
|---|---:|---|
| Non-pointed, fusion-simple integral modular | 31 | [Uniform ordinary ledger](data/ordinary_release.json), R01 |
| Same, nonsquare dimension | 33 | [Rank 31](verify_nonsquare.py), [rank 32](verify_nonsquare32.py), R02–R03 |
| Same, odd dimension / odd nonsquare dimension | 57 / 75 | [Paired recursion](code/odd_frontier.py), R06 |
| Divisor-minimal slightly degenerate obstruction | 34 (full rank) | [Reduced-table check](verify_scoped.py), R05 |
| Weakly integral braided non-WGT, any / odd dimension | 9 / 41 | [Formal row test](code/formal_check.py), R08, and the odd twist-counting argument |
| Weakly integral non-WGT, any / odd dimension | 7 / 13 | [Small ranks](code/unbraided_small.py), [odd rank eleven](code/odd_unbraided.py), R09–R10 |

Strict weak integrality gives bounds 8 without braiding and 11 with braiding. Further fixed-dimension conclusions include 101 at the modular boundary 29,160,000 and 553 at the odd braided boundary 5⁶·7⁶·23⁴.

## Ordinary rank ledger

```sh
python3 ordinary_release.py
python3 ordinary_release.py --replay --ranks 22 23 24 25 26 27 --jobs 4
python3 ordinary_release.py --gmp --ranks 22 23 24 25 26 27
python3 ordinary_release.py --replay --ranks 28 29 30 --jobs 4
```

The same native and GMP kernels generate every rank from 6 to 30. Smaller ranks cannot satisfy the minimum-codegree inequality. At ranks 6–27 the manifest divides the search by its first denominator. Higher ranks use finer partitions. Rank 30 also uses the certified independent nonsquare partition and the explicit theorem leaves. The ordinary no-nine theorem discards the square-class-one root with first denominator 3 before recursion. Each surviving list is checked against the same [terminal screen](code/current_screen.py). The [ledger](data/ordinary_release.json) gives complete first-exclusion counts and the SHA-256 digest of the sorted lists at each rank. The separate GMP option compares a second, unpartitioned traversal with every row of the published small-rank return.

For reduced tables, paired odd ordinary tables, and non-simple categories use the applicable suite:

```sh
python3 verify_scoped.py reduced --replay
python3 verify_scoped.py paired-odd
python3 verify_scoped.py formal --replay
python3 verify_nonsquare.py --ranks 31 --replay
python3 verify_nonsquare32.py --replay
python3 square_slices.py
python3 code/q10_crosscheck.py
python3 code/unbraided_small.py
python3 code/odd_unbraided.py --verify data/odd_unbraided.json
```

## Corroboration and continuation

An earlier rank-28 generator returned 9,458 lists. Filtering these by the stronger prefix conditions gives exactly the 1,940 lists of the present enumeration, compared row for row in [verify_revision.py](verify_revision.py). The earlier rank-specific caps at codegree nine and the signed-Egyptian-sum check are retained for byte-for-byte regression of historical receipts. The uniform ordinary proof uses the no-nine theorem at the initial root.

The structural proof excluding fourteen denominator fours through rank 31 has a separate arithmetic check. Its 251 degree lists are excluded by generator activity (133 at three, 118 at five). The code, expected outputs and GMP replay are in [rank31_fourteen.py](rank31_fourteen.py). These corroborations do not enter the final lower bounds.

The [full continuation](frontier31.py) has 1,234 jobs. The seed contains 1,090, the fourteen-fours theorem closes 46, and 98 are unresolved. Each remaining job has square dimension, largest simple degree one quarter of the square root, and at most thirteen objects of this degree. Check these counts in a fresh directory with `python3 frontier31.py --prepare-only --work rank-work/seed-audit`. Start or resume with `python3 frontier31.py --jobs 4`. A completed run writes `rank-work/31/rank31_completed.zip`. The command `python3 frontier31.py --verify rank-work/31/rank31_completed.zip` checks this return.

The [229-page comparison index](supplement/BDP-statement-index-229p.tsv) is historical. The [statement index for this release](../certificates/statement-index.tsv) is generated directly from the final cross-references.
