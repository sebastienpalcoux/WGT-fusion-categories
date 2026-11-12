# Rank bounds and proof calculations

## The general rank bound

Every non-pointed simple integral modular category has rank at least 32. The [full article](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/main/papers/Manuscript.pdf) and the [17-page self-contained paper](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/main/papers/Rank32.pdf) use the same degree equation, local character restrictions and finite verification. The full article also enables its proved dimension-family exclusions within that search. The self-contained paper supplies a complete proof using the same program with those additional filters disabled.

There is one [optimized verifier](../papers/rank-bound-supplement/verify.py), one pair of rank-parameterized generators, and one [certificate](../papers/rank-bound-supplement/certificate-families.json.gz) for the entire range proved in the full article. Ranks 1–9 are excluded analytically. The default computation also checks rank 9 and exhausts every rank through 31. From the repository root, check the published record with

```sh
python3 papers/rank-bound-supplement/verify.py --families --check papers/rank-bound-supplement/certificate-families.json.gz
```

To regenerate every list with four workers and a fresh checkpoint directory, run

```sh
python3 papers/rank-bound-supplement/verify.py --families --jobs 4 --work rank-results
```

The full-article run has 3,012 tasks, 411,325 square-dimension profiles and 122 nonsquare-dimension profiles, with no survivors. The [supplement guide](../papers/rank-bound-supplement/README.md) explains the records, resuming a run, and choosing another rank. For example,

```sh
python3 papers/rank-bound-supplement/verify.py --families --rank 32 --jobs 4 --work rank-32
```

The same recursions accept any positive rank. Bounds whose proof has an upper-rank or support hypothesis are applied only under that hypothesis. Higher-rank computations are exploratory, and the published certificate covers the theorem's range.

## Further rank bounds from the full article

The following computations have additional hypotheses. Their [manifest entries](../certificates/manifest.json) specify the input files, commands and outcomes. The [statement index](../certificates/statement-index.tsv) identifies their locations in the article.

| Category | Lower bound | Calculation |
|---|---:|---|
| Non-pointed simple integral modular, nonsquare dimension | 33 | Unified calculation R01 through 31; [nonsquare rank-32 check](verify_nonsquare32.py), R03 |
| Non-pointed simple integral modular, odd dimension / odd nonsquare dimension | 57 / 75 | [Paired recursion](code/odd_frontier.py), R06 |
| Square dimension, largest degree smaller than √D/4 / square root not divisible by four / largest degree smaller than √D/5 | 36 / 39 / 46 | [Conditional square searches](square_slices.py), R04 |
| Divisor-minimal slightly degenerate obstruction | 34 (full rank) | [Reduced-table check](verify_scoped.py), R05 |
| Weakly integral braided non-WGT, any / odd dimension | 9 / 41 | [Formal row test](code/formal_check.py), R08, and the odd twist-counting argument |
| Weakly integral non-WGT, any / odd dimension | 7 / 13 | [Small ranks](code/unbraided_small.py), [odd rank eleven](code/odd_unbraided.py), R09–R10 |

Strict weak integrality gives bounds 8 without braiding and 11 with braiding. The fixed-dimension conclusions include rank at least 101 at the modular boundary 29,160,000 and at least 553 at the odd braided boundary 5⁶·7⁶·23⁴.

The shared code under `current/` and `code/` supplies the conditional square, fermion-quotient and odd-dimension calculations. Earlier general searches and the former continuation records are available in the [publication snapshot](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/independent-referee-2026-11-10). The [release guide](../certificates/README.md) describes the current inputs and checks.
