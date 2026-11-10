# Verification of the rank bounds

These programs and data accompany *Prime factors of fusion-category dimensions: weak group-theoreticity and decomposition*. The enumeration constructs the degrees from the rank and the dimension equation. The published rank-15 list in the parent repository is an independent comparison.

| Category | Rank at least | Verification |
| --- | ---: | --- |
| Non-pointed simple integral modular | 31 | [Ranks through 29](verify.py), [complete rank 30](verify_rank30.py) |
| Same, nonsquare dimension | 33 | [Ranks through 31](verify_nonsquare.py), [rank 32](verify_nonsquare32.py) |
| Same, odd dimension / odd nonsquare dimension | 57 / 75 | [Odd-dimensional enumeration](code/odd_frontier.py) and [certificates](verify_revision.py) |
| Divisor-minimal slightly degenerate obstruction, full rank | 34 | [Fermionic quotient](verify.py) |
| Weakly integral braided non-WGT, any / odd dimension | 9 / 41 | [Small-rank check](code/formal_check.py), [proof by counting odd twists](../papers/Manuscript.pdf#page=127) |
| Weakly integral non-WGT, any / odd dimension | 7 / 13 | [Small ranks](code/unbraided_small.py), [odd rank eleven](code/odd_unbraided.py) |

For strictly weakly integral non-WGT categories the bounds are 8 without braiding and 11 with braiding. The proof uses the integral neutral component of the dimension grading. [Theorem 20.1](../papers/Manuscript.pdf#page=122)

At dimension 29,160,000 the modular obstruction has rank between 101 and 46,774. At dimension 5⁶·7⁶·23⁴ every braided non-WGT category has rank at least 553. All proper divisors are covered, and the forced twist orders give the bound at this dimension. [Theorem 26.1](../papers/Manuscript.pdf#page=218) · [Corollary 24.11](../papers/Manuscript.pdf#page=202)

## Checking completed ranks

Use Python 3.10 or later, a C++17 compiler with unsigned 128-bit integers and GNU MP development headers and libraries. Run from this directory.

```sh
python3 verify.py
python3 verify_after30.py
python3 code/unbraided_small.py
python3 code/odd_unbraided.py --verify data/odd_unbraided.json
```

The first two commands check the rank partitions and stored lists, repeat the exclusion tests and regenerate short component calculations. The [rank-30 return](data/revision/rank30_completed.zip) contains 96,119 different degree lists. The [rank-32 nonsquare manifest](data/revision/nonsquare/rank32.json) describes 14 complete roots and 3,081 lists. The [conditional square-dimension data](data/rank31/square_slices.json) give the bounds 36, 39 and 46 for the three hypotheses in the article.

To regenerate searches from their arithmetic inputs, run the following commands separately. They have different computational costs.

```sh
python3 verify.py --full
python3 verify_revision.py --replay --jobs 12
python3 verify_rank30.py --replay --jobs 12
python3 verify_nonsquare32.py --replay
python3 rank31_fourteen.py --replay
python3 square_slices.py
```

The [completed rank-28 archive](data/ordinary/rank28_completed.zip) and [its projection onto the optimized conditions](code/compare28.py) supply a separate comparison. The native and arbitrary-precision kernels can be compared on the same problems with [the kernel regression test](optimized/test_engine.py).

## Continuing at rank 31

```sh
python3 frontier31.py --jobs 12
```

The [bundled seed](data/rank31/rank31_seed.zip) covers 1,090 jobs. The [fourteen-fours certificate](rank31_fourteen.py) closes another 46. The remaining 98 have square dimension and largest degree one quarter of its square root. The same command resumes an interrupted run. A complete run writes `rank-work/31/rank31_completed.zip`. Check it with `python3 frontier31.py --verify rank-work/31/rank31_completed.zip`.

## Mathematical inputs

The [prime-family catalogue](code/catalogue.py), [complementary twist fibers](code/fiber_packets.py), [generator activity](code/generator_activity.py), [indicator divisibility](code/regular_indicator.py), and [combined exclusion test](code/current_screen.py) implement necessary conditions. The [historical statement index](supplement/BDP-statement-index.tsv) records labels, theorem numbers and pages of the [229-page coauthor edition](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/f2ef6eedb58ade6adcdb5ed9caeba1836b4ba8ce/papers/Manuscript.pdf). Its finite word and four-prime tests are available in the [main verification directory](../code/README.md#prime-factor-families). Earlier generator versions are retained where archived receipts use them.
