# Arithmetic verification

Run the commands below from this directory. The [manifest](../certificates/manifest.json) gives the hypotheses and expected output of the 27 dimension, fixed-dimension rank and comparison checks, with their permanent labels. The [statement index](../certificates/statement-index.tsv) supplies current printed numbers and pages in the [full article](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.pdf). The two historical comparisons link to the original coauthor edition. The [general rank guide](../rank/README.md) gives the separate commands for unrestricted and odd-dimensional ranks. The programs in this directory use Python 3.10 or later and its standard library. The largest support count uses a C++20 compiler with unsigned 128-bit integer support.

## Dimensions and ranks

| Command | What it checks |
| --- | --- |
| `python3 check_dimension_bounds.py` | The 16 entries in the proof of the four-dimension bound and their four column minima. See [the dimension comparison](#dimension-comparison) below. |
| `python3 check_frontier_rank.py` | The two boundary recursions for ranks 95 through 100 in the one-neutral case at dimension 29160000 with two and five separated. They leave five rank-100 profiles, excluded by the tame-phase argument in the article. |
| `python3 verify_rank_5400.py` | The character-block exclusions at rank 10 in the commutative case and rank 12 in the noncommutative case for simple perfect categories under the dimension-5400 hypotheses of the unbraided boundary theorem. |
| `python3 verify_rank15_published.py` | The smallest formal codegree for each of the 2341 [rank-15 dimension types](../data/published-rank15-perfect.txt) from the [pinned classification](../data/README.md). |
| `python3 verify_small_ranks.py` | The eight rank-six dimensions at or above 5400 and the reciprocal tuples at odd ranks nine and eleven in [Corollary 18.2 of the 229-page coauthor edition](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/f2ef6eedb58ade6adcdb5ed9caeba1836b4ba8ce/papers/Manuscript.pdf#page=112). |
| `python3 verify_odd13.py` | The 520 odd reciprocal tuples at rank 13, the 17 passing the cutoff and codegree test, and the two requiring a character calculation modulo seven. |
| `python3 verify_witnesses.py` | The individual witnesses in [certificates.json](../certificates/certificates.json), including both disjoint fibers and their norm vectors. |

### Dimension comparison

[The dimension checker](check_dimension_bounds.py) implements the sufficient exponent bounds used in the proof of the four-dimension bound. Dong’s theorem covers at most two repeated primes. For each support size from three through six, the program tests the smallest primes in increasing order, first including 2 and then using only odd primes. Replacing a support by the smallest primes of the same parity preserves the exponent tests, including the position of 2, and can only decrease the product. With six or more primes, the product of the first six prime powers already gives the stated lower bound.

The priority queue visits exponent vectors in order of increasing product. For the braided bounds it starts at exponent two at every prime. For the bounds without braiding it starts at exponent one in the category dimension and doubles the exponents before applying the same tests to its center. The program checks every entry of the four-row table against the proof and exits with an error on any discrepancy. It prints the unrestricted and odd column minima on each line, first with braiding and then without it.

```text
29160000 6087156075
5400 121275
```

## Prime-factor families

The [four-prime guide](four_primes/README.md) gives the commands for the finite trace windows and the periodic projections with an unbounded opposite exponent. Surviving projection classes satisfy necessary conditions only.

For small exposed difference sets and the exceptional four-position support, run

```sh
python3 exposed_support.py
python3 -c 'from exposed_support import words; print(words(29,4,3))'
python3 -c 'from exposed_support import words; print([(s,words(s,6,2)) for s in (83,89,97,101,103,107)])'
python3 word_closure.py check 23 4 --ell 3 --classes
```

For the residual tests and the larger normalized-support counts, run

```sh
python3 short_words.py
c++ -O3 -std=c++20 strip_extra.cpp -o strip_extra
./strip_extra
```
