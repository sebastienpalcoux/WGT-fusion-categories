# Arithmetic verification

The programs use Python 3.10 or later and its standard library. The largest support count uses a C++20 compiler with unsigned 128-bit integer support. Run the commands below from this directory. The [manifest](../certificates/manifest.json) records the hypotheses, manuscript labels and expected output for each of the 24 checks.

## Dimensions and ranks

| Command | What it checks |
| --- | --- |
| `python3 check_dimension_bounds.py` | The four ordered-exponent minima. The article also derives them by hand. |
| `python3 check_frontier_rank.py` | The two rank recursions at dimension 29160000 for ranks 95 through 100. They leave five profiles at rank 100 for the tame-phase argument. |
| `python3 verify_rank_5400.py` | The commutative and noncommutative character profiles at the two lower rank endpoints for simple perfect categories of dimension 5400. |
| `python3 verify_rank15_published.py` | The tally of the 2341 [rank-15 dimension types](../data/published-rank15-perfect.txt) from the [pinned classification](../data/README.md). |
| `python3 verify_odd13.py` | The 520 odd reciprocal tuples at rank 13, the 17 passing the cutoff and codegree test, and the two requiring a character calculation modulo seven. |
| `python3 verify_witnesses.py` | The individual witnesses in [certificates.json](../certificates/certificates.json), including both disjoint fibers and their norm vectors. |

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
