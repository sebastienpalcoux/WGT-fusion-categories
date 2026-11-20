# Structural inputs at rank 44

The covering census uses finite bounds on the odd-prime twists and dimension supports of degree lists with many copies of denominator four. The three calculations below give 39 possible mixed signatures, each with a nontrivial twist at both two and three. The second calculation eliminates the two further signatures with activity at two and five.

Run these short proofs from the supplement directory.

```sh
python3 arithmetic/isolation44.py
python3 arithmetic/two_five44.py
python3 arithmetic/high_four44.py
```

The [structural records](verification/structural/) contain bounded rank-specific outputs. The complete rank-44 census, including all-primary branches and both dimension types, is audited by [the Table 1 replay](census/recheck.py).
