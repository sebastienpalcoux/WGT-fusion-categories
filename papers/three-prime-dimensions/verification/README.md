# Bounded arithmetic checks

The short paper has a self-contained, computer-free proof. These independently
written arithmetic checks accompanied the referee report and are supplied only
as supporting checks outside that proof. They do not certify the categorical
arguments or arbitrary-conductor statements.

Requires Python 3.10 or later and SymPy. From this directory, run:

```sh
python3 -B verify_three_prime_review.py --output /tmp/three-prime-review.json
```

The output must equal `review_checks.json` after JSON parsing. The checks cover
the scalar endpoint inequalities, bounded cyclotomic words, all four-position
supports containing zero at prime level 23, dyadic identities, and small modular
examples. They do not classify all modular data or construct a category at either
endpoint. The calculation was rerun on 17 November 2026 and exactly reproduced
the supplied output.

The original report's result-by-result audit register is an assessment of the
source it reviewed, not an executable proof certificate for the revised paper.
It is therefore not used as evidence of a full proof certification here.
