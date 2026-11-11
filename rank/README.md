# Rank bounds and proof calculations

## The general rank bound

Every non-pointed simple integral modular category has rank at least 32. The 17-page paper gives a self-contained proof. Its supplement provides one verifier and one certificate for the entire range. [PDF](../papers/Rank32.pdf) · [TeX](../papers/Rank32.tex) · [Verification supplement](../papers/rank-bound-supplement/README.md)

## Further rank bounds from the full article

For each set of hypotheses, the table gives the lower bound and the computation used in the 252-page article. The [manifest](../certificates/manifest.json) records inputs, commands and terminal counts. The [statement index](../certificates/statement-index.tsv) gives the corresponding locations in the article.

| Category | Lower bound | Calculation |
|---|---:|---|
| Non-pointed simple integral modular, nonsquare dimension | 33 | [Rank 31](verify_nonsquare.py), [rank 32](verify_nonsquare32.py), R02–R03 |
| Non-pointed simple integral modular, odd dimension / odd nonsquare dimension | 57 / 75 | [Paired recursion](code/odd_frontier.py), R06 |
| Divisor-minimal slightly degenerate obstruction | 34 (full rank) | [Reduced-table check](verify_scoped.py), R05 |
| Weakly integral braided non-WGT, any / odd dimension | 9 / 41 | [Formal row test](code/formal_check.py), R08, and the odd twist-counting argument |
| Weakly integral non-WGT, any / odd dimension | 7 / 13 | [Small ranks](code/unbraided_small.py), [odd rank eleven](code/odd_unbraided.py), R09–R10 |

Strict weak integrality gives bounds 8 without braiding and 11 with braiding. Further fixed-dimension conclusions include 101 at the modular boundary 29,160,000 and 553 at the odd braided boundary 5⁶·7⁶·23⁴.

## Reproducing the full article

The [rank guide in the publication snapshot](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/rank/README.md) gives the original commands for R01–R12, including the ordinary rank ledger, the reduced and paired tables, and the historical cross-checks. The [release guide](../certificates/README.md) describes how to reproduce that edition.
