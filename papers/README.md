# Manuscripts

## Short, self-contained proofs

These two papers present specific conclusions of the larger project with the hypotheses, arguments and arithmetic needed to read each proof on its own. They make the results easier to study and check without working through the full article.

- **A rank bound for simple integral modular categories** (21 pages). Every non-pointed simple integral modular category has rank at least 32. [PDF](Rank32.pdf) · [TeX](Rank32.tex) · [Unified verification supplement](rank-bound-supplement/README.md)
- **Weakly integral braided fusion categories below dimension $5400^2$** (29 pages). Such categories are weakly group-theoretical, and the Drinfeld-center argument gives a bound of 5400 without braiding. [PDF](Below_5400_squared.pdf) · [TeX](Below_5400_squared.tex)

For the rank paper, run `python3 papers/rank-bound-supplement/verify.py --jobs 4` from the repository root. The [rank calculations](../rank/README.md) and [full article](Manuscript.pdf) address the broader results. The [earlier 17-page version](Rank31.pdf) proved rank at least 31.
