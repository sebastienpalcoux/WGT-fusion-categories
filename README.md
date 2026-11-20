<a id="wgt-for-fusion-categories"></a>
# Weak group-theoreticity and rank bounds

Research on fusion categories by **Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux**. This repository brings together dimension theorems, rank bounds, broader research manuscripts, and the code and data accompanying each work.

**[All manuscripts](papers/README.md)** · [Dimension bounds](papers/dimension-bounds.md) · [Rank bounds](rank/README.md) · [Verification](VERIFICATION.md) · [Versions and citations](papers/VERSIONS.md)

## Where to start

The project asks two related questions: **how small can a weakly integral fusion category be if it is not weakly group-theoretical, and how few simple objects can a non-pointed simple integral modular category have?** Here dimension means Frobenius–Perron dimension, and rank counts simple objects up to isomorphism. All categories are over ℂ.

### Dimension bounds: begin with a uniform theorem

Start with the **[29-page dimension paper](papers/Below_5400_squared.pdf)**. It proves that every weakly integral braided fusion category of dimension **strictly below 5400² = 29,160,000** is weakly group-theoretical. Passing to the Drinfeld center gives **dimension below 5400 without assuming a braiding**. The proof is computer-free, with the arithmetic included in the text.

Continue with the **[47-page three-prime paper](papers/three-prime-dimensions/three_prime_dimensions.pdf)** to study how prime factorization constrains a possible counterexample. For distinct primes $p,q,r$, it proves weak group-theoreticity for braided categories of dimension a proper divisor of $p^6q^6r^4$; within this family, **5400² is the sole possible even-dimensional exception**. In the odd-dimensional part of the same family, it proves solvability strictly below $5^6\cdot7^6\cdot23^4$. This proof is also computer-free. Both endpoints remain open.

### Rank bounds: choose the scope you need

The rank bounds below apply to **non-pointed simple integral modular fusion categories**. The existence of such a category remains open.

| Your aim | Recommended reading | What it adds |
| --- | --- | --- |
| Understand the general bound in a focused paper | **[24 pages](papers/rank-bounds/general/Rank_bound.pdf)** | Proves **rank at least 45**, using theoretical restrictions followed by an exact finite enumeration. |
| Study both the general and odd-dimensional cases | **[37 pages](papers/rank-bounds/general-and-odd/Rank_bounds.pdf)** | Includes the general bound **45** and proves **217 in odd dimension**. Start here directly if you want both cases. |
| Study the wider rank theory and the classification results | **[Revised complete manuscript, 98 pages](papers/rank-bounds/full/Rank-modular-data.pdf)** | Includes both bounds, **219 in odd nonsquare dimension**, and the classification of **nontrivial perfect weakly integral modular categories through rank 44** in Section 10. |

For the shorter earlier argument, the **[original 17-page paper](papers/Rank32.pdf)** proves rank at least **32**. The new papers strengthen this to **45**. The classification in the complete manuscript concerns *perfect* categories, whose only invertible simple object is the unit; its hypotheses differ from those of the simple-category rank theorem.

### Broader sources and open problems

- **[252-page project manuscript](papers/Manuscript.pdf):** consult this for the broader prime-exponent families, dimension bounds, earlier rank results and decomposition theorems. It is the work recorded as **Hexagon 2611.00233** and supplies the numbered earlier results cited by the later papers.
- **[423-page three-prime research draft](papers/research-drafts/p6q6r4/p6q6r4_manuscript.pdf):** continue here from the 47-page paper for the detailed investigation of $p^6q^6r^4$ and the unresolved dimension **5400²**. It separates proved theoretical restrictions from stronger numerical conclusions requiring additional computational coverage.

### Checking a proof or following a reference

The **[verification guide](VERIFICATION.md)** identifies each paper's companion and explains which commands check stored results and which regenerate a census. The rank papers use separate editions of their programs and records. The **[versions and citations guide](papers/VERSIONS.md)** identifies the Hexagon edition and fixed code snapshots used in their references.

<a id="short-self-contained-proofs"></a>
## Dimension bounds

| Manuscript | Result and scope | Files |
| --- | --- | --- |
| **Weakly integral braided fusion categories below dimension 5400²** · 29 pages | Weak group-theoreticity for dimension below $5400^2=29,160,000$. Passing to the Drinfeld center gives the bound 5400 without braiding. | [PDF](papers/Below_5400_squared.pdf) · [TeX](papers/Below_5400_squared.tex) |
| **Braided fusion categories with dimension dividing p⁶q⁶r⁴** · 47 pages | Weak group-theoreticity for proper divisors of the family dimension, with restrictions on a possible obstruction at the full dimension. | [PDF](papers/three-prime-dimensions/three_prime_dimensions.pdf) · [TeX](papers/three-prime-dimensions/three_prime_dimensions.tex) · [Guide](papers/three-prime-dimensions/README.md) |

The [dimension guide](papers/dimension-bounds.md) explains the relation between these focused papers, the broader project, and the open $5400^2$ endpoint.

## Rank bounds

**The three new manuscripts are major improvements of the 17-page rank-32 paper: the general lower bound rises from 32 to 45.** The two versions treating odd dimension prove a lower bound of **217**, improving the bound of 57 in the 252-page Hexagon manuscript.

These category bounds concern **non-pointed simple integral modular fusion categories over the complex numbers**. The new papers prove more general algebraic statements under the precise Drinfeld-duality and modular-data hypotheses given in their theorems.

| Manuscript | Result and scope | Files |
| --- | --- | --- |
| **A rank bound for simple integral modular data** · 24 pages | The focused proof of rank at least **45**. | [PDF](papers/rank-bounds/general/Rank_bound.pdf) · [TeX](papers/rank-bounds/general/Rank_bound.tex) · [Guide](papers/rank-bounds/general/README.md) |
| **Rank bounds for simple integral modular data** · 37 pages | The focused proof of rank at least **45**, and **217** in odd dimension. | [PDF](papers/rank-bounds/general-and-odd/Rank_bounds.pdf) · [TeX](papers/rank-bounds/general-and-odd/Rank_bounds.tex) · [Guide](papers/rank-bounds/general-and-odd/README.md) |
| **The rank of simple integral modular data** · revised complete version, 98 pages | The bounds **45 / 217 / 219**, the last for odd nonsquare dimension; a classification of nontrivial perfect weakly integral modular categories through rank 44; further theoretical restrictions. | [PDF](papers/rank-bounds/full/Rank-modular-data.pdf) · [TeX](papers/rank-bounds/full/Rank-modular-data.tex) · [Guide](papers/rank-bounds/full/README.md) |
| **A rank bound for simple integral modular categories** · 17 pages | The original proof of rank at least **32**, retained as a separate manuscript with its original companion. | [PDF](papers/Rank32.pdf) · [TeX](papers/Rank32.tex) · [Guide](papers/rank-bound-supplement/README.md) |

[Compare the three new manuscripts](papers/rank-bounds/README.md) · [Rank results and verification](rank/README.md)

## Broader research manuscripts

| Manuscript | Role | Files |
| --- | --- | --- |
| **Prime factors of fusion-category dimensions: weak group-theoreticity, rank and decomposition** · 252 pages | The broad project manuscript: dimension families, the earlier rank bounds, and decomposition results. This is the work recorded as **Hexagon 2611.00233**. | [PDF](papers/Manuscript.pdf) · [TeX](papers/Manuscript.tex) · [Hexagon](https://hexagonmath.org/2611.00233) · [Verification](certificates/README.md) |
| **Arithmetic restrictions on braided fusion categories of dimension p⁶q⁶r⁴** · 423 pages | The detailed three-prime investigation, with proved restrictions and an explicit account of which stronger numerical conclusions remain conditional. | [PDF](papers/research-drafts/p6q6r4/p6q6r4_manuscript.pdf) · [TeX](papers/research-drafts/p6q6r4/p6q6r4_manuscript.tex) · [Guide](papers/research-drafts/p6q6r4/README.md) |

<a id="three-prime-arithmetic-frontier"></a>
## Research approach and open frontiers

For the rank-bound investigation, as for the dimension-bound investigation, we tried to push the AI model to its limits through repeated attempts to strengthen the theory and improve the bounds. The focused papers make individual results easier to study; the larger manuscripts record the broader arguments and the frontiers reached.

The dimension $5400^2$ remains unresolved. The new rank work also leaves rank **45**, and odd rank **217**, open. Surviving arithmetic configurations are necessary data for further investigation; they do not establish the existence of a category. The [research guide](papers/research-drafts/README.md) and the individual companions describe the scope of the completed calculations.

## Code, data and references

Start with the **[verification guide](VERIFICATION.md)** to select the companion for the paper you are reading. The [original arithmetic programs](code/README.md), [rank-15 input data](data/README.md), and [252-page manuscript's manifest](certificates/manifest.json) retain their own documented scope. The [new rank companions](papers/rank-bounds/README.md#verification) are separate editions.

The **[version and citation guide](papers/VERSIONS.md)** distinguishes the Hexagon manuscript, the focused papers, and the fixed code snapshots used in existing references.
