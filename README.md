<a id="wgt-for-fusion-categories"></a>
# Weak group-theoreticity and rank bounds

Research on fusion categories by **Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux**. This repository brings together dimension theorems, rank bounds, broader research manuscripts, and the code and data accompanying each work.

**[All manuscripts](papers/README.md)** · [Dimension bounds](papers/dimension-bounds.md) · [Rank bounds](rank/README.md) · [Verification](VERIFICATION.md) · [Versions and citations](papers/VERSIONS.md)

## Where to start

| Question | Reading route |
| --- | --- |
| Below which dimensions is weak group-theoreticity proved? | The [29-page dimension theorem](papers/Below_5400_squared.pdf), followed by the [47-page three-prime paper](papers/three-prime-dimensions/three_prime_dimensions.pdf). Both proofs are computer-free. |
| How large must a non-pointed simple integral modular category be? | The [24-page general rank proof](papers/rank-bounds/general/Rank_bound.pdf), or the [37-page general and odd-dimensional proof](papers/rank-bounds/general-and-odd/Rank_bounds.pdf). |
| What are the wider results and remaining questions? | The [complete rank manuscript](papers/rank-bounds/full/Rank-modular-data.pdf), the [252-page project manuscript](papers/Manuscript.pdf), and the [three-prime research draft](papers/research-drafts/p6q6r4/p6q6r4_manuscript.pdf). |
| Which computations accompany a particular paper? | The [verification guide](VERIFICATION.md) maps each manuscript to its own programs, records and instructions. |

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
