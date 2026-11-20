# Dimension bounds

[Overview](../README.md) · [All manuscripts](README.md) · [Rank bounds](../rank/README.md) · [Verification](../VERIFICATION.md) · [Versions and citations](VERSIONS.md)

The dimension investigation asks how small a weakly integral fusion category can be if it is not weakly group-theoretical. The repository offers two focused computer-free papers and two broader sources.

## Focused proofs

### The dimension bound 5400²

**Weakly integral braided fusion categories below dimension 5400²** — 29 pages. Every weakly integral braided fusion category of Frobenius–Perron dimension strictly below $5400^2=29,160,000$ is weakly group-theoretical. Applying the theorem to the Drinfeld center gives the bound **5400** without a braiding assumption.

[PDF](Below_5400_squared.pdf) · [TeX source](Below_5400_squared.tex)

### The three-prime family

**Braided fusion categories with dimension dividing p⁶q⁶r⁴** — 47 pages. This paper proves weak group-theoreticity at proper divisors and gives necessary conditions for an obstruction at the full dimension. In the even case the only possible endpoint within this family is $2^6 3^6 5^4=5400^2$. The full-dimension endpoint remains unresolved. Within the same three-prime family, odd-dimensional categories below $5^6 7^6 23^4=514,422,090,765,625$ are solvable; that endpoint also remains unresolved.

[PDF](three-prime-dimensions/three_prime_dimensions.pdf) · [TeX source](three-prime-dimensions/three_prime_dimensions.tex) · [Paper guide](three-prime-dimensions/README.md) · [Supporting checks](three-prime-dimensions/verification/README.md)

Both focused proofs are computer-free. The bounded arithmetic checks attached to the three-prime paper are supplementary checks outside its proof.

## Broader sources

| Source | How it complements the focused papers |
| --- | --- |
| [252-page project manuscript](Manuscript.pdf) · [TeX](Manuscript.tex) · [Hexagon](https://hexagonmath.org/2611.00233) | A broader catalogue of prime-exponent families, square-full dimension bounds, rank results and decomposition theorems. Its original computations have a [dedicated guide](../certificates/README.md). |
| [423-page three-prime research draft](research-drafts/p6q6r4/p6q6r4_manuscript.pdf) · [TeX](research-drafts/p6q6r4/p6q6r4_manuscript.tex) | The detailed investigation of $p^6q^6r^4$, including theoretical restrictions at the unresolved endpoint and numerical claims whose additional coverage requirements are stated explicitly. Its [compact companion](research-drafts/p6q6r4/README.md) reproduces specified bounded checks. |

The 47-page paper presents results from the joint Burciu–Dong–Palcoux project in a shorter, self-contained form. The [research guide](research-drafts/README.md) explains the role of the longer texts.

## The frontier

The investigation repeatedly sought stronger theoretical restrictions, with the aim of pushing the AI model as far as possible on the dimension problem. The current endpoint is still **5400²**: none of the linked manuscripts claims to settle it or to construct a non-WGT category there. The later [rank investigation](../rank/README.md) follows the same approach and substantially improves the earlier rank bounds.

For exact theorem and page references, consult the [version guide](VERSIONS.md). The [verification map](../VERIFICATION.md) separates the original project's computations from the compact three-prime checks.
