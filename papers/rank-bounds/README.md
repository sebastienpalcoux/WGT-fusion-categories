# Improved rank bounds

[Overview](../../README.md) · [All manuscripts](../README.md) · [Dimension bounds](../dimension-bounds.md) · [Rank guide](../../rank/README.md) · [Verification](../../VERIFICATION.md) · [Versions and citations](../VERSIONS.md)

These three manuscripts by **Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux** substantially improve the original [17-page rank-32 paper](../Rank32.pdf). They offer three scopes of the new rank work, with a separate matching companion for each.

## Choose a manuscript

| Edition | Pages | What it contains | Files |
| --- | ---: | --- | --- |
| **General bound** — *A rank bound for simple integral modular data* | 24 | A focused proof of rank at least **45**. | [PDF](general/Rank_bound.pdf) · [TeX](general/Rank_bound.tex) · [Guide](general/README.md) |
| **General and odd bounds** — *Rank bounds for simple integral modular data* | 37 | A focused proof of rank at least **45**, and at least **217** in odd dimension. | [PDF](general-and-odd/Rank_bounds.pdf) · [TeX](general-and-odd/Rank_bounds.tex) · [Guide](general-and-odd/README.md) |
| **Revised complete manuscript** — *The rank of simple integral modular data* | 98 | Bounds **45 / 217 / 219**, the last in odd nonsquare dimension; broader theoretical restrictions; perfect-category classification through rank 44. | [PDF](full/Rank-modular-data.pdf) · [TeX](full/Rank-modular-data.tex) · [Guide](full/README.md) |

The category bounds concern **non-pointed simple integral modular fusion categories over ℂ**. The main algebraic theorems apply to non-pointed based-simple commutative integral fusion rings with a symmetric algebraic Drinfeld map and rotation-admissible cyclotomic modular data, under the precise definitions in each manuscript.

The 98-page file is the revised complete edition shared on **20 November 2026**. The 24- and 37-page files are the extracts shared on **19 November 2026**. Their manuscript bytes and companion code are retained as supplied.

## Improvement over the earlier work

| Hypothesis | Earlier result | New result |
| --- | --- | --- |
| General | **32**, in the [17-page paper](../Rank32.pdf) and the [252-page project manuscript](../Manuscript.pdf) | **45**, in all three new papers |
| Odd dimension | **57**, in the 252-page project manuscript | **217**, in the 37-page and complete papers |
| Odd nonsquare dimension | **75**, in the 252-page project manuscript | **219**, in the complete paper |

For rank bounds, as for dimension bounds, the investigation deliberately tried to push the AI model to its limits by repeatedly seeking stronger theoretical restrictions. The current work still leaves **rank 45** and **odd rank 217** undecided. Its surviving arithmetic configurations are not constructed categories.

The original [17-page manuscript](../Rank32.pdf), [source](../Rank32.tex), and [verification package](../rank-bound-supplement/README.md) remain a separate reading route.

## The additional classification

The complete paper classifies the **nontrivial perfect weakly integral modular categories through rank 44**, up to braided equivalence. It gives 108 classes in ranks 22, 32, 37 and 44. At rank 44, the 18 deleted-module $A_5$-gaugings are weakly group-theoretical and not group-theoretical. Thus the paper identifies 44 as the least rank for a perfect weakly integral modular category that is not group-theoretical.

This classification concerns perfect categories; the lower bound 45 concerns non-pointed **simple integral** modular categories. See Section 10 of the [complete manuscript](full/Rank-modular-data.pdf) for the precise statements and proofs.

## Verification

| Edition | Its own companion |
| --- | --- |
| 24 pages | [Instructions](general/README.md#verification) · [Verifier](general/supplement/verify.py) · [Reference census](general/supplement/census/reference-census.zip) |
| 37 pages | [Instructions](general-and-odd/README.md#verification) · [Verifier](general-and-odd/supplement/verify.py) · [Reference census](general-and-odd/supplement/census/reference-census.zip) |
| Revised 98 pages | [Supplement guide](full/rank-bound-supplement/README.md) · [Proof map](full/rank-bound-supplement/proof-map.md) · [Verifier](full/rank-bound-supplement/verify.py) |

The same filenames can contain different implementations in different editions. Run the companion supplied with the selected paper. The complete paper's PDF links to its adjacent `rank-bound-supplement/`; the extracts retain their respective `supplement/` directories.

The complete manuscript distinguishes a bounded audit, checks of stored execution records, full odd replays, and fresh general-census generation. These are separate operations, described in its [guide](full/rank-bound-supplement/README.md). Lean formalization of the 24-page paper was still in progress when these editions were shared.

## References

All three papers cite **[Hexagon 2611.00233v3](https://hexagonmath.org/2611.00233v3)** for specified earlier results. That identifier refers to the 252-page project manuscript, not to these additions. See [Versions and citations](../VERSIONS.md) for the preserved paths and fixed computational snapshots.
