# Versions and citations

[Overview](../README.md) · [All manuscripts](README.md) · [Dimension bounds](dimension-bounds.md) · [Rank bounds](../rank/README.md) · [Verification](../VERIFICATION.md)

Use the **title and edition of the document actually consulted**. The repository contains related manuscripts with different scopes, theorem numbers and verification packages. Its catalogue groups them for reading; it does not identify them as interchangeable versions of one paper.

## Manuscript identities

| Document | Identity and relationship |
| --- | --- |
| [252-page project manuscript](Manuscript.pdf) | *Prime factors of fusion-category dimensions: weak group-theoreticity, rank and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. This work is recorded as **Hexagon 2611.00233**. [Record](https://hexagonmath.org/2611.00233) · [Version 3](https://hexagonmath.org/2611.00233v3). |
| [17-page rank-32 paper](Rank32.pdf) | *A rank bound for simple integral modular categories*. A separate focused proof of the general lower bound 32, with its original [verification companion](rank-bound-supplement/README.md). |
| [29-page dimension paper](Below_5400_squared.pdf) | *Weakly integral braided fusion categories below dimension 5400²*. A separate computer-free proof of the dimension theorem. |
| [47-page family paper](three-prime-dimensions/three_prime_dimensions.pdf) | *Braided fusion categories with dimension dividing p⁶q⁶r⁴*. A self-contained extraction attributing its results and methods to the joint project. |
| [423-page family draft](research-drafts/p6q6r4/p6q6r4_manuscript.pdf) | *Arithmetic restrictions on braided fusion categories of dimension p⁶q⁶r⁴*. A separate detailed investigation, including references to Hexagon version 3. |
| [24-page general rank paper](rank-bounds/general/Rank_bound.pdf) | *A rank bound for simple integral modular data*. The general rank-45 extract shared on 19 November 2026. |
| [37-page general and odd rank paper](rank-bounds/general-and-odd/Rank_bounds.pdf) | *Rank bounds for simple integral modular data*. The rank-45/217 extract shared on 19 November 2026. |
| [98-page complete rank paper](rank-bounds/full/Rank-modular-data.pdf) | *The rank of simple integral modular data*. The revised complete version shared on 20 November 2026. Only this revised complete edition is supplied in the new rank collection. |

All three new rank manuscripts cite **[Hexagon 2611.00233v3](https://hexagonmath.org/2611.00233v3)** under the bibliography key `BDPcore`, including numbered results from that edition. That reference belongs to the 252-page project manuscript. Adding the new papers here does not create a new version of the Hexagon record or change those references.

The original 17-page paper proves the **general bound 32**. The earlier **odd-dimensional bound 57** belongs to the 252-page project manuscript. The new general and odd bounds are **45 and 217**; the complete rank paper also gives **219 in odd nonsquare dimension**.

## Hexagon history

The [Hexagon version history](https://hexagonmath.org/2611.00233v3) distinguishes a **229-page version 1**, a **252-page version 2**, and a **252-page version 3**. Version 3 records the improved earlier general rank bound 32 and its verification through rank 31. A page or theorem reference to a previous edition should continue to use that edition.

## Fixed code and manuscript snapshots

| Reference | Purpose |
| --- | --- |
| [433d1c45795202991f1693e00d9770e4cccf05f6](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/433d1c45795202991f1693e00d9770e4cccf05f6) | The exact code commit cited throughout the 252-page manuscript and by the 17-page rank-32 paper. |
| [f2ef6eedb58ade6adcdb5ed9caeba1836b4ba8ce](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/f2ef6eedb58ade6adcdb5ed9caeba1836b4ba8ce) | The historical 229-page coauthor edition. Existing data/code guides use its specific page and statement numbers. |
| [independent-referee-2026-11-10](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/independent-referee-2026-11-10) | The retained publication-snapshot branch. Its recorded immutable commit is [8b44b13483d82d204194549b638cbfe2693437af](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/8b44b13483d82d204194549b638cbfe2693437af). |

Existing manuscript paths and these fixed references are preserved. For a reproducible reference to one of the new papers or its companion, use GitHub's permanent link to the relevant file at the commit consulted. Keep the paper and its supplement at the same edition.

## Matching a paper to its computations

The [verification guide](../VERIFICATION.md) is the entry point. In particular, the original `papers/rank-bound-supplement/` supplies the finite verification for the rank-32 result; the new editions have separate companions under `papers/rank-bounds/`. Their same-named programs can differ, so their files and commands should not be mixed.

The complete rank PDF uses local links to its adjacent `rank-bound-supplement/` folder. Keep that folder beside the PDF when working from a checkout or downloaded repository. Each extract similarly retains its own `supplement/` directory. The manuscript PDFs and TeX sources are published as supplied.
