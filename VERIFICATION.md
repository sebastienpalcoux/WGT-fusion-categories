# Verification guide

[Overview](README.md) · [All manuscripts](papers/README.md) · [Dimension bounds](papers/dimension-bounds.md) · [Rank bounds](rank/README.md) · [Versions and citations](papers/VERSIONS.md)

Choose the companion belonging to the **paper and edition you are reading**. The repository contains several generations of rank calculations and two different dimension investigations. Their theorem numbers, inputs and commands have distinct scopes.

## Manuscript-to-companion map

| Manuscript | Companion | What it supplies |
| --- | --- | --- |
| 29-page dimension theorem | [Manuscript](papers/Below_5400_squared.pdf) | A computer-free proof. No computational certificate is needed for this paper. |
| 47-page three-prime theorem | [Bounded arithmetic checks](papers/three-prime-dimensions/verification/README.md) | Supporting exact checks outside the computer-free proof. |
| 24-page general rank theorem | [Edition guide](papers/rank-bounds/general/README.md) · [Verifier](papers/rank-bounds/general/supplement/verify.py) | The finite arithmetic and archived general census accompanying the rank-45 extract. |
| 37-page general and odd rank theorems | [Edition guide](papers/rank-bounds/general-and-odd/README.md) · [Verifier](papers/rank-bounds/general-and-odd/supplement/verify.py) | The extract's general census and odd-dimensional calculations for the bounds 45 and 217. |
| Revised 98-page complete rank manuscript | [Supplement guide](papers/rank-bounds/full/rank-bound-supplement/README.md) · [Proof map](papers/rank-bounds/full/rank-bound-supplement/proof-map.md) | The general covering corpus, terminal exclusions, odd calculations, independent checks, and finite arithmetic used with the perfect-category classification. |
| Original 17-page rank-32 paper | [Original rank companion](papers/rank-bound-supplement/README.md) | The original independent-mode enumeration and certificate through rank 31. |
| 252-page project manuscript | [Release guide](certificates/README.md) · [Manifest](certificates/manifest.json) · [Statement index](certificates/statement-index.tsv) | Original dimension, fixed-dimension, rank and support calculations. The general rank calculation also enables the manuscript's dimension-family filters. |
| 423-page three-prime draft | [Compact companion](papers/research-drafts/p6q6r4/verification/README.md) | The specified bounded checks and preserved source records. It does not establish stronger numerical conclusions that require omitted exhaustive inventories. |

## Running checks

Follow the commands in the selected edition's guide. The original rank-32 companion, the two new extracts, and the revised complete rank manuscript have different drivers. In particular, `--full` does not have a universal meaning across these packages.

The revised complete rank manuscript documents both a bounded audit and full production replays. A terminal replay checks the supplied lists; regenerating an exhaustive census is a different operation. The guides state which operation a command performs and what assumptions its conclusion uses.

The [original arithmetic programs](code/README.md) and [rank-15 input](data/README.md) retain their historical scope. Earlier exploratory programs under `rank/continuation/` are described in their own [guide](rank/continuation/README.md); they are not the companion for the new rank-45 theorem.

## Files and proof status

Checksums identify the exact files in an edition. The mathematical soundness of its reductions, the completeness of enumeration, and the validity of terminal exclusions are separate parts of the argument. The supplied documentation distinguishes these parts.

The three new manuscripts and their companions are supplied as separate editions. Lean formalization of the 24-page paper was still in progress when these manuscripts were shared; this collection does not assert a completed Lean certificate.

For stable code citations and the numbered references into Hexagon, see [Versions and citations](papers/VERSIONS.md).
