# WGT for fusion categories

*Prime factors of fusion-category dimensions: weak group-theoreticity, rank and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. [Full article (PDF, 252 pages)](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.pdf) · [TeX source](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.tex)

The [rank bounds](rank/README.md) give rank at least 31 for a non-pointed simple integral modular category, 33 in nonsquare dimension, 57 in odd dimension and 75 in odd nonsquare dimension. The revised paper distinguishes the ordinary, reduced and arbitrary-center cases.

The introduction gives the rank bounds by hypothesis and a catalogue of the prime families. The [release guide](certificates/README.md) links each proof computation to its exact hypotheses and the [statement index](certificates/statement-index.tsv).

A shorter paper proves that every weakly integral braided fusion category of dimension below 5400² is weakly group-theoretical. Applying this result to the Drinfeld center gives the bound 5400 without a braiding assumption. The proof develops the shared principal-block and local Weyl arguments and gives all finite arithmetic in a form that can be checked by hand. [Short paper (PDF, 29 pages)](papers/Below_5400_squared.pdf) · [TeX source](papers/Below_5400_squared.tex)

| Material | Contents |
| --- | --- |
| [Rank bounds](rank/README.md) | Degree enumeration, mathematical exclusion tests, certificates and a resumable rank-31 search. |
| [Verification code](code/README.md) | Commands for the finite arithmetic used in the paper. |
| [Input data](data/README.md) | The pinned list of 2341 perfect rank-15 dimension types, with provenance and checksum. |
| [Verification manifest](certificates/manifest.json) | Hypotheses, commands, expected outputs, and locations in the full article for the dimension, character-table and rank calculations. |
| [Explicit certificates](certificates/certificates.json) | The terminal trace residues for seven dyadic cases and the two disjoint fibers for each of five rank-100 profiles. |

The [dimension checker](code/README.md#dimension-comparison) reproduces all 16 entries in the dimension theorem and prints the four dimension bounds. Run it from the repository root with

```sh
python3 code/check_dimension_bounds.py
```

A short check of the explicit certificates is available from the repository root.

```sh
python3 code/verify_witnesses.py
```

It tests the recorded residue and norm conditions. The enumeration programs in the manifest recover the complete lists under the hypotheses of the corresponding proofs. The boundary-rank proof excludes the five remaining rank-100 profiles by showing that two quartets of different degrees would have to have equal degree.
