# WGT for fusion categories

*Prime factors of fusion-category dimensions: weak group-theoreticity, rank and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. [Full article (PDF, 252 pages)](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.pdf) · [TeX source](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.tex)

The [rank bounds](rank/README.md) in the full article give rank at least 31 for a non-pointed simple integral modular category, 33 in nonsquare dimension, 57 in odd dimension and 75 in odd nonsquare dimension. The short rank paper below extends the general bound to 32.

The introduction gives the rank bounds by hypothesis and a catalogue of the prime families. The [release guide](certificates/README.md) links each proof computation to its exact hypotheses and the [statement index](certificates/statement-index.tsv).

## Short, self-contained proofs

These two manuscripts extract individual results from the larger project. They give readers a shorter route to the statements and their proofs, with the required arithmetic and verification material available separately.

- **Rank at least 32** for a non-pointed simple integral modular category (21 pages). [PDF](papers/Rank32.pdf) · [TeX](papers/Rank32.tex) · [Unified verification supplement](papers/rank-bound-supplement/README.md)
- **Weak group-theoreticity below dimension $5400^2$** for weakly integral braided fusion categories (29 pages). The Drinfeld-center argument also gives the bound 5400 without a braiding assumption. [PDF](papers/Below_5400_squared.pdf) · [TeX](papers/Below_5400_squared.tex)

| Material | Contents |
| --- | --- |
| [Short manuscripts](papers/README.md) | Both focused proofs and their source files. |
| [Rank bounds](rank/README.md) | Degree enumeration, mathematical exclusion tests and the unified verification of ranks 2–31. |
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
