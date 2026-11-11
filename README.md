# WGT for fusion categories

*Prime factors of fusion-category dimensions: weak group-theoreticity, rank and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. [Full article (PDF, 252 pages)](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.pdf) · [TeX source](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.tex)

The full article gives [rank bounds](rank/README.md) under several hypotheses and a catalogue of prime decompositions. The [release guide](certificates/README.md) links its computations to their hypotheses and [statements](certificates/statement-index.tsv).

## Short, self-contained proofs

Each paper extracts one theorem from the larger project and includes the arguments needed to read and check it independently.

- **Rank at least 32** for every non-pointed simple integral modular category (17 pages). A single program verifies the finite arithmetic in the last step. [PDF](papers/Rank32.pdf) · [TeX](papers/Rank32.tex) · [Verification supplement](papers/rank-bound-supplement/README.md)
- **Weak group-theoreticity below dimension $5400^2$** for weakly integral braided fusion categories (29 pages). The proof is computer-free. Applying it to the Drinfeld center gives the bound 5400 without a braiding assumption. [PDF](papers/Below_5400_squared.pdf) · [TeX](papers/Below_5400_squared.tex)

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
