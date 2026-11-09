# WGT for fusion categories

*Prime factors of fusion-category dimensions: weak group-theoreticity and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. [Full article (PDF, 229 pages)](papers/Manuscript.pdf) · [TeX source](papers/Manuscript.tex)

A shorter paper proves that every weakly integral braided fusion category of dimension below 5400² is weakly group-theoretical. Applying this result to the Drinfeld center gives the bound 5400 without a braiding assumption. The proof develops the shared principal-block and local Weyl arguments and gives all finite arithmetic in a form that can be checked by hand. [Short paper (PDF, 29 pages)](papers/Below_5400_squared.pdf) · [TeX source](papers/Below_5400_squared.tex)

| Material | Contents |
| --- | --- |
| [Verification code](code/README.md) | Commands for the finite arithmetic used in the paper. |
| [Input data](data/README.md) | The pinned list of 2341 perfect rank-15 dimension types, with provenance and checksum. |
| [Verification manifest](certificates/manifest.json) | Hypotheses, commands, expected outputs, and locations in the full article for 27 calculations and checks. |
| [Explicit certificates](certificates/certificates.json) | The terminal trace residues for seven dyadic cases and the two disjoint fibers for each of five rank-100 profiles. |

The [dimension checker](code/README.md#dimension-comparison) reproduces all 16 entries in the proof of Theorem 17.1 and prints the four dimension bounds. Run it from the repository root with

```sh
python3 code/check_dimension_bounds.py
```

A short check of the explicit certificates is available from the repository root.

```sh
python3 code/verify_witnesses.py
```

It tests the recorded residue and norm conditions. The enumeration programs in the manifest recover the complete lists under the hypotheses of the corresponding proofs. Section 24.4 of the full article excludes the five remaining rank-100 profiles by showing that two quartets of different degrees would have to have equal degree.
