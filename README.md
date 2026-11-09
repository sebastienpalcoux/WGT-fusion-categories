# WGT for fusion categories

*Prime factors of fusion-category dimensions: weak group-theoreticity and decomposition*, by Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux. [Full article (PDF)](papers/Manuscript.pdf) · [TeX source](papers/Manuscript.tex)

A shorter paper gives a standalone proof that every weakly integral braided fusion category of dimension below 5400² is weakly group-theoretical. Its finite arithmetic can be checked by hand. [Short paper (PDF, 29 pages)](papers/Below_5400_squared.pdf) · [TeX source](papers/Below_5400_squared.tex)

| Material | Contents |
| --- | --- |
| [Verification code](code/README.md) | Commands for the finite arithmetic used in the paper. |
| [Input data](data/README.md) | The pinned list of 2341 perfect rank-15 dimension types, with provenance and checksum. |
| [Verification manifest](certificates/manifest.json) | Hypotheses, commands, expected outputs, and locations in the full article for 24 calculations and checks. |
| [Explicit certificates](certificates/certificates.json) | The terminal trace residues for seven dyadic cases and the two disjoint fibers for each of five rank-100 profiles. |

A short check of the explicit certificates is available from the repository root.

```sh
python3 code/verify_witnesses.py
```

It tests the recorded residue and norm conditions. The completeness of both lists is checked separately by the enumeration programs in the manifest. The final rank-100 obstruction uses the mathematical argument on tame twists.
