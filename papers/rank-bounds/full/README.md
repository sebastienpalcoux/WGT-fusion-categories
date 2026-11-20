# The rank of simple integral modular data

**Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux · revised complete version · 98 pages**

[PDF](Rank-modular-data.pdf) · [TeX source](Rank-modular-data.tex) · [Computational supplement](rank-bound-supplement/README.md) · [Compare the three editions](../README.md) · [All manuscripts](../../README.md)

This is the revised complete manuscript shared on **20 November 2026**. Theorem A and its category corollary give rank at least **45** in general, **217** in odd dimension, and **219** in odd nonsquare dimension. The category statement concerns non-pointed simple integral modular fusion categories over ℂ; the algebraic theorem uses the precise Drinfeld-duality, modular-data and rotation hypotheses stated in the paper.

The general bound substantially strengthens the original [17-page rank-32 paper](../../Rank32.pdf). The [24-page extract](../general/README.md) isolates the general proof, and the [37-page extract](../general-and-odd/README.md) treats the general and odd cases.

## Wider results and remaining ranks

Section 10 classifies the nontrivial perfect weakly integral modular categories through rank **44**, up to braided equivalence, giving **108 classes**. At rank 44 it includes 24 twisted doubles of $A_6$ and 18 $A_5$-gaugings of the pointed category on the deleted binary permutation module. The latter are weakly group-theoretical and not group-theoretical.

The first open simple-category ranks remain **45** in general and **217** in odd dimension. The manuscript records 204 general starting tuples with 462 refinements at rank 45, and 1,628 odd graph reservations at rank 217. These are configurations requiring further restrictions, not examples of categories.

## Compilation and file layout

With a standard LaTeX installation, run from this directory:

```sh
latexmk -pdf Rank-modular-data.tex
```

Keep `rank-bound-supplement/` beside the PDF and TeX. The PDF's local links use that relative directory. The manuscript and all companion programs and data are supplied as in the revised package.

## Verification

The [supplement guide](rank-bound-supplement/README.md) distinguishes the fixed covering corpus, terminal tests, completed odd-run records, fresh enumeration, and independent checks. The [proof map](rank-bound-supplement/proof-map.md) links the mathematical implications to their implementations.

A bounded audit, including independent finite checks and native-prefix comparisons, can be run from the supplement directory:

```sh
cd rank-bound-supplement
python3 verify.py --quick --independent --native-prefix
```

This requires Python 3.10 or newer, SymPy, and a C++17 compiler with undefined-behavior sanitizer support. It does not regenerate the general denominator census or replay the two full odd continuation recursions. The guide gives the commands for those separate operations.

The [edition manifest](MANIFEST.sha256) records the published manuscript, guide, and complete companion. Run `sha256sum -c MANIFEST.sha256` from this directory to check file integrity.

## References and edition

The bibliography retains its numbered references to **[Hexagon 2611.00233v3](https://hexagonmath.org/2611.00233v3)**, the separate 252-page project manuscript. This publication supplies the revised complete rank manuscript without changing the earlier project manuscript, its rank-32 extract, or the Hexagon record. See [Versions and citations](../../VERSIONS.md).
