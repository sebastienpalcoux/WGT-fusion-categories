# A rank bound for simple integral modular data

**Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux · 24 pages**

[PDF](Rank_bound.pdf) · [TeX source](Rank_bound.tex) · [Compare the three editions](../README.md) · [All manuscripts](../../README.md) · [Versions and citations](../../VERSIONS.md)

This is the focused general-case proof shared on 19 November 2026. It proves rank at least **45** under the algebraic hypotheses of Theorem A, and hence for every non-pointed simple integral modular fusion category over ℂ. It improves the original [17-page rank-32 result](../../Rank32.pdf).

For the odd-dimensional bound 217, use the [37-page paper](../general-and-odd/README.md). For the wider theory and perfect-category classification, use the [revised 98-page complete manuscript](../full/README.md).

## Compilation

With a standard LaTeX installation, run from this directory:

```sh
latexmk -pdf Rank_bound.tex
```

## Verification

The complete supplied companion is in [supplement/](supplement/). Use Python 3.10 or newer. From this directory:

```sh
python3 supplement/verify.py
```

The default command checks the archived degree lists and structural arithmetic. It does not regenerate the complete general census. Its [reference archive](supplement/census/reference-census.zip) contains the original producer sources and completed records.

The [driver](supplement/verify.py) also documents `--partitions` for native task-partition regeneration and `--full --jobs N` for fresh census generation through rank 44. Those options additionally require Linux and a C++17 compiler; a full enumeration is a separate, larger computation.

The supplied [companion checksums](supplement/SHA256SUMS) are checked by running `sha256sum -c SHA256SUMS` **inside `supplement/`**. The [edition manifest](MANIFEST.sha256) checks the manuscript, guide and companion from this directory.

The companion is preserved for this extract. Use its own files rather than substituting the similarly named programs from another edition. Lean formalization was still in progress when this manuscript was shared.

## References

The bibliography cites [Hexagon 2611.00233v3](https://hexagonmath.org/2611.00233v3) for specified earlier results. The new paper is a separate document, and its original PDF and TeX are supplied without manuscript edits.
