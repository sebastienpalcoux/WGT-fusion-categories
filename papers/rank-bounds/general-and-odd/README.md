# Rank bounds for simple integral modular data

**Sebastian Burciu, Jingcheng Dong and Sébastien Palcoux · 37 pages**

[PDF](Rank_bounds.pdf) · [TeX source](Rank_bounds.tex) · [Compare the three editions](../README.md) · [All manuscripts](../../README.md) · [Versions and citations](../../VERSIONS.md)

This is the focused proof of the general and odd-dimensional statements of Theorem A, shared on 19 November 2026. It proves rank at least **45**, and at least **217 in odd dimension**, under the theorem's algebraic hypotheses. These bounds apply to non-pointed simple integral modular fusion categories over ℂ.

The general result improves the [17-page rank-32 theorem](../../Rank32.pdf); the odd result improves the rank-57 theorem in the [252-page project manuscript](../../Manuscript.pdf). The [24-page extract](../general/README.md) treats only the general bound. The [revised 98-page complete paper](../full/README.md) includes the wider theory and perfect-category classification.

## Compilation

With a standard LaTeX installation, run from this directory:

```sh
latexmk -pdf Rank_bounds.tex
```

## Verification

The exact companion for this edition is in [supplement/](supplement/). Use Python 3.10 or newer. From this directory:

```sh
python3 supplement/verify.py --jobs 3
```

The default command replays the terminal tests and **both complete odd censuses**, along with structural arithmetic. It does not regenerate the general denominator census. The [driver](supplement/verify.py) documents `--partitions` for native task-partition regeneration and `--full --jobs N` for fresh general-census generation. These additional native options require Linux and a C++17 compiler.

The [reference archive](supplement/census/reference-census.zip) contains the producer sources and completed general records. For a bounded audit of the revised complete edition instead, use that edition's own [instructions](../full/rank-bound-supplement/README.md#commands-and-scope).

The supplied [companion checksums](supplement/SHA256SUMS) are checked by running `sha256sum -c SHA256SUMS` **inside `supplement/`**. The [edition manifest](MANIFEST.sha256) checks the manuscript, guide and companion from this directory.

## References

The bibliography cites [Hexagon 2611.00233v3](https://hexagonmath.org/2611.00233v3) for specified earlier results. This paper retains its own original manuscript files and matching companion; the complete paper's similarly named programs are a different edition.
