# Rank bounds and proof calculations

[Overview](../README.md) · [All manuscripts](../papers/README.md) · [Dimension bounds](../papers/dimension-bounds.md) · [Verification](../VERIFICATION.md) · [Versions and citations](../papers/VERSIONS.md)

## The general rank bound

The new manuscripts strengthen the general bound for **non-pointed simple integral modular fusion categories over ℂ** from **32 to 45**. The versions treating odd dimension give **217**; the revised complete paper also gives **219 in odd nonsquare dimension**.

| Reading route | Scope | Paper and matching companion |
| --- | --- | --- |
| Focused general proof · 24 pages | Rank at least 45 | [PDF](../papers/rank-bounds/general/Rank_bound.pdf) · [Guide](../papers/rank-bounds/general/README.md) |
| General and odd proof · 37 pages | Rank at least 45, and 217 in odd dimension | [PDF](../papers/rank-bounds/general-and-odd/Rank_bounds.pdf) · [Guide](../papers/rank-bounds/general-and-odd/README.md) |
| Revised complete rank manuscript · 98 pages | Bounds 45 / 217 / 219, wider theory, perfect-category classification through rank 44 | [PDF](../papers/rank-bounds/full/Rank-modular-data.pdf) · [Guide](../papers/rank-bounds/full/README.md) |
| Original general proof · 17 pages | Rank at least 32 | [PDF](../papers/Rank32.pdf) · [Original companion](../papers/rank-bound-supplement/README.md) |

The [new rank collection](../papers/rank-bounds/README.md) compares the precise hypotheses and additional results. The complete paper's classification concerns **perfect weakly integral** modular categories; its simple-category rank theorem concerns **non-pointed simple integral** modular categories.

For this investigation, as for the [dimension investigation](../papers/dimension-bounds.md), we tried to push the AI model to its limits through repeated theoretical improvements. Rank 45 and odd rank 217 remain undecided.

## Verification for the new manuscripts

Follow the instructions for the selected edition. The [complete rank supplement](../papers/rank-bounds/full/rank-bound-supplement/README.md) distinguishes bounded audits, full odd replays, and fresh general enumeration. The [24-page](../papers/rank-bounds/general/README.md#verification) and [37-page](../papers/rank-bounds/general-and-odd/README.md#verification) extracts have their own drivers and source records.

Their companions are separate from the original rank-32 programs described below. The [verification map](../VERIFICATION.md) and [version guide](../papers/VERSIONS.md) identify the corresponding manuscripts and fixed references.

## Original rank-32 calculation

The original general theorem gives rank at least 32. The [252-page project manuscript](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/main/papers/Manuscript.pdf) and the [17-page self-contained paper](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/main/papers/Rank32.pdf) use the same degree equation, local character restrictions and finite verification. The 252-page manuscript also enables its proved dimension-family exclusions within that search. The self-contained paper supplies a complete proof using the same program with those additional filters disabled.

There is one [optimized verifier](../papers/rank-bound-supplement/verify.py), one pair of rank-parameterized generators, and one [certificate](../papers/rank-bound-supplement/certificate-families.json.gz) for the entire range proved in the full article. Ranks 1–9 are excluded analytically. The default computation also checks rank 9 and exhausts every rank through 31. From the repository root, check the published record with

```sh
python3 papers/rank-bound-supplement/verify.py --families --check papers/rank-bound-supplement/certificate-families.json.gz
```

To regenerate every list with four workers and a fresh checkpoint directory, run

```sh
python3 papers/rank-bound-supplement/verify.py --families --jobs 4 --work rank-results
```

The full-article run has 3,012 tasks, 411,325 square-dimension profiles and 122 nonsquare-dimension profiles, with no survivors. The [supplement guide](../papers/rank-bound-supplement/README.md) explains the records, resuming a run, and choosing another rank. For example,

```sh
python3 papers/rank-bound-supplement/verify.py --families --rank 32 --jobs 4 --work rank-32
```

The same recursions accept any positive rank. Bounds whose proof has an upper-rank or support hypothesis are applied only under that hypothesis. Higher-rank computations are exploratory, and the published certificate covers the theorem's range.

## Further rank bounds from the full article

This section records the additional bounds in the **252-page project manuscript**, with its original numbering and computations. Its general and odd bounds precede the improvements above. The following computations have additional hypotheses. Their [manifest entries](../certificates/manifest.json) specify the input files, commands and outcomes. The [statement index](../certificates/statement-index.tsv) identifies their locations in the article.

| Category | Lower bound | Calculation |
|---|---:|---|
| Non-pointed simple integral modular, nonsquare dimension | 33 | Unified calculation R01 through 31; [nonsquare rank-32 check](verify_nonsquare32.py), R03 |
| Non-pointed simple integral modular, odd dimension / odd nonsquare dimension | 57 / 75 | [Paired recursion](code/odd_frontier.py), R06 |
| Square dimension, largest degree smaller than √D/4 / square root not divisible by four / largest degree smaller than √D/5 | 36 / 39 / 46 | [Conditional square searches](square_slices.py), R04 |
| Divisor-minimal slightly degenerate obstruction | 34 (full rank) | [Reduced-table check](verify_scoped.py), R05 |
| Weakly integral braided non-WGT, any / odd dimension | 9 / 41 | [Formal row test](code/formal_check.py), R08, and the odd twist-counting argument |
| Weakly integral non-WGT, any / odd dimension | 7 / 13 | [Small ranks](code/unbraided_small.py), [odd rank eleven](code/odd_unbraided.py), R09–R10 |

Strict weak integrality gives bounds 8 without braiding and 11 with braiding. The fixed-dimension conclusions include rank at least 101 at the modular boundary 29,160,000 and at least 553 at the odd braided boundary 5⁶·7⁶·23⁴.

The shared code under `current/` and `code/` supplies the conditional square, fermion-quotient and odd-dimension calculations. Earlier general searches and the former continuation records are available in the [publication snapshot](https://github.com/sebastienpalcoux/WGT-fusion-categories/tree/independent-referee-2026-11-10). The [release guide](../certificates/README.md) describes the current inputs and checks.
