# Manuscript and compact exact checks

The complete theoretical discussion is in `p6q6r4_manuscript.tex` and
`p6q6r4_manuscript.pdf`. There is no theoretical appendix or supplement.
The [revision response](revision_letter.pdf) explains the referee corrections
and the precise status of the computational conclusions.

## Build

Run `bash build.sh` with a working pdfLaTeX installation and the standard
packages named in the source. Three passes resolve cross-references.

## Verification scope

See [verification/README.md](verification/README.md) for the compact exact
checks and their inputs, outputs and limitations. The files reproduce
bounded arithmetic and specific obstructions; they do not certify the
large field-and-row inventories. Statements depending on those inventories
are explicitly conditional in the manuscript. Retained source code is not,
by itself, evidence of a completed exhaustive computation.

The family theorem and the general theoretical arguments are distinguished
from those conditional numerical conclusions. The dimension $5400^2$
remains unresolved. Necessary field-and-duality states are not constructed
fusion rings or categories.

## Integrity and size

`SHA256SUMS.txt` authenticates the delivered package files except itself.
Run `sha256sum -c SHA256SUMS.txt` in this directory. Hashes check file integrity,
not mathematical correctness.

Apart from the two manuscripts' TeX/PDF files, the combined size of all
accompanying material is below 10,000,000 bytes. No large certificate archive or split archive is required by the
compact checks. Historical large replay commands are not supported by this
release; use only the commands in the current verification guide.
