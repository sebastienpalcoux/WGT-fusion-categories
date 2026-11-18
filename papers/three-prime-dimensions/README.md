# Braided fusion categories with three-prime dimensions

[Manuscript (PDF, 47 pages)](three_prime_dimensions.pdf) ·
[LaTeX source](three_prime_dimensions.tex) ·
[Revision response (PDF)](revision_letter.pdf) ·
[Response source](revision_letter.tex)

This self-contained, computer-free paper treats braided fusion categories with
dimension dividing $p^6q^6r^4$. The main theorem is unchanged by the revision.
The endpoint dimensions remain unresolved; necessary arithmetic data do not
construct a fusion ring or a category.

The [bounded arithmetic checks](verification/README.md) are supporting checks
outside the proof. They do not certify categorical arguments or arbitrary-conductor
statements. The [423-page research draft](../research-drafts/p6q6r4/p6q6r4_manuscript.pdf)
records the broader investigation and distinguishes theoretical results from
numerical conclusions conditional on complete computational coverage.

## Build and integrity

Compile `three_prime_dimensions.tex` with pdfLaTeX, repeating until references
are resolved. Compile `revision_letter.tex` separately for the response.
Run `sha256sum -c SHA256SUMS.txt` in this directory to check file integrity.
Hashes are not a certificate of mathematical correctness.

[All manuscripts](../README.md) · [Revision notes](../../REVISION_NOTES.md) ·
[Accompanying-material size report](../../SIZE_REPORT.json)
