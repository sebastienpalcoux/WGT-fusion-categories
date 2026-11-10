# Rank-15 dimension types

[published-rank15-perfect.txt](published-rank15-perfect.txt) contains one sorted dimension type per line. The exhaustive list of 2341 remaining perfect rank-15 types is an input from Section 4.3 of Dong and Palcoux, *A new criterion for integral modular categorification*. [arXiv:2510.10311v1](https://arxiv.org/abs/2510.10311v1) · [Pinned upstream file](https://github.com/sebastienpalcoux/Fusion-Categories/blob/06f3a37b086c0591e086a3aed9a9458f4270165f/NewCritPaper/Rank15/RemainingPerfect.txt)

The [rank-15 check](../code/verify_rank15_published.py) verifies the checksum and the number of distinct types, then computes the smallest formal codegree of each type. In [Corollary 18.2 of the 229-page coauthor edition](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/f2ef6eedb58ade6adcdb5ed9caeba1836b4ba8ce/papers/Manuscript.pdf#page=112) the respective counts for codegrees 2, 3, 4 and 5 are 459, 24, 1857 and 1. The sole type with smallest codegree 5 has dimension 720. The merged manuscript generates its rank profiles independently of this classification. [Current rank certificates](../rank/README.md)

SHA-256: `04a25855339d84fed23fb0e67cdef00b7783e2f1816d83fe679e3e143ac75dbc`.
