# Reproducing the arithmetic

Use the commit named in the article’s **Code, data and certificates** section. GitHub’s **Download ZIP** at that commit contains all sources and inputs. Commands below run from the repository root unless indicated. Use Python 3.10 or later with assertions enabled. The compiled searches require C++17 with 128-bit integers and, for the full ordinary searches, GNU MP and its C++ bindings. The seven-position search requires C++20.

The [corresponding manuscript](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.pdf) and [TeX source](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/papers/Manuscript.tex) are in the named publication snapshot. Its paper cites this computation commit. The code commit was created first so its hash could be printed in the paper.

The [manifest](manifest.json) identifies each calculation by its hypotheses and permanent TeX labels. The [statement index](statement-index.tsv) maps those labels to the printed numbers and pages in the revised article. The index also gives each hyperref destination and the title where available. Two historical comparisons in the manifest give their original edition and pages explicitly. For R01, the manifest maps every recorded terminal reason to its proof under `ordinary_recorded_terminal_reasons`, describes each pruning test under `ordinary_filter_map`, and gives the ranks 28–30 leaf partitions under `rank_leaf_partitions`.

The input and program digests cover the code, data and receipts used at this commit. A digest identifies bytes. Partition coverage and the mathematical validity of each rejection require separate checks.

```sh
python3 certificates/check_release.py
python3 certificates/check_release.py --cases D00 D01 D26
python3 certificates/check_release.py --cases all
```

The first command checks the file inventory. The second regenerates the dimension minima and the boundary profiles, then checks the seven dyadic trace witnesses and the two disjoint fibers for each of the five rank-100 profiles. The third executes all 27 dimension, fixed-dimension, support-window and comparison commands, including the largest C++ support search. Only an elapsed-time suffix is ignored in the seven odd-window outputs.

## Which object is checked?

| Domain | Candidate generation | Verification and conclusion |
|---|---|---|
| Ordinary integral, non-pointed, fusion-simple modular data | Start with the finite pairs `(q,B)` with `q` square-free and `q B² < r`. The prefix lcm, residual reciprocal-square mass, number of unused labels and denominator order determine a finite next interval. At two remaining labels solve the integral square equation. `(q,B)=(1,3)` is excluded by the ordinary no-codegree-nine theorem. | [Ordinary release](../rank/ordinary_release.py) checks the unique unit, degree-square divisibility, absence of codegree nine, uniqueness and digest of each list, and recomputes the terminal exclusions. At ranks 28–30 it checks the complete manifest against the leaf receipts, including theoretical shortcuts and splits. `--replay` also regenerates every list. Through rank 27 a whole GMP traversal supplies a separate comparison. R01 lists the exact counts. |
| Nonsquare/paired ordinary data | Use the stated square-free kernels, and count dual pairs only in odd dimension. Count and order comparisons use integers; the odd mass upper bound rounds each term upwards. | R02–R04 and R06 compare the generated partitions or the per-root lists and the terminal screen under their separate hypotheses. |
| Positive quotient of a divisor-minimal fermionic obstruction | Use the fermion orbit basis, half-rank and reduced dimension. Apply dimension families to **twice** the table dimension. | R05 uses only proved reduced-table conditions. Its terminal list is independent of the ordinary codegree-nine test. |
| General weakly integral categories | Generate integer formal-codegree tuples from the reciprocal identity under the relevant commutativity and grading reductions. | R08–R10 give the rank-eight formal and arbitrary-center checks. Some terminal tuples require a categorical descent or a character argument in the paper. |
| Boundary dimensions | In the rank-100 calculation fix a twist configuration, enumerate possible degrees and multiplicities, then solve the two disjoint-fiber budgets. At 5400 start with the commutative or matrix-block dimension equations. | D01–D02 regenerate the boundary lists. D26 verifies explicit arithmetic witnesses. The argument excluding the last five modular lists compares two quartet degrees and uses the tame-phase theorem in the paper. |
| Finite windows and all-height projections | Enumerate axis multiplicities, signed Gram entries, trace residues and inertia floors for the indicated prime positions and windows. The all-height run uses multiplicative periods after the size cap is reached. | D06–D19 record every tested endpoint and output. D16 and D17 implement the same projection separately. A surviving signature certifies only compatibility with these particular tests. |
| Seven-position calculation | Normalize the affine support and grow the opposite-pair mask. Scalar bounds reduce the problem to explicitly listed primes; separate four-pair and seven-position enumerations occur. | D24 and D25 give full input and output, including 7,648,354 normalized seeds at 127 and 33,726,047 at 167. D23 covers the distinct four-position calculation at 23. |

## Invariants at the two finite compressions

In the support calculation bit `i` represents an opposite pair `{x,−x}` in the multiplicative quotient by sign. Rotating a mask dilates the represented set. Canonicalization by the least rotation preserves all later tests, which range over all relative dilations. A mask in the closure queue is a necessary envelope for a constituent. After an intersection the allowed opposite differences are represented by the pairs of `x`, `2x`, and `x±y`. A size or overlap contradiction discards a seed. The seven-position program skips intersections of more than six pairs. It enqueues a new difference mask only if twice the size of the full signed set (including zero) is less than the prime. When the zero-based queue index exceeds 40, it returns “potentially feasible”. Omission of a consequence or reaching the work cap is inconclusive. Each seed reported as excluded has already incurred a contradiction. The Python four-position verifier reconstructs stored closure steps with explicit residues. Its certificate can be generated and rechecked separately:

```sh
python3 code/word_closure.py check 23 4 --ell 3 --write /tmp/wgt-23.json
python3 code/word_closure.py verify /tmp/wgt-23.json
```

In the four-prime computation an axis signature retains its norm cost and trace residue. A common-degree signature retains both endpoint coordinates and hence its two squared costs and signed Gram contribution. When types with the same signature have different mandatory multiplicities, retaining the **smallest** floor enlarges the allowed choices. Every genuine table therefore maps to a feasible state. In the unbounded exponent calculation the residue of the opposite prime power is periodic modulo the product of the two omitted primes. The size entering the inertia floor is kept exactly until it reaches its cap. The remaining tail is exhausted by one full period. A finite height box alone does not establish this all-height result. [Independent projection](../code/four_primes/independent_projection.py) · [Gram and trace recursion](../code/four_primes/unbounded_four.py)

## Work and limits

The dimension entries record observed wall times on one shared worker. Times depend on compiler and hardware. The rank entries give machine-independent work measures where available. A clean rank-30 stored-list check, including compilation and the terminal screens, took 113 seconds and peaked at 392 MiB on a shared Linux worker. Candidate verification reads compact archives into memory. Full traversal is a CPU computation and can take substantially longer, particularly for rank 30 (about 1.83 billion nodes). The number of workers is configurable. Checked native overflow transfers to the GMP kernel, or returns an incomplete run when the selected stand-alone program has no fallback. An interrupted job, unknown phase-mass test, or exhausted work cap never proves exclusion.

The [rank guide for the full article](https://github.com/sebastienpalcoux/WGT-fusion-categories/blob/independent-referee-2026-11-10/rank/README.md) records the status of each computation at publication. For the general bound of 32, the [self-contained rank paper](../papers/Rank32.pdf) and its [verification supplement](../papers/rank-bound-supplement/README.md) give the proof and the complete calculation.
