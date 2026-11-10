# Independent rank shards for the ordinary modular case

The search takes a rank between 6 and 48 and returns a finite superset of candidate degree vectors for a nontrivial simple perfect integral modular category. Several proved WGT exponent boxes are removed at the two-term stage. The remaining vectors require the manuscript's subsequent categorical sieves. The upper accepted rank is a software parameter and carries no assertion about running time or completed proofs.

From the package root, use `python3 continue_ranks.py --ranks 24 --jobs 12`.
The archive to return is `rank-work/rank24/rank24_result.zip`. The default
split depth in that launcher is twelve. The exact partition can be prepared
without a full run using `python3 continuation/run_shards.py 24 --depth 12
--work rank-work/rank24 --prepare-only` from the package root. The coordinator
recompiles when necessary, regenerates the manifest, checks the per-shard
receipts and resumes missing work. Use a fresh work directory when changing
rank, depth or executable source. The command `python3
continuation/verify_result.py result.zip --replay` repeats the manifest,
all screening and the entire unsharded rank. Without `--replay` it checks
the returned data and reruns the screens. The source is selected by the
hash in the certificate. `rank20_search.cpp` preserves the code used on
the original laptop run.

## Search rules

Write the dimension as `q M²`, where `q` is square-free. A nonunit degree `x` contributes `1/b²` to `q`, with `b=M/x`. The vacuum contributes `1/M²`. Denominators of nonunits are sorted in ascending order with repetitions. The state records their lcm `L` and positive residual numerator `E` over `L²`. With `h` nonunit positions left, the next denominator lies above `sqrt(L²/E)` and at most `sqrt((h+1/36)L²/E)`. The upper bound uses the minimum nonunit degree six. After selecting `b`, put `g=gcd(L,b)`, `k=b/g`, and replace `(L,E)` by `(k L, k² E−(L/g)²)`. Formal codegrees below eight and prime codegrees are excluded as in the manuscript.

For each odd prime `p` of the global dimension, write `m=(p−1)/2`. Galois covariance and twist divisibility require a run of at least `m` equal denominators with `p|q b`. At even valuation, i.e. `p∤q`, either an eligible run has at least `2m` entries or there must be room for three whole packets of `m` in the eligible runs. Indeed a twist of order divisible by `p²` already provides at least `p m` equal labels. Otherwise the rational local Gauss sum equates the positive masses at all nonzero residues modulo `p`. At least one packet occurs on each square class. If only two packets are available, the two common squared dimensions must coincide. In particular a new prime outside `qL` must satisfy `p−1 ≤ h`. The prefix bound on further packets is `floor((h+r)/m)`, where `r` is the current eligible run length modulo `m` (zero if the final run is ineligible).

With two remaining nonunits and `k=M/L`, both degrees `d ≤ e` divide `L` and are prime to `k`. Indeed every prime of `k` must miss two nonunit degrees, while it divides all prefix degrees. Thus `E k²=1+d²+e²`. A prime 2 or 5 in `k` is impossible by reduction modulo 4 or 5. A prime at least 5 outside `qL` cannot enter with only two nonunit labels left. If possible prime support has size at most two, the two-prime WGT theorem applies. On a support of size three the heights at 2 and 5 are already fixed, and the ordered three-prime bounds from BDP can be applied. This also covers smaller actual support. For a short denominator interval the program enumerates the larger tail degree using `b=M/e`. As deleting the last nonunit leaves the lcm unchanged, `k=b/gcd(b,L)` and `e=L/gcd(b,L)`. For a long interval it pairs divisors of `L` by `e² ≡ −1−d² (mod E)`, then tests the square quotient exactly. All denominator, root, lcm, residual, exponent-test, square and counter arithmetic with unbounded inputs uses GMP integers; bit masks and packet counters are bounded by the explicitly checked rank range.

For a small independent partition check, run `run_shards.py 16 --depth 4 --work check16 --jobs 2`. The closed rank-16 example produces zero vectors after the built-in filters. The equivalent monolithic command is `check16/search --all 16`.

## Small formal codegrees

The continuing search applies the proved multiplicity caps of two at eight and four at nine for the largest degree. An even largest degree at codegree eight occurs exactly twice. A degree at codegree nine is prime to three. In either case at least four distinct dimension primes are necessary. If the prefix lcm is L, any new denominator prime outside qL is at most h+1 with h positions left. These prefix tests are enabled in full at rank 24 and above. See the manuscript section on large multiplicities for proofs.
