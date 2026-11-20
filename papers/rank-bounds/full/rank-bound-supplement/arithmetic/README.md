# Arithmetic in the proofs

Run `python3 verify.py --quick` from the supplement directory for the bounded rank audit. It checks receipts, selected exact inequalities, the terminal census and the necessary-condition screens without replaying the final odd continuations. The default `python3 verify.py` also replays both complete production odd certificates; `--full` adds the auxiliary arithmetic replays. Denominator generation is separate. All proof calculations use integers or exact fractions. The [proof map](../proof-map.md) states each calculation's hypotheses, including the distinct rank scopes in the combined high-four proposition and the actual-row, full prime-free support, unit-weight conditions for the full-spectrum projection.

| Source | Calculation |
|---|---|
| [Block traces](../census/arithmetic/block_trace.py), [comparison replay](../census/recheck.py) | Simultaneous principal-block trace congruences on every list in the complete reference census, with conservative resource guards |
| [Neutral-orbit checker](../verification/neutral_orbits/check.py) | Odd isolated-prime active costs, disjoint moved-neutral reservations, and mixed-component neutral orbits; see the [focused audit](../verification/neutral_orbits/README.md) for the saved 201 and exploratory 209 inputs |
| [three_prime_cutoff.py](three_prime_cutoff.py) | Three-prime square support using nine cutoff profiles and 18,095 recursion nodes |
| [support_five_moments.py](support_five_moments.py) | Five-prime support by pairwise Weyl moments and exact mass intervals |
| [odd_edge_mass.py](odd_edge_mass.py) | Exclusion of the odd edge of product 15 |
| [gram_packets.py](gram_packets.py) | All 1,440 matrices on the joint five- and seven-packets |
| [three_prime_rank.py](three_prime_rank.py), [three_primes.py](three_primes.py) and [full_census.py](full_census.py) | Three-prime square support through rank 44: row quotas, inverse-residue bounds, and exact rational intervals |
| [support_five.py](support_five.py) | Separated five-prime square supports through rank 44: simultaneous proper supports and exact mass/Gauss bounds; six-prime bound through rank 47 |
| [large_four.py](large_four.py) | Support choices and omission bounds when denominator four occurs at least nine times. The [signature generator](frontier_cases.py) covers ranks through 43 |
| [phase_deficits.py](phase_deficits.py) | Positive Gauss deficits on isolated primary components |
| [gauss_divisors.py](gauss_divisors.py) | Signed divisor identities and the exact degree-square check at the two possible roots |
| [odd_nonsquare.py](odd_nonsquare.py) | Multiplicity bounds, effective Weyl columns and degree constraints in odd nonsquare dimension |
| [odd_square.py](odd_square.py), [edges.py](edges.py) and [edge_budget.py](edge_budget.py) | Odd square dimensions through rank 201, including the prime-row projection, inertia and the two phase moments |
| [full_blocks.py](full_blocks.py), [dense_budget.py](dense_budget.py), [dense_boundary.py](dense_boundary.py) and [small_phases.py](small_phases.py) | Prime-free multiplicities, ordinary mass on isolated primes and occupancies at denominator 15 |
| [isolated_gauss.py](isolated_gauss.py) | Separated Gauss sums, denominator divisibility and primitive-phase parity |
| [joint_gauss.py](joint_gauss.py), [nine_gram.py](nine_gram.py) and [ordinary_tail.py](ordinary_tail.py) | Simultaneous primitive-phase differences, the Gram bound at nine and ordinary mass at a signed Gauss root |
| [perfect_counts.py](perfect_counts.py) | Perfect classification through rank 44, including the prime-power and neutral-orbit lists, cohomology orbit counts, and the sector ranks and global dimension of the deleted permutation gaugings |

In the isolated-component recursion, two equal dyadic charges \(1/d\) with \(d>32\) carry to one charge \(1/(d/2)\). The replacement uses one fewer label and its underlying degree denominator divides the previous one, so sorting after all carries preserves the attainable dimension roots. The verifier allows repeated charges at the smallest level and retains unrestricted multiplicities when a separate quota is imposed at denominator eight. [Gauss deficit recursion](phase_deficits.py)

The source-bound auxiliary odd-square check covers rank 201: 13,013 prime profiles lead to 750 isolated and 174 mixed reservations, all excluded. Its receipt binds 15 production source files and records every indexed outcome. [odd_square.py](odd_square.py) is the complete auxiliary entrypoint. The optional comparison in `verification/odd_uniform/` applies additional inertia, denominator-floor and norm conditions to these source-bound inputs. The rank bound uses the single rank-209 ceiling in [the square certificate](../odd_square/certificate.json). The exact state and charge disjunction, reflected phase moments and outward-rounded bounds are proved in the manuscript. These arithmetic reservations are not candidate fusion rings or modular data.

From the supplement directory, use the exact native continuation with a C++17 compiler and GMP development files:

```sh
python3 build_odd_native.py
ODD_NATIVE=1 ODD_NATIVE_NODES=0 python3 arithmetic/odd_square.py --jobs 10
```

Set `ODD_NATIVE=0` for the slower Python recursion. Reaching a resource guard is unresolved and never certifies an exclusion. The [receipt audit](../verify_odd_receipts.py) regenerates the reservations, checks every source fingerprint and validates every indexed outcome without repeating the final search.
