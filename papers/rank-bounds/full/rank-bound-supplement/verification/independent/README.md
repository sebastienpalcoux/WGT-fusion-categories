# Independent finite verifier

Run from the supplement root, with Python 3.10 or newer and SymPy:

```sh
python3 verification/independent/check.py
```

The checker imports no production modules. It reconstructs the five terminal
conditions for all 18,199 lists in `census/reference-census.zip`; its defect-one
test enumerates every permissible tiling instead of using the production
greedy recurrence. It also recalculates the odd starting counts, the last odd
nonsquare equation and finite group/permutation arithmetic. The cohomology
orbit calculations assume the groups and actions established in Section 10.

## Files and witness schema

- `check.py`: separate arithmetic implementation.
- `certificate.json`: deterministic receipt binding this verifier and the
  covering archive by SHA-256, with all table cells and finite-check counts.
- `terminal-witnesses.jsonl.gz`: one JSON record for every covering-census row.

Each record has `rank`, `task`, `row_sha256`, `stage` and `witness` fields.
`row_sha256` hashes the exact uncompressed input JSON line, without its final
newline. `stage` is the first failed condition in Table 1. Witnesses contain:

| Stage | Recorded evidence |
|---|---|
| `primary_mass` | Prime, minimum degree, eligible mass and the two integer sides of the violated inequality. |
| `prime_square_gcd` | Prime, squarefree kernel, row degree, all prime-free degrees and their gcd. |
| `defect_one` | Prime and prime-free denominator multiplicities; the separate exhaustive tiling recursion finds no tiling. |
| `principal_block` | Prime, complete eligible-degree set and its gcd; eligibility is recomputed from every required inequality. |
| `ternary_trace` | Modulus, eligible signed degrees, all attainable subset residues and the unattainable required residue. |

Running the checker regenerates these witnesses in memory and compares every
byte with the ledger. `--write` regenerates the receipt and ledger after an
intentional, reviewed change. `--archive PATH` selects another compatible
covering archive; `--fresh PATH` additionally compares emitted lists from a
fresh run with the covering corpus. The latter is an inclusion comparison,
not a verification of the fresh run's coverage.

## Verification boundary

This is an exact finite arithmetic audit, not a proof-assistant certificate.
It does not regenerate the native task partition or issue a witness for every
discarded prefix. The general rank conclusion also depends on the necessity
of the manuscript's structural reductions and on completeness of the
source-bound producer run. `census/audit.py` verifies the declared task and
receipt accounting separately. The odd-square continuation recursions are
replayed by `odd_square/check.py`, not by this checker.

The implementation is adapted from the independently supplied referee
verifier (original SHA-256
`91c7624bb72571acf986ce95fc74efd219021efd0aba88dd89a74f874b0bd6e9`).
The release adaptation adds input validation, explicit terminal witnesses,
deterministic source/input bindings and a portable command-line interface.
Its arithmetic algorithms and their scope were inspected before integration.
