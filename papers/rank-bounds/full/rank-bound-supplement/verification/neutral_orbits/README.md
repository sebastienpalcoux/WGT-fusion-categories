# Neutral-orbit checks for saved odd-square reservations

From the supplement root, run:

```sh
python3 verification/neutral_orbits/check.py
```

This short, exact check reads saved arithmetic reservations, applies the
neutral-orbit necessary conditions, and compares the result with the two
deterministic receipts in this folder. It uses only the Python standard
library. It runs neither profile regeneration nor a final Gauss search.

| Rank | Saved reservations | Active-cost exclusions | Neutral-field exclusions | Mixed-neutral exclusions | Retained |
|---:|---:|---:|---:|---:|---:|
| 201 | 924 | 32 | 54 | 4 | 834 |
| 209 | 1,401 | 62 | 69 | 2 | 1,268 |

These rows count arithmetic reservations. A retained reservation is not a
fusion ring or modular datum. In particular, the second row does not certify
a complete modular-data census or an exclusion at rank 209. These auxiliary screens do not replace the two complete rank-209
continuations used for the odd rank bound.

## Inputs and provenance

Rank 201 uses the source-bound
[`arithmetic/odd-square-cases.json`](../../arithmetic/odd-square-cases.json).
Its entries and ordering are bound to the saved profile output. The
rank-201 execution receipt records the complete final calculation; this
neutral-orbit screen does not repeat that calculation.

[`rank209-input.json`](rank209-input.json) is the saved exploratory support and
multiplicity output at rank 209. It is bound byte-for-byte to the saved execution output. The 16,512 recorded prime profiles yielded 1,401 saved reservations.
No final `odd_square.check_case` or Gauss verification at rank 209 is asserted.

The saved run imported the rank-201 module, set `odd_square.RANK = 209`, and
called `odd_square.run()`. Its import-time `PRIMES` list still ended at 101,
as initialized at rank 201. This does not omit a rank-209 prime-support
profile: the only additional prime in the updated range is 103, and a
support containing it requires at least three distinct odd primes. The
support bound gives rank at least `2*(3+5+103)-1 = 221`, exceeding 209.
Consequently both prime ranges give the same arithmetic profile output
at this rank.

[`provenance.json`](provenance.json) binds both inputs, the rank-201 receipt,
and its fifteen production source files by SHA-256. The two derived receipts
also bind this checker and the provenance file. An ordinary invocation only
reads and checks these files. `--rank 201` or `--rank 209` selects one input.
`--supplement-root PATH` selects an explicitly located supplement.
After an intentional checker change, `--write-receipts` regenerates only
this folder's derived receipts; it does not alter the source-bound execution certificates.

### Optional profile reproduction

The following explicit command updates both the rank and the prime list,
regenerates only the arithmetic support and multiplicity profiles, and checks
their serialization against the saved input. It writes a temporary file,
leaving the saved input and all receipts unchanged. It performs no final
Gauss search. This is separate from the default short receipt check above.

```sh
python3 - <<'PY'
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path("arithmetic").resolve()))
import odd_square

odd_square.RANK = 209
odd_square.PRIMES = tuple(
    p for p in range(3, (odd_square.RANK + 1) // 2 + 1, 2)
    if all(p % d for d in range(3, isqrt(p) + 1, 2))
)
count, rows = odd_square.run()
payload = json.dumps({"count": count, "rows": rows}, indent=1).encode("utf-8")
saved = Path("verification/neutral_orbits/rank209-input.json").read_bytes()
if payload != saved:
    raise RuntimeError("Regenerated profiles differ from the saved rank-209 input")
with tempfile.NamedTemporaryFile(prefix="odd-rank209-profiles-",
                                 suffix=".json", delete=False) as temporary:
    temporary.write(payload)
    temporary_path = temporary.name
print("Prime profiles:", count)
print("Saved reservations:", len(rows))
print("Temporary output:", temporary_path)
print("SHA-256:", hashlib.sha256(payload).hexdigest())
PY
```

The serialization intentionally has no final newline, matching the saved
input's SHA-256 fingerprint.

## Mathematical conditions and counting

The script uses quadratic-field incidence, neutral inertia, the fourfold
active-orbit obstruction, and the mixed neutral-orbit bound proved in the
manuscript. At an isolated absent prime, writing
`d = 2^v_2(p-1)`, it reserves at least `d*(p-1)` active labels. An optimistic
mass test determines when the active count must be below `p*(p-1)`; only then
does it reserve the additional moved or pure neutral labels supplied by the
sparse-activity argument.

For a present valuation-one prime denominator, the complete marked class is
already booked. The script subtracts its neutral members before reserving
any additional moved neutral labels. The field reservations for different
isolated primes are disjoint by field localization.

Two mass tests are deliberately separate: the first combines the base
reservations with moved neutral fields; the second combines the base with
pure neutral labels and mixed neutral labels. They never add the moved-field
and pure-neutral counts together. Pure neutral labels have denominator
supported on one component; mixed neutral labels meet at least two components.
The 80-label ternary-edge reservation consists of labels active at 3 or 5,
and is therefore disjoint from the all-neutral reservations.

All residual denominators have floor 15 because the full prime-denominator
classes and the full denominator-nine class have already been booked. Every
mass comparison uses exact fractions. Unhandled twist graphs are retained.
Malformed inputs or altered source bindings cause an error, not an exclusion.

The fourfold obstruction is available in the checker but is not invoked by
any marked class in these two saved input lists. Its presence therefore adds
no exclusions to the counts above.
