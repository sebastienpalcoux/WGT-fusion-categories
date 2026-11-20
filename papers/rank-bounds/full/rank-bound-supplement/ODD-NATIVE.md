# Exact odd-dimensional continuation

The native continuation implements the same exact integer recursion as `joint_gauss.py` on supports (3,5,7) and (3,5,7,11) with the marked-eleven class. It uses GMP integers throughout. Python generates the exact admissible states and evaluates every terminal ordinary-tail condition. The bridge partitions large state families into bounded batches; all admissible states and Gauss charges are retained.

From the supplement directory, build with a C++17 compiler and GMP development files:

```sh
python3 build_odd_native.py
ODD_NATIVE=1 ODD_NATIVE_NODES=0 python3 arithmetic/odd_square.py --jobs 10
```

Alternatively, `python3 verify.py --full --jobs 10 --native-odd` rebuilds the helper and runs the complete proof arithmetic. Omit `--native-odd`, or set `ODD_NATIVE=0` for the direct script, to use the slower Python recursion.

`ODD_NATIVE_NODES` may be set to a positive node limit per batch. Reaching that guard raises an unresolved-case error; it never proves exclusion. A missing executable, malformed protocol or interrupted continuation is likewise unresolved.

The archived rank-201 replay excluded all 924 reservations in 365.20 seconds with three workers. The receipt hashes the C++ source, Python bridge, build helper and remaining dependencies, as well as every indexed outcome. The generated executable is platform-specific and is not distributed. Run `python3 verify_odd_receipts.py` to check the receipt and regenerate the reservation list without repeating the final searches.
