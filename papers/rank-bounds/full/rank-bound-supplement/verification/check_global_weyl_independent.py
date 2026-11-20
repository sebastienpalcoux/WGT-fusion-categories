#!/usr/bin/env python3
"""Independent finite exponent checks; no rank search or profile generation."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
from time import process_time

HERE = Path(__file__).resolve().parent
ROOT = HERE/'odd_uniform'
sys.path.insert(0, str(ROOT))
import global_weyl


def actual_signs(j7, j11, powers3, powers5):
    result = []
    for a, b in product(powers3, powers5):
        M = 3**a*5**b*7*11
        residues = j7*(M//7) % 7, j11*(M//11) % 11
        if residues in ((1, 1), (6, 10)):
            result.append((a, b))
    return result


def main():
    start = process_time()
    targets = heights = solutions = 0
    for j7, j11 in product((2, 4, 6), (2, 4, 6, 8, 10)):
        marked = ((7, j7), (11, j11))
        actual = actual_signs(j7, j11, range(1, 31), range(1, 31))
        assert global_weyl.possible_targets((3, 5, 7, 11), marked) == bool(actual)
        targets += 1
        for low, high in ((1, 1), (1, 3), (2, 2), (2, 5), (5, 8)):
            actual = actual_signs(j7, j11, range(low, high+1), range(1, 31))
            assert global_weyl.possible_targets((3, 5, 7, 11), marked, 3, high, low) == bool(actual)
            targets += 1
        row = [[3, 5, 7, 11], [[0, 0, 0], [1, 6*j7+2, 6*j7],
               [1, 10*j11+2, 10*j11]], [], 44]
        actual = actual_signs(j7, j11, range(2, 32), range(2, 32))
        failure, bounds = global_weyl.global_constraints(row)
        assert (failure is not None) == (not actual)
        if actual:
            for q, column in ((3, 0), (5, 1)):
                assert bounds.get(q, 2) <= min(pair[column] for pair in actual)
                heights += 1
            solutions += len(actual)
    output = {'status': 'passed', 'scope': 'Independent small finite exponent groups only.',
        'target_comparisons': targets, 'height_lower_bound_comparisons': heights,
        'actual_shared_sign_solutions_checked': solutions,
        'module_sha256': sha256((ROOT/'global_weyl.py').read_bytes()).hexdigest(),
        'audit_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'cpu_seconds': process_time()-start}
    Path(__file__).with_suffix('.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
