#!/usr/bin/env python3
"""Independent exact audit of local Cayley--Hamilton coefficient tests.

Only four previously excluded control pools are read.  No released module is
imported.  Field products use explicit cyclotomic polynomials; field traces
sum distinct Galois embeddings rather than using the explorer's trace code.
This script does not claim completeness of a row pool or a new exclusion.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json


POLYNOMIALS = {
    1: (-1, 1),
    3: (1, 1, 1),
    4: (1, 0, 1),
    9: (1, 0, 0, 1, 0, 0, 1),
    15: (1, -1, 0, 1, -1, 1, 0, -1, 1),
}


class Field:
    def __init__(self, meta):
        self.c = meta['conductor']
        self.p = POLYNOMIALS[self.c]
        self.n = len(self.p) - 1
        self.h = frozenset(meta['subgroup'])
        self.degree = meta['degree']
        self.one = (1,) + (0,) * (self.n - 1)
        self.zero = (0,) * self.n
        self.powers = tuple(self.reduce((0,) * i + (1,)) for i in range(self.c))
        units = {i for i in range(self.c) if gcd(i, self.c) == 1}
        assert self.h <= units and (1 % self.c) in self.h
        assert all((a * b) % self.c in self.h for a in self.h for b in self.h)
        reps = []
        while units:
            a = min(units)
            reps.append(a)
            units -= {(a * h) % self.c for h in self.h}
        self.reps = tuple(reps)
        assert len(reps) == self.degree

    def reduce(self, values):
        v = list(values) + [0] * max(0, self.n - len(values))
        for i in range(len(v) - 1, self.n - 1, -1):
            a = v[i]
            for j in range(self.n):
                v[i - self.n + j] -= a * self.p[j]
        return tuple(v[:self.n])

    @lru_cache(None)
    def mul(self, a, b):
        out = [0] * (2 * self.n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return self.reduce(out)

    @lru_cache(None)
    def auto(self, a, u):
        return tuple(sum(x * self.powers[(i * u) % self.c][j]
                         for i, x in enumerate(a)) for j in range(self.n))

    def star(self, a):
        return self.auto(a, (self.c - 1) % self.c)

    @lru_cache(None)
    def trace(self, a):
        assert all(self.auto(a, h) == a for h in self.h)
        total = tuple(sum(self.auto(a, u)[j] for u in self.reps)
                      for j in range(self.n))
        assert total[1:] == (0,) * (self.n - 1), total
        return total[0]


def pack(r):
    return (r['dimension'], r['matrix_trace'], r['matrix_norm'],
            r['selfdual'], tuple(tuple(v) for v in r['values']))


def audit_pool(path):
    data = json.loads(path.read_text())
    target = data['target']
    assert target['kind'] == 'symplectic'
    D, F = target['dimension'], target['matrix_codegree']
    fields = [Field(o) for o in data['orbits']]
    codegrees = [o['formal_codegree'] for o in data['orbits']]
    assert all(D % f == 0 for f in codegrees)
    assert Q(1, D) + sum(Q(f.degree, c) for f, c in zip(fields, codegrees)) + Q(2, F) == 1
    assert 1 + sum(f.degree for f in fields) == target['rank'] - 4
    unit = pack(data['unit'])
    assert unit == (1, 2, 2, True, tuple(f.one for f in fields))

    @lru_cache(None)
    def star(x):
        return x[:4] + (tuple(f.star(v) for f, v in zip(fields, x[4])),)

    @lru_cache(None)
    def scalar(x, y):
        # Bilinear rather than Hermitian pairing; all embeddings are present.
        return Q(x[0] * y[0], D) + sum(
            Q(f.trace(f.mul(a, b)), c)
            for f, c, a, b in zip(fields, codegrees, x[4], y[4]))

    @lru_cache(None)
    def cubic(x, y, z):
        return Q(x[0] * y[0] * z[0], D) + sum(
            Q(f.trace(f.mul(f.mul(a, b), v)), c)
            for f, c, a, b, v in zip(fields, codegrees, x[4], y[4], z[4]))

    formula_checks = 0

    @lru_cache(None)
    def square(x, y, same=False, is_unit=False):
        nonlocal formula_checks
        det = Q(x[2], 2)
        # Direct matrix CH plus coefficient-trace orthogonality.
        direct = cubic(x, x, star(y)) + x[1] * (int(same) - scalar(x, star(y))) - det * Q(y[1], F)
        # Independently evaluate the centralized polynomial used by attachment.
        perron = (x[0] ** 2 - x[1] * x[0] + det) * y[0]
        central = Q(perron, D)
        for f, c, a, b in zip(fields, codegrees, x[4], y[4]):
            v = f.mul(a, a)
            v = tuple(vj - x[1] * aj + det * ej for vj, aj, ej in zip(v, a, f.one))
            central += Q(f.trace(f.mul(v, f.star(b))), c)
        central += x[1] * int(same) - det * int(is_unit)
        assert direct == central
        formula_checks += 1
        return direct

    @lru_cache(None)
    def norm(x, y):
        return cubic(x, star(x), star(y)) + Q(x[2] * y[1], 2 * F)

    def fixed_pair(x, y):
        return scalar(x, y) + scalar(x, star(y)) + Q(x[1] * y[1], F)

    raw = {r['id']: r for r in data['rows']}
    groups = defaultdict(dict)
    determinants = 0
    for r in data['rows']:
        x = pack(r)
        assert r['dual_id'] in raw
        rd = raw[r['dual_id']]
        assert rd['dual_id'] == r['id'] and pack(rd) == star(x)
        if x[3]:
            assert r['dual_id'] == r['id']
        assert scalar(x, unit) + Q(x[1], F) == 0
        assert scalar(x, star(x)) + Q(x[2], F) == 1
        assert fixed_pair(x, x) == (2 if x[3] else 1)
        recovered = (x[1] ** 2 - F * (int(x[3]) - scalar(x, x))) / 2
        assert recovered == Q(x[2], 2) and recovered.denominator == 1
        assert square(x, unit, False, True) == int(x[3])
        determinants += 1
        if not x[3] and r['dual_id'] < r['id']:
            continue
        fixed = tuple(tuple(a + b for a, b in zip(v, f.star(v)))
                      for f, v in zip(fields, x[4]))
        key = (x[3], x[0], x[1], fixed)
        groups[key][(x[2], x[4])] = x
    classes = [list(v.values()) for _, v in sorted(groups.items())]

    def valid(z, x, y, stronger=False):
        if z.denominator != 1 or z < 0 or z * y[0] > x[0] ** 2:
            return False
        if stronger and (z > y[0] or z * y[0] > x[0] ** 2 - int(x[3])):
            return False
        return True

    def unary(x, stronger=False):
        u = square(x, x, True)
        v = square(x, star(x), x[3])
        assert norm(x, x) == u == norm(x, star(x))
        if not (valid(u, x, x, stronger) and valid(v, x, star(x), stronger)):
            return False
        return (not stronger or x[0] * (u + (0 if x[3] else v)) <= x[0] ** 2 - int(x[3]))

    old = [[x for x in c if unary(x)] for c in classes]
    strong = [[x for x in c if unary(x, True)] for c in classes]
    baseline = [(i, j) for i, j in combinations(range(len(classes)), 2)
                if fixed_pair(classes[i][0], classes[j][0]) == 0]

    def orientations(x):
        # A nonselfdual pair contributes TWO labels even if values coincide.
        return [x] if x[3] else [x, star(x)]

    def pair_good(x, y, stronger=False, include_norm=False):
        for a, b in ((x, y), (y, x)):
            for aa in orientations(a):
                for bb in orientations(b):
                    if not valid(square(aa, bb), aa, bb, stronger):
                        return False
                    if include_norm and not valid(norm(aa, bb), aa, bb, stronger):
                        return False
                if stronger:
                    own = a[0] * (square(aa, aa, True) +
                        (0 if a[3] else square(aa, star(aa), False)))
                    other = b[0] * sum(square(aa, bb) for bb in orientations(b))
                    if own + other > a[0] ** 2 - int(a[3]):
                        return False
        return True

    old_unary_edges = [(i, j) for i, j in baseline if old[i] and old[j]]
    old_edges = [(i, j) for i, j in old_unary_edges
                 if any(pair_good(x, y) for x, y in product(old[i], old[j]))]
    both_edges = [(i, j) for i, j in old_unary_edges
                  if any(pair_good(x, y, include_norm=True) for x, y in product(old[i], old[j]))]
    strong_unary_edges = [(i, j) for i, j in baseline if strong[i] and strong[j]]
    strong_edges = [(i, j) for i, j in strong_unary_edges
                    if any(pair_good(x, y, True) for x, y in product(strong[i], strong[j]))]
    expected = {
        1284815: (90, 114, 3, 323, 1, 1),
        461210: (167, 239, 1, 3393, 0, 0),
        47412: (285, 319, 31, 5522, 79, 62),
        836343: (76, 82, 1, 788, 0, 0),
    }
    got = (len(classes), sum(map(len, classes)), sum(bool(c) for c in old),
           len(baseline), len(old_unary_edges), len(old_edges))
    assert got == expected[target['state_id']], (got, expected[target['state_id']])
    return dict(state_id=target['state_id'], pool_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                matrix_codegree=F, raw_rows=len(raw), determinant_checks=determinants,
                fixed_classes=len(classes), variants=sum(map(len, classes)),
                square_unary_classes=sum(bool(c) for c in old),
                square_unary_variants=sum(map(len, old)), baseline_edges=len(baseline),
                edges_after_unary=len(old_unary_edges), edges_after_square_pairs=len(old_edges),
                edges_after_square_and_norm_pairs=len(both_edges),
                stronger_unary_classes=sum(bool(c) for c in strong),
                stronger_unary_variants=sum(map(len, strong)),
                stronger_edges_after_unary=len(strong_unary_edges),
                stronger_edges_after_pairs=len(strong_edges),
                formula_equivalence_checks=formula_checks,
                status='PASS_INDEPENDENT_BOUNDED_CONTROL_AUDIT',
                additional_target_exclusion=False)


def det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
          - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
          + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def check_group(name, elements, mul, inv, linear, trace, F):
    D = len(elements)
    def A(x, y):
        return sum(Q(a*b, D) for a, b in zip(linear(x), linear(y)))
    def H(x, y):
        return 2 * F * (int(x == inv(y)) - A(x, y)) - trace(x) * trace(y)
    count = 0
    for x, y, z in product(elements, repeat=3):
        w = inv(z)
        L = sum(Q(a*b*c, D) for a, b, c in zip(linear(x), linear(y), linear(w)))
        J = (2*L + trace(x) * (int(y == z) - A(y, w))
             + trace(y) * (int(x == z) - A(x, w))
             + trace(z) * (int(x == inv(y)) - A(x, y))
             - Q(trace(x) * trace(y) * trace(z), F))
        triple = (x, y, w)
        delta = -det3([[H(a, b) for b in triple] for a in triple]) / (4*F*F)
        p, q = int(mul(x, y) == z), int(mul(y, x) == z)
        assert J == p + q and delta == (p - q)**2
        assert delta.denominator == 1 and delta >= 0
        root = isqrt(delta.numerator)
        assert root*root == delta and sorted(((J+root)/2, (J-root)/2)) == sorted((p, q))
        assert J >= root and (J - root) % 2 == 0
        count += 1
    return dict(group=name, exact_jordan_and_commutator_checks=count, status='PASS')


def group_checks():
    s3 = tuple(permutations(range(3)))
    def smul(p, q): return tuple(p[q[i]] for i in range(3))
    def sinv(p): return tuple(p.index(i) for i in range(3))
    def sign(p): return (-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
    q8 = tuple(product((1, -1), range(4)))
    def qmul(x, y):
        a, i = x; b, j = y
        if i == 0: return a*b, j
        if j == 0: return a*b, i
        if i == j: return -a*b, 0
        return a*b*(1 if (i,j) in ((1,2),(2,3),(3,1)) else -1), 6-i-j
    def qinv(x): return x if x[1] == 0 else (-x[0], x[1])
    def qlin(x):
        return tuple(1 if x[1] == 0 else a if x[1] == 1 else b if x[1] == 2 else a*b
                     for a,b in product((1,-1), repeat=2))
    return [check_group('S3', s3, smul, sinv, lambda p:(1,sign(p)),
                        lambda p:sum(p[i] == i for i in range(3))-1, 3),
            check_group('Q8', q8, qmul, qinv, qlin,
                        lambda x:2*x[0] if x[1] == 0 else 0, 4)]


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--explorer', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    results = [audit_pool(p) for p in sorted((args.explorer/'data/known_results').glob('state_*/pool.json'))]
    assert {r['state_id'] for r in results} == {1284815,461210,47412,836343}
    report = dict(status='PASS', arithmetic='Exact rational arithmetic and explicit Galois embeddings',
                  imports_released_modules=False, formal_proof_assistant_verification=False,
                  scope='Four already-excluded controls; no pending target pools or new target exclusions',
                  controls=results, general_coefficient_identity_controls=group_checks())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
