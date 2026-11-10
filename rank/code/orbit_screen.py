"""Check that each degree multiplicity is a sum of admissible Galois orbit sizes."""
from collections import Counter
from itertools import combinations

from trace_rows import nonunit_trace_row


def possible_orbits(degrees, x, m, f, principal):
    return [h for h in range(1, min(f-1, m) + 1)
            if nonunit_trace_row(degrees, x, specified_h=h, principal=principal)]


def first_obstruction(degrees, principal=True):
    total = sum(d*d for d in degrees)
    for x, m in sorted(Counter(degrees).items(), reverse=True):
        f = total // (x*x)
        if x == 1 or f > 120:
            break
        sizes = possible_orbits(degrees, x, m, f, principal)
        reachable = {0}
        for j in range(1, m+1):
            if any(j-h in reachable for h in sizes):
                reachable.add(j)
        if m not in reachable:
            return x, f, m, sizes
    return None


def small_codegree_obstruction(degrees):
    """Apply the categorical bounds on every degree of codegree eight or nine."""
    total = sum(d*d for d in degrees)
    for degree, multiplicity in Counter(degrees).items():
        f = total // (degree*degree)
        if f == 9 and (degree % 3 == 0 or multiplicity > 4):
            return f
        if f == 8:
            if multiplicity > 4 or (degree == max(degrees) and multiplicity > 2):
                return f
            if degree % 2 == 0:
                if multiplicity % 2 or (degree == max(degrees) and multiplicity != 2):
                    return f
                eligible = [d for d in degrees if d > 1 and (degree//2) % d == 0]
                if not any(abs(u + v - w) == 1 for u, v, w in combinations(sorted(eligible), 3)):
                    return f
    return None


def screen_survivors(rows, principal=True, local_bounds=True):
    excluded = Counter()
    survivors = []
    for row in rows:
        small = small_codegree_obstruction(row) if local_bounds else None
        if small is not None:
            excluded[str(small)] += 1
            continue
        rejection = first_obstruction(row, principal)
        if rejection is None:
            survivors.append(row)
        else:
            excluded[str(rejection[1])] += 1
    return {'excluded_at_codegree': dict(sorted(excluded.items())), 'survivors': survivors}
