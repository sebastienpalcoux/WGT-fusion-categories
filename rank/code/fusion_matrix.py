"""Necessary multiplication-matrix constraints for simple perfect integral tables.

Input is a JSON list of degree types, a JSON record with a "survivors" list,
or JSONL with one type per line. Pure Python; no solver or floating point.
"""
from bisect import bisect_right
from collections import Counter
from functools import reduce
from itertools import combinations_with_replacement
from math import gcd


def weighted_sums(terms, target):
    """Bit j records reachability of j using one choice from each finite domain."""
    reachable = 1
    mask = (1 << (target + 1)) - 1
    for weight, domain in terms:
        if len(domain) == 1:
            continue  # The value is zero after the lower-bound shift.
        if domain[-1] == len(domain) - 1:
            remaining = domain[-1]
            step = 1
            while remaining:
                part = min(step, remaining)
                reachable |= reachable << (part * weight)
                reachable &= mask
                step *= 2
                remaining -= part
        else:
            old = reachable
            reachable = 0
            for value in domain:
                reachable |= old << (value * weight)
            reachable &= mask
    return reachable


def project(domains, weights, target, budget_limit, size_limit):
    """Bounded exact coordinate projections for a single weighted sum."""
    lower = [values[0] for values in domains]
    remaining = target - sum(w * lo for w, lo in zip(weights, lower))
    if remaining < 0:
        return None
    free = [i for i, values in enumerate(domains) if len(values) > 1]
    if not free:
        return domains if remaining == 0 else None
    divisor = reduce(gcd, (weights[i] for i in free))
    if remaining % divisor:
        return None
    remaining //= divisor
    if remaining > budget_limit:
        return domains
    terms = {}
    for i in free:
        w = weights[i] // divisor
        stop = bisect_right(domains[i], lower[i] + remaining // w)
        terms[i] = w, [c - lower[i] for c in domains[i][:stop]]
    if not (weighted_sums(terms.values(), remaining) >> remaining & 1):
        return None
    result = list(domains)
    for i, (w, values) in terms.items():
        if len(values) > size_limit:
            continue
        other = weighted_sums((term for j, term in terms.items() if j != i), remaining)
        result[i] = [k + lower[i] for k in values if other >> (remaining - w*k) & 1]
        if not result[i]:
            return None
    return result


def sum_domain(domain, copies, limit):
    """Possible sums of 'copies' elements of a domain, with repetitions."""
    if copies == 0:
        return {0}
    if len(domain) == domain[-1] - domain[0] + 1:
        return set(range(copies*domain[0], min(limit, copies*domain[-1])+1))
    values = {0}
    for _ in range(copies):
        values = {a+b for a in values for b in domain if a+b <= limit}
    return values


def check_generator(mult, x, depth=0, budget_limit=7000, size_limit=100, passes=5):
    ds = tuple(sorted(mult))
    targets = {a: x*a*mult[a] - (a == x) for a in ds}

    def key(a, b):
        return (a,b) if a <= b else (b,a)

    def bound(a, b):
        return min(x*a//b, x*b//a, a*b//x)

    cells = {}
    for a,b in combinations_with_replacement(ds, 2):
        hi = min(targets[a]//b, targets[b]//a, mult[a]*mult[b]*bound(a,b))
        cells[a,b] = range(hi+1)

    def propagate(cells, levels):
        for _ in range(passes):
            changed = False
            for a in ds:
                line = [cells[key(a,b)] for b in ds]
                result = project(line, ds, targets[a], budget_limit, size_limit)
                if result is None:
                    return "degree-class sum"
                if mult[a] > 1:
                    local = [range(min(values[-1], mult[b]*bound(a,b))+1)
                             for b,values in zip(ds,result)]
                    ordinary = project(local, ds, x*a, budget_limit, size_limit)
                    if ordinary is None:
                        return "single-label sum"
                    if a == x:
                        adjoint = project(local, ds, x*a-1, budget_limit, size_limit)
                        if adjoint is None:
                            return "adjoint sum"
                    for i, values in enumerate(result):
                        if len(values) > size_limit:
                            continue
                        possible = sum_domain(ordinary[i], mult[a] - int(a == x), values[-1])
                        if a == x:
                            possible = {v+c for v in possible for c in adjoint[i] if v+c <= values[-1]}
                        result[i] = sorted(possible.intersection(values))
                        if not result[i]:
                            return "single-label lift"
                for b,values,new in zip(ds,line,result):
                    if len(new) < len(values):
                        cells[key(a,b)] = new
                        changed = True
            if not changed:
                break

        # The unit is connected to exactly the class x. Every label must be
        # generated by X and hence reachable in this relaxed support graph.
        reached = {x}
        work = [x]
        while work:
            a = work.pop()
            for b in ds:
                if b not in reached and cells[key(a,b)][-1] > 0:
                    reached.add(b)
                    work.append(b)
        if len(reached) < len(ds):
            return "disconnected support"

        if levels:
            choices = [(2*(a//x + b//x) + len(values), a, b)
                       for (a,b), values in cells.items() if 1 < len(values) <= 8]
            if choices:
                _,a,b = min(choices)
                for value in cells[a,b]:
                    branch = cells.copy()
                    branch[a,b] = (value,)
                    if propagate(branch, levels-1) is None:
                        break
                else:
                    return "all integer branches fail"
        return None

    return propagate(cells, depth)


def first_obstruction(degrees, max_degree=80, depth=4):
    mult = Counter(degrees)
    if any(d < 1 for d in mult) or mult.pop(1,0) != 1:
        raise ValueError("The degrees must be positive and the unit unique.")
    generators = [x for x in sorted(mult) if x <= max_degree]
    for d in (0, depth) if depth else (0,):
        for x in generators:
            reason = check_generator(mult, x, d)
            if reason:
                return x, reason
    return None


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--depth', type=int, default=4)
    args = parser.parse_args()
    content = args.input.read_text()
    try:
        data = json.loads(content)
    except json.JSONDecodeError:
        data = [json.loads(row) for row in content.splitlines() if row.strip()]
    if isinstance(data, dict):
        data = data['survivors']
    if data and isinstance(data[0], int):
        data = [data]
    reasons = [first_obstruction(ds, depth=args.depth) for ds in data]
    print(json.dumps(dict(rows=len(data), excluded=sum(t is not None for t in reasons),
                          verdicts=reasons), indent=2))
