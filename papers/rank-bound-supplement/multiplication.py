"""Necessary integral constraints for multiplication by a smallest nonunit."""
from collections import Counter


def multiplication_obstruction(degrees):
    multiplicity = Counter(degrees[1:])
    values = sorted(multiplicity)
    x = values[0]
    targets = {a: x*a*multiplicity[a] - int(a == x) for a in values}

    def entry_bound(a, b):
        return min(x*a//b, x*b//a, a*b//x)

    domains = {}
    for i, a in enumerate(values):
        for b in values[i:]:
            upper = min(targets[a]//b, targets[b]//a,
                        multiplicity[a]*multiplicity[b]*entry_bound(a, b))
            domains[a, b] = range(upper+1)

    def reachable(row, target, omitted=None):
        """Bit n is set if the included degree classes can sum to n."""
        bits = 1
        mask = (1 << (target+1))-1
        for b, options in zip(values, row):
            if b == omitted:
                continue
            if isinstance(options, range):
                # An interval of multiplicities is a bounded number of copies.
                # Binary blocks build exactly the same shifts in logarithmic time.
                bits = (bits << (b*options.start)) & mask
                remaining, block = len(options)-1, 1
                while remaining:
                    take = min(block, remaining)
                    bits = (bits | bits << (b*take)) & mask
                    remaining -= take
                    block *= 2
            else:
                updated = 0
                for k in options:
                    updated |= bits << (k*b)
                bits = updated & mask
        return bits

    for _ in range(5):
        changed = False
        for a in values:
            target = targets[a]
            if target > 3000:
                continue
            cells = [(min(a, b), max(a, b)) for b in values]
            row = [domains[cell] for cell in cells]
            if not (reachable(row, target) >> target):
                return True
            for b, cell in zip(values, cells):
                if len(domains[cell]) < 2:
                    continue
                others = reachable(row, target, omitted=b)
                allowed = [k for k in domains[cell] if (others >> (target-b*k)) & 1]
                if len(allowed) < len(domains[cell]):
                    if allowed[-1]-allowed[0]+1 == len(allowed):
                        domains[cell] = range(allowed[0], allowed[-1]+1)
                    else:
                        domains[cell] = allowed
                    changed = True
            if multiplicity[a] > 1:
                caps = [min(domains[cell][-1], multiplicity[b]*entry_bound(a, b))
                        for b, cell in zip(values, cells)]
                one_label = reachable([range(k+1) for k in caps], x*a)
                if not (one_label >> (x*a)) or (a == x and not ((one_label >> (x*a-1)) & 1)):
                    return True
        if not changed:
            break
    return False
