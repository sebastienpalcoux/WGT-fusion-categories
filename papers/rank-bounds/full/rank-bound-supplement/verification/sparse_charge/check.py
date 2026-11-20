"""Replay the fixed-rank parity and leading-Gauss-coefficient mass bounds.

The input is the saved list of exact-rank 209 reservations surviving the
earlier component envelope. Arithmetic is integral or rational throughout.
Run this module from anywhere in a complete supplement tree.
"""
from functools import lru_cache
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
import argparse
import json


def rounded(n, step):
    return (n + step - 1) // step * step


def twopart(n):
    return n & -n


@lru_cache(None)
def power_residues(q, modulus, height, parity):
    """A specified power, or the entire specified-parity coset of <q^2>."""
    if height:
        return (pow(q, height, modulus),)
    start = pow(q, parity, modulus)
    step = pow(q, 2, modulus)
    orbit = [start]
    value = start * step % modulus
    while value != start:
        orbit.append(value)
        value = value * step % modulus
    return tuple(orbit)


@lru_cache(None)
def cofactor_residues(support, p, modulus, fixed_items, parities):
    fixed = dict(fixed_items)
    residues = {1}
    for q, parity in zip(support, parities):
        if q != p:
            allowed = power_residues(q, modulus, fixed.get(q, 0), parity)
            residues = {x * y % modulus for x in residues for y in allowed}
    return frozenset(residues)


def branch_mass(support, marked, absent, joint, nines, rank, parities,
                height_two_at_three, sparse, top_counts):
    """Disjoint populations, followed by a denominator-15 free tail."""
    heights = {}
    for p, bit in zip(support, parities):
        minimum = 1 if p in marked or joint and p == 5 else 2
        if p == 3 and not height_two_at_three:
            minimum = 3
        if p in sparse or len(support) == 3 and p != 3 and bit:
            minimum = max(minimum, 3)
        heights[p] = minimum + (bit - minimum) % 2

    # The three-free reservation and its possible overlap with the edge
    # reservation are the proved input bounds in the component census.
    free_three = (6 if height_two_at_three else 8) if nines >= 42 else 0
    overlap = min(8, free_three) if joint else 0
    joint_count = 80 if joint else 0
    booked = nines + joint_count + free_three - overlap
    population = rounded(booked, 8)
    mass = Fraction(nines, 81) + Fraction(joint_count-overlap, 225)
    mass += Fraction(free_three, max(15, 3**heights[3])**2)
    mass += Fraction(population-booked, 225)

    for p, (_, total, _) in marked.items():
        field_size = rounded(total, max(8, 2*twopart(p-1)))
        population += field_size
        mass += Fraction(total, p*p) + Fraction(field_size-total, (3*p)**2)

    for p in absent:
        d = twopart(p-1)
        if p in sparse:
            top = top_counts[p]
            field_size = rounded(max(top, d*(p+1)), max(8, 2*d))
            high_denominator = p ** ((heights[p] + 2) // 2)
            mass += Fraction(field_size-top, p**4) + Fraction(top, high_denominator**2)
        else:
            # In this branch only the mandatory active labels are known
            # to have pure denominator >=p^2.
            active = 2*p*(p-1)
            field_size = rounded(active, max(8, 2*d))
            mass += Fraction(active, p**4) + Fraction(field_size-active, max(15, p)**2)
        population += field_size

    remaining = rank - population - 1
    if remaining < 0:
        return None
    vacuum_denominator = prod(p**heights[p] for p in support)
    return mass + Fraction(remaining, 225) + Fraction(1, vacuum_denominator**2)


def feasible(row, rank):
    primes, classes, edges, nines = row
    support = tuple(primes)
    if not support or support[0] != 3 or not nines or edges not in ([], [[3, 5]]):
        return True, None
    marked = {p: c for p, c in zip(support[1:], classes) if c[0]}
    if any(c[0] != 1 for c in marked.values()) or edges and 5 in marked:
        return True, None
    if len(support) == 3 and marked:
        return False, None  # Prime-free trace theorem for an actual prime-square row.
    absent = [p for p in support if p != 3 and p not in marked and not (edges and p == 5)]
    unknown = [p for p in support if p not in marked]
    best = None
    for bits in product(range(2), repeat=len(unknown)):
        h = {p: 1 for p in marked}
        h.update(zip(unknown, bits))
        parities = tuple(h[p] for p in support)
        odd_part = prod(p for p in support if h[p])
        if odd_part**2 % 16 != rank % 16:
            continue
        high_three = h[3] or len(support) == 3 and nines >= 9
        for small_three in ((False,) if high_three else (True, False)):
            fixed = {p: 1 for p in marked}
            if small_three:
                fixed[3] = 2
            fixed_items = tuple(sorted(fixed.items()))

            for sign in (1, -1):
                if any((-sign*pow(c[2]//(p-1), -1, p)) % p not in
                       cofactor_residues(support, p, p, fixed_items, parities)
                       for p, c in marked.items()):
                    continue
                top_counts = {}
                for p in absent:
                    if not h[p]:
                        continue  # The sparse branch would give odd height.
                    modulus = p*p
                    residues = cofactor_residues(support, p, modulus, fixed_items, parities)
                    for c in range(2, 2*p, 2):
                        if -sign*pow(c, -1, modulus) % modulus in residues:
                            top_counts[p] = (p-1)*c
                            break
                # No admissible c allows only the dense branch at that prime.
                for modes in product((False, True), repeat=len(top_counts)):
                    sparse = {p for p, mode in zip(top_counts, modes) if mode}
                    upper = branch_mass(support, marked, absent, bool(edges), nines,
                                        rank, parities, small_three, sparse, top_counts)
                    if upper is not None and (best is None or upper > best):
                        best = upper
                    if best is not None and best >= 1:
                        return True, best
    return False, best


def audit(supplement):
    prior_path = supplement/'verification/odd_uniform/current-inertia-and-component-receipt.json'
    screen = next(x for x in json.loads(prior_path.read_text())['screens'] if x['rank'] == 209)
    input_file = supplement/screen['input_path']
    input_bytes = input_file.read_bytes()
    if sha256(input_bytes).hexdigest() != screen['input_sha256']:
        raise ValueError('Input checksum mismatch')
    rows = json.loads(input_bytes)['rows']
    retained = screen['retained_indices']
    if len(set(retained)) != len(retained) or len(rows) != screen['input_reservations']:
        raise ValueError('Invalid reservation indices or dimensions')
    rejections = []
    for index in retained:
        keep, upper = feasible(rows[index], 209)
        if not keep:
            rejections.append({'index': index, 'upper': None if upper is None else str(upper)})
    return {'rank': 209, 'input_sha256': screen['input_sha256'],
            'previous_receipt_sha256': sha256(prior_path.read_bytes()).hexdigest(),
            'tested': len(retained), 'retained': len(retained)-len(rejections),
            'rejected': rejections}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplement', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = audit(args.supplement)
    path = Path(__file__).with_name('receipt.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2)+'\n')
    elif result != json.loads(path.read_text()):
        raise ValueError('Independent result and saved receipt differ')
    print(f"{result['tested']} tested; {len(result['rejected'])} excluded; {result['retained']} retained.")
