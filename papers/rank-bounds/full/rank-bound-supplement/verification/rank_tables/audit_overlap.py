"""Independently check each exclusion witness, without regenerating degree multisets."""
from collections import Counter
from itertools import combinations_with_replacement
from math import gcd, lcm, prod
from pathlib import Path
from zipfile import ZipFile
import argparse
import gzip
import json

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'census'))
from principal_overlap import phase_failure


def valuation(x, p):
    e = 0
    while x % p == 0:
        e += 1
        x //= p
    return e


def verify(row, D, witness, primes):
    p = witness["prime"]
    if p not in primes:
        raise ValueError("invalid prime")
    e = valuation(D, p)
    multiplicity = Counter(row)
    unit_degrees = {a for a in row if a > 1 and a % p}
    allowed = {}
    principal = set()
    for c in set(row) - {1}:
        if D // (c * c) % (p * p):
            continue
        energy = 1
        for a in unit_degrees:
            common = gcd(c, a)
            if common == 1 or lcm(c, a)**2 >= D:
                break
            n = multiplicity[a]
            residue = common % p
            q = max(n, -(-(n * residue * (p - residue)) // (p - 1)))
            while (q - n * common**2) % p:
                q += 1
            energy += q * (a // common)**2
        else:
            if energy <= D // (c * c):
                principal.add(c)
    for a in set(row):
        allowed[a] = {
            c for c in principal
            if a == 1 or (gcd(a, c) > 1 and lcm(a, c)**2 < D)
        }

    if witness["test"] == "principal_product":
        more = tuple(witness["additional_degrees"])
        if (set(witness["p_unit_degrees"]) != unit_degrees
                or any(a not in multiplicity for a in more)
                or sum(valuation(a, p) for a in more) >= e):
            return False
        possible = principal.intersection(*(allowed[a] for a in (*unit_degrees, *more)))
        return (possible == set(witness["possible_columns"])
                and all(2*valuation(c, p) > sum(min(valuation(c, p), valuation(a, p))
                                                for a in more) for c in possible))

    if witness["test"] != "primary_phase" or set(witness["possible_columns"]) != principal:
        return False
    # Explicit set merging is independent of the production union-find.
    parts = [{a} for a in multiplicity]
    constant = {a for a in multiplicity if multiplicity[a] == 1}
    for a, b in combinations_with_replacement(multiplicity, 2):
        if a == b and multiplicity[a] == 1:
            continue
        products = [c*c // (gcd(c, a) * gcd(c, b)) for c in allowed[a] & allowed[b]]
        if gcd(*products) == 1:
            continue
        joined = set().union(*(part for part in parts if a in part or b in part))
        parts = [part for part in parts if part.isdisjoint(joined)] + [joined]
        constant |= joined
    labels = phases = mass = 0
    h = max(1, (p-1)//2)
    for group in parts:
        if 1 in group or any((D//(a*a)) % p for a in group):
            continue
        if group <= constant:
            if p >= 5:
                continue
            phases += 1
            labels += sum(multiplicity[a] for a in group)
            mass += sum(multiplicity[a]*a*a for a in group)
        else:
            if len(group) != 1:
                raise ValueError("invalid phase partition")
            a = next(iter(group))
            slots = multiplicity[a] - multiplicity[a] % h
            labels += slots
            phases += slots
            mass += slots*a*a
    if [labels, phases, mass] != [witness["max_active"], witness["max_phases"], witness["max_mass"]]:
        return False
    phase_min = label_min = h
    if p > 2 and e % 2 == 0:
        phase_min = p-1
        kernel = [t for t in primes if t > 2 and valuation(D,t) % 2]
        base = prod((t+1)//2 for t in kernel)
        choices = [base * (2 if valuation(D,2) % 2 else 1)]
        choices += [base // ((t+1)//2) * ((t+3)//2) for t in kernel]
        label_min = max(choices)*phase_min
    return ((phase_min, label_min) == (witness["min_phases"], witness["min_active"])
            and (phases < phase_min or labels < label_min
                 or 2*p*min(set(row)-{1})*mass < (p-1)*D))


def main(archive):
    tally = Counter()
    # All primes up to 2*44+1. Verify they exhaust each input dimension.
    candidates = [p for p in range(2,90) if all(p%q for q in range(2,p))]
    with ZipFile(archive) as z:
        for r in range(1,45):
            for line in gzip.decompress(z.read(f'ranks/{r:02}/raw-census.jsonl.gz')).splitlines():
                row = json.loads(line)[-1]
                D = sum(a*a for a in row)
                factors = [(p, valuation(D, p)) for p in candidates if D%p == 0]
                if prod(p**e for p,e in factors) != D:
                    raise ValueError("prime support mismatch")
                witness = phase_failure(row, D, factors)
                if witness is None:
                    tally["admitted"] += 1
                    continue
                if not verify(row, D, witness, [p for p,e in factors]):
                    raise ValueError(f"invalid witness at rank {r}: {row} {witness}")
                tally[witness["test"]] += 1
    print(json.dumps(dict(tally), sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('archive', type=Path)
    main(parser.parse_args().archive)
