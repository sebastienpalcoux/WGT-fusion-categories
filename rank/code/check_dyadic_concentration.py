"""Audit overlapping concentration tests on a completed frontier return."""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from zipfile import ZipFile


@lru_cache(maxsize=8192)
def prime_divisors(number):
    result = []
    p = 2
    while p * p <= number:
        if number % p == 0:
            result.append(p)
            while number % p == 0:
                number //= p
        p = 3 if p == 2 else p + 2
    if number > 1:
        result.append(number)
    return result


def profile_tests(degrees):
    """Return exact mass and norm tests on the largest codegree-sixteen class."""
    counts = Counter(degrees)
    x = max(counts)
    multiplicity = counts[x]
    dimension = sum(m * d * d for d, m in counts.items())
    if dimension != 16*x*x or multiplicity <= 8:
        return set()
    hits = set()
    smallest = min(d for d in counts if d > 1)
    odd = {d: m for d, m in counts.items() if d % 2}
    for y in odd:
        cost = sum(m * (d // gcd(d, y))**2 for d, m in odd.items())
        if x % 2 == 0:
            cost += (2*multiplicity - 16) * (x // gcd(x, y))**2
        if cost * y*y > dimension:
            hits.add('odd_row_norm')
        if y > 1 and any(d > 1 and gcd(d, y) == 1 for d in odd):
            hits.add('coprime_odd_pair')
    for p in prime_divisors(2*x):
        packet = 1 if p == 2 else (p - 1) // 2
        capacity = 0
        even_capacity = 0
        for d, m in counts.items():
            if d == 1 or dimension // (d*d) % p:
                continue
            weight = (m - m % packet)*d*d
            capacity += weight
            if d % 2 == 0:
                even_capacity += weight
        threshold = (p - 1)*dimension
        if 2*p*smallest*capacity < threshold:
            hits.add('original_activity')
        if p != 2 and 2*p*smallest*even_capacity < threshold:
            hits.add('odd_prime_even_capacity')
    if x % 2 and multiplicity > 8:
        hits.add('odd_codegree16_cap')
    return hits


def audit(archive):
    result = Counter()
    new_tests = {'odd_row_norm', 'coprime_odd_pair', 'odd_prime_even_capacity', 'odd_codegree16_cap'}
    with ZipFile(archive) as source:
        for line in source.open('receipts.jsonl'):
            for row in json.loads(line)['rows']:
                result['raw_profiles'] += 1
                hits = profile_tests(row)
                if hits & new_tests:
                    result['union_new_tests'] += 1
                    if 'original_activity' not in hits:
                        result['new_union_past_original_activity'] += 1
                if 'odd_prime_even_capacity' in hits and 'original_activity' not in hits:
                    result['even_capacity_past_original_activity'] += 1
                result.update(hits)
    return {'archive_sha256': sha256(archive.read_bytes()).hexdigest(), 'counts': dict(sorted(result.items()))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--expect', type=Path, help='compare to stored audit')
    args = parser.parse_args()
    report = audit(args.archive)
    if args.expect and json.loads(args.expect.read_text()) != report:
        raise SystemExit('audit differs from stored result')
    print(json.dumps(report, indent=2))
