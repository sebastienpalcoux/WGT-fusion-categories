"""Necessary conditions for an ordinary simple perfect integral modular table."""
from collections import Counter
from fiber_packets import inspect
from generator_activity import mass_obstruction
from regular_indicator import obstructions
from mixed_row import test as mixed_row_test
from rank_audit import fast_verdict, validate
from class_mass import verdict as balance
from families import factors


def screen(rows, rank):
    counts = Counter()
    survivors = []
    for row in rows:
        row = list(row)
        dimension = validate(row, rank)
        counted = inspect(row, dyadic=True)
        prime = next((p for p, r in counted['primes'].items()
                      if r['available'] < r['required']), None)
        if prime:
            reason = f'complementary fiber count at {prime}'
        elif hit := mass_obstruction(row):
            reason = f'generator activity at {hit[0]}'
        elif hits := obstructions(row):
            reason = f'regular indicator at {hits[0][0]}'
        else:
            reason = fast_verdict(row, dimension)
            if reason == 'open':
                check = mixed_row_test(row)
                if check and check != 'open':
                    reason = 'mixed ternary row: '+check
            if reason == 'open':
                valuation = factors(dimension)
                odd = tuple(p for p,e in valuation if p > 2 and e % 2)
                for p,e in valuation:
                    if p > 2 and e % 2 == 0:
                        result = balance(row,p,odd_primes=odd)
                        if result != 'passes':
                            reason = f'balanced prime phases at {p}: {result}'
                            break
        counts[reason] += 1
        if reason == 'open':
            survivors.append(row)
    return {'rank':rank, 'profiles':len(rows), 'first_exclusions':dict(sorted(counts.items())),
            'survivors':survivors}
