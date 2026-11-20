"""Gauss-sum check when no twist is active at two distinct primes."""
from fractions import Fraction
from functools import lru_cache
from math import ceil, gcd, lcm, prod

import argparse
import hashlib
import json
import sys
from pathlib import Path


def twopart(n):
    return n & -n


def rounded(n, m):
    return -(-n // m) * m


class Incomplete(Exception):
    """Work limit reached with a possible completion outstanding."""



def separate(row, rank, budget, witnesses):
    """Ordered active-denominator recursion on an edgeless reservation."""
    support, classes, edges, nines = row
    if (edges or not nines or not support or support[0] != 3
            or len(classes) != len(support) - 1 or len(support)<3
            or any(v not in (0,1) for v, m, a in classes)):
        return "unsupported", 0
    marked = {p:(v,m,a) for p,(v,m,a) in zip(support[1:],classes) if v}
    if any(p < 11 for p in marked):
        return "unsupported", 0
    multiplicities = {p:v[1] for p,v in marked.items()}
    active = {p:v[2] for p,v in marked.items()}
    other = tuple(p for p in support if p not in marked)
    per_field = {3:nines, **multiplicities}
    for p in other:
        if p != 3:
            per_field[p] = twopart(p-1)*(p+1)
    for p in per_field:
        per_field[p] = rounded(per_field[p], max(8, 2*twopart(p-1)))
    marked_population = sum(per_field[p] for p in marked)
    caps = {p:rank-1-sum(per_field.values())+per_field[p] for p in support}
    sparse = {p for p in other if caps[p] < 2*p*(p-1)}
    free_three = 8 if ((len(support)==3 and nines>=42)
                       or (len(support)==4 and nines>=72)) else 0
    base = prod(p if p in marked else p**3 if p in sparse or (p==3 and free_three)
                else p*p for p in support)
    proper = 10 if len(support)==4 and len(marked)>=2 else 6 if marked else 0
    neutral = []
    if len(support)==3:
        neutral = [(witnesses(p,nines),(9 if nines>=9 else 3)*p*p)
                   for p in support if p!=3]
    if not proper:
        neutral.append((max(0,6-sum(n for n,b in neutral)),15))
    room = rank-1-nines-sum(multiplicities.values())-proper-sum(n for n,b in neutral)
    charge = 1-sum((Fraction(a,p*(p-1)) for p,a in active.items()),Fraction())
    quota = tuple((4 if p%4==1 else 2)*(p-1) for p in other)
    spare = rank-1-nines-sum(multiplicities.values())
    residual_mass = 1-Fraction(nines,81)-sum(Fraction(m,p*p) for p,m in multiplicities.items())

    def allowed(b):
        rem = b
        for p in support:
            while rem%p==0:
                rem //= p
        if rem!=1 or b in support or any(b%(p*p)==0 for p in marked):
            return False
        for p in support[1:]:
            k = b
            while k%p==0:
                k //= p
            u = b//k
            if k==1 or u==1:
                continue
            a = pow(k,-1,u)
            pair = 2*Fraction(min(a,u-a)**2,u*u)
            if p in marked:
                j = active[p]//(p-1)
                norm_budget = 2*p*j-j*j-multiplicities[p]
                if pair*p*p+4 > norm_budget:
                    return False
            else:
                x = min(Fraction(spare,p**4),residual_mass,Fraction(p*p-1,2*p*p))
                cap = Fraction(p+1,p-1)*x-Fraction(p*p,(p-1)**2)*x*x
                if pair > cap:
                    return False
        return True

    free_floor = next((b for b in range(15,101,2) if allowed(b)),101)
    nodes = 0
    # Retain an upper bound on vacuum mass in integer units. Rounding
    # every positive term upward can only permit additional completions.
    scale = 1 << 75

    @lru_cache(None)
    def upper_weight(b):
        return (scale + b * b - 1) // (b * b)

    free_weight = upper_weight(free_floor)
    initial_mass = (upper_weight(base) + nines * upper_weight(9)
                    + sum(m * upper_weight(p) for p, m in multiplicities.items())
                    + proper * upper_weight(33)
                    + sum(n*upper_weight(b) for n,b in neutral) + room * free_weight)
    @lru_cache(None)
    def next_levels(last):
        result = []
        for p in other:
            b = p * p
            while b <= last:
                b *= p
            result.append(b)
        return tuple(result)

    def visit(last, slots, num, den, mass, multiple, credits, counts):
        nonlocal nodes
        nodes += 1
        if nodes > budget:
            raise Incomplete()

        # num/den is the exact signed residual; reduction is unnecessary.
        if slots < 0 or mass < scale or num == 0:
            return False
        missing = [max(0, q - c) for q, c in zip(quota, credits)]
        total_missing = sum(missing)
        # Each closed p-class and the complete class of nines are already paid.
        # Pure active counts below are incremental beyond the recorded classes.
        population=[]
        for i,p in enumerate(other):
            d=(p-1)&-(p-1)
            if p==3:
                low=nines+max(counts[i]+missing[i],free_three)
            else:
                low=max(d*(p+1), counts[i]+missing[i])
                if caps[p]<p*(p-1):
                    low=max(low,counts[i]+missing[i]+2*d)
            m=max(8,2*d)
            population.append((low+m-1)//m*m)
        if sum(population)+marked_population>rank-1: return False
        if total_missing > slots:
            return False
        levels = next_levels(last)
        multiple=lcm(multiple,*(b*b//p if p in sparse else b
                         for p,b,n in zip(other,levels,missing) if n))
        if num*multiple < -den: return False
        # Bound the positive Gauss tail per required active sector.
        # Upper weights rounded upwards; exact comparison in integers.
        if num>0:
            weights=[(p*scale+(p-1)*b*b-1)//((p-1)*b*b) for p,b in zip(other,levels)]
            ceiling=sum(n*v for n,v in zip(missing,weights))+(slots-total_missing)*max(weights)
            if (num*multiple-den)*scale>den*multiple*ceiling: return False

        loss = sum(n * (free_weight - upper_weight(b)) for n, b in zip(missing, levels))
        if mass - loss < scale:
            return False
        if total_missing == 0 and den % abs(num) == 0:
            M = den // abs(num)
            if M % multiple == 0:
                return True

        for i, p in enumerate(other):
            limit = slots // (p - 1)
            if limit == 0:
                continue
            b = levels[i]
            while True:
                required = lcm(multiple, b)
                left = num * required
                square = b * b
                # -1/required <= num/den <= 1/required + 3*slots/(2*b^2).
                if left < -den or 2 * (left - den) * square > 3 * slots * den * required:
                    break

                # The current coefficient c decreases the residual by c/h.
                # At the following level, apply the same two inequalities
                # with the old slot count, which is an upper bound.
                h = square // p
                common = lcm(den, h)
                numerator = num * (common // den)
                unit = common // h
                next_square = min(next_levels(b))**2
                lower_num = (2 * (left - den) * next_square - 3 * slots * den * required) * h
                lower_den = 2 * den * required * next_square
                first = max(1, -(-lower_num // lower_den))
                last_c = min(limit, (left + den) * h // (den * required))

                # A largest primitive-square-class imbalance would need at
                # least p^2 packets. Higher phase orders contribute even c.
                step = 2 if limit < p*p else 1
                first += (-first) % step
                for c in range(first, last_c + 1, step):
                    used = (p - 1) * c
                    next_multiple = required
                    if limit < p*p:
                        next_multiple = lcm(required, h // gcd(h, c))
                    new_credits = list(credits)
                    new_credits[i] += p * c
                    new_mass = mass - used * (free_weight - upper_weight(b))
                    if visit(b, slots - used, numerator - c * unit, common,
                             new_mass, next_multiple, new_credits, counts[:i]+(counts[i]+used,)+counts[i+1:]):
                        return True
                b *= p
        return False

    try:
        for c in range(nines//2+1):
            residual = charge-Fraction(c,27)
            credits = [3*c if p==3 else 0 for p in other]
            if visit(9, room, residual.numerator, residual.denominator, initial_mass,
                     base, credits, (0,)*len(other)):
                return "retained", nodes
    except Incomplete:
        return "limit", nodes
    return "excluded", nodes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supplement',type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--budget',type=int,default=1500000)
    parser.add_argument('--receipt',type=Path)
    args=parser.parse_args()
    base=args.supplement.resolve()
    sys.path.insert(0,str(base/'arithmetic'))
    from odd_square import witnesses
    old=base/'verification/odd_uniform/current-inertia-and-component-receipt.json'
    saved=json.loads(old.read_text())
    sc=next(screen for screen in saved['screens'] if screen['rank']==209)
    source=base/sc['input_path']
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    if digest!=sc['input_sha256']:
        raise ValueError('Input digest does not match the source receipt')
    rows=json.loads(source.read_text())['rows']
    outcomes={i:separate(rows[i],209,args.budget,witnesses)
              for i in sc['retained_indices'] if not rows[i][2]}
    report={'rank':209,'input_sha256':digest,'scope':'edgeless saved reservations',
            'parent_screen_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
            'budget_per_row':args.budget,'outcomes':outcomes,
            'total_integer_nodes':sum(n for _,n in outcomes.values()),
            'retained':sum(x!='excluded' for x,_ in outcomes.values())}
    if args.receipt:
        args.receipt.write_text(json.dumps(report,indent=1,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='outcomes'}))
    if not outcomes or report['retained']:
        raise SystemExit(1)

if __name__=='__main__':
    main()
