"""Enumerate short sums of reciprocal prime powers with signed error 1/M."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt, lcm


def denominator_options(p, least, ceiling):
    """Effective reciprocal denominator and its underlying degree denominator."""
    d = 32 if p == 2 else p**3
    while d <= ceiling:
        if d >= least:
            b = isqrt(d)
            if b*b != d:
                b = isqrt(p*d)
            yield d, b
        d *= p if p == 2 else p*p


def sparse_budget(residual, multiple, odd_primes, budget, max_eights=None, protected=()):
    """Test all term counts within the label budget in one sorted recursion."""
    primes = (2,)+tuple(odd_primes)
    def prime_free_part(n, p):
        while n % p == 0:
            n //= p
        return n
    def visit(value, required, least, left, missing, eights):
        if not value:
            # With a protected prime, the terminal 1/M either has this prime
            # absent from the remaining denominators, or is smaller than its
            # entire positive contribution. Otherwise the caller uses a prime
            # in the target denominator absent from every remaining term.
            assert protected
            return
        if not missing and abs(value.numerator)==1 and value.denominator % required==0:
            yield value
        max_terms = left-sum(p-2 for p in missing)
        if max_terms <= 0 or abs(value)>F(max_terms,least)+F(1,required):
            return
        if value<0 and -value >= F(1,required):
            return
        # An underlying pure p-power denominator b divides M, as does
        # the prime-to-p part of the lcm accumulated so far.
        for p in primes:
            cost=max(1,p-1)
            rest_missing=missing-{p}
            if left-cost < sum(q-1 for q in rest_missing):
                continue
            coprime = prime_free_part(required,p)
            error = abs(value)
            ceiling = max(2*max_terms/error, 4/(coprime*error)**2)
            if error > F(1,required):
                ceiling = min(ceiling,max_terms/(error-F(1,required)))
            for d,b in denominator_options(p,least,ceiling):
                k = eights + (p == 2 and b == 8)
                if max_eights is not None and k > max_eights:
                    continue
                # Pairs of equal dyadic digits carry to the preceding power.
                # Carrying only weakens the label and divisor requirements.
                # At 1/32 repetition preserves the lower denomination bound.
                next_least=d+1 if p==2 and d>32 and max_eights is None else d
                yield from visit(value-F(1,d),lcm(required,d if p in protected else b),next_least,left-cost,rest_missing,k)
    yield from visit(residual,multiple,1,budget,set(odd_primes),0)

@lru_cache(None)
def pure_check(rank,rows,support,modes,neutral):
    prime=[p for p,z in zip(support,modes) if z==2]
    odd=tuple(p for p in support if p not in prime)
    M0=8
    for p,z in zip(support,modes):
        M0*= p if z==2 else p**(3 if p>=5 else 2)
    from itertools import product
    out=set()
    for js in product(*((1,) if p==5 else range(1,p) for p in prime)):
        ms=tuple(p if j==1 else (p-1)*j for p,j in zip(prime,js))
        budget=rank-1-rows-sum(ms)-neutral
        if budget<sum(p-1 for p in odd):continue
        rem=1-F(rows,16)-sum(F(n,p*p) for n,p in zip(ms,prime))
        protected={p for p in odd if budget-sum(t-1 for t in odd if t!=p)<p*(p-1)}
        for a in range(-rows,rows+1):
            val=1-F(rows-a,16)-sum(F(j,p) for j,p in zip(js,prime))
            if val+F(1,M0)<=0 or val>2*rem+F(1,M0):continue
            for e in sparse_budget(val,M0,odd,budget,protected=protected):
                M=e.denominator
                if M%M0:continue
                rest=M
                for p in (2,)+tuple(support):rest=prime_free(rest,p)
                if rest==1:out.add(M)
    return tuple(sorted(out))
def prime_free(M,p):
    while M%p==0:M//=p
    return M
