"""Certified sufficient tests for BDP's three-repeated-prime divisor families.

Input to covered_three is a factored square-full profile. The caller supplies
distinct increasing primes. A string cites the first certificate; None means
no certificate from this function. All arithmetic is exact.
"""
from functools import lru_cache
from itertools import permutations
from math import gcd, isqrt


def g(h):
    return max(1, (1 + isqrt(4*h-3)) // 2)


def support(s, k, odd=True):
    """Scalar exclusion at k positions. odd=False uses mixed-parity exposure."""
    if k <= 1:
        return True
    e = k*(k-1)//2
    ell = 3 if odd else 2
    if ell*(1+2*e) >= (ell-1)*s:
        return False
    # On each interval of constant g(h), its right endpoint maximizes the
    # quotient, and maxima at complete endpoints increase strictly.
    j = g(e)
    return all((s-3)*g(h) > 2*h*(h-1) for h in (e, j*(j-1)) if h)


@lru_cache(maxsize=4096)
def factors(n):
    out = []
    t = 2
    while t*t <= n:
        a = 0
        while n % t == 0:
            n //= t
            a += 1
        if a:
            out.append((t, a))
        t = 3 if t == 2 else t + 2
    if n > 1:
        out.append((n, 1))
    return tuple(out)


@lru_cache(maxsize=4096)
def divisors(n):
    ds = [1]
    for p, a in factors(n):
        ds = [d * p**i for d in ds for i in range(a+1)]
    return tuple(sorted(ds))


def stabilizer(s, bound, real=False):
    return max(d//gcd(d, 2) if real else d for d in divisors(s-1) if d <= bound)


# Finite word certificates printed in BDP 8.9, 20.16 and 20.19, and after
# 20.11.  Thresholds here apply to odd ordinary data only.
POSITION_THRESHOLD = (0, 0, 5, 11, 29, 47, 83, 127)
SPECIAL = {577:10, 1783:14, 180547:64, 3877:18, 17981:30, 26449:34,
           29453:35, 66889:46}


def equality_possible(H, h, j):
    if j >= 2 and h == j*(j+1):
        return False
    for ell, a in factors(h-j):
        if a % 2 == 0:
            continue
        for b, _ in factors(H):
            if b == ell:
                continue
            if b == 2:
                return False
            # In a cyclic group of even size 2^a u, the order of ell
            # is even precisely when ell^u is nontrivial.
            u = b-1
            while u % 2 == 0:
                u //= 2
            if pow(ell, u, b) != 1:
                return False
    return True


@lru_cache(maxsize=8192)
def word(s, k, quotient_budget=0):
    """Odd ordinary word exclusion. Unknown words receive False.

    quotient_budget=0 uses only the order-two quotient. Positive values
    bound the total number of additional occupancy nodes examined per word.
    """
    if quotient_budget < 0:
        raise ValueError('quotient_budget must be nonnegative')
    if k <= 1 or support(s,k) or (k <= 7 and s >= POSITION_THRESHOLD[k]):
        return True
    if k <= SPECIAL.get(s, 0):
        return True
    e = k*(k-1)//2
    if 3*(1+2*e) >= 2*s:
        return False
    H = (s-1)//2
    # Check 8.38 and the order-two quotient of 20.16. No subset search.
    budget = [quotient_budget]
    jmax = g(e)
    for j in range(jmax, 0, -1):
        lo, hi = j*(j-1)+1, min(e, j*(j+1))
        if hi*(hi-1) < (H-1)*j:
            # A terminal partial interval can be weaker than its predecessor.
            if j != jmax:
                break
            continue
        for h in range(hi, lo-1, -1):
            delta = h*(h-1) - (H-1)*j
            if delta < 0:
                break
            if delta == 0 and not equality_possible(H, h, j):
                continue
            if H % 2 == 0 and delta < 2*(h-1):
                high = isqrt(h-j+delta)
                high -= (high-h) % 2
                if high < 0 or high*high < h-j-delta:
                    continue
            if budget[0] and delta < 2*(h-1):
                if any(not quotient_possible(H,b,h,j,budget) for b in divisors(H)
                       if 2 < b and budget[0]):
                    continue
            return False
    return True


def quotient_possible(H,b,h,j,budget):
    """Necessary test from BDP 20.16. Budget exhaustion means possible."""
    cap = H//b
    # Only nonempty residue classes need be stored. A translation sends
    # an occupied class to zero; the equivalent shifted population is tested.
    pop = []
    def search(left, start):
        if not budget[0]:
            return True
        budget[0] -= 1
        if left == 0:
            corr = {}
            for i,x in pop:
                for k,y in pop:
                    d = (k-i)%b
                    corr[d] = corr.get(d,0) + x*y
            return (len(corr) == b and min(corr.values()) >= cap*j
                    and corr[0] >= cap*j + h-j)
        if start >= b or left > cap*(b-start):
            return False
        for i in range(start, b-(left+cap-1)//cap+1):
            if not pop and i != 0:
                break
            for x in range(min(left,cap),0,-1):
                pop.append((i,x))
                allowed = search(left-x,i+1)
                pop.pop()
                if allowed:
                    return True
        return False
    # Prevent recursive depth overflow on large optional jobs. No exclusion
    # is drawn from a skipped population problem.
    if min(h,b)>500:
        return True
    return search(h,0)


def is_prime(n):
    return n >= 2 and (n == 2 or (n % 2 and all(n % t for t in range(3,isqrt(n)+1,2))))


def all_words_above(n, k, quotient_budget=0):
    # Every integer s>=upper meets the scalar inequalities. If W is unknown
    # at any prime in the residual finite interval, return no certificate.
    if support(n+1,k) or (k <= 7 and n+1 >= POSITION_THRESHOLD[k]):
        return True
    e = k*(k-1)//2
    j = g(e)
    H = max(h*(h-1)//g(h) for h in (e,j*(j-1)) if h)
    upper = max(2*H+6, 3*(1+2*e)//2+2)
    if k <= 7:
        upper = min(upper, POSITION_THRESHOLD[k])
    for s in range(n+1, upper):
        if is_prime(s) and not word(s,k,quotient_budget):
            return False
    return True


def separated(t, a, u, b, s, odd):
    M = t**(a//2) * u**(b//2)
    h = max(t,u)
    if s*h*t*u > M*M:
        return True
    lhs = (M//h)**2-1
    for scale in ((1,2) if odd else (1,)):
        v = (-pow(scale*s*t, -1, u)) % u
        w = (-pow(scale*s*u, -1, t)) % t
        if lhs < scale*s*(t*v+u*w):
            return True
    return False


def covered_three(primes, exponents, quotient_budget=0):
    """First BDP divisor-family certificate, or None; inputs are factored."""
    p,q,r = primes
    a,b,c = exponents
    if not (2 <= p < q < r and min(exponents) >= 2):
        raise ValueError("supply three increasing primes and exponents >= 2")
    return _covered(p,q,r,a,b,c,quotient_budget)


@lru_cache(maxsize=8192)
def _covered(p,q,r,a,b,c,quotient_budget):
    def W(s,k):
        return word(s,k,quotient_budget)

    if a <= 5 or b <= 4 or c <= 3:
        return '4.12/7.3/8.3'
    if all(x <= cap for x, cap in zip(sorted((a,b,c)), (4,5,7))):
        return '8.1'
    if p == 2 and a <= 6 and (q,r) != (3,5):
        return '21.2'
    h = a//2
    if q > p**h:
        return '7.11 (gap)'
    odd = p != 2
    if q > p**(h-1):
        if b <= 5 or (odd and (p-1)*(r-1) >= stabilizer(r,p)*p**(a-2)):
            return '7.11 (companion)'
        if b%2 or c%2:
            cert = _covered(p,q,r,a,b-b%2,c-c%2,quotient_budget)
            if cert:
                return '7.11 (parity) + ' + cert
    for (t,x),(u,y),(s,_) in permutations(zip((p,q,r),(a,b,c))):
        if t < u and separated(t,x,u,y,s,odd):
            return '4.14'
        if s > 2 and x >= 4 and support(s,t**((x-2)//2),odd):
            return '8.10'
    k = 1+r//q
    if c <= 4:
        if odd and k < p and W(r,k):
            return '8.2/8.7/8.9 (word)'
        if q >= 5 and (odd or (r-1)*(r-2)%q) and b%2 and support(r,k,odd):
            cert = _covered(p,q,r,a,b-1,c,quotient_budget)
            if cert:
                return '8.11 + ' + cert
        for t,x,u in ((p,a,q),(q,b,p)) if odd else ((p,a,q),):
            G = t**(x//2)
            if r > G+1 and (W(r,1+G//u) if odd else support(r,1+G//u,False)):
                return 'after 8.11 (finite largest prime)'
    if not odd:
        return None
    if b <= 5:
        j = q//p
        if 3*(1+j*(j+1)) < 2*q and (c <= 4 or (r-q)*(p-1)>q):
            return '8.19'
        if a >= 6 and c <= 4 and q > p**(h-2):
            return '8.20'
    if c <= 4:
        M = p**h//q
        if M+1 < p and q >= (M+1)**2:
            return '8.12'
    # For results stated with even caps take CEIL, including odd input a,b.
    H, K = (a+1)//2, (b+1)//2
    U, V = p**(H-1), q**(K-1)
    if c <= 4:
        m = p**H//q
        if (p-1)*q > U and (q >= (m+1)**2 or H >= 6 or
                              all(all_words_above(j*q,j+1,quotient_budget) for j in range(1,m+1))):
            return '20.1/20.2'
        if r > max(U,V) and W(r,1+max(U//q,V//p)):
            return '20.3'
    if H >= 3 and r>U and W(r,1+U//q):
        if all(W(r,p**t) or (p-1)*(r-1)>=p**(2*(H-t-1))*stabilizer(r,p**t)
               for t in range(1,H-1)):
            return '20.6'
    if a <= 6:
        square = p*p
        if (c <= 4 or r > square) and W(r,1+min(r,square)//q):
            return '20.10'
        if p >= 5 and 2*q > square:
            return '20.14'
        if p >= 5 and r < square:
            kq, kr = 1+square//r, 1+square//q
            middle = 2*r > square+2 and support(r,kr,False)
            general = (q >= 5*(kq+1) and support(q,kq,False)
                       and support(r,kr,False))
            if middle or general:
                cost = sum((s-1)//stabilizer(s,p,True) for s in (q,r))
                if (p-1)*cost >= square:
                    return '20.13 (rectangular support)'
    if b <= 6 and r > q*q and W(r,1+q*q//p):
        return '20.15'
    return None
