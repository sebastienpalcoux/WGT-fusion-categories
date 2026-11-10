from functools import lru_cache
from itertools import combinations
from math import comb, isqrt

def least_overlap(h):
    j = (isqrt(1 + 4*h) - 1) // 2
    return j + (j*(j+1) < h)

def words(s, k, ell=3):
    """Exclude all full-weight s-primary words on <= k positions at a nonunit row.
    Returns (status, count of normalized upper bounds, first surviving upper bound).
    s must be an odd prime. Use ell=2 in either parity, ell=3 in odd ordinary data.
    """
    if s < 3 or s % 2 == 0 or any(s%d == 0 for d in range(3,isqrt(s)+1,2)):
        raise ValueError("expected an odd prime")
    if not 2 <= k <= s or ell not in (2,3):
        raise ValueError("invalid position or order bound")
    H = (s-1)//2
    E = comb(k,2)
    if ell*(1+2*E) >= (ell-1)*s:
        return "exposure undecided", 0, None
    h0 = next((h for h in range(1,E+1)
               if (H-1)*least_overlap(h) <= h*(h-1)), None)
    if h0 is None:
        return "excluded by scalar bound", 0, None
    need = least_overlap(h0)
    b = comb(k-1,2)
    if max(E*E, s*b*b+(E-b)**2) >= s*H*need:
        return "higher conductor undecided", 0, None

    for gen in range(2,s):
        units = [1]
        for _ in range(1,H):
            units.append(units[-1]*gen % s)
        if len({min(x,s-x) for x in units}) == H:
            break
    bits = [0]*s
    for j,x in enumerate(units):
        bits[x] = bits[-x] = 1 << j
    full = (1 << H)-1

    def rotate(a,j):
        return (a << j) & full | a >> (H-j)

    def canonical(a):
        return min(rotate(a,j) for j in range(H))

    @lru_cache(maxsize=131072)
    def difference(a):
        # a represents the symmetric set {0} union all encoded +/- pairs.
        xs = [x for i,x in enumerate(units) if a & (1<<i)]
        out = a
        for x in xs:
            out |= bits[2*x % s]
        for x,y in combinations(xs,2):
            out |= bits[(x-y) % s] | bits[(x+y) % s]
        return out

    @lru_cache(maxsize=None)
    def contradiction(a):
        known, seen = [a], {a}
        for i,e in enumerate(known):
            if k == 6 and i == 60:
                break
            for f in known[:i+1]:
                for j in range(H):
                    inter = e & rotate(f,j)
                    if inter.bit_count() < need:
                        return True
                    nxt = difference(inter)
                    h = nxt.bit_count()
                    if h < h0:
                        return True
                    if (1+2*h)*ell >= (ell-1)*s or (k == 6 and h >= 19):
                        continue
                    nxt = canonical(nxt)
                    if nxt not in seen:
                        known.append(nxt)
                        seen.add(nxt)
        return False

    total = comb(s-2,k-2)
    first = None
    survivors = 0
    for tail in combinations(range(2,s), k-2):
        A = (0,1,*tail)
        d = 0
        for x,y in combinations(A,2):
            d |= bits[(x-y) % s]
        if d.bit_count() < h0:
            continue
        if any((d & rotate(d,j)).bit_count() < need for j in range(1,H//2+1)):
            continue
        if not contradiction(canonical(d)):
            survivors += 1
            if first is None:
                first = A
    return ("excluded" if survivors == 0 else f"{survivors} undecided"), total, first

if __name__ == "__main__":
    for s,k in ((23,4),(31,4),(47,5),(53,5),(59,5),(61,5)):
        print(s,k,words(s,k,3))
    # The six-position checks are longer. Each has zero survivors.
    # for s in (83,89,97,101,103,107):
    #     print(s,6,words(s,6,2))
