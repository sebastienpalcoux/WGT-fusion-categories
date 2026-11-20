"""Conservative counts for (3,t)-active slots in odd square dimension."""
from itertools import product
from functools import lru_cache
from math import isqrt

@lru_cache(None)
def axis(t,k,j,h,v,primitive=None, primitive_small=None):
    """Feasibility of a nontrivial p-phase with t the complementary prime.
    k joint pairs, j at minimum B=3t; h solely-p pairs, v at B.
    """
    if not k:return True
    if not h:return False
    small=j+v
    other=k+h-small
    square=25 if t==3 and other<3 else 9
    if other < (3 if t%4==3 else 2):
        # Higher t-adic height than B cannot then occur. The remaining
        # denominators are supported on 3 and the other prime, first height t.
        q2=t*t
        u=3 if t!=3 else 5  # For p=7, 7 generates the same squares modulo 9.
        squares={pow(u,2*a,q2) for a in range(q2)}
        residues={small%q2}
        for _ in range(other):
            residues={(a+b)%q2 for a in residues for b in squares}
        if 0 not in residues:return False
    # t divides each coefficient; dual-pair count after cancelling pure powers.
    # A minimum denominator bearing nonzero coefficient costs its square in units of 1/B^2.
    n=k+h
    # The primitive mass is that at any single primitive phase. Its
    # residue modulo eight is the number of dual pairs in that phase.
    if primitive is None:
        choices = range(k//(t-1)+1)
        primitive = tuple(c for c in choices if k >= t*(t-1)//2 or
                           (k-(t-1)*c)%8 == 0)
    else:
        primitive = (primitive,)
    if all((k-c)%2 for c in primitive):
        return False
    def only_primitive(c):
        least=1+(c-1)%8 if c else 0
        return k < c+(t-2)*least+t*(t-1)//2
    primitive = tuple(c for c in primitive if (
        not only_primitive(c) or h >= 4 or (c-(t+1)*h)%8==0))
    if not primitive:
        return False
    # Weighted first and second moments bound the small coefficients jointly.
    # Bounds below are in units 1/(square B^2) for the smallest denominator B.
    joint_max = (square-1)*j+k
    neutral_max = (square-1)*v+h
    budget = (joint_max+neutral_max)//square
    bounded = [a for a in range(-isqrt(budget),isqrt(budget)+1) if a%t==0]
    # A pair with even coefficient has a disjoint partner with identical
    # degree, phase, and coefficient. At a fixed denominator these are tiles
    # of size two; odd coefficients use tiles of size one.
    small_costs = [set() for _ in range(v+1)]
    small_costs[0].add((0, 0, 0))
    for used in range(v+1):
        for a in bounded:
            size = (t-1) & -(t-1) if a%2 == 0 else 1
            if used+size > v:
                continue
            small_costs[used+size].update(
                (u+size*a, w+size*a*a, e+(size if a%2 == 0 else 0))
                for u, w, e in small_costs[used] if w+size*a*a <= budget)
    small = small_costs[v]
    high = h-v
    big_costs = [set() for _ in range(high+1)]
    big_costs[0].add((0, 0))
    for used in range(high+1):
        for a in range(8):
            size = (t-1) & -(t-1) if a%2 == 0 else 1
            if used+size <= high:
                big_costs[used+size].update(
                    ((u+size*a)%8, (w+size*a*a)%8)
                    for u,w in big_costs[used])
    large = big_costs[high]
    for c in primitive:
        if (k-c)%2:
            continue
        for uc in ((primitive_small,) if primitive_small is not None else range(c+1)):
            c_low = square*uc
            c_high = (square-1)*uc+c
            for A,Q,even_small in small:
                if h<4 and only_primitive(c) and t==5 and c==2 and uc==1:
                    # Three reciprocal squares cannot sum to (1+d^-2)/6.
                    continue
                if only_primitive(c):
                    parity_ok = False
                    dyadic = (t-1) & -(t-1)
                    for even in range(0,h+1,8 if t%4==1 else 4):
                        # In the fixed complementary phase, take invariants
                        # under its local sign. Every even coefficient moves
                        # in p-inertia orbits of length a multiple of dyadic.
                        # The vacuum supplies Steinberg in this same space.
                        if (2*dyadic*k < (t-1)*(2*dyadic+even*(dyadic-1))):
                            continue
                        if not even_small<=even<=even_small+h-v:
                            continue
                        if even==0:
                            if uc==c and c and ((c-(t+1)*h)%(8*((t+1)&-(t+1)))):
                                continue
                            if (c-(t+1)*h)%8 or A!=-t*v or Q!=t*t*v:
                                continue
                            gap_upper=(t+1)*((square-1)*v+h)-square*uc
                            gap_lower=(square-1)*uc+c-square*(t+1)*v
                            if gap_upper<0 or gap_lower<0 or (gap_upper==0 and c>uc) or (gap_lower==0 and h>v):
                                continue
                        # Integer-lattice mean at the adjacent residue classes.
                        # The sum of (a+t)^2 is t^2 or more on even a/t,
                        # and nonnegative on odd a/t. We retain exact small sums.
                        big_even=even-even_small
                        coefficient=(t+1)**2
                        # Tail identity plus the two residue-class minima;
                        # exact small coefficients are not relaxed.
                        upper=(coefficient*(h-v)-t*t*big_even
                               +square*((2*t+1)*v-Q-2*t*A))
                        if upper < (t+1)*square*uc:
                            continue
                        # Reflect the same lattice inequality about +t.
                        # Large even coordinates have a negative coefficient;
                        # dropping their positive weights gives an upper bound.
                        odd_large=h-v-big_even
                        reflected=((3*t-1)*c_high
                                   +square*(2*t*A-Q-(2*t-1)*v)
                                   +(t-1)**2*odd_large)
                        if reflected<0:
                            continue
                        parity_ok=True
                        break
                    if not parity_ok:
                        continue
                if ((h-c-A)%8,(n-Q)%8) not in large:
                    continue
                delta = max(0,square*(v-A)-c_high,c_low-square*(v-A))
                upper = (h-v)*(joint_max-square*v+2*c_high+2*square*A-square*Q)
                if delta*delta<=upper and (h>v or Q):
                    return True
    return False

def blocks(cells,weight,cap,max_high=None):
    # One complete equal-mass phase family (primitive 3t, or one order-9 class).
    yield ((0,0),)*cells
    for residue in range(2, min(8,cap//(weight*cells))+1,2):
        options=range(residue,cap//weight+1,8)
        for ns in product(options,repeat=cells):
            if sum(ns)*weight>cap:continue
            for js in product(*(range(max(0,n-(max_high//weight if max_high is not None else n)),n+1) for n in ns)):
                if any(j%2 or (n-j)%2 for n,j in zip(ns,js)):
                    continue  # Independent primary signs on joint labels.
                if max_high is not None and weight*sum(n-j for n,j in zip(ns,js))>max_high:continue
                # Common mass times B^2, since the next denominator is at least 3B.
                lower=max(js);upper=min(9*j+n-j for n,j in zip(ns,js))
                if any(n>j and j==lower for n,j in zip(ns,js)):
                    if 9*lower>=upper:continue
                elif 9*lower>upper:continue
                yield tuple(zip(ns,js))

@lru_cache(None)
def shapes(t,slots,max_high_labels=None):
    # All active counts are in dual pairs; each pair counted once.
    cap=slots//2; w=(t-1)//2
    max_high=cap if max_high_labels is None else min(cap,max_high_labels//2)

    L=[(4,w),(2,3*w),(2,3*w)]
    choices=[tuple(blocks(c,h,cap-2-(t-1),max_high)) for c,h in L]
    out={}
    for idx in product(*(range(len(c)) for c in choices)):
        A,B,C=(choices[i][j] for i,j in enumerate(idx))
        rows=[A[:2],A[2:],B,C];ww=(1,1,3,3)
        z=sum(w*v*(k+j) for ((k,_),(j,_)),v in zip(rows,ww))
        if not z or z+2+t-1>cap:continue
        # Enumerate exact solely-active mixed counts on every occupied phase
        # orbit. No mixed solely-active label can have an unoccupied phase:
        # J=C=0 would force its Weyl coefficient to be1, contrary to p|a.
        # The number above the least denominator remains only a lower bound.
        row_fib=[(t,v,w*sum(n for n,j in row),w*sum(j for n,j in row),row[0][0],row[0][1]) for v,row in zip(ww,rows)]
        col_fib=[(3,w,sum(v*row[c][0] for v,row in zip(ww,rows)),sum(v*row[c][1] for v,row in zip(ww,rows)),A[c][0],A[c][1]) for c in (0,1)]
        large_z=sum(w*v*sum(n-j for n,j in row) for v,row in zip(ww,rows))
        if large_z>max_high:continue
        # store all nondominated cost triples (non-minimal joint,solely-3,solely-t)
        opts={(0,0,large_z,0,0)}
        for side,fs in enumerate([row_fib,col_fib]):
            for q,v,k,j,c,uc in fs:
                if not k:continue
                o=set()
                for h in range(1,(cap-z)//v+1):
                    for a in reversed(range(h+1)):
                        if (h-a)*v<=max_high-large_z and axis(q,k,j,h,a,c,uc):
                            o.add((h*v,(h-a)*v));break
                opts={(x+(n if side==0 else 0),y+(n if side==1 else 0),zc,xc+(cost if side==0 else 0),yc+(cost if side==1 else 0)) for x,y,zc,xc,yc in opts for n,cost in o if z+x+y+n<=cap and zc+xc+yc+cost<=max_high}
        for x,y,zc,xc,yc in opts:
            key=(z,x,y,zc)
            out.setdefault(key,set()).add((xc,yc))
    shapes=[]
    for (z,x,y,zc),cs in out.items():
        for xc,yc in cs:
            if not any(u<=xc and v<=yc and (u,v)!=(xc,yc) for u,v in cs):
                shapes.append(tuple(2*a for a in (z,x,y,zc,xc,yc)))
    for high_t,multiple,q,min_pair,extra in ((True,t*w,3,2,2),(False,9,t,t-1,t-1)):
        for k in range(min_pair,cap//multiple+1):
            for h in range(1,(cap-extra)//multiple-k+1):
                for j in range(k+1):
                    for v in reversed(range(h+1)):
                        if not axis(q,k,j,h,v):continue
                        rec=(k*multiple,extra,h*multiple,(k-j)*multiple,0,(h-v)*multiple)
                        if not high_t:rec=(rec[0],rec[2],rec[1],rec[3],rec[5],rec[4])
                        if sum(rec[3:])<=max_high:shapes.append(tuple(2*a for a in rec))
                        break
    return shapes

