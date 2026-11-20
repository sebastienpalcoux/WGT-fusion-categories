"""Exact Gauss charges and phase imbalance on an odd 3--5 twist component."""
from bisect import bisect_left, bisect_right
import dense_budget
import ordinary_tail
import nine_gram
from small_phases import small_counts
from full_blocks import count as free_count, neutral_pair_floors
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd, isqrt, lcm, prod


def powers(p, cap):
    value = 1
    while value <= cap:
        yield value
        value *= p


@lru_cache(None)
def levels(cap, isolated):
    mixed = {a*b for a in powers(3,cap) for b in powers(5,cap) if 15 <= a*b <= cap}
    return sorted(mixed | {b for p in isolated for b in powers(p,cap) if b >= p*p})


@lru_cache(None)
def next_mixed(b, divisor=15, floor=15):
    """Least allowed 3,5-smooth multiple of divisor above the floor."""
    target=max(b,floor)
    values=levels(1<<(3*target).bit_length(),())
    return next(v for v in values[bisect_left(values,target):] if v%divisor==0)


@lru_cache(None)
def next_single(b,p):
    value=p*p
    while value<b:
        value*=p
    return value


def next_level_cap(residual,scale,slots,vacuum_factor):
    """Exact positive root of the remaining-charge bound with M >= factor*b."""
    value=abs(residual)
    if not value:
        raise ValueError("zero charge residual")
    return (scale+isqrt(scale*scale+6*slots*vacuum_factor**2*value*scale))//(2*vacuum_factor*value)


def upper_units(scale,denom):
    """Round up when an anticipated floor does not divide the node scale."""
    square=denom*denom
    return (scale+square-1)//square


@lru_cache(None)
def single_packets(b, slots, joint_order):
    """Group the single-prime packets by charge, retaining their primitive counts."""
    def packets(p):
        if b % p:
            return [(0,0,0)]
        return [((p-1)*(i+p*j),p*i+p*(p-1)*j,i)
                for j in range(slots//(p*(p-1))+1)
                for i in range(slots//(p-1)-p*j+1)
                if b%15 or (not i if p==3 and joint_order==45 else not j)]
    groups = {}
    for a,c,i in packets(3):
        for f,d,j in packets(5):
            if a+f <= slots:
                groups.setdefault((a,f,c+d),[]).append((i,j))
    return tuple((*key,tuple(values)) for key,values in groups.items())



def zero_moments(p,V,Z,lo,hi,lattice=None):
    """Necessary fiber moments, with zero coefficients on mass Z and J in [lo,hi]."""
    if V<=Z:
        return False
    step=p if lattice is None else lattice
    q=p-1
    J=max(lo,min(hi,((p*p-1)*V-q*q*Z)//2))
    if V*((p*p-1)*J-q*q*Z)<J*J+q*q*Z*J:
        return False
    T=V-Z
    centers={(q*V-lo)//(q*step*T),(q*V-hi)//(q*step*T),-1,1}-{0}
    for h in centers:
        J=hi if q+step*(2*h+1)>=0 else lo
        if q*(V+J)-step*(2*h+1)*(q*V-J)+q*step*step*h*(h+1)*T<0:
            return False
    return True


def triple_mass(rank,nines,shape,w5,wq,q,root):
    full_q=max(wq,free_count(q,nines))
    z,x,y,hard_z,hard_x,hard_y=shape
    activity=2*(q-1)
    spare = rank-1-nines-activity-full_q-z-x-y
    for a in range(0,min(x,w5)+1,2):
        for b in range(0,min(x-a,8)+1,2):
            other=w5-a+8-b
            if other>spare:
                continue
            large=max(0,hard_x-a-b)
            small=x-a-b-large
            mass=(Q(1,root**2)+Q(nines,81)+Q(activity,q**4)+Q(wq,(9*q**3)**2)+Q(full_q-wq,q**6)
                  +Q(w5,225**2)+Q(b,135**2)+Q(8-b,27**2)
                  +Q(z-hard_z+y-hard_y+small+spare-other,225)
                  +Q(hard_z+hard_y+large,45**2))
            if mass>=1:
                return True
    return False


def possible(rank,support,nines,shape,witnesses,marked=None,profile_shapes=None):
    """Test a whole shape. Joint denominators above fifteen are chosen in order."""
    marked={} if marked is None else marked
    support=tuple(support)
    z,x,y,big_z,big_x,big_y=shape
    if profile_shapes is None:
        profile_shapes=(shape,)
    if any(s[0]!=z or s[3]!=big_z for s in profile_shapes):
        raise ValueError("joint counts differ within a profile family")
    states=tuple(sorted({(s[1],s[2],s[4],s[5]) for s in profile_shapes}))
    x,y,big_x,big_y=(min(t[i] for t in states) for i in range(4))
    shape=(z,x,y,big_z,big_x,big_y)
    small_five_rigid=all(s[2]<72 for s in profile_shapes)
    minimal_five_allowed=any(s[2]==40 for s in profile_shapes)
    if any(v!=(1,2*p,2*(p-1)) for p,v in marked.items()):
        raise ValueError("unexpected prime-denominator class",marked)
    isolated=tuple(p for p in support if p not in (3,5) and p not in marked)
    initial_counts=tuple(2*(p-1) for p in isolated)
    occupied=1+nines+sum(v[1] for v in marked.values())
    used=Q(nines,81)+sum(Q(v[1],p*p) for p,v in marked.items())
    deficit=1-sum(Q(2,p) for p in marked)
    if any(dense_budget.possible(rank,support,marked,nines,s,p,witnesses)
           for p in isolated for s in profile_shapes):
        raise ArithmeticError("unresolved dense isolated prime",support,nines,shape)
    if len(support)==3 and not marked:
        q,=isolated
        if nines<9:
            raise ValueError("nine-row support requires at least nine labels")
        w5,wq=(witnesses(t,nines) for t in (5,q))
        root=27*25*q**3
        if not triple_mass(rank,nines,shape,w5,wq,q,root):
            return False
        neutral=max(wq,free_count(q,nines))
        omit3=8
    elif len(support) in (4,5):
        parts={3:9,5:5,**{p:p for p in marked},**{p:p**3 for p in isolated}}
        pair_floors=neutral_pair_floors(support,marked,parts,nines)
        neutral=2*len(pair_floors)
        root=45*prod(p**3 for p in isolated)*prod(marked)
        w5=omit3=0
    else:
        raise ValueError("unresolved joint support",support,marked)
    # A sparse isolated prime forces a neutral pure pair outside the
    # p-free reservations when their maximum possible number is too small.
    pure_quota=rank-occupied-z-x-y-sum(initial_counts)
    for p in isolated:
        if pure_quota>=2*p*p*(1-Q(1,p-1)-Q(1,6*p)):
            raise ArithmeticError("unresolved pure neutral quota",p,pure_quota)
    neutral+=2*len(isolated)
    qcaps=tuple(min((p-1)**2,2*(p-1)*((rank-occupied-neutral-z-x-y
                        -sum(initial_counts)+2*(p-1))//(2*(p-1)))) for p in isolated)
    if any(cap<2*(p-1) for p,cap in zip(isolated,qcaps)):
        return False
    qexact=tuple(cap==2*(p-1) for p,cap in zip(isolated,qcaps))
    if z not in (32,48,64) or z < big_z:
        raise ValueError("joint phase order is unresolved",shape)
    # Across all embeddings the means at orders 15 and 45 are 1/8 and 0.
    joint_charge=Q(1) if z==48 else Q(7,8)
    packet=12 if z==48 else 4
    small_joint=z-big_z
    slots=rank-occupied-neutral-small_joint
    used+=Q(small_joint,225)
    deficit-=joint_charge*Q(small_joint,225)

    caps=small_counts(z,big_z,rank-occupied-neutral-z-sum(initial_counts))

    # At each node all weights have the common denominator 8 S².
    # S also contains the fixed floors; only L is asserted to divide M.
    scale_root=lcm(root,675)
    root_scale=8*scale_root**2
    def exact_units(value):
        result=value*root_scale
        if result.denominator!=1:
            raise ArithmeticError("incomplete common denominator",value)
        return result.numerator
    initial_mass=exact_units(used)
    initial_charge=exact_units(deficit)
    rate_scale=lcm(8,*(p-1 for p in isolated))
    joint_numerator=8 if z==48 else 7
    joint_order=45 if z==48 else 15
    @lru_cache(None)
    def primary_parts(L):
        current={}
        for p in support:
            power=1
            rest=L
            while rest%p==0:
                power*=p
                rest//=p
            current[p]=power
        return current

    @lru_cache(None)
    def vacuum_factor(L):
        current=primary_parts(L)
        return min([L//(current[3]*current[5])]
                   +[p*(L//current[p]) for p in isolated])

    @lru_cache(None)
    def neutral_floors(L,nine):
        current=primary_parts(L)
        if len(support)==3:
            q,=isolated
            if nine:
                return (9*current[q],)*(neutral//2-1)
            return ((9*current[q],)*(wq//2)
                    +(current[q],)*((neutral-2-wq)//2))
        return tuple(neutral_pair_floors(support,marked,current,nines,nine=nine))

    def neutral_mass(L,scale,nine=False):
        return 2*sum(upper_units(scale,b) for b in neutral_floors(L,nine))

    def limits(k,x,y,q,o3,o5,unit,scale):
        # Integer floors suffice because the two imbalances have integral coordinates.
        three=((k-y-q-o3-o5)*unit+o3*min(unit,scale//729)+o5*min(unit,scale//50625))//2
        five=(k-q-max(x,o3+o5))*unit//4
        return three,five

    @lru_cache(None)
    def summary(states):
        return (*tuple(min(s[i] for s in states) for i in range(4)),
                max(s[0] for s in states),max(s[1] for s in states),
                max(s[0]+s[1] for s in states))

    def visit(last,k,U,mass,L,S,states,zleft,qs,qmax,qfirst,o3,o5,balances,eligible,cosine,cosine_square,mixed_three,mixed_five,mixed_joint,small_three,small_five,count_three,count_five,three45,five45):
        if not states:
            return False
        scale=8*S*S
        outside=sum(qs)+zleft
        next_mix=upper_units(scale,next_mixed(last+1))
        charge_room=U+scale//L
        if charge_room<0:
            return False
        future_three=charge_room if z==48 else (2*charge_room+2)//3
        future_five=(4*charge_room+4)//5
        neutral_upper=(neutral_mass(L,scale)+scale//(L*L)
                       +2*sum(upper_units(scale,b or next_single(last+1,p))
                              for p,b,exact in zip(isolated,qfirst,qexact) if not (b and exact)))
        # Each exact state must satisfy its own mass and slot bounds before
        # coordinate envelopes are used to share the remaining search.
        retained=[]
        for sx,sy,shx,shy in states:
            if k<max(sx,shx,o3+o5)+max(sy,shy)+outside:
                continue
            free=k-sx-sy-outside
            high_mix=upper_units(scale,next_mixed(last+1,floor=45))
            Aupper=(sx-shx)*next_mix+shx*high_mix
            Bupper=(sy-shy)*next_mix+shy*high_mix
            pure_floor=min(next_single(last+1,p) for p in (3,5)+isolated)
            excess_upper=((Aupper+1)//2 if z!=48 else 0)+(Bupper+3)//4
            excess_upper+=(free*upper_units(scale,pure_floor)+1)//2
            excess_upper+=sum((n*upper_units(scale,next_single(last+1,p))+p-2)//(p-1)
                              for p,n in zip(isolated,qs))
            if U-scale//L+mass-scale>excess_upper:
                continue
            V=mixed_three+min(sx*next_mix,future_three)
            if 24*V<mixed_joint or (z!=48 and 8*(mixed_five+min(sy*next_mix,future_five))<mixed_joint):
                continue
            lo=mixed_joint
            hi=lo+zleft*upper_units(scale,next_mixed(last+1,floor=45))
            Z=small_three*(scale//225)
            lattice=5
            if z==32 and big_z==16:
                total_x=sx+count_three
                total_y=sy+count_five
                if total_x<=20:
                    odd=(total_x-16)*upper_units(scale,45)
                    if 22*V<3*lo-50*odd or 18*mixed_three>7*hi+50*odd:
                        continue
                    Z+=max(0,three45-(total_x-16))*(scale//2025)
                    if total_x==16:
                        lattice=10
                W=mixed_five+min(sy*next_mix,future_five)
                odd=(total_y-32)*upper_units(scale,45)
                if 7*W<2*lo-9*odd or 5*mixed_five>4*hi+9*odd:
                    continue
                Z5=small_five*(scale//225)
                if total_y<=44:
                    Z5+=max(0,five45-(total_y-32))*(scale//2025)
                if not zero_moments(3,W,Z5,lo,hi):
                    continue
                if total_x<=20:
                    # Compare final charge excess with ordinary mass. The
                    # two reflected windows give E >= 9J/616-25O3/22-9O5/28.
                    odd3=(total_x-16)*upper_units(scale,45)
                    excess=(86*lo-308*mixed_three-154*mixed_five-700*odd3-198*odd)
                    free=k-sx-sy-outside
                    room=mass+charge_room+neutral_upper+free*upper_units(scale,15)-scale
                    if 616*room<excess:
                        continue
            if z in (32,48) or (z==64 and big_z):
                # Common phase mass makes every small coefficient vanish.
                # Joint mass may still range over the indicated interval.
                if not zero_moments(5,V,Z,lo,hi,lattice):
                    continue
            retained.append((sx,sy,shx,shy))
        states=tuple(retained)
        if not states:
            return False
        x,y,hx,hy,max_x,max_y,max_sum=summary(states)
        left=max(x,hx,o3+o5)
        right=max(y,hy)
        spare=k-left-right-outside
        bound=scale//L
        if hx>x or hy>y or spare<0 or not U or U < -bound:
            return False
        # In each primitive-five fiber, subtract the two moment identities.
        # A nonzero coefficient at fifteen saturates the budget, so its
        # solely-five mass is exactly 5/225 (in pair units).  With fewer
        # than thirteen pairs in a fiber this forces all five pairs to be at15.
        # The remaining fibers have zero coefficient at every small label.
        if z==32 and big_z==0 and small_five_rigid and last>=15:
            good=False
            totalZ=small_five*(scale//225)
            maximum=mixed_five+min(max_y*next_mix,future_five)
            for saturated in range(3):
                sat=20*saturated*(scale//225)
                if totalZ<sat or mixed_five<sat:
                    continue
                J=(32-16*saturated)*(scale//225)
                Z=totalZ-sat
                V=min(maximum-sat,5*J//4)
                if not J:
                    if minimal_five_allowed and small_five==40:
                        good=True;break
                    continue
                if V<Z or V*(8*J-4*Z)<J*J+4*Z*J:
                    continue
                T=V-Z
                if T<=0:
                    continue
                A2=2*V-J
                center=A2//(6*T)
                if any(2*(V+J)-3*(2*h+1)*A2+18*h*(h+1)*T<0
                       for h in {center,-1,1}-{0}):
                    continue
                good=True;break
            if not good:
                return False
        # Couple the phase windows with the entire remaining ordinary mass.
        # Pure activity uses only positions outside the exact mixed quotas.
        unassigned=max(0,k-x-y-outside)
        pure_three=unassigned*upper_units(scale,next_single(last+1,3))
        if z==48:
            balances={(d,e) for d,e in balances if 2*d<=pure_three}
        else:
            room3=9*(scale-mass)+9*mixed_joint-16*mixed_three+16*pure_three
            balances={(d,e) for d,e in balances if 50*d<=room3}
        if not balances:
            return False
        # Average the real-part Gram inequality after all already chosen
        # levels. The remaining mixed three-only positions have floor 45.
        demand=max(x,hx)
        tail=(neutral_mass(L,scale,True)+scale//(L*L)
              +min((k-outside-right-demand)*(scale//729)
                   +demand*min(scale//2025,scale//(last*last)),
                   (k-outside-right-o5)*(scale//729)+o5*(scale//50625)))
        H=nines*(scale//81)+eligible+tail
        # u=4 times the unknown remaining real part; the minimizer of the
        # squared Gram inequality is the projection of H-4*cosine onto its interval.
        u=min(4*tail,max(-2*tail,H-4*cosine))
        if (4*cosine+u)**2>2*H*(8*cosine_square+4*tail+u):
            return False
        if ((0,0,0,0) in states and not outside and (0,0) in balances
                and scale%abs(U)==0 and (scale//abs(U))%L==0
                and ordinary_tail.possible(scale//abs(U),support,marked,L,nines,k,neutral,o3,o5,
                                           mass,scale,witnesses,qfirst,qexact)):
            return True
        base_book=mass+neutral_mass(L,scale)+scale//L**2
        pure_mass=2*sum(upper_units(scale,b or next_single(last+1,p)) for p,b,exact in zip(isolated,qfirst,qexact) if not (b and exact))
        book=base_book+pure_mass
        free=k-max(x,hx)-right-outside
        if (book+max(0,U+bound)+(8-joint_numerator)*zleft*(scale//(8*45**2))
                +free*(scale//225)<scale):
            return False
        if 3*mass+2*(U-bound)>3*scale:
            return False
        factor=vacuum_factor(L)
        # Solve the exact quadratic bound from M >= factor*b.
        cap=next_level_cap(U,scale,k,factor)
        candidates=levels(1<<cap.bit_length(),isolated)
        for b in candidates[bisect_right(candidates,last):bisect_right(candidates,cap)]:
            b2=b*b
            if 2*(U-bound)*b2>3*k*scale:
                break
            prime=next((p for p in isolated if b%p==0),None)
            multiple=lcm(L,b2//prime if prime else b)
            S1=lcm(S,multiple)
            up=(S1//S)**2
            scale1=scale*up
            unit=scale1//b2
            book1=base_book*up+2*sum(upper_units(scale1,level or next_single(b,p)) for p,level,exact in zip(isolated,qfirst,qexact) if not (level and exact))
            # In each reserved twist sector the future levels belong to
            # a multiplicative semigroup. Replace a shared numerical floor
            # by its next member, rounding upper masses outwards.
            mix=upper_units(scale1,next_mixed(b))
            high_mix=upper_units(scale1,next_mixed(b,floor=45))
            prime_masses=[n*upper_units(scale1,next_single(b,p)) for p,n in zip(isolated,qs)]
            isolated_mass=sum(prime_masses)
            three_mass=(max(x,hx)-hx)*mix+hx*high_mix
            five_mass=(right-hy)*mix+hy*high_mix
            ordinary=three_mass+five_mass+zleft*high_mix+isolated_mass
            room=book1+free*(scale1//225)+ordinary-scale1
            if room<0:
                break
            # One unassigned position can be neutral at 15, or active from b
            # upwards. A positive multiple of charge makes their upper bounds equal.
            iso_charge=sum((p*rate_scale//(p-1))*m for p,m in zip(isolated,prime_masses))
            required_charge=((rate_scale if z==48 else 3*rate_scale//2)*three_mass
                            +(5*rate_scale//4)*five_mass+iso_charge
                            +(joint_numerator*rate_scale//8)*zleft*high_mix)
            if 675*rate_scale*room<2*(b2-225)*(rate_scale*(U-bound)*up-required_charge):
                break
            if len(support)==3:
                # The omission sets are disjoint from one another. Their
                # intersections with a mixed three-sector have larger floors.
                best=0
                overlap_three=upper_units(scale1,next_mixed(b,135,135))
                overlap_five=upper_units(scale1,next_mixed(b,225,225))
                for five_overlap in range(0,min(max(x,hx),o5)+1,2):
                    for three_overlap in range(0,min(max(x,hx)-five_overlap,o3)+1,2):
                        excess=max(0,hx-five_overlap-three_overlap)
                        ordinary_mixed=max(x,hx)-five_overlap-three_overlap-excess
                        available=k-outside-right-max(x,hx)-o3-o5+five_overlap+three_overlap
                        if available<0:
                            continue
                        upper=(book1+right*mix+isolated_mass+zleft*high_mix
                            +excess*high_mix+ordinary_mixed*mix
                            +three_overlap*overlap_three+five_overlap*overlap_five
                            +(o3-three_overlap)*(scale1//729)+(o5-five_overlap)*(scale1//50625)
                            +available*(scale1//225))
                        best=max(best,upper)
                if best<scale1:
                    break
            elif (book1+o3*(scale1//729)+o5*(scale1//50625)+(k-o3-o5-right-outside)*(scale1//225)
                    +five_mass+isolated_mass+zleft*high_mix<scale1):
                break
            max3,max5=limits(k,x,right,outside,o3,o5,unit,scale1)
            feasible=[(d*up,e*up) for d,e in balances if d*up<=max3 and e*up<=max5]
            if not feasible:
                break
            cap_charge=U*up+scale1//multiple
            if prime:
                choices=((0,0,0,q,8*prime*(q//(prime-1)),((0,0),))
                         for q in range(2*(prime-1),min(k,qmax[isolated.index(prime)])+1,2*(prime-1)))
            else:
                joint_choices=range(0,min(k,zleft)+1,packet) if b%15==0 and b>15 else (0,)
                if z==32 and big_z==16 and b%15==0 and b>15:
                    joint_choices=(0,zleft) if zleft and zleft<=k else (0,)
                choices=((a,f,j,0,8*c+joint_numerator*j,pairs) for j in joint_choices
                         for a,f,c,pairs in single_packets(gcd(b,15),min(k-j,max_sum) if b%15==0 else k-j,joint_order)
                         if j+a+f and (b%15 or (a<=max_x and f<=max_y)))
            for a,f,j,q,ncharge,packets in choices:
                if mass*up+(a+f+j+q)*unit>=scale1:
                    continue
                charge=ncharge*(unit//8)
                if charge>cap_charge:
                    continue
                next_three=mixed_three*up+(a*unit if b%15==0 else 0)
                next_five=mixed_five*up+(f*unit if b%15==0 else 0)
                next_joint=mixed_joint*up+j*unit
                joint_upper=next_joint+(zleft-j)*unit
                if (16*next_three>9*joint_upper
                        or (2*next_five>joint_upper if z==48 else 4*next_five>5*joint_upper)):
                    continue
                k1=k-a-f-j-q
                mixed=b%15==0
                at_first=(last==9)
                if at_first and caps is not None:
                    a15,f15=(a,f) if b==15 else (0,0)
                    if a15 not in caps[0] or f15 not in caps[1]:
                        continue
                remaining=set()
                for sx,sy,shx,shy in states:
                    if at_first and caps is not None:
                        if sx not in caps[0][a15] or sy not in caps[1][f15]:
                            continue
                    if mixed and (a>sx or f>sy):
                        continue
                    nx=sx-a if mixed else sx
                    ny=sy-f if mixed else sy
                    nhx=max(0,shx-a) if mixed and b>15 else shx
                    nhy=max(0,shy-f) if mixed and b>15 else shy
                    if nhx<=nx and nhy<=ny:
                        # Every subsequent mixed level is strictly above15.
                        remaining.add((nx,ny,nx,ny))
                if not remaining:
                    continue
                next_states=tuple(sorted(remaining))
                x1,y1,hx1,hy1,*_=summary(next_states)
                qs1=tuple(max(0,n-q) if t==prime else n for t,n in zip(isolated,qs))
                qfirst1=tuple(b if t==prime and not level else level for t,level in zip(isolated,qfirst))
                qmax1=tuple(n-q if t==prime else n for t,n in zip(isolated,qmax))
                o3new=max(0,o3-a) if b%27==0 else o3
                o5new=max(0,o5-a) if b%225==0 else o5
                other=sum(qs1)+zleft-j
                if k1<max(x1,hx1,o3new+o5new)+max(y1,hy1)+other:
                    continue
                m3,m5=limits(k1,x1,max(y1,hy1),other,o3new,o5new,unit,scale1)
                next_balances=set()
                for i,t in packets:
                    for d,e in feasible:
                        ds=[abs(d+v*unit) for v in range(-i,i+1,2) if abs(d+v*unit)<=m3]
                        if ds:
                            es=[abs(e+v*unit) for v in range(-t,t+1,2) if abs(e+v*unit)<=m5]
                            next_balances.update((dd,ee) for dd in ds for ee in es)
                exact_pair=2 if prime and qexact[isolated.index(prime)] else 0
                active= a if b%9==0 else 0
                first=packets[0][0] if active else 0
                if next_balances and visit(b,k1,U*up-charge,mass*up+(a+f+j+q+exact_pair)*unit,multiple,S1,
                        next_states,zleft-j,qs1,qmax1,qfirst1,o3new,o5new,next_balances,eligible*up+active*unit,cosine*up-first*unit,
                        cosine_square*up+(active-first)*(unit//2),next_three,next_five,next_joint,small_three+(a if b==15 else 0),small_five+(f if b==15 else 0),count_three+(a if mixed else 0),count_five+(f if mixed else 0),three45+(a if b==45 else 0),five45+(f if b==45 else 0)):
                    return True
        return False

    tail=nine_gram.ceiling(rank,support,nines,shape,marked,witnesses)
    import native_gauss
    if native_gauss.applicable(locals()):
        starts=[(c,[i for i in range(c%2,c+1,2) if nine_gram.possible(nines,c,i,tail)])
                for c in range(nines//2+1)]
        return native_gauss.possible(locals(),[(c,i) for c,i in starts if i])
    for c in range(nines//2+1):
        balances={(i*(root_scale//81),0) for i in range(c%2,c+1,2)
                  if nine_gram.possible(nines,c,i,tail)}
        if balances and visit(9,slots,initial_charge-c*(root_scale//27),initial_mass,root,
                              scale_root,states,big_z,initial_counts,qcaps,tuple(0 for p in isolated),omit3,w5,balances,0,
                              (nines-3*c)*(root_scale//81),0,0,0,small_joint*(root_scale//225),0,0,0,0,0,0):
            return True
    return False
