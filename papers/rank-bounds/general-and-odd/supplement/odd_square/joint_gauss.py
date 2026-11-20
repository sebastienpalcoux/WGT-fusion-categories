"""Exact Gauss charges and phase imbalance on an odd 3--5 twist component."""
from bisect import bisect_left, bisect_right
import dense_budget
import ordinary_tail
import nine_gram
from small_phases import small_counts
from profile_window import phase_window
from pure_window import power_brackets
from full_blocks import omission_count as free_count, neutral_pair_floors
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import gcd, isqrt, lcm, prod


def dyadic_pure_costs():
    # Minimal slots for a prescribed mod-eight charge and primitive
    # three/five imbalance parities, allowing every pure-power phase.
    costs={(0,0,0):0}
    for _ in range(32):
        for (r,a,b),v in tuple(costs.items()):
            for t,c,x,y in ((2,3,1,0),(6,6,0,0),(4,5,0,1),(20,20,0,0)):
                k=((r+c)%8,a^x,b^y)
                costs[k]=min(costs.get(k,2000),v+t)
    return {(r,i,j):costs.get((r,i,j),2000) for r in range(8) for i in range(2) for j in range(2)}

PURE_DYADIC=dyadic_pure_costs()


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


@lru_cache(None)
def three_layout_prefix(z,population,last):
    """Conditional Egyptian fractions for a rigid component at fifteen."""
    from itertools import product
    @lru_cache(None)
    def fiber(left, gap, floor):
        if not left:
            return ((),) if not gap else ()
        if gap<=0:
            return ()
        cap=isqrt(left*gap.denominator//gap.numerator)
        patterns=[]
        for b in (i for i in levels_small(cap) if i>=floor):
            sq=b*b
            for count in range(1,min(left,gap.numerator*sq//gap.denominator)+1):
                residual=gap-Q(count,sq)
                patterns.extend(((b,count),)+part for part in
                                fiber(left-count,residual,b+1))
        return patterns
    def levels_small(cap):
        return sorted(a*b for a in powers(3,cap) for b in powers(5,cap//a))
    if population in (12,18):
        sizes=(3,3) if z==32 else (3,3,3)
    elif z==32:
        sizes=(3,11)
    else:
        sizes=(11,11,11)
    # The three ninth-root phases are in one Galois orbit, with
    # identical denominator multiplicities.
    fs={h:fiber(h,Q(1,3),3) for h in set(sizes)}
    choices=(product(*(fs[h] for h in sizes)) if z==32 else
             ((a,a,a) for a in fs[sizes[0]]))
    result=set()
    for selection in choices:
        items=[(b,2*n) for f in selection for b,n in f if 15*b<=last]
        result.add((sum(n for b,n in items),sum(n for b,n in items if b==3),
                    sum((Q(n,b*b) for b,n in items),Q())))
    return result


def shapes_root_admissible(z,big,states,caps):
    if z not in (32,48) or big or caps is None:
        return True
    a,b=caps
    return any(x in a.get(0,()) for x,_,_,_ in states)


def possible(rank,support,nines,shape,witnesses,marked=None,profile_shapes=None,seven=None):
    """Test a whole shape. Joint denominators above fifteen are chosen in order."""
    marked={} if marked is None else marked
    support=tuple(support)
    z,x,y,big_z,big_x,big_y=shape
    if z in (32,48) and 0<big_z<z:return False
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
    initial_counts=tuple((4 if p%4==1 else 2)*(p-1) for p in isolated)
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
        # Vacuum orthogonality of an isolated pure q-row requires both
        # a truly neutral three-free and a truly neutral five-free pair.
        extra_neutral=4
        neutral+=extra_neutral
        omit3=6
        w5=max(0,w5-2)
    elif len(support) in (4,5):
        extra_neutral=0
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
    if len(support)==3 and nines>=9:
        # Three- and five-free witnesses are disjoint and neither can be
        # pure at the third prime. Only the mixed solely-three count x
        # can intersect these built-in reservations.
        overlap=x
        if not big_z and z in (32,48) and all(
                row[1]==(12 if z==32 else 18) for row in profile_shapes):
            overlap=0  # all denominators there are 45
        pure_quota-=max(0,8+w5-overlap)
    for p in isolated:
        if pure_quota>=2*p*p*(1-Q(1,p-1)-Q(1,6*p)):
            raise ArithmeticError("unresolved pure neutral quota",p,pure_quota)
    neutral+=2*len(isolated)
    qcaps=tuple(min((p-1)**2,2*(p-1)*((rank-occupied-neutral-z-x-y
                        -sum(initial_counts)+floor)//(2*(p-1)))) for p,floor in zip(isolated,initial_counts))
    if any(cap<floor for floor,cap in zip(initial_counts,qcaps)):
        return False
    qexact=tuple(cap==2*(p-1) for p,cap in zip(isolated,qcaps))
    initial_first=tuple(0 for p in isolated)
    prebooked=0
    if seven is not None:
        if isolated!=(7,):raise ValueError("seven activity split outside domain")
        total,square,cube=seven
        if not initial_counts[0]<=total<=qcaps[0]:return False
        assert total%12==square%12==cube%12==0 and square+cube<=total
        initial_counts=(total,)
        prebooked=square+cube
        qexact=(total==12,)
        bfirst=49 if square else 343 if cube else 0
        initial_first=(bfirst,)
        # The largest active level is known, or bounded below, throughout
        # the entire branch. Its charge digit lies in {2,4,6}, a 7-unit.
        biggest=2401 if prebooked<total else 343 if cube else 49
        root=lcm(root,biggest**2//7)
        shell=(2 if qexact[0] and bfirst else 0)
        used+=Q(square+(shell if square else 0),49**2)+Q(cube+(shell if cube and not square else 0),343**2)
        deficit-=Q(7*square,6*49**2)+Q(7*cube,6*343**2)
        qcaps=(total-prebooked,)

    if z not in (32,48,64) or z < big_z:
        raise ValueError("joint phase order is unresolved",shape)
    # Across all embeddings the means at orders 15 and 45 are 1/8 and 0.
    joint_charge=Q(1) if z==48 else Q(7,8)
    packet=12 if z==48 else 8
    small_joint=z-big_z
    slots=rank-occupied-neutral-small_joint-prebooked
    used+=Q(small_joint,225)
    deficit-=joint_charge*Q(small_joint,225)

    caps=small_counts(z,big_z,rank-occupied-neutral-z-sum(initial_counts))
    if not shapes_root_admissible(z,big_z,states,caps):return False

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
    exact_sign=-1 if 11 in marked and support==(3,5,7,11) else 0
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

    limit_multiplier=100000
    @lru_cache(None)
    def root_multipliers(fixed7):
        return sorted(a*b*c for a in powers(3,limit_multiplier)
                      for b in powers(5,limit_multiplier//a)
                      for c in ((1,) if fixed7 else powers(49,limit_multiplier//(a*b))))
    @lru_cache(None)
    def minimum_root(L):
        if isolated!=(7,) or len(marked)>1 or any(p!=11 for p in marked):return L
        known=primary_parts(L)
        fixed7=False
        digits=(2,4,6)
        if seven is not None:
            total,sq,cu=seven
            unplaced=total-sq-cu
            fixed7=not unplaced
            if fixed7:digits=((cu or sq)//6,)
            else:digits=tuple(range(2,unplaced//6+1,2))
        for multiplier in root_multipliers(fixed7):
            value=L*multiplier
            qroot=known[7]
            other=value//qroot
            while other%7==0:other//=7
            for sign in (-1,1):
                if 11 in marked and (2*(value//11)+sign)%11:continue
                if any((other*digit+sign)%49==0 for digit in digits):
                    return value
        return L*limit_multiplier
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
            new=(current[3]*q*q,current[5]*q*q)
            if nine:
                return ((9*current[q],)*((neutral-extra_neutral)//2-1)
                        +tuple(lcm(9,b) for b in new))
            return ((9*current[q],)*(wq//2)
                    +(current[q],)*((neutral-extra_neutral-2-wq)//2)+new)
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

    def visit(last,k,U,mass,L,S,states,zleft,qs,qmax,qfirst,o3,o5,balances,eligible,cosine,cosine_square,mixed_three,mixed_five,mixed_joint,small_three,small_five,count_three,count_five,three45,five45,lock_joint=False,three75=0,five75=0,three135=0,five135=0):
        global remaining_visits
        remaining_visits -= 1
        if remaining_visits < 0:
            raise RuntimeError("budget")
        if not states:
            return False
        scale=8*S*S
        outside=sum(qs)+zleft
        next_mix=upper_units(scale,next_mixed(last+1))
        root_floor=minimum_root(L)
        charge_room=U+((scale+root_floor-1)//root_floor)
        if charge_room<0:
            return False
        future_three=charge_room if z==48 else (2*charge_room+2)//3
        future_five=(4*charge_room+4)//5
        neutral_upper=(neutral_mass(L,scale)+upper_units(scale,root_floor)
                       +2*sum(upper_units(scale,b or next_single(last+1,p))
                              for p,b,exact in zip(isolated,qfirst,qexact) if not (b and exact)))
        # Each exact state must satisfy its own mass and slot bounds before
        # coordinate envelopes are used to share the remaining search.
        retained=[]
        residual_bounds=[]
        bound=((scale+root_floor-1)//root_floor)
        pure3=upper_units(scale,next_single(last+1,3))
        pure5=upper_units(scale,next_single(last+1,5))
        pure_charge=(max(6*pure3,5*pure5)+3)//4
        iso_units=[upper_units(scale,next_single(last+1,p)) for p in isolated]
        pure_charge=max([pure_charge]+[(m*p+p-2)//(p-1) for m,p in zip(iso_units,isolated)])
        q_reserve_mass=sum(n*m for n,m in zip(qs,iso_units))
        q_reserve_charge=sum((n*p*m+p-2)//(p-1) for n,p,m in zip(qs,isolated,iso_units))
        for sx,sy,shx,shy in states:
            if k<max(sx,shx,o3+o5)+max(sy,shy)+outside:
                continue
            free=k-sx-sy-outside
            base_charge=((7*zleft//8 if z!=48 else zleft)+(3*sx//2 if z!=48 else sx)+5*sy//4)
            compatible=False
            for qt in product(*(range(n,c+1,2*(p-1)) for p,n,c in zip(isolated,qs,qmax))):
                pure_slots=free+sum(qs)-sum(qt)
                if pure_slots<0:continue
                fixed=base_charge+sum(p*(n//(p-1)) for p,n in zip(isolated,qt))
                for d,e in balances:
                    dp=((d//8)^(sx//2 if z==64 else 0))%2
                    ep=((e//8)^(sy//4))%2
                    if any(PURE_DYADIC[((U//8-fixed-sign)%8,dp,ep)]<=pure_slots for sign in (-3,-1,1,3)):
                        compatible=True;break
                if compatible:break
            if not compatible:continue
            high_mix=upper_units(scale,next_mixed(last+1,floor=45))
            Aupper=(sx-shx)*next_mix+shx*high_mix
            Bupper=(sy-shy)*next_mix+shy*high_mix
            pure_floor=min(next_single(last+1,p) for p in (3,5)+isolated)
            excess_upper=((Aupper+1)//2 if z!=48 else 0)+(Bupper+3)//4
            excess_upper+=(free*upper_units(scale,pure_floor)+1)//2
            excess_upper+=sum((n*upper_units(scale,next_single(last+1,p))+p-2)//(p-1)
                              for p,n in zip(isolated,qs))
            if U-((scale+root_floor-1)//root_floor)+mass-scale>excess_upper:
                continue
            if not big_z and last>=15 and z in (32,48,64):
                parts=tuple(next_single(last+1,p) for p in (3,5)+isolated)
                pure_mass_unit=max(upper_units(scale,level) for level in parts)
                remaining5=max(0,o5-sx)
                remaining3=max(0,o3-max(0,sx-o5))
                neutral_loss=(remaining5*(scale//225-upper_units(scale,225))+
                              remaining3*(scale//225-upper_units(scale,27)))

                lo3=mixed_three
                hi3=lo3+min(sx*next_mix,future_three)
                lo5=mixed_five
                hi5=lo5+min(sy*next_mix,future_five)
                def windows(side,total,v,w,low,high):
                    w75=(three75 if side=='three' else five75) if last>=75 else None
                    w135=(three135 if side=='three' else five135) if last>=135 else None
                    for a,b in phase_window(z,side,total,v,w if last>=45 else None,w75,w135):
                        l=max(low,(a.numerator*scale)//a.denominator)
                        u=min(high,(b.numerator*scale+b.denominator-1)//b.denominator)
                        if l<=u:yield (l-low,u-low)
                one=windows('three',sx+count_three,small_three,three45,lo3,hi3)
                two=list(windows('five',sy+count_five,small_five,five45,lo5,hi5))
                surviving=False; bound_abs=None
                room_pure=free*pure_charge+q_reserve_charge
                ordinary_other=neutral_upper+free*(scale//225)+q_reserve_mass
                lx=4 if z==48 else 6
                for u,v in one:
                    for a,b in two:
                        excess=mass+v+b+ordinary_other-scale
                        if excess<neutral_loss or mass+u+a>scale:continue
                        slots_active=free-remaining5
                        cost=scale//225-pure_mass_unit
                        if cost>0:slots_active=min(slots_active,excess//cost)
                        if slots_active<0:continue
                        bracket=(power_brackets(parts,slots_active,qs,qmax,isolated)
                                 if last<=2500 and slots_active<=44 else None)
                        # Bounds on four times the further pure charge plus signed error.
                        gl=4*U-lx*v-5*b
                        gu=4*U-lx*u-5*a
                        if gl>4*(bound+slots_active*pure_charge+q_reserve_charge) or gu < -4*bound:continue
                        if bracket:
                            bs,ends,starts=bracket
                            lower=gl if exact_sign==-1 else gl-4*bound
                            upper=gu+4*bound
                            pos=bisect_left(ends,-((-lower*bs)//(4*scale)))
                            if pos==len(ends) or starts[pos]*4*scale>upper*bs:continue
                        if (45*(mass+v+b+ordinary_other-scale)
                            +16*(lx*v+5*b-4*U+4*bound+4*q_reserve_charge)<0):continue
                        # For p=3,5 the largest possible pure activity at 25 or
                        # above has combined mass+64/45*charge <=1/225 per slot.
                        least=max(0,gl,-gu)
                        bound_abs=least if bound_abs is None else min(bound_abs,least)
                        surviving=True
                if not surviving:continue
                residual_bounds.append((bound_abs,k-sx-sy))
            V=mixed_three+min(sx*next_mix,future_three)
            if 24*V<mixed_joint or (z!=48 and 8*(mixed_five+min(sy*next_mix,future_five))<mixed_joint):
                continue
            lo=mixed_joint
            hi=lo+zleft*upper_units(scale,next_mixed(last+1,floor=45))
            if z in (32,48):
                if 24*mixed_three>hi or 24*V<lo:
                    continue
                if small_three:
                    continue
                if not zleft:
                    # All joint mass is known. Every remaining nonzero
                    # solely-3 charge is at most the remaining exact mass.
                    gap=lo-24*mixed_three
                    if (gap==0)!=(sx==0):
                        continue
                    if not big_z and sx+count_three in ((12,28) if z==32 else (18,66)):
                        population=sx+count_three
                        prefix=three_layout_prefix(z,population,last)
                        num=(mixed_three*225)
                        if not any(total==count_three and at45==three45 and 
                                   num*weight.denominator==weight.numerator*scale
                                   for total,at45,weight in prefix):
                            continue
                    if not big_z and (sx+count_three==(12 if z==32 else 18)):
                        n45=sx+count_three
                        if ((last<45 and count_three!=0) or
                            (last>=45 and three45!=n45)):
                            continue


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
        bound=((scale+root_floor-1)//root_floor)
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
        # In a primitive third-root test on a fixed fifth-root phase
        # there are four active pairs and precisely four even coefficients.
        # A nonzero coefficient at B=15 saturates the whole fiber, forcing
        # five solely-active pairs at B. A nonzero even at 45 would occur
        # twice with scaled cost at least 2*6*5/9>6.
        if z==32 and big_z==0 and last>=15:
            feasible_fibers=False
            for sat in range(3):
                booked=20*sat*(scale//225)
                if booked>mixed_five or small_five<20*sat:continue
                j=(32-16*sat)*(scale//225)
                if j==0:
                    if (small_five==40 and mixed_five==booked and
                        any(sy==0 for _,sy,_,_ in states)):
                        feasible_fibers=True
                    continue
                odd_max=max(sy+count_five-20*sat-(32-16*sat)
                            for _,sy,_,_ in states)
                if odd_max<0:continue
                zero=(small_five-20*sat)*(scale//225)
                zero+=max(0,five45-odd_max)*(scale//2025)
                upper=mixed_five+min(max_y*next_mix,future_five)-booked
                if zero_moments(3,upper,zero,j,j):
                    feasible_fibers=True
                    break
            if not feasible_fibers:return False
        # Couple the phase windows with the entire remaining ordinary mass.
        # Pure activity uses only positions outside the exact mixed quotas.
        unassigned=max(0,k-x-y-outside)
        pure_three=unassigned*upper_units(scale,next_single(last+1,3))
        if z in (32,48):
            balances={(d,e) for d,e in balances if 2*d<=pure_three}
        else:
            room3=9*(scale-mass)+9*mixed_joint-16*mixed_three+16*pure_three
            balances={(d,e) for d,e in balances if 50*d<=room3}
        if not balances:
            return False
        # Average the real-part Gram inequality after all already chosen
        # levels. The remaining mixed three-only positions have floor 45.
        demand=max(x,hx)
        tail=(neutral_mass(L,scale,True)+upper_units(scale,root_floor)
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
                and ordinary_tail.possible(scale//abs(U),support,marked,L,nines,k+extra_neutral,neutral-extra_neutral,o3,o5,
                                           mass,scale,witnesses,qfirst,qexact)):
            return True
        base_book=mass+neutral_mass(L,scale)+upper_units(scale,root_floor)
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
        if z in (32,48) and not zleft:
            diff=mixed_joint-24*mixed_three
            # Cancel all remaining solely-three mass and charge symbolically.
            # Common denominator 48 clears its fixed 3/2 or 1 charge rate.
            rest=48*U-(3 if z==32 else 2)*diff
            if rest:
                cap=min(cap,next_level_cap(rest,48*scale,k-x,factor))

        if residual_bounds and all(v for v,k0 in residual_bounds):
            cap=min(cap,max(next_level_cap(v,4*scale,k0,factor) for v,k0 in residual_bounds))

        candidates=levels(1<<cap.bit_length(),isolated)
        for b in candidates[bisect_right(candidates,last):bisect_right(candidates,cap)]:
            b2=b*b
            if 2*(U-bound)*b2>3*k*scale:
                break
            prime=next((p for p in isolated if b%p==0),None)
            if seven is not None and prime==7 and b<=343:
                continue
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
            if z in (32,48) and not zleft:
                # All remaining mixed three slots together have exact mass;
                # their next-level upper mass may always be lowered to it.
                fixed_remaining=(mixed_joint-24*mixed_three)*up
                limited=(fixed_remaining+23)//24
                three_mass=min(three_mass,limited)
                ordinary=three_mass+five_mass+isolated_mass
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
            cap_charge=U*up+((scale1+minimum_root(multiple)-1)//minimum_root(multiple))
            if prime:
                choices=((0,0,0,q,8*prime*(q//(prime-1)),((0,0),))
                         for q in range(2*(prime-1),min(k,qmax[isolated.index(prime)])+1,2*(prime-1)))
            else:
                joint_choices=range(0,min(k,zleft)+1,packet) if b%15==0 and b>15 else (0,)
                if z in (32,48) and b%15==0 and b>15:
                    joint_choices=(0,zleft) if zleft and zleft<=k else (0,)
                if z==64 and big_z and zleft==64 and b>15 and b%15==0:
                    after=next_mixed(b+1)
                    joint_choices=[j for j in joint_choices if (j<=32 or j==64)
                                   and (after**2<=2*b*b or j in (0,32,64))]
                if z==64 and lock_joint and b>15 and b%15==0:
                    joint_choices=[j for j in joint_choices if j in (0,zleft)]
                choices=((a,f,j,0,8*c+joint_numerator*j,pairs) for j in joint_choices
                         for a,f,c,pairs in single_packets(gcd(b,15),min(k-j,max_sum) if b%15==0 else k-j,joint_order)
                         if j+a+f and (b%15 or (a<=max_x and f<=max_y)))
            for a,f,j,q,ncharge,packets in choices:
                if z in (32,48) and not big_z and a and b%15==0 and b!=45:
                    if all(sx+count_three==(12 if z==32 else 18) for sx,_,_,_ in states):
                        continue
                if mass*up+(a+f+j+q)*unit>=scale1:
                    continue
                charge=ncharge*(unit//8)
                if charge>cap_charge:
                    continue
                next_three=mixed_three*up+(a*unit if b%15==0 else 0)
                next_five=mixed_five*up+(f*unit if b%15==0 else 0)
                next_joint=mixed_joint*up+j*unit
                joint_upper=next_joint+(zleft-j)*unit
                if z in (32,48) and (24*next_three>joint_upper or (b==15 and a)):
                    continue
                if not big_z and z in (32,48):
                    rem=next_joint-24*next_three
                    if rem<0:continue
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
                    if mixed and z==32:
                        i=0  # the complete mixed third-root imbalance vanishes
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
                        cosine_square*up+(active-first)*(unit//2),next_three,next_five,next_joint,small_three+(a if b==15 else 0),small_five+(f if b==15 else 0),count_three+(a if mixed else 0),count_five+(f if mixed else 0),three45+(a if b==45 else 0),five45+(f if b==45 else 0),lock_joint or (z==64 and zleft==64 and j==32),three75+(a if b==75 else 0),five75+(f if b==75 else 0),three135+(a if b==135 else 0),five135+(f if b==135 else 0)):
                    return True
        return False

    tail=nine_gram.ceiling(rank,support,nines,shape,marked,witnesses)
    for c in range(nines//2+1):
        balances={(i*(root_scale//81),0) for i in range(c%2,c+1,2)
                  if nine_gram.possible(nines,c,i,tail)}
        if balances and visit(9,slots,initial_charge-c*(root_scale//27),initial_mass,root,
                              scale_root,states,big_z,tuple(n-prebooked for n in initial_counts),qcaps,initial_first,omit3,w5,balances,0,
                              (nines-3*c)*(root_scale//81),0,0,0,small_joint*(root_scale//225),0,0,0,0,0,0):
            return True
    return False
