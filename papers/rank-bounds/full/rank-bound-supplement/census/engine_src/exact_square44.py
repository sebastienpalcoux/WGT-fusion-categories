"""Unbounded rank-44 square continuation with existing necessary prefix rules.

All arithmetic is integral. This module has no search cutoff. The guards
are ports of the already certified square.cpp tests and the completed
rank-44 high-four structural calculation. A bounded external pilot may
interrupt a run but must report it as unresolved.
"""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from math import gcd,isqrt,lcm
from pathlib import Path
import json

from denominators import smooth_interval,smooth_over
from terminal_fast import terminal
from primes import primes_up_to
from square_exact import packet_room
from prime_cut_prefix import prime_cut_prefix_failure

BASE=Path(__file__).resolve().parents[1]

@lru_cache(None)
def high_four_cases():
    data=json.loads((BASE/'arithmetic/high-four44-mixed.json').read_text())
    assert data['rank']==44 and data['remaining']==39
    cases=[]
    for row in data['signatures']:
        support=frozenset((2,*row['support']))
        primary=frozenset(p for p,z in zip(row['support'],row['modes']) if z==2)
        cases.append((row['rows'],support,primary))
    return tuple(cases)

def _part(n,p):
    power=1
    while n%p==0:n//=p;power*=p
    return power

def _mixed_need(bins,phases=2):
    full=sum(n//4 for n in bins)
    if full>=phases:return 0
    left=phases-full;missing=4*left
    residues=Counter(n%4 for n in bins)
    for r in (3,2,1):
        used=min(left,residues[r]);missing-=used*r;left-=used
    return missing

def _separated_three(prefix,counts,last,h):
    # Exact translation of square.cpp separated_three; its optional lcm
    # ceiling is retained, so this port introduces no stronger cut.
    bins=[0]*8;eligible=[];W=1
    for b,n in sorted(counts.items()):
        two=_part(b,2);odd=b//two
        rest=odd
        while rest%3==0:rest//=3
        if two>1 and odd>1 and rest==1:
            bins[min(7,two.bit_length()-2)]+=n
        if b%3==0 and rest==1:
            q=b//gcd(W,b)
            if W>(1<<58)//q:break
            W*=q;eligible.append([b,W,n])
    # The native bins include every prefix mixed-three label, even when
    # the lcm loop above stops. Recompute those initial bins separately.
    bins=[0]*8
    for b,n in counts.items():
        two=_part(b,2);odd=b//two;rest=odd
        while rest%3==0:rest//=3
        if two>1 and odd>1 and rest==1:
            bins[min(7,two.bit_length()-2)]+=n
    begin=0
    for j,(bj,bound,_) in enumerate(eligible):
        if bound>=last:break
        unit=256*bound*bound
        tail=h*((unit+last*last-1)//(last*last))
        for b in prefix:
            if b<=bj or b%3:continue
            rest=b//_part(b,2)
            while rest%3==0:rest//=3
            if rest==1:tail+=(unit+b*b-1)//(b*b)
            if tail>=256:break
        if tail>=256:continue
        changed=True
        while changed:
            changed=False
            for p in (2,3):
                high=0;idx=0;number=0
                for k in range(begin,j+1):
                    b,_,n=eligible[k]
                    if not n:continue
                    part=_part(b,p)
                    if part>high:high=part;idx=k;number=0
                    if part==high:number+=n
                if number==1:
                    b=eligible[idx][0];eligible[idx][2]=0;changed=True
                    if b%2==0:bins[min(7,_part(b,2).bit_length()-2)]-=1
        begin=j+1
    return _mixed_need(bins)

def _width_room(prefix,counts,L,h,support):
    # Native width_room's omission/pure-packet clauses, reapplied after
    # every exact child on a high-four, three-mixed support.
    last=prefix[-1]
    active=0
    two_height=max(8,L&-L)
    have_odds=sum(b%two_height==0 for b in prefix)
    odds_other=sum(b%two_height==0 and b%5!=0 for b in prefix)
    proper5=sum(b>5 and b%5==0 for b in prefix)
    inert=3-min(3,have_odds)
    if counts[5] and last>5:
        if proper5>4:return False
        inert=max(inert,4-proper5+2-min(2,odds_other))
    pure_counts={p:sum(n for b,n in counts.items() if _part(b,p)==b) for p in support-{2,3}}
    if sum(p-min(p,n) for p,n in pure_counts.items())+inert>h:return False
    active=_separated_three(prefix,counts,last,h)
    pure_slots=active;prime_inert=0
    for p in sorted(support-{2,3}):
        pure_count=pure_counts[p]
        pure_slots+=p-min(p,pure_count)
        power_in_L=_part(L,p);vp=0;remaining=power_in_L
        while remaining>1:remaining//=p;vp+=1
        present=bool(counts[p])
        if present and vp>1:return False
        power=p**((vp+2)//2)
        paid=False
        while power<last:
            if counts[power]>=p-1:paid=True;break
            power*=p
        if not paid:
            count=counts[last] if power==last else 0
            active+=p-1-min(p-1,count)
        if not present and last>p:
            eventual=last+1 if power>last else power*power//p
            paid_free=sum(b%eventual==0 for b in prefix) if eventual<=last else 0
            needed=3 if p==5 else 2
            prime_inert=max(prime_inert,needed-min(needed,paid_free))
    return pure_slots+inert<=h and active+max(inert,prime_inert)<=h

NORM_SCALE=1<<48
NORM_ROOT=1<<24

def _upper_mass(b):
    return 1 if b>=NORM_ROOT else (NORM_SCALE+b*b-1)//(b*b)

@lru_cache(maxsize=8192)
def _component_contribution(b,p):
    if b%p:return 0,0
    power=_part(b,p);cofactor=b//power
    if cofactor==1:return _upper_mass(b),0
    a=pow(cofactor,-1,power)
    value=NORM_ROOT*min(a,power-a)//power
    return 0,value*value

def complete(item,task,rank,diagnostics=None):
    if rank!=44 or task['q']!=1:
        raise ValueError('this continuation is restricted to square rank44')
    forced=task.get('prefix',[task['B']])
    prefix=list(item['exact'])
    if (not prefix or prefix[:len(forced)]!=forced or len(prefix)>rank-3 or
        any(type(z) is not int or z<=0 for z in prefix) or prefix!=sorted(prefix)):
        raise ValueError('invalid square arbitrary-precision prefix')
    all_primes=frozenset(primes_up_to(rank-1))
    if item.get('primes') or not all(smooth_over(z,all_primes) for z in prefix):
        raise ValueError('invalid square handoff support')
    L,E=1,1
    for z in prefix:
        g=gcd(L,z);scale=z//g
        L,E=L*scale,E*scale*scale-(L//g)**2
        if E<=0:return
    counts=Counter(prefix)
    component_primes=tuple(p for p in all_primes if p>=5)
    pure_mass={p:0 for p in component_primes};mixed_cost={p:0 for p in component_primes}
    for b in prefix:
        for p in component_primes:
            pure,mixed=_component_contribution(b,p)
            pure_mass[p]+=pure;mixed_cost[p]+=mixed
    info=diagnostics if diagnostics is not None else {}
    def note(key):info[key]=info.get(key,0)+1

    def visit(L,E,h):
        note('nodes')
        last=prefix[-1];fours=counts[4];fives=counts[5]
        # Native component_room, with outward-rounded accumulated masses.
        isolated_five=bool(fives or fours>=9 or (counts[7] and fours+counts[6]+27>rank))
        for p in component_primes:
            if not mixed_cost[p] or (p==5 and not isolated_five):continue
            floor=p
            while floor<last and floor<NORM_ROOT:floor*=p
            pure=pure_mass[p]+h*_upper_mass(floor)
            if (p-1)*mixed_cost[p]>(p+1)*pure:
                note('component_norm');return
        if fours>13 or fives>5:
            note('class_count');return
        if fives and fives!=5 and (last>5 or fives+h<5):
            note('class_count');return
        if fours==12 and fives==5 and (counts[6] or counts[8]>2):
            note('twelve_four');return
        if fours==11 and fives==5 and counts[6]>=4:
            note('eleven_four');return
        seen={p for p in all_primes if L%p==0}
        if len(seen)>4:
            note('support_count');return
        allowed=seen|{2}|{p for p in all_primes if p<=h+1}
        if len(seen)==4:allowed=seen
        if fours>=9:
            present={p for p in all_primes if p>=5 and counts[p]}
            compatible=[]
            for exact_fours,support,primary in high_four_cases():
                if (exact_fours<fours or (last>4 and exact_fours!=fours)
                    or not seen<=support or not present<=primary
                    or any(p<last and not counts[p] for p in primary)):
                    continue
                compatible.append(support)
            if not compatible:
                note('high_four');return
            allowed &= set().union(*compatible)
            common=set.intersection(*(set(s) for s in compatible)) | seen
            if len(common)==4 and not _width_room(prefix,counts,L,h,common):
                note('width_room');return
        if len(allowed)<4 or not packet_room(counts,last,seen,h):
            note('packet_support');return
        if fives and L%25==0:
            note('prime_height');return
        for p in seen:
            if p<7 or not counts[p]:continue
            if L%(p*p)==0 or (counts[p]<p-1 and (last>p or counts[p]+h<p-1)):
                note('prime_height');return

        missing5=0
        if fives==5 and last>5:
            proper=[z for z in prefix if z>5 and z%5==0]
            number=len(proper)
            if number>4 or any(z%25 not in (5,20) for z in proper):
                note('five_shape');return
            if fours==12 and any(z<80 or z%4 for z in proper):
                note('five_shape');return
            missing5=4-number
            if (missing5>h or len(set(proper))+missing5<3
                or (fours>=4 and (16-fours)*(sum(z%4==0 for z in proper)+missing5)<fours)):
                note('five_count');return
            trace=sum((1 if z%25==5 else -1)*(L//z) for z in proper)
            magnitude=abs(trace)
            if missing5==0:
                if magnitude!=1 or lcm(*proper)!=L:
                    note('five_trace');return
            elif 6*magnitude*last>(6*missing5+1)*L:
                note('five_trace');return

        lo=max(last,isqrt(L*L//E)+1)
        hi=isqrt((36*h+1)*L*L//(36*E))
        if lo>hi:
            note('mass');return
        if prime_cut_prefix_failure(prefix,h,tuple(sorted(allowed)),lo,
                Fraction(E,L*L),required_primes=seen,
                vacuum_upper=Fraction(1,max(L*L,36*lo*lo))) is not None:
            note('prime_cut_mass');return
        if h==2:
            note('terminal_calls')
            for z,k,d,e in terminal(L,E,lo,hi,allowed,fixed_support=False):
                M=k*L
                suffix=(z,M//d)
                if fives==5:
                    future=[b for b in suffix if b>5 and b%5==0]
                    if (len(future)!=missing5 or any(b%25 not in (5,20) for b in future)
                        or (fours==12 and any(b<80 or b%4 for b in future))):
                        continue
                    allproper=proper+future
                    if (len(set(allproper))<3 or lcm(*allproper)!=M
                        or abs(sum((1 if b%25==5 else -1)*(M//b) for b in allproper))!=1):
                        continue
                actual=tuple(p for p in all_primes if M%p==0)
                if len(actual)!=4:continue
                row=[1,d,e]+[M//b for b in reversed(prefix)]
                if any(sum(x%p==0 for p in actual)<2 for x in row[1:]):continue
                if any(sum(x%p!=0 for x in row[1:])<(3 if p in (2,5) else 2) for p in actual):continue
                note('rows');yield row
            return

        # If every remaining label must be a proper five-free witness,
        # use its required v5(b)=1 already in the smooth iterator.
        if missing5==h:
            choices=(5*z for z in smooth_interval((lo+4)//5,hi//5,allowed-{5})
                     if z%5 in (1,4) and (fours!=12 or z%4==0))
        else:
            choices=smooth_interval(lo,hi,allowed)
        for z in choices:
            if fives==5 and z%5==0 and z>5 and (z%25 not in (5,20) or
                (fours==12 and (z<80 or z%4))):continue
            g=gcd(L,z);scale=z//g
            prefix.append(z);counts[z]+=1
            contributions=[_component_contribution(z,p) for p in component_primes]
            for p,(pure,mixed) in zip(component_primes,contributions):
                pure_mass[p]+=pure;mixed_cost[p]+=mixed
            yield from visit(L*scale,E*scale*scale-(L//g)**2,h-1)
            for p,(pure,mixed) in zip(component_primes,contributions):
                pure_mass[p]-=pure;mixed_cost[p]-=mixed
            counts[z]-=1
            if not counts[z]:del counts[z]
            prefix.pop()
    yield from visit(L,E,rank-1-len(prefix))
