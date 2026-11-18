#!/usr/bin/env python3
"""Independent raw-minor audit of the orthogonal f=6 catalogue.

Does not import the producer.  Enumerates every bounded multiset without a
determinant pruning rule, computes all Pluecker coordinates by the Leibniz
formula, and contracts the alternating three-form directly.  The positivity
check uses separating lines in the fixed plane, rather than kernel circuits.
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
P3=list(permutations(range(3)))
SGN={p:(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)) for p in P3}

def determinant(columns):
    return sum(SGN[p]*columns[0][p[0]]*columns[1][p[1]]*columns[2][p[2]] for p in P3)

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def involutions(n):
    # Enumerate disjoint transposition sets; the unused labels are fixed.
    edges=list(combinations(range(n),2))
    def visit(start,used,p):
        yield tuple(p)
        for h in range(start,len(edges)):
            i,j=edges[h]
            if i in used or j in used:continue
            p[i],p[j]=j,i
            yield from visit(h+1,used|{i,j},p)
            p[i],p[j]=i,j
    yield from visit(0,set(),list(range(n)))

def positive(columns,duality):
    v=[tuple(columns[i][k]+(columns[j][k] if i!=j else 0) for k in range(3))
       for i,j in enumerate(duality) if i<=j]
    normal=next((cross(a,b) for a,b in combinations(v,2) if any(cross(a,b))),None)
    assert normal is not None,'fixed pure space must have rank two'
    assert all(dot(normal,a)==0 for a in v)
    for a in v:
        if not any(a):continue
        side=cross(normal,a);values=[dot(side,b) for b in v]
        if (all(x>=0 for x in values) or all(x<=0 for x in values)) and any(values):return False
    return True

def signature(g,d,s,signed=True):
    out=[]
    for i,j in enumerate(d):
        if i>j:continue
        t=s[i] if signed else abs(s[i])
        if i==j:out.append((True,g[i][i],0,t))
        else:
            assert (g[i][i]+g[i][j])%2==0
            out.append((False,(g[i][i]+g[i][j])//2,(g[i][i]-g[i][j])//2,t))
    return tuple(sorted(out))

def key(columns,d,s):return (tuple(map(tuple,columns)),tuple(d),tuple(s))

def run():
    frozen=json.loads((HERE/'orthogonal_six_catalogue.json').read_text())
    states={key(r['columns'],r['duality'],r['traces']):r for r in frozen['states']}
    assert len(states)==len(frozen['states'])
    vec=[]
    for a,b,c in product(range(-3,4),repeat=3):
        v=(a,b,c)
        if not 0<a*a+b*b+c*c<=11:continue
        if next(x for x in v if x)>0:vec.append(v)
    vec.sort();counts=Counter(oriented_columns=len(vec));frames=[];found=set()
    families=Counter();signed=set();positive_signed=set();opposites=Counter()
    # Raw bounded multiset inventory: unlike the producer recursion, every
    # multiset of the given size is visited and then its norm is tested.
    for extra_count in range(5):
        for ids in combinations_with_replacement(range(len(vec)),extra_count):
            extra=tuple(vec[i] for i in ids)
            if sum(dot(a,a) for a in extra)>11:continue
            counts['normalized_inputs']+=1
            columns=((1,0,0),(0,1,0),(0,0,1))+extra;n=len(columns)
            coordinates={abc:determinant([columns[i] for i in abc]) for abc in combinations(range(n),3)}
            if sum(x*x for x in coordinates.values())!=12:continue
            support=[abc for abc,x in coordinates.items() if x]
            if any(all(i in abc for abc in support) for i in range(n)):continue
            frames.append(columns)
            for signs in product((-1,1),repeat=extra_count):
                signed_columns=columns[:3]+tuple(tuple(sign*a for a in v) for v,sign in zip(extra,signs))
                def alternating(i,j,k):
                    if len({i,j,k})<3:return 0
                    p=(i,j,k);sign=(-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3))
                    ids=tuple(sorted(p));value=coordinates[ids]
                    for h in ids:
                        if h>=3:value*=signs[h-3]
                    return sign*value
                g=tuple(tuple(sum(alternating(i,a,b)*alternating(j,a,b) for a,b in combinations(range(n),2))
                              for j in range(n)) for i in range(n))
                assert sum(g[i][i] for i in range(n))==36
                for d in involutions(n):
                    if sum(g[i][d[i]] for i in range(n))!=12:continue
                    if any(g[d[i]][d[j]]!=g[i][j] for i in range(n) for j in range(n)):continue
                    for tail in product(range(-2,3),repeat=extra_count):
                        head=tuple(-sum(signed_columns[j+3][i]*tail[j] for j in range(extra_count)) for i in range(3))
                        if any(abs(x)>2 for x in head):continue
                        s=head+tail
                        if any(s[d[i]]!=s[i] for i in range(n)):continue
                        if 8-sum(x*x for x in s) not in (0,4,8):continue
                        if any((s[i]*s[i]-g[i][d[i]])%4 for i in range(n)):continue
                        if any((s[i]*s[i]+g[i][i])%2 or s[i]*s[i]+g[i][i]>=12 for i in range(n)):continue
                        counts['trace_states']+=1
                        modulus=2 if any(x%2 for x in s) else 4
                        if any((g[i][j]+s[i]*s[j])%modulus for i in range(n) for j in range(n)):continue
                        k=key(signed_columns,d,s);assert k in states and k not in found;found.add(k)
                        r=states[k];assert tuple(map(tuple,r['gram']))==g
                        counts['integral_states']+=1;families[r['family']]+=1
                        signed.add(signature(g,d,s))
                        pos=positive(signed_columns,d)
                        assert pos==r['positive']
                        if pos:
                            counts['positive_states']+=1;counts['positive_'+r['family']]+=1
                            positive_signed.add(signature(g,d,s))
                            # Verify the advertised positive rational witness independently.
                            from fractions import Fraction
                            w=list(map(Fraction,r['positive_witness']))
                            assert all(x>0 for x in w) and all(w[i]==w[d[i]] for i in range(n))
                            assert all(sum(g[i][j]*w[j] for j in range(n))==0 for i in range(n))
                            for i,j in combinations(range(n),2):
                                if s[i]+s[j] or g[i][i]+g[j][j]+2*g[i][j]:continue
                                opposites[r['family']]+=1
                                assert not(i==d[i] and j==d[j]),'opposite selfdual active images'
                                assert j!=d[i],'an active image is the negative of its dual'
    assert found==set(states)
    assert sorted(frames)==sorted(tuple(map(tuple,f)) for f in frozen['frames'])
    counts['norm12_frames']=len(frames)
    counts['integral_signed_signatures']=len(signed)
    counts['positive_signed_signatures']=len(positive_signed)
    for name,n in families.items():counts['integral_'+name]=n
    assert dict(counts)==frozen['counts'],(counts,frozen['counts'])
    result={'passed':True,'independent_methods':['Leibniz determinants','direct three-form contractions',
             'raw multiset enumeration','disjoint-transposition involutions','separating-line positive-cone test'],
            'counts':dict(sorted(counts.items())),
            'opposite_images_control':dict(opposites),
            'opposite_images_scope':'No opposite pair among selfdual active labels, and no active X with image(X*)=-image(X); arbitrary distinct nonselfdual labels may be opposite.'}
    (HERE/'orthogonal_six_independent_audit.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
