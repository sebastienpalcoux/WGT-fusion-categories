#!/usr/bin/env python3
"""Exact, bounded checks for the referee review of three_prime_dimensions.tex.

Run: python verify_three_prime_review.py [--output review_checks.json]
Requires Python 3.10+ and SymPy. No network, private files, or category database.
These checks are NOT a proof or categorification of modular data.
"""
import argparse
from collections import Counter
from itertools import combinations, combinations_with_replacement
import json
import math
from pathlib import Path
import sympy as sp

x = sp.Symbol('x')

def cyclotomic_vectors(level, modulus=None):
    polynomial = sp.Poly(sp.cyclotomic_poly(level, x), x, modulus=modulus)
    if modulus:
        polynomial = polynomial.sqf_part()
    degree = polynomial.degree()
    result = []
    for exponent in range(level):
        remainder = sp.rem(sp.Poly(x**exponent, x, modulus=modulus), polynomial)
        result.append(tuple(int(remainder.nth(i)) % modulus if modulus
                            else int(remainder.nth(i)) for i in range(degree)))
    return result

def sum_vectors(vectors, indices):
    return tuple(sum(vectors[j][i] for j in indices) for i in range(len(vectors[0])))

def word_search(level, mass, residues, allowed):
    """Positive words on mu_level also encode signed words at half-even level.

    Squarefree parts of Phi_level mod p test vanishing at EVERY prime above p.
    Equality is checked in the exact cyclotomic power basis over Z.
    """
    vectors = cyclotomic_vectors(level)
    modular = {p: cyclotomic_vectors(level, p) for p in residues}
    hits = set()
    tested = 0
    for length in range(mass+1):
        for indices in combinations_with_replacement(range(level), length):
            tested += 1
            good = True
            for p, residue in residues.items():
                vector = sum_vectors(modular[p], indices)
                if any((value-(residue if i == 0 else 0)) % p
                       for i, value in enumerate(vector)):
                    good = False
                    break
            if good:
                hits.add(sum_vectors(vectors, indices))
    assert hits <= allowed, (level, mass, 'unexpected word', hits-allowed)
    return dict(level=level, maximum_mass=mass, residues=residues,
                words_tested=tested, distinct_values=len(hits),
                permitted_values=len(allowed), all_values_classified=True,
                all_permitted_values_attained=(hits == allowed))

def scalar(h):
    return next(k for k in range(1, h+1) if k*(k+1) >= h)

def support_checks():
    histogram, passes, survivors = Counter(), Counter(), []
    for triple in combinations(range(1, 23), 3):
        support = {0, *triple}
        differences = {(a-b) % 23 for a in support for b in support}
        h = (len(differences)-1)//2
        histogram[h] += 1
        intersections = [differences & {(b*a) % 23 for a in differences}
                         for b in range(1, 23)]
        if min((len(v)-1)//2 for v in intersections) < 2:
            continue
        passes[h] += 1
        if all(len({(a-b) % 23 for a in intersection for b in intersection})
               >= len(differences) for intersection in intersections):
            survivors.append(triple)
    dilation_classes = sorted({min(tuple(sorted((b*a) % 23 for a in triple))
                                   for b in range(1, 23)) for triple in survivors})
    def affine_canonical(triple):
        support = {0, *triple}
        return min(tuple(sorted((b*(a-origin)) % 23 for a in support))
                   for origin in support for b in range(1, 23))
    affine_classes = sorted({affine_canonical(t) for t in survivors})
    assert histogram == {3:44, 4:264, 5:616, 6:616}
    assert len(survivors) == 176 and len(dilation_classes) == 8
    assert len(affine_classes) == 2
    A = {1,3,4,5,9}
    overlaps = [len(A & {(a+t) % 11 for a in A}) for t in range(11)]
    assert overlaps == [5]+[2]*10
    return dict(level=23, four_position_supports_containing_zero=1540,
                difference_pair_histogram=dict(sorted(histogram.items())),
                supports_passing_two_pair_overlap=dict(sorted(passes.items())),
                supports_passing_one_difference_descent=len(survivors),
                dilation_class_representatives=dilation_classes,
                unweighted_affine_class_representatives=affine_classes,
                eleven_mask_overlaps=overlaps,
                scope='Prime level only; one necessary descent test; no category existence claim.')

def modular_sanity():
    # Standard integral data of Rep(D(S_3)); these entries are rational and real.
    S = sp.Matrix([[1,1,2,3,3,2,2,2], [1,1,2,-3,-3,2,2,2],
                   [2,2,4,0,0,-2,-2,-2], [3,-3,0,3,-3,0,0,0],
                   [3,-3,0,-3,3,0,0,0], [2,2,-2,0,0,4,-2,-2],
                   [2,2,-2,0,0,-2,-2,4], [2,2,-2,0,0,-2,4,-2]])
    dimensions = [1,1,2,3,3,2,2,2]
    assert S*S.T == 36*sp.eye(8)
    omega=(-1+sp.sqrt(-3))/2
    T=sp.diag(1,1,1,1,-1,1,omega,sp.conjugate(omega))
    assert (S*T*S-6*T.inv()*S*T.inv()).applyfunc(sp.simplify) == sp.zeros(8)
    characters = S*sp.diag(*[sp.Rational(1,d) for d in dimensions])
    assert all(a.is_Integer for a in characters)
    fusion = []
    for a in range(8):
        for b in range(8):
            for c in range(8):
                n = sum(S[a,j]*S[b,j]*S[c,j]/(36*dimensions[j]) for j in range(8))
                assert n.is_Integer and n >= 0
                fusion.append(int(n))
    # Primary-twist equality classes: twists are (1,1,1,1,-1,1,w,w^2).
    phases = {2:[0,0,0,0,1,0,0,0], 3:[0,0,0,0,0,0,1,2]}
    blocks = {}
    for p in (2,3):
        grouped = {}
        for j in range(8):
            signature = tuple(int(characters[i,j]) % p for i in range(8))
            grouped.setdefault(signature, []).append(j)
        blocks[p] = list(grouped.values())
        for block in blocks[p]:
            for i in range(8):
                for j in range(8):
                    if phases[p][i] != phases[p][j]:
                        assert sum(S[z,i]*S[z,j] for z in block) == 0
    return dict(category='Rep(D(S_3))', rank=8, dimension=36,
                fusion_coefficients_verified=len(fusion), residue_blocks=blocks,
                block_orthogonality=True, modular_STS_relation=True,
                caveat='A sanity check for unrestricted preliminary identities, not a simple perfect example.')

def main():
    result = {}
    result['dimensions'] = dict(even=2**6*3**6*5**4, odd=5**6*7**6*23**4,
                                odd_square_root=5**3*7**3*23**2)
    assert result['dimensions'] == {'even':29160000, 'odd':514422090765625,
                                     'odd_square_root':22680875}
    result['scalar_endpoint'] = [dict(h=h, u=scalar(h), left=10*scalar(h), right=h*(h-1),
                                     passes=10*scalar(h)<=h*(h-1)) for h in range(1,7)]
    W=sp.Matrix([[0,-1],[1,0]]);U=sp.Matrix([[1,1],[0,1]]);P=sp.Matrix([[0,1],[1,0]])
    h=U*W*U*W.inv();k=U*W*U.inv()*W.inv();R=P*(sp.eye(2)+2*h)
    assert h**3 == -sp.eye(2)
    assert k**3 % 8 == 5*sp.eye(2)
    assert (W*k*W.inv()-k.inv()) % 8 == sp.zeros(2)
    assert (R*h*R.inv()-h.inv()) == sp.zeros(2)
    assert int(R.det()) % 8 == 1
    result['dyadic_matrices'] = dict(h=h.tolist(),k=k.tolist(),h_cubed_mod8=(h**3%8).tolist(),
                                     k_cubed_mod8=(k**3%8).tolist(),conjugator_determinant=int(R.det()))
    assert all(L%4 == 0 for L in range(1,33) if (32-2*L)**2%32 == 0)
    assert all(L%4 == 0 for L in range(1,65) if (64-2*L)**2%64 == 0)
    assert all((u+v)%4 == 0 for u in range(16) for v in range(16) if (u*u+v*v)%8 == 0)
    result['dyadic_congruences'] = 'All residue classes tested; divisibility assertions hold.'
    primes=list(sp.primerange(3,200))
    checked=0
    for p,q,r in combinations(primes,3):
        m=min(r,p*p)//q
        if m:
            assert m<=p-2 and p*(m*(m+1)+1)<(p-1)*r
            checked+=1
        if 2*q<p*p:
            j=q//p
            assert 3*(1+j*(j+1))<2*q
    result['bounded_prime_inequalities'] = dict(prime_bound=200,exposure_instances=checked)
    result['supports'] = support_checks()
    result['modular_sanity'] = modular_sanity()
    phi=(1+sp.sqrt(5))/2;phi_other=(1-sp.sqrt(5))/2
    A5=sp.Matrix([[1,1,1,1,1],[3,-1,0,phi,phi_other],
                  [3,-1,0,phi_other,phi],[4,0,1,-1,-1],[5,1,-1,0,0]])
    assert (A5*sp.diag(1,15,20,12,12)*A5.T/60).applyfunc(sp.simplify)==sp.eye(5)
    result['hypothesis_sanity'] = dict(category='Rep(A_5) with symmetric braiding',
        degrees=[1,3,3,4,5],dimension=60,character_orthogonality=True,
        warning='Not modular. It only disproves an extension obtained by deleting modularity.')
    words=[]
    for level,r in ((30,5),(42,7)):
        v=cyclotomic_vectors(level)
        allowed={tuple(4*a for a in v[0])}
        for t in range(1,r):
            allowed.add(tuple(a+3*b for a,b in zip(v[0],v[(level//r)*t])))
        if r==5:
            # The exceptional signed branch: two 1's, one omega, one omega^2.
            for a in range(4):
                for b in range(4):
                    if a==b: continue
                    exponents=[(level//r)*(t+1)+(level//3)*(1 if t==a else 2 if t==b else 0)
                               for t in range(4)]
                    allowed.add(tuple(-z for z in sum_vectors(v,[e%level for e in exponents])))
        words.append(word_search(level,4,{3:4,r:4},allowed))
    v=cyclotomic_vectors(30)
    allowed={sum_vectors(v,[6*a,6*b]) for a,b in combinations_with_replacement(range(5),2)}
    words.append(word_search(30,2,{5:2},allowed))
    v=cyclotomic_vectors(60);allowed=set()
    for dyadic in range(4):
        shift=15*dyadic
        allowed.add(tuple(3*a for a in v[shift]))
        allowed.add(tuple(a-b for a,b in zip(v[(shift+20)%60],v[(shift+40)%60])))
    words.append(word_search(60,3,{2:1,3:0},allowed))
    result['bounded_cyclotomic_word_checks']=words
    result['limitations']=['No classification or existence of a fusion category is computed.',
                          'No arbitrary-conductor theorem follows from the finite searches.',
                          'The categorical indicator and fermion-sector proofs require mathematical verification.']
    args=argparse.ArgumentParser();args.add_argument('--output',default='review_checks.json');opts=args.parse_args()
    Path(opts.output).write_text(json.dumps(result,indent=2,default=int)+'\n')
    print(json.dumps(result,indent=2,default=int))

if __name__=='__main__':main()
