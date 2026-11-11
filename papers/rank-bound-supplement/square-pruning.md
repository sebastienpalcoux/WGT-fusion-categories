# Square-dimension exclusions for the rank-31 computation

The references below are to *A rank bound for simple integral modular categories*. Write the global dimension as \(D=M^2\), the degree of a simple label as \(d_Z\), its denominator as \(b_Z=M/d_Z\), and its formal codegree as \(f_Z=b_Z^2\). A set of labels has normalized mass \(\sum_Z d_Z^2/D=\sum_Z b_Z^{-2}\). The vacuum has denominator \(M\) and is included in all complements. Write \(s=S/M\), \(T=\operatorname{diag}(\theta_Z)\), and let \(\xi\) be the multiplicative central charge. The vacuum modular relation is \(sT\mathbf d=\xi T^{-1}\mathbf d\). Consequently, any set of rows supported on a column set has mass at most the mass of those columns, by unitarity (equation (4.1)).

## Three-prime heights and packet counts

Lemma 1.2 supplies two nonunit degrees prime to any given prime \(p\mid M\), and every nonunit degree has at least two distinct prime factors (Lemma 1.1). Thus at least three primes divide \(M\). With precisely three primes, the sets of nonunit labels prime to the different primes are pairwise disjoint. Let \(u_p\) be the number of nonunit labels of degree prime to \(p\). Reducing \(D=1+\sum_{Z\ne0}d_Z^2\) modulo four, three and five respectively gives \(u_2\equiv3\pmod4\), \(u_3\equiv2\pmod3\) and \(u_5\ge3\). For the last inequality, one or two residues from the nonzero squares \(1,4\) cannot sum to \(-1\pmod5\) with two terms; the one-term possibility is excluded by Lemma 1.2.

If four or more denominators are four, the manuscript's “Prime powers in a prefix” argument shows their rows are supported on multiples of four and their nonunit hit degrees have gcd one. In three-prime support this forces \(8\mid M\). Similarly if a denominator is five, Lemma 3.3 shows that the nonunit degrees with nonzero absolute field trace in that row have gcd one and divide \((M/5)/\gcd(M/5,5)\). Their common multiple must contain all three primes. Otherwise, since every degree has at least two prime factors, all those degrees share a prime. Hence \(125\mid M\). For each \(p^a\Vert M\) at least one four-row hit degree is prime to \(p\), and at least one five-row trace degree is prime to \(p\) whenever the respective hypotheses apply. Their denominators are multiples of \(\operatorname{lcm}(4,p^a)\) and \(\operatorname{lcm}(25,p^a)\), respectively. The two witnesses may be different.

In an ordered prefix of denominators write \(L\) for the lcm, \(E=L^2(1-\sum_{\rm prefix}b^{-2})\), and \(h\) for the remaining nonunit slots. At an odd prime \(p\mid L\), the active labels occur in packets of \(j=(p-1)/2\) equal-degree labels with \(p\mid b\). A closed eligible run of multiplicity \(m\) can therefore contribute at most \(j\lfloor m/j\rfloor\) active positions. Count every position in the eligible open run, since it might grow. Call the sum of these upper bounds \(C_p\). For the subset of known primes with positive deficit, Proposition 2.4 gives a union of at least \(\sum(p-1)\) active positions, while at most \(\sum C_p\) of them are in the prefix. Hence \(H=\sum_{p\mid L,\,p\text{ odd}}\max(0,p-1-C_p)\le h\). Adding any unseen odd support prime to that subset shows \(p-1\le h-H\), because it has no active prefix label. This bounds the available prime set. The argument does not assume that the active sets themselves are disjoint.

If precisely three primes are available, all occur. Write \(M=kL\). For each prime \(p\) consider whether \(p\mid k\), giving eight cases in total. If so, no known degree is \(p\)-free and its full denominator height is at least \(pL_p\), where \(L_p=p^{v_p(L)}\). Otherwise the known \(p\)-free positions are exactly the prefix denominators divisible by \(L_p\), and the height is exactly \(L_p\). Impose the lower heights just proved, discarding a case if its exact height is too small. Known free sets must be disjoint. For each prime reserve the least number of suffix slots satisfying the lower bound and the indicated residue for \(u_p\), and at least one slot when a four-row or five-row witness for that prime is missing. If the total is \(v>h\), discard the case. Otherwise each reserved slot needs a denominator multiple of the height bound. At least one slot per prime must satisfy the largest applicable bound amongst the height and its lcm with 4 or 25. Round each such bound up to its least multiple at least the last prefix denominator; denote the resulting lower bounds by \(a_i\). For the strongest bound this rounding is performed *after* the maximum is taken, which remains valid for the witness supplying that maximum.

The reserved mass is at most \(W/L^2\), where \(W=\sum_{i=1}^v\lceil L^2/a_i^2\rceil\). The vacuum contributes at most \(1/L^2\). For the next denominator \(b\), each of the other \(h-v\) terms contributes at most \(1/b^2\). Thus if \(E>W+1\) one may impose \(b^2\le (h-v)L^2/(E-W-1)\). The search takes the largest bound among feasible cases. When two slots remain, \(k\) is prime to ten by the two deletion witnesses and the sum of squares modulo four and five. When three remain and \(k\) is even, all three remaining degrees must be odd, giving \(Ek^2\equiv4\pmod8\), so \(v_2(k)=1\) and \(E\) is odd.

## Mixed classes at denominators four and five

Write \(I=\{b=4\}\), \(J=\{b=5\}\) with multiplicities \(i,j\), and let \(A,B\) be the respective sums of twists. Their twists are dyadic and quinary roots, respectively, by Lemma 2.2. If their complement has mass \(t/400\), the Gauss sum, of absolute value \(M\) at every embedding, gives

\[
 |25\sigma A+16\tau B|\le t+400/M .
\]

Here \(\sigma,\tau\) can vary independently. The variance of the conjugates of the dyadic integer \(A\) is an integer. Indeed the integral power basis of a two-power cyclotomic field is orthonormal for the mean squared norm, and only its constant basis vector has nonzero mean. A disk of radius less than one for \(\sigma A\) therefore forces \(A\in\mathbb Z\). Square Galois covariance puts \(B\) in the quadratic subfield \(\mathbb Q(\sqrt5)\), since the multiset of twists on its degree class is invariant under all squares. Distinct conjugates of a quadratic integer in this field differ by at least \(\sqrt5\). For rational \(A,B\), reduction at the respective primes and the mean of a five-power root give \(A\equiv i\pmod2\), \(B\equiv j\pmod5\) and \(-j/4\le B\le j\).

### Twelve fours and six fives, or ten fours and nine fives

The two configurations leave respectively \(t=4,6\); these multiplicities are exact by mass. Here \(20\mid M\) and another prime divides \(M\). If \(M\le60\), one has \(M=60\) and the complement after removing the vacuum has squared mass 35 or 53. All nonunit degrees are at least six, so neither is a sum of their squares.

For \(M>60\) the displayed disk has radius less than 13. Hence \(A\) is rational and two conjugates of \(B\) differ by less than \(26/16<\sqrt5\), so \(B\) is rational. For \((12,6)\), the allowed \(B\) are 1 and 6; for \((10,9)\), they are \(-1,4,9\). In both cases \(A\) is even, and the only allowed absolute value of \(25A+16B\) is \(t\). Their contribution is therefore \(\varepsilon tD/400\), \(\varepsilon=\pm1\). Taking the real part of \(\varepsilon\) times the Gauss sum and averaging conjugates yields

\[
 \sum_{Z\notin I\cup J}d_Z^2\bigl(1+\operatorname{mean}(\varepsilon\theta_Z)\bigr)\le M .
\]

For a prime \(p\mid M\) different from 2 and 5, all \(p\)-active labels are in this complement. At such a label \(-\varepsilon\theta_Z\ne1\), and a nonidentity root of unity has mean at most \(1/2\). Thus their total squared mass \(V_p\le2M\). Lemma 2.5 gives \(2pgV_p\ge (p-1)D\) for the least nonunit degree \(g\), hence \(g\ge M/6\). Every nonunit denominator now divides 60 since it is an integer at most six. Multiplying the identity \(1-M^{-2}=\sum_{Z\ne0}b_Z^{-2}\) by \(60^2\) forces \(M\mid60\), a contradiction.

### Thirteen fours and four fives at rank 31

The complementary mass is \(\delta=11/400\). In the case of three support primes the preceding height argument gives \(125\mid M\), so \(M\ge1500\). With four or more support primes \(M\ge420\), and \(M\ge2100\) if \(25\mid M\). Thus \(M\ge420\) throughout, and \(M>1000\) whenever \(25\mid M\).

The disk now has radius less than 12, so both twist sums are rational as above. All four quinary twists have order dividing five: a higher phase would require a square orbit of at least ten labels. The conditions \(A\) odd, \(-1\le B\le4\), \(B\equiv4\pmod5\) leave just \((1,-1)\) and \((-3,4)\).

For \((-3,4)\) the contribution is \(-\delta D\), so the argument of the preceding exclusion gives \(g\ge M/6\). There are thirteen remaining nonunit labels at rank 31, all with denominator at least seven because \(\delta<1/36\). This is a contradiction.

For \((1,-1)\) the four twists in \(J\) are the four primitive fifth roots, since \(1+B=0\) is a vanishing sum of five fifth roots. Let \(H\) consist of the complement labels with five-power twist, including the vacuum, and write \(\eta\) for its mass. Now the two classes contribute \(9D/400\). Taking the real part of the Gauss sum and averaging gives

\[
 \sum_{Z\notin I\cup J}d_Z^2(1+\operatorname{mean}(\theta_Z))
 \le D/200+M,\qquad \eta\le \tfrac43(1/200+1/M),
\]

since a five-power root has mean at least \(-1/4\). Each \(J\)-row is supported on the principal five-block \(P\), sees only columns of five-power twist, and misses \(I\), whose codegree is prime to five (Lemmas 1.4, 1.5 and 2.2). If \(X\in J\setminus P\), symmetry excludes every \(J\)-column from its support. The mass bound (4.1) would give \(1/25\le\delta\), impossible. Hence \(J\subset P\).

Put \(y=M/5\), \(\zeta=\zeta_5\), and write \(\chi_X(Z)=S_{XZ}/y\). By Lemma 1.3 the character field of \(X\in J\) lies in \(\mathbb Q(\zeta_{25})\); its degree is at most four since the orbit lies in \(J\). Cyclicity of the Galois group then puts it in \(\mathbb Q(\zeta)\). In the modular relation for that row all supported twists have five-power order. Hence the two-power root \(\xi\) belongs to a five-power cyclotomic field, giving \(\xi=\pm1\). Define

\[
 R_X=5\xi\theta_X^{-1}-\sum_{Z\in J}\theta_Z\chi_X(Z)
     ={25\over M}\sum_{Z\in H}s_{XZ}\theta_Zd_Z\in\mathbb Z[\zeta].
\]

Then \(\sum_{X\in J}\operatorname{mean}(|R_X|^2)\le625\eta\), by unitarity of \(s\) in the last expression and after any conjugation (Galois permutes the rows of \(s\) in square dimension). For local five-inertia \(\sigma\) one also has \(R_{\sigma X}=\sigma^2R_X\). Indeed replace \(Z\) by \(\sigma Z\) in the sum over \(J\), using twist covariance and Galois symmetry in each of the two character indices.

If \(5\mid y\), the principal-block congruence \(\chi_X(Z)\equiv y\pmod{1-\zeta}\) for \(Z\in J\) puts every \(R_X\) in \((1-\zeta)\mathbb Z[\zeta]\). Each nonzero element there has mean squared norm at least \(5/2\). To check this, present it by five integer coefficients on the fifth roots. Their sum is divisible by five, so subtract a common integer to arrange coefficient sum zero. Now the mean squared norm equals \(5/4\) times the sum of the coefficient squares, which is at least two. Covariance and the four primitive twists force at least two nonzero values once any is nonzero. Yet \(625\eta\le (2500/3)(1/200+1/M)<5\) since \(M>1000\). If all values vanish, project \(sT(\mathbf d 1_J)\) onto \(J\). Its entries there are \(\xi\theta_X^{-1}y\), so this projection already has the norm of its input. Its vacuum coordinate must vanish, although it equals \(y^2B/M\ne0\).

If \(5\nmid y\), each \(J\)-row is principal and meets every five-free degree by its character congruence. Hence all five-free degrees have five-power twist. Every five-active label has a five-free degree, since activation requires \(5\mid b_Z\) and \(v_5(M)=1\). On any nontrivial five-phase eigenspace, every other local upper unipotent is therefore scalar. Its local group preserves that space and is generated by conjugates of the upper unipotent, so its Weyl fixes each coordinate line there. By Lemma 2.2 the codegree of any five-active label has no other prime factor and must be 25. Thus all twists on \(H\) are one. The last bound improves to \(\eta\le(1/200+1/M)/2\). On the other hand \(R_X\equiv -4y\ne0\pmod{1-\zeta}\). Each of the four nonzero integers has mean squared norm at least one by its algebraic norm, whereas \(625\eta\le (625/2)(1/200+1/420)<4\).

### Thirteen fours, three fives and a six

More generally, suppose that \(j\) equals one or three and the mass of \(R=\{Z:5\mid b_Z,\ b_Z\ne5\}\) is less than \(1/25\). This includes the vacuum, because \(5\mid M\) and \(M>5\). These conditions are impossible.

If \(25\mid M\), the \(j\) character rows in \(J\) form Galois orbits of size dividing \([\mathbb Q(\zeta_{25}):\mathbb Q]=20\) and at most three. As \(j\) is odd, at least one row is rational. Each nonzero entry in that row has nonzero absolute field trace, so its nonunit support degrees divide \(M/25\) by Lemma 3.3. The same divisibility holds for the vacuum. Thus the whole support is in \(R\), contradicting (4.1).

If \(25\nmid M\), all \(J\)-rows are supported on the principal five-block, whose labels satisfy \(5\mid b\) (Lemmas 1.4–1.5). Any nonprincipal \(J\)-row misses all its own degree class by symmetry, leaving it supported on \(R\). This again contradicts (4.1), so all \(J\)-rows are principal. The local argument in the preceding paragraph now places every five-active label in \(J\); Proposition 2.4 requires at least four such labels.

In the stated prefix \(13/16+3/25+1/36=3457/3600\); these labels are disjoint from \(R\). Its mass, including the vacuum, is therefore at most \(143/3600<1/25\). A fourth five cannot fit. The generator uses only this fixed prefix as a constant-time test.
