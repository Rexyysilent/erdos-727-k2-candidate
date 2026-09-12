# A quadratic-family candidate for the k = 2 case of Erdős problem 727

**Date:** 7 September 2026  
**Status:** Candidate proof developed in this conversation. It has been checked algebraically and against exact numerical calculations, but has not been independently reviewed or formally verified. It does **not** resolve the assertion for every fixed k >= 2. The accompanying Python script verifies individual examples, not the asymptotic proof.

## 1. The problem and the proposed result

Erdős problem 727 asks whether, for every fixed integer \(k\geq2\), infinitely many positive integers \(n\) satisfy
\[
 ((n+k)!)^2\mid(2n)!.
\]
The original paper [1, p. 90] states this question, and the Erdős Problems page [2] still lists even the case \(k=2\) as open.

Write
\[
 F(t)=6t^2+11t+3,\qquad t\geq0.
\]
The candidate result is as follows.

**Proposed theorem.** Put
\[
 \rho=1-\frac12\log3
       -4\sum_{J=4}^{\infty}2^{1-J}\log\frac{J+1}{J}.
\]
Then
\[
 \liminf_{T\to\infty}
 \frac1T\#\{t\in[T,2T]\cap\mathbb Z:
             ((F(t)+2)!)^2\mid(2F(t))!\}\geq\rho.
 \tag{1}
\]
Here \(\rho=0.25922348\ldots\). In particular,
\[
 \rho\geq\frac34-\frac12\log3>\frac15>0.
 \tag{2}
\]
Thus, if the argument below withstands independent checking, it proves infinitude in the first open case \(k=2\).

All logarithms in this note are natural. Endpoints in integer intervals contribute only \(O(1)\) and have no effect on the limits.

## 2. The exact divisibility criterion

Let
\[
 B(n)=\binom{2n}{n},\qquad D(n)=(n+1)(n+2),\qquad
 C_p(n)=v_p(B(n)).
\]
Cancelling \((n!)^2\) gives the exact equivalence
\[
 ((n+2)!)^2\mid(2n)!
 \quad\Longleftrightarrow\quad D(n)^2\mid B(n).
 \tag{3}
\]
Consequently the required inequalities are
\[
 C_p(n)\geq2v_p(D(n))\quad\text{for every prime }p\mid D(n).
 \tag{4}
\]
No check is needed at the other primes, since \(B(n)\) is an integer.

By Legendre's formula,
\[
 C_p(n)=\sum_{j\geq1}
 \left(\left\lfloor\frac{2n}{p^j}\right\rfloor
       -2\left\lfloor\frac n{p^j}\right\rfloor\right)
 =\sum_{j\geq1}\mathbf1_{\{\{n/p^j\}\geq1/2\}}.
 \tag{5}
\]
These are the carries when doubling \(n\) in base \(p\); see also [3]. Each term is zero or one.

If \(p\geq7\), \(i\in\{1,2\}\), and \(p^e\mid n+i\), then for \(1\leq j\leq e\),
\[
 n\equiv-i\pmod{p^j},\qquad
 \{n/p^j\}=1-\frac{i}{p^j}>\frac12.
 \tag{6}
\]
There are therefore at least \(e\) initial carries. When \(v_p(D(n))=1\), condition (4) fails exactly when there is no additional carry at any level \(p^2,p^3,\ldots\).

## 3. The factorization that makes the construction useful

Set \(m=6t+5\), so
\[
 n=F(t)=\frac{m(m+1)}6-2.
\]
Then
\[
 n+1=\frac{(m-2)(m+3)}6=(2t+1)(3t+4),
 \qquad
 n+2=\frac{m(m+1)}6=(t+1)(6t+5).
 \tag{7}
\]
Work with
\[
 \mathcal M_X=\{m\in[X,2X]\cap\mathbb Z:m\equiv5\pmod6\},
 \qquad N_X=\#\mathcal M_X=X/6+O(1).
\]
The four shifts of \(m\) are indexed by
\[
 \mathcal C=\{0,1,-2,3\}.
\]
Their pairwise differences have absolute value at most five. Hence a prime \(p\geq7\) divides at most one of these four factors, and
\[
 v_p(D(n))=v_p(m+c)
\]
for that unique \(c\).

For \(p\geq7\) dividing \(m+c\), write \(m+c=ap\). Put
\[
 d_c=1-2c,\qquad
 i_c=\begin{cases}2,&c=0,1,\\1,&c=-2,3.\end{cases}
\]
An exact calculation gives
\[
 n=\frac{a^2p^2+d_cap}{6}-i_c,
 \qquad
 \frac n{p^2}=\frac{a^2}{6}+\frac{d_ca}{6p}-\frac{i_c}{p^2}.
 \tag{8}
\]
For fixed \(p,c\), the condition \(m\equiv5\pmod6\) puts \(a\) in one residue class modulo six. Its admissible interval has length \(X/p\), so the number of choices is
\[
 L_{p,c}=\frac{X}{6p}+O(1).
 \tag{9}
\]
The fractional parts of \(a^2/6\) are independent of which unit residue class \(p\) occupies modulo six:

| c | d_c | {a²/6} |
|---:|---:|---:|
| 0 | 1 | 1/6 |
| 1 | -1 | 0 |
| -2 | 5 | 1/2 |
| 3 | -5 | 2/3 |

The signs of the small correction in (8) matter at both zero and one half.

## 4. Fixed primes cause a density-zero exceptional set

**Lemma 1.** For each fixed prime \(p\),
\[
 \#\{t\in[T,2T]:C_p(F(t))<2v_p(D(F(t)))\}=o_p(T).
 \tag{10}
\]

**Proof.** First, \(v_p(D(F(t)))\) is bounded in probability. Indeed, by (7),
\[
 D(F(t))=(2t+1)(3t+4)(t+1)(6t+5).
\]
Each linear factor is primitive. For any fixed prime, it either has no roots modulo \(p^h\), or exactly one. The event that the product has valuation at least \(4h\) is therefore contained in at most four residue classes modulo \(p^h\), and its limiting density is at most \(4p^{-h}\).

It remains to show that \(C_p(F(t))\) tends to infinity in probability. Fix a carry bound \(K\) and consider residues modulo \(p^L\). Every base-\(p\) digit at least \(\lceil p/2\rceil\) forces a carry when doubling, regardless of an incoming carry. The number of residues with at most \(K\) carries among the first \(L\) levels is consequently
\[
 O_{p,K}\bigl(L^K\lceil p/2\rceil^L\bigr)=o(p^L).
 \tag{11}
\]
For \(p=2,3\), the polynomial \(F\) is a permutation modulo \(p^L\) for every \(L\): it is a permutation modulo \(p\), and \(F'(t)=12t+11\) is a unit everywhere. Thus (11) pulls back to a density tending to zero.

For \(p\geq5\), use
\[
 24F(t)+49=(12t+11)^2.
\]
Fix \(H\), and discard the residue class \(p^H\mid12t+11\), of density \(p^{-H}\). For \(L>2H\), the map \(t\mapsto F(t)\pmod{p^L}\), restricted to the remaining residues, has fibers of size at most \(2p^H\). To see this, an equation \(y^2\equiv z\pmod{p^L}\) with \(v_p(y)=s<H\) has at most \(2p^s\) solutions modulo \(p^L\); both 12 and 24 are units. Hence the pullback of the set in (11) still has density tending to zero as \(L\to\infty\), for fixed \(H\). Let \(H\to\infty\).

Taking first a bound on \(v_p(D)\), then letting the number of available digit levels grow, proves (10). This argument does not assume independence of \(v_p(D)\) and \(C_p(F(t))\). ∎

## 5. Higher powers of nonfixed primes can be discarded cheaply

Fix \(P\geq7\). A prime square \(p^2\mid m+c\), with \(p>P\), excludes at most
\[
 \frac{X}{6p^2}+O(1)
\]
parameters in \(\mathcal M_X\). Necessarily \(p\leq\sqrt{2X+3}\). Summing over the four shifts and these primes gives
\[
 \#\{m\in\mathcal M_X:\exists p>P,c\in\mathcal C,
          p^2\mid m+c\}=O(X/P)+O(\sqrt X).
 \tag{12}
\]
Thus its normalized contribution is \(O(1/P)+o(1)\).

Primes \(p\leq P\) are dealt with by Lemma 1. After (12), every remaining troublesome prime has \(v_p(D)=1\), so only one additional carry is needed. The limit \(P\to\infty\) will be taken at the end.

## 6. A uniform estimate for the small growing primes

**Lemma 2.** There is a function \(\varepsilon(\delta)\to0\) as \(\delta\downarrow0\) such that the normalized number of parameters having a prime \(7\leq p\leq X^\delta\) with \(v_p(D)=1\) and \(C_p(n)=1\) is at most
\[
 \varepsilon(\delta)+o(1).
 \tag{13}
\]

**Proof.** Use the variable \(t\); its interval has length \(X/6+O(1)\). At the four roots modulo \(p\) of the factors in (7), the derivative \(F'(t)\) is respectively \(5,-5,-1,1\). It is therefore a unit for every \(p\geq7\).

On each root class \(t\pmod p\), the map \(t\mapsto F(t)\pmod{p^L}\) is a bijection onto the corresponding class \(n\equiv-i_c\pmod p\). This follows inductively: at each step the next digit is multiplied by the unit \(F'(t)\pmod p\).

Take \(L=\lfloor\log_p(X/12)\rfloor\). The first carry is forced. If \(C_p(n)=1\), there is no carry at levels \(2,\ldots,L\). At every lift, at most \((p+1)/2\) choices for the next digit avoid a carry. There are at most \(((p+1)/2)^{L-1}\) eligible residues modulo \(p^L\) in the specified root class. Each occurs at most \(O(X/p^L)\) times in the interval, because \(p^L\leq X/12\). Thus, per root, the number of bad parameters is at most
\[
 \frac{CX}{p}\left(\frac{p+1}{2p}\right)^{L-1}
 \leq\frac{C'X}{p}
       \exp\left(-c\frac{\log X}{\log p}\right),
 \tag{14}
\]
with absolute positive constants. One can use \((p+1)/(2p)\leq4/7\).

Partial summation using \(\pi(y)\ll y/\log y\) shows
\[
 \limsup_{X\to\infty}
 \sum_{7\leq p\leq X^\delta}\frac1p
       \exp\left(-c\frac{\log X}{\log p}\right)
 \ll\int_{1/\delta}^{\infty}\frac{e^{-cu}}u\,du.
 \tag{15}
\]
Changing the positive absolute constants does not affect the assertion. The integral tends to zero. Summing (14) over the four roots proves (13). ∎

## 7. The uniform equidistribution lemma

This is the main analytic step. It uses the standard second-derivative estimate [4]: if a real quadratic phase has constant second derivative of size \(\lambda>0\), then on an interval of \(L\) consecutive integers,
\[
 \left|\sum e(f(b))\right|\ll L\sqrt\lambda+\lambda^{-1/2},
 \qquad e(x)=e^{2\pi i x}.
 \tag{16}
\]
The bound is independent of the linear coefficient and of the starting point of the interval.

**Lemma 3.** Fix \(\delta,\eta>0\) and an integer \(J\geq4\). Suppose
\[
 X^\delta\leq p,\qquad
 p\leq X^{1/2-\eta},\qquad p\leq X^{2/J-\eta}.
 \tag{17}
\]
For each \(c\in\mathcal C\), as \(a\) runs through the admissible progression in (9), the vector
\[
 \left(\{n/p^2\},\ldots,\{n/p^J\}\right)
 \tag{18}
\]
is equidistributed in \([0,1)^{J-1}\), uniformly in the primes in (17).

**Proof.** Write \(a=a_0+6b\). The length of the interval of \(b\)'s is \(L_{p,c}\asymp X/p\). Apply the multidimensional Weyl criterion to any fixed nonzero integer frequency vector \((h_2,\ldots,h_J)\).

If only \(h_2\) is nonzero, remove the integer-valued variation of \(h_2a^2/6\) in (8). The remaining phase is linear in \(b\), with slope \(h_2d_c/p\). For large \(X\), the prime \(p\) exceeds the fixed nonzero integer \(|h_2d_c|\). A geometric-series estimate gives an exponential sum \(O_h(p)\), and hence a normalized bound
\[
 O_h(p^2/X)=o(1)
 \tag{19}
\]
uniformly by (17).

Otherwise, let \(j_0\geq3\) be the smallest index with \(h_{j_0}\ne0\). Again remove the integer-valued quadratic part of the \(h_2\) term, if present. The resulting quadratic phase in \(b\) has second derivative
\[
 12\sum_{j=3}^J h_jp^{2-j}\asymp_h p^{2-j_0}.
 \tag{20}
\]
There is no asymptotic cancellation: the first nonzero term dominates because the frequencies are fixed and \(p\to\infty\).

Applying (16) and dividing by \(L_{p,c}\asymp X/p\) gives
\[
 O_h\left(p^{-(j_0-2)/2}+\frac{p^{j_0/2}}X\right)
 \leq O_h\left(X^{-\delta/2}+\frac{p^{J/2}}X\right)=o(1).
 \tag{21}
\]
This is uniform in (17). Finite trigonometric approximations to boxes, or the Weyl criterion followed by such approximations, prove uniform equidistribution. In particular,
\[
 \#\{a:\{n/p^j\}<1/2\ \text{for }2\leq j\leq J\}
  =\bigl(2^{1-J}+o(1)\bigr)L_{p,c},
 \tag{22}
\]
uniformly in these prime ranges. ∎

A one-dimensional version will also be used: for
\[
 X^{1/2+\eta}\leq p\leq X^{2/3-\eta},
 \tag{23}
\]
\(\{n/p^3\}\) is uniformly distributed over the same progression. Its quadratic phase has second derivative of size \(1/p\), and its normalized exponential sums are
\[
 O_h\left(p^{-1/2}+p^{3/2}/X\right)=o(1).
 \tag{24}
\]
In particular, its probability of being below one half is \(1/2+o(1)\).

The proof of this lemma does not invoke a heuristic that different digit levels are independent. Equation (22) is the conclusion of the exponential-sum calculation.

## 8. Primes above the square-root parameter scale

Fix \(\eta>0\). If \(p\geq X^{1/2+\eta}\), then uniformly for \(m\in\mathcal M_X\),
\[
 a/p=(m+c)/p^2=o(1).
\]
Use the exact formula (8), together with the table in Section 3. For all sufficiently large \(X\),
\[
 \{n/p^2\}\begin{cases}
 <1/2,&c=0,\\
 >1/2,&c=1,-2,3.
 \end{cases}
 \tag{25}
\]
For \(c=1\), the correction to an integer is negative, so the fractional part is near one, not zero. For \(c=-2\), the correction to one half is positive, since \(5ap>6\). For \(c=3\), the fractional part is close to two thirds.

Therefore, above this scale, the three factors \(m+1,m-2,m+3\) always provide the needed second carry. Only the factor \(m\) can be bad.

For \(p>X^{2/3-\eta}\), it is enough to bound all incidences \(p\mid m\), without considering carries. By the reciprocal-prime asymptotic
\[
 \sum_{X^\alpha<p\leq X^\beta}\frac1p
  =\log(\beta/\alpha)+o(1),\qquad 0<\alpha<\beta,
 \tag{26}
\]
their normalized count has limsup at most
\[
 \log\frac1{2/3-\eta}.
 \tag{27}
\]
The upper endpoint is actually \(2X+3\); its logarithmic exponent tends to one, so it has the same limit.

For the intermediate range (23), only \(c=0\) matters, and being bad requires \(\{n/p^3\}<1/2\). By (24), its normalized contribution is at most
\[
 \frac12\log\frac{2/3-\eta}{1/2+\eta}+o(1).
 \tag{28}
\]
Sending \(\eta\downarrow0\), the combined contributions (27)--(28) tend to
\[
 \log\frac32+\frac12\log\frac43=\frac12\log3.
 \tag{29}
\]

## 9. Primes below that scale

Fix \(\delta>0\). Apart from narrow bands around the finitely many exponent boundaries, every prime in \([X^\delta,X^{1/2}]\) lies in a range
\[
 X^{2/(J+1)+\eta}<p<X^{2/J-\eta},\qquad J\geq4,
 \tag{30}
\]
for one of finitely many \(J\)'s (bounded in terms of \(\delta\)). For \(J=4\), the upper inequality also enforces the needed margin below \(X^{1/2}\); for \(J>4\) that margin is automatic, after reducing \(\eta\) if necessary.

A bad prime has no carry at any of the levels \(p^2,\ldots,p^J\), so Lemma 3 bounds its conditional frequency by \(2^{1-J}+o(1)\). This is an upper bound even if there are additional digit levels above \(J\).

For each of the four linear factors, (9), (22), and (26) give the normalized contribution
\[
 \leq 2^{1-J}\log\frac{J+1}{J}+o(1)
\]
as the boundary margin tends to zero. Summing over the four factors and the finitely many bands gives at most the corresponding partial sum of
\[
 4\sum_{J=4}^{\infty}2^{1-J}\log\frac{J+1}{J}.
 \tag{31}
\]

For completeness, the narrow excluded bands introduce no difficulty. A band \(X^{\alpha-\eta}<p<X^{\alpha+\eta}\), with fixed \(\alpha>0\), contributes at most the count of all divisibility incidences, which after normalization is
\[
 O\left(\log\frac{\alpha+\eta}{\alpha-\eta}\right)+o(1)=O_\alpha(\eta)+o(1).
\]
There are only finitely many such boundaries at fixed \(\delta\). Their total tends to zero with \(\eta\). The small interval next to the cutoff \(X^\delta\) can be handled in exactly the same way, or absorbed into Lemma 2 with a slightly larger cutoff.

All occurrences of the \(O(1)\) error in (9) are harmless in prime summation: there are only \(O(X/\log X)\) relevant primes. The uniform \(o(1)\) errors in (22) and (24), multiplied by \(X/p\), sum to \(o_\delta(X)\), since \(\sum_{X^\delta<p\leq2X+3}1/p=O_\delta(1)\).

## 10. Combining the estimates and keeping the limits in order

Let \(\mathcal B_X\) be the parameters for which (3) fails. Fix \(P\), then a small \(\delta\), then boundary margins \(\eta\). The preceding estimates cover every possible prime obstruction:

- primes \(p\leq P\): density zero by Lemma 1;
- repeated factors \(p^2\mid m+c\) for \(p>P\): normalized cost \(O(1/P)+o(1)\) by (12);
- remaining primes below \(X^\delta\): normalized cost at most \(\varepsilon(\delta)+o(1)\) by Lemma 2;
- all remaining primes: Sections 8--9 and their boundary bands.

Take \(\limsup_{X\to\infty}\) first. Next let the boundary margins tend to zero, then let \(\delta\downarrow0\) and \(P\to\infty\). The union bound gives
\[
 \limsup_{X\to\infty}\frac{\#\mathcal B_X}{N_X}
 \leq \frac12\log3
       +4\sum_{J=4}^{\infty}2^{1-J}\log\frac{J+1}{J}
 =1-\rho.
 \tag{32}
\]
The bound is strictly less than one without relying on numerical approximation. Indeed,
\[
 \log\frac{J+1}{J}\leq\frac1J\leq\frac14\qquad(J\geq4),
\]
so
\[
 4\sum_{J=4}^{\infty}2^{1-J}\log\frac{J+1}{J}
 \leq\sum_{J=4}^{\infty}2^{1-J}=\frac14.
\]
This proves (2) and yields a positive lower proportion of good parameters.

Finally, the affine substitution \(m=6t+5\) takes the \(m\)-intervals used above into the corresponding \(t\)-intervals, up to \(O(1)\) changes at their endpoints. Thus (32) is exactly the proposed lower bound (1). The polynomial \(F(t)\) is strictly increasing for \(t\geq0\), so infinitely many good parameters yield infinitely many distinct \(n\). ∎

## 11. Exact checks and limits of what has been proved here

At \(t=5\), equivalently \(m=35\), the construction gives \(n=208\). Then
\[
 D(n)=209\cdot210=2\cdot3\cdot5\cdot7\cdot11\cdot19.
\]
The corresponding valuations of \(\binom{416}{208}\) are
\[
 (v_2,v_3,v_5,v_7,v_{11},v_{19})=(3,3,2,2,2,2).
\]
Every required valuation is two, proving \((210!)^2\mid416!\).

An exact scan of \(0\leq t<1000\) gives 120 valid parameters. The first values of \(n\) in this family are
\[
208,\ 3430,\ 8175,\ 17440,\ 19435,\ 34578,\ 46374,\ 47435,\ldots
\]
These are numerical checks only. The lower bound in (1) is asymptotic and does not assert that every finite sample has at least that fraction of successes.

The argument does not control \(n+3,\ldots,n+k\). For general \(k\), condition (3) becomes
\[
 \left(\prod_{i=1}^k(n+i)\right)^2\mid\binom{2n}{n},
\]
and the four-linear-factor argument no longer covers its denominator primes. No implication from the proposed \(k=2\) result to all \(k\) is asserted.

**The main checkpoints for an independent mathematical review are:** the uniform exponential-sum bounds (19)--(24); the simple-root digit count (14); and the aggregation over primes with the order of limits in Section 10. The exact computations do not verify those analytic steps.

## References

[1] P. Erdős, R. L. Graham, I. Z. Ruzsa, and E. G. Straus, “On the prime factors of the central binomial coefficient,” *Mathematics of Computation* 29 (1975), 83--92. The factorial question appears on p. 90. Original paper: `https://renyi.hu/~p_erdos/1975-27.pdf`.

[2] Thomas Bloom, *Erdős Problems*, problem 727 and its discussion thread; accessed 7 September 2026. `https://www.erdosproblems.com/727` and `https://www.erdosproblems.com/forum/thread/727`.

[3] Kevin Ford and Sergei Konyagin, “Divisibility of the central binomial coefficient,” *Transactions of the American Mathematical Society* 374 (2021), 923--953; arXiv:1909.03903. This is background for the carry-count formulation; its main theorem is not used to infer the proposed quadratic-family result.

[4] Olivier Robert, “On van der Corput's k-th derivative test for exponential sums,” *Indagationes Mathematicae* 27 (2016), 559--589, DOI: 10.1016/j.indag.2015.11.009, Theorem 1 (second-derivative test). Only this standard estimate is used in Section 7.

[5] Leo Goldmakher, “A quick proof of Mertens' theorem,” Williams College, two-page mathematical note. `https://web.williams.edu/Mathematics/lg5/mertens.pdf`.

The additional standard analytic input (26) is Mertens' reciprocal-prime asymptotic [5], which also follows from the prime number theorem and partial summation. The proposed construction and estimates specific to it are the argument developed in this note, not claims attributed to the references.
