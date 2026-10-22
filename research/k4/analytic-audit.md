# Analytic audit and the quadratic obstruction at k = 4

Date: 22 October 2026. Base revision: `8c1e3e50040224b5d3394f91f8c491be8307cd47`.

This is an independent agent audit of the mathematical argument, not human peer
review or a proof-assistant certificate. The audit found no specific mathematical
error in the main k=2 and k=3 analytic estimates after checking their hypotheses.
It does find a rigorous obstruction to the most direct k=4 extension: a quadratic
cannot have all four consecutive translates split over the rationals.

**Final status: no k=4 infinitude proof was obtained.** The main continuation is
the partial-splitting family `m=(y^2-1681)/552`, `n=m(m+1)/210-3`. The audited
estimates bound its failure contribution from primes below the parameter height
by 0.5625676641884855292. This is a partial prime-range bound, not a bound on all
failures and not a lower bound for successful factorial divisibility. The
remaining quadratic prime divisors above that height have not been controlled.

Sections 2--3 audit the existing argument and its quadratic obstruction.
Sections 4--5 develop the quartic carry estimates and the main 552 candidate.
Sections 6--7 record supplementary constructions and bounded searches. The
non-even alternatives improve some lower-prime costs but have worse large-prime
residue behavior; they do not complete the missing step. No probabilistic
independence heuristic is used as a theorem in this note.

## 1. What was read and independently recomputed

- `erdos_727_k2_candidate.md`, especially Sections 4--10.
- `papers/k3-extension-2026-09-17.tex`, all mathematical sections.
- The k=3 signed residue table was recomputed with exact modular arithmetic.
- Robert's cited second-derivative estimate and the fixed-modulus reciprocal-prime
  asymptotic were checked against their sources.

The exact k=3 table is reproduced by the following pairs, where each row is
`(unit square class, six residues, number of bad shifts)`:

```text
(1,   [15, 196,   1, 0, 175,  36], 3)
(79,  [135,154,  79, 0, 175, 114], 1)
(109, [165,154, 109, 0, 175, 144], 0)
(121, [135,196, 121, 0, 175, 156], 0)
(151, [165,196, 151, 0, 175, 186], 0)
(169, [15, 154, 169, 0, 175, 204], 1)
```

Every square class has eight unit preimages modulo 210. The exact average bad
count is therefore `(3+1+0+0+0+1)/6 = 5/6`.

## 2. Audit of the existing analytic proof

### 2.1 Exact reduction and initial carries

For any positive integer k,

\[
((n+k)!)^2\mid(2n)!
\quad\Longleftrightarrow\quad
\left(\prod_{i=1}^k(n+i)\right)^2\mid\binom{2n}{n}.
\]

Legendre's formula gives

\[
C_p(n)=v_p\binom{2n}{n}
 =\sum_{j\geq1}\mathbf1_{\{\{n/p^j\}\geq1/2\}}.
\]

If p>2k and p^e divides n+i, then all first e summands equal one. In the
existing quadratic families, sufficiently large p divides at most one of the
fixed linear factors. After removing repeated prime factors, the required
valuation is two and a failure is exactly the absence of every carry after the
first. No independence assumption is needed in this reduction.

### 2.2 Fixed primes

The product of primitive linear factors has a valuation bounded in probability:
if there are R factors, a product valuation at least Rh implies that one factor
has valuation at least h. Its upper limiting density is at most R p^{-h}.

Among residues modulo p^L, at most

\[
O_{p,K}\bigl(L^K\lceil p/2\rceil^L\bigr)=o(p^L)
\]

have at most K carries. Indeed, a digit at least ceil(p/2) forces a carry,
independently of the incoming carry. Hence only K digit positions may contain
such digits.

For the denominator primes of the two existing families, the polynomial is
linear with unit slope modulo p and has unit derivative. Hensel lifting makes it
a permutation modulo every p^L. For the other primes, completion of the square
works as written in the manuscript. Outside the event p^h divides the affine
square-root variable, the relevant fibers have size at most 2p^h. First let
L tend to infinity, then h tend to infinity. This is a valid argument for
`C_p(F(t)) -> infinity in probability`; it does not assume independence from
the denominator valuation.

### 2.3 Repeated nonfixed primes

For a fixed shift c, m+c=O(X). Thus p^2 | m+c implies p=O(sqrt(X)). The count
of the relevant residue class is X/(Qp^2)+O(1), where Q is fixed. Summation gives
O(X/P)+O(sqrt(X)), and after normalization the cost is O(1/P)+o(1). The fixed
cutoff P is removed only after the height limit. This step is sound.

### 2.4 Small growing primes

For p outside the finite exceptional set, each linear factor gives a simple
root of the polynomial F(t)+i modulo p. On that root class the map

\[
t\pmod{p^L}\longmapsto F(t)\pmod{p^L}
\]

is bijective onto the residue class `n = -i (mod p)`. This follows digit by
digit from the unit derivative. There are at most `((p+1)/2)^(L-1)` possible
lifts that lack every later carry.

Taking L=floor(log_p(X/(2Q))) ensures that each lift occurs O(X/p^L) times in
the interval. Uniformly for p<=X^delta the conditional failure count is bounded
by

\[
\frac{A_Q X}{p}\exp\left(-c\frac{\log X}{\log p}\right).
\]

Partial summation using pi(y)=O(y/log y) bounds its normalized sum by a constant
times

\[
\int_{1/\delta}^{\infty}\frac{e^{-cu}}{u}\,du+o(1).
\]

This tends to zero with delta. The constants may depend on the fixed family;
the proof never needs them to be uniform in Q.

### 2.5 The Fourier step: the important cancellation is valid

On a fixed large-prime root class, write m+c=ap and a=a0+Qb. The exact formula is

\[
n=\frac{a^2p^2+d_c ap}{Q}-i_c,
\qquad d_c=1-2c\ne0.
\]

For the h2 term in a fixed Fourier frequency vector,

\[
\frac{a^2}{Q}
 =\frac{a_0^2}{Q}+2a_0b+Qb^2.
\]

Its variation in b is integer-valued. Removing that integer polynomial from
the phase is legitimate. If no higher frequency is present, the remaining
slope is h2*d_c/p. The geometric sum is O_h(p), since p eventually exceeds the
fixed nonzero integer |h2*d_c|. With progression length N~X/p the normalized
bound is O_h(p^2/X).

If j0>=3 is the first nonzero higher frequency, the residual second derivative
is

\[
2Q\sum_{j=3}^J h_jp^{2-j}
 =2Qh_{j_0}p^{2-j_0}(1+O_h(p^{-1})).
\]

There is no asymptotic cancellation because the Fourier frequencies are fixed
before p tends to infinity. Robert's Theorem 1 gives

\[
\frac{|S|}{N}
 \ll_h p^{-(j_0-2)/2}+\frac{p^{j_0/2}}X.
\]

This tends to zero uniformly if

\[
X^\delta\leq p\leq X^{1/2-\eta},\qquad
p\leq X^{2/J-\eta}.
\]

For each fixed box approximation there are only finitely many frequencies, so
the uniform estimates imply the asserted uniform box count. The exact digit
relations have frequencies depending on p; they do not contradict this fixed
frequency equidistribution conclusion.

### 2.6 Large-prime signs and prime-class weighting

For p>=X^(1/2+eta), a/p=O(X^(-2eta)). Thus n/p^2 is a fixed residue s/Q plus a
uniformly small correction. At s=0 the bad sign is d_c>0; at s=Q/2 it is
d_c<0. In particular, a negative correction at zero produces a fractional
part near one. The k=3 zero column is correctly marked good.

For intermediate primes p<X^(2/3-eta), the single coordinate n/p^3 has second
derivative of size 1/p. Its normalized sum bound is
O(p^(-1/2)+p^(3/2)/X)=o(1), uniformly in the relevant root classes.

Keliher--Lee equation (6) is the unconditional fixed-modulus estimate

\[
\sum_{p\leq y,\ p\equiv v\ (Q)}\frac1p
 =\frac{\log\log y}{\varphi(Q)}+M(Q,v)+O_Q(1/\log y).
\]

Subtracting its values at X^alpha and X^beta gives precisely the needed
`log(beta/alpha)/phi(Q)+o(1)`. The GRH-dependent results later in that paper are
not being used. All Q here are fixed.

### 2.7 Errors and coverage

Because every relevant large prime in a linear factor is O(X), the total of the
O(1) progression errors is O(pi(CX))=o(X). Uniform Fourier errors multiply
X/p; their sum is o(X) on every fixed exponent range bounded away from zero.

At fixed delta there are finitely many exponent bands. Their omitted boundary
strips have total reciprocal-prime mass tending to zero when eta tends to zero.
The order `X -> infinity`, then `eta -> 0`, then `delta -> 0`, then `P -> infinity`
is valid. No infinite sum of nonuniform errors is required.

### Audit conclusion

The above checks support the internal validity of the existing k=2/k=3
argument. They are not a substitute for outside mathematical review. Minor
presentation issue: the k=3 TeX bibliography refers to two PDF paths not present
in the pinned repository. The original Markdown manuscript supplies the
small-prime and Fourier arguments needed to fill those references.

## 3. General quadratic template and its exact limitation

Suppose, hypothetically, that a fixed Q and residue r give

\[
n=\frac{m(m+1)}Q-k,\qquad m=Qt+r,
\]

and each n+i splits into two rational linear factors in m. Write the shifts as
c and 1-c, where

\[
c(c-1)=Q(k-i).
\]

After absorbing fixed denominator primes, there are 2k simple roots. Define the
same boundary-aware bad indicator and its mean b over units modulo Q. The proof
template gives

\[
\liminf \text{(successful parameter proportion)}
 \geq 1-\frac{k}{2}S-\frac b2\log3,
\qquad
S=4\sum_{J\geq4}2^{1-J}\log\frac{J+1}{J}.
\]

For k=4 this would require b<1.1233430797845023..., since
S=0.1914703720396162.... This numerical criterion is moot for fully split
quadratics: such a family cannot exist.

**Proposition (barrier).** There is no quadratic F(t)=At^2+Bt+C over Q with
A!=0 for which F(t)+1, F(t)+2, F(t)+3 and F(t)+4 all split over Q.

**Proof.** The discriminants

\[
\Delta_i=B^2-4A(C+i),\qquad i=1,2,3,4,
\]

would all be rational squares. They form an arithmetic progression with
nonzero common difference -4A. Clearing their square denominators gives four
distinct integer squares in arithmetic progression, contradicting Fermat's
four-squares theorem. QED.

This blocks every quadratic with eight rational linear factors, not merely a
bounded search in Q. It does not disprove Erdős 727 for k=4, or block higher
degree or nonpolynomial constructions.

## 4. Quartic continuation: what the analytic audit can establish

The construction agent supplied

\[
F(x)=\frac{x^4+419x^2+43260}{210},\qquad x\equiv14\pmod{210}.
\]

It has the eight quadratic factors

\[
\begin{aligned}
210(F+1)&=(x^2+189)(x^2+230),\\
210(F+2)&=(x^2+195)(x^2+224),\\
210(F+3)&=(x^2+209)(x^2+210),\\
210(F+4)&=(x^2-x+210)(x^2+x+210).
\end{aligned}
\]

### 4.1 Why a direct reuse of the old prime sum is invalid

For x~X these factors are O(X^2). A root count modulo p gives O(X/p)+O(1), but
now primes extend to O(X^2). The accumulated O(1) error can exceed the entire
parameter interval. Consequently a residue table and Mertens' theorem alone
do not control p>X. The old argument's O(X) bound on prime divisors was essential.

### 4.2 Repeated factors and fixed primes remain manageable

For p outside finitely many discriminant and resultant primes, a square divisor
p^2 of any quadratic factor specifies at most two root classes modulo p^2.
Since p^2<=O(X^2), one has p=O(X), and

\[
\sum_{P<p\ll X}\left(O(X/p^2)+O(1)\right)
 =O(X/P)+O(X/\log X).
\]

Hence repeated primes can still be discarded at cost O(1/P)+o(1).

For any fixed p, valuation tails of the fixed nonzero polynomial factors tend
to zero. Carry scarcity also has density zero under a nonconstant polynomial
map: discard a small p-adic neighborhood of its finitely many critical roots,
then use a bounded-fiber lifting argument. This can be written directly for
the explicit quartic. It does not solve the growing-prime ranges.

### 4.3 Uniform Fourier estimates below X^(4/5)

Condition on a simple root r modulo p of one of F+i, and impose x=14 (mod 210).
Then x=r+210pb after choosing a suitable representative r. The b-interval has
length N~X/(210p).

For a fixed Fourier vector in levels 2,...,J, remove all integer-valued leading
terms coming from levels j<=4. If a frequency at j>=5 is present, let j0 be the
smallest such j. The fourth derivative is a fixed nonzero constant times

\[
\sum_{j=5}^J h_jp^{4-j}\asymp_h p^{4-j_0}.
\]

The explicit fourth-derivative theorem of Arias de Reyna gives, with
lambda~p^(4-j0),

\[
|S|/N\ll N^{-1/8}+\lambda^{1/14}+(\lambda N^4)^{-1/8}.
\]

This tends to zero if p grows and p<=X^(4/J-eta), because
lambda N^4~X^4/p^j0.

There is also a shorter classical Weyl proof for this particular polynomial
phase. Its fourth-degree leading coefficient alpha is asymptotic to a fixed
nonzero constant times p^(4-j0). Put t=1/|alpha|, q=round(t), and a=sign(alpha).
Then `(a,q)=1` and

\[
|\alpha-a/q|=\frac{|q-t|}{tq}\leq\frac1{2tq}\leq q^{-2}
\]

for all sufficiently large t. Weyl's inequality bounds the normalized sum by
`O(N^epsilon (q^(-1)+N^(-1)+q/N^4)^(1/8))`. Here
q asymptotic to p^(j0-4), so q/N^4 is O(p^j0/X^4). On every fixed exponent
band with positive margins, all three terms save a fixed power of X. Choosing
epsilon after the margins absorbs N^epsilon. This argument does not require
alpha to have an exact denominator q. The same reasoning verifies the
construction agent's cubic and quadratic small-leading-coefficient estimates.

If all frequencies at j>=5 vanish, the highest remaining level j*=2,3,4
produces a phase of degree j*-1. Its leading coefficient has exact denominator
p outside finitely many exceptional primes. The lower coefficients may have
larger p-power denominators, so a complete-period argument modulo p would be
incorrect. The classical Weyl inequality is uniform in those lower
coefficients and applies to the rational leading coefficient.

- Pure level 2: a geometric sum is O(p), giving cancellation for p<X^(1/2-eta).
- Highest level 3: a quadratic Weyl bound gives cancellation for
  `N^2 >> p`, hence p<X^(2/3-eta).
- Highest level 4: a cubic Weyl bound gives cancellation for
  `N^3 >> p`, hence p<X^(3/4-eta).

For the degree-two and degree-three leading coefficients, the required
nonvanishing follows from fixed resultants. Put
P_i(x)=x^4+419x^2+43260+210i. None shares a root over Qbar with P_i', P_i'', or
P_i'''. For the second derivative the obstruction is a divisor of
`679555+7560i`; for the third derivative it is a divisor of `43260+210i`.
The first derivative gives the already fixed discriminant primes. All are
absorbed in P before limits are taken.

The resulting usable coordinate sets are:

| Prime exponent alpha, p=X^alpha | Equidistributed coordinates sufficient here |
|---|---|
| `alpha<1/2`, `p<X^(4/J-eta)`, `J>=8` | levels 2,...,J |
| `1/2<alpha<2/3`, `p<X^(4/J-eta)`, `J=6,7` | levels 3,...,J |
| `2/3<alpha<3/4` | levels 4,5 |
| `3/4<alpha<4/5` | level 5 |

Every one of the eight irreducible quadratics has `1+chi_D(p)` roots for
nonexceptional p. A fixed-modulus reciprocal-prime sum therefore has average
one root per quadratic, for a total average root count eight. Weighting the
uniform conditional bounds by these counts gives the provisional partial
failure bound

\[
\begin{split}
C_{<4/5}={}&8\sum_{J\geq8}2^{1-J}\log\frac{J+1}{J}
 +8\sum_{J=6}^7 2^{2-J}\log\frac{J+1}{J}\\
&+2\log(9/8)+4\log(16/15)\\
={}&0.6175856642659768\ldots.
\end{split}
\]

Boundary strips are removed in the same ordered-limit manner as before. This
bound only handles primes below X^(4/5) and is not a positive-density theorem
for successful parameters.

### 4.4 A deterministic large-prime fact for the two new factors

For g_sigma=x^2+sigma*x+210=ap, sigma in {-1,1}, the exact identity is

\[
\frac{F(x)}{p^2}
 =\frac{a^2}{210}-\frac{\sigma ax}{105p}-\frac4{p^2}.
\]

For p>X^(3/2+eta), the correction has magnitude O(X^(-2eta)). When x=14 mod210,
the sigma=1 factor gives a^2=0 mod210 and a negative correction, hence a
fractional part near one. The sigma=-1 factor gives a^2 equal to 154 or 196
mod210 for every unit prime class, both above 105, and a positive correction.
Thus both new factors supply a second carry at all sufficiently large primes
in this range. This removes their very largest prime divisors, but leaves the
range X^(4/5)<p<=X^(3/2+eta) unhandled.

### Quartic conclusion

The quartic construction controls factor sizes, permits a real partial Fourier
analysis, and has a useful deterministic large-prime property. It does not yet
provide the missing count for the remaining prime ranges. A claimed k=4 proof
that substitutes the old Mertens union bound above p=X would contain a gap.

## 5. Stronger averaged carry estimates from quadratic-root equidistribution

This section improves the elementary partial bound in Section 4. It uses a
deep published theorem, rather than assuming that quadratic roots have random
digits. The improvement still stops below the parameter height X.

### 5.1 The precise external input

Toth's Theorem 6.1.3 states that, for a fixed integral quadratic of nonsquare
discriminant, its normalized roots r/p are equidistributed in [0,1] as p runs
through primes in a fixed arithmetic progression that supplies infinitely many
roots. This is an arithmetic-progression version of the Duke--Friedlander--
Iwaniec quadratic-root theorem; Toth treats both signs of the discriminant.

The theorem is reproduced in the author's 2018 dissertation, Chapter 6,
printed page 96 (PDF page 98). No assertion about roots modulo p^2 is used.
Partial summation gives the corresponding reciprocal-prime-weighted statement
on every fixed exponent interval X^alpha<p<X^beta, 0<alpha<beta. In particular,
if E is a fixed finite union of intervals, its weighted root mass is
`|E|` times the total weighted root mass, plus o(1).

### 5.2 A finite-differencing lemma

Let g be such an irreducible quadratic. For each root r modulo p, suppose a
polynomial phase on an interval of N consecutive integers has degree d and
leading coefficient

\[
\theta_{p,r}=A\frac rp+\frac Bp\pmod1,
\qquad A\in\mathbb Z\setminus\{0\},\quad B\in\mathbb Z,
\]

where A and B are fixed. The lower coefficients and the interval's starting
point may depend arbitrarily on p and r. Assume N tends to infinity uniformly
over the considered roots, for example N asymptotic to X/(Mp) with fixed M and
p<=X^(1-eta). Then

\[
\sum_{X^\alpha<p<X^\beta}\frac1p
 \sum_{g(r)=0\ (p)}
 \left|\frac1N\sum_b e(f_{p,r}(b))\right|=o(1)
\qquad(0<\alpha<\beta<1).
\]

The same assertion holds with primes restricted to any fixed progression.

Here is a proof that avoids a potentially misleading use of an N^epsilon
Weyl bound. Fix a differencing cutoff H. Repeated van der Corput differencing,
d-1 times with shifts at most H, reduces the normalized sum to finitely many
linear sums whose slopes are

\[
d!h_1\cdots h_{d-1}\theta_{p,r},
\qquad1\leq h_j\leq H.
\]

The other error terms tend to zero as H tends to infinity, after N tends to
infinity. For fixed H and epsilon>0, exclude the roots r/p in the finite union
of arcs on which one of

\[
\bigl\|d!A h_1\cdots h_{d-1}(r/p)\bigr\|<2\epsilon.
\]

The length of this union is O_H(epsilon). Outside it, B/p is harmless for
large p and every resulting normalized linear sum is O_H(1/(N epsilon)).
All bounds are uniform in the lower coefficients. Toth's theorem bounds the
weighted mass of the excluded roots by O_H(epsilon)+o(1). First let X tend
to infinity, then epsilon to zero, then H to infinity. For d=1 the same
argument is just the geometric-sum estimate.

This is an averaged estimate: individual prime roots can be exceptional. It
uses only fixed arcs, not a shrinking-interval equidistribution theorem.

### 5.3 Application to a centered quartic and a quadratic factor

Write

\[
F(x)=\frac{x^4+Bx^2+Lx+C}{D},\qquad
g(x)=x^2+ux+v\mid D(F(x)+i),
\]

with all coefficients fixed. Refine the admissible parameter progression if
necessary to a fixed modulus M divisible by D. On a root class write
x=R+Mpb, with R reducing to r modulo p. Removing integer-valued polynomial
terms, the leading coefficients of carry levels 2, 3, 4 are respectively

\[
\begin{array}{c|c|c}
\text{level}&\text{degree in }b&\text{leading coefficient modulo one}\\ \hline
2&1&\displaystyle\frac{M}{D}\frac{(4u^2-4v+2B)r+4uv+L}{p}\\[1ex]
3&2&\displaystyle\frac{M^2}{D}\frac{-6ur+B-6v}{p}\\[1ex]
4&3&\displaystyle\frac{4M^3}{D}\frac rp.
\end{array}
\]

These formulas follow by reducing P'(r), P''(r)/2 and P'''(r)/6 modulo
g(r), for P(x)=x^4+Bx^2+Lx+C. In the quartic constructions considered here,
`4u^2-4v+2B` is nonzero for every irreducible quadratic factor. Therefore:

- Levels 2 and 4 are jointly equidistributed in the weighted prime-root
  average throughout p<X^(1-eta).
- If u is nonzero, levels 2, 3, 4 are jointly equidistributed throughout the
  same range. This additional level 3 is available for both new odd factors.
- If u=0, the pure level-3 quadratic coefficient is a fixed nonzero constant
  divided by p. Its second-derivative bound remains uniform only up to
  p<X^(2/3-eta). No extension of that coordinate is claimed here.

For a joint Fourier vector, take the highest nonzero level among 2,3,4.
The preceding finite-differencing lemma applies to its leading coefficient;
all lower coefficients may be ignored for this purpose. If a level j>=5 is
present, the fourth-derivative estimate of Section 4.3 applies uniformly when
p<X^(4/J-eta). Thus these are joint box estimates, not marginal estimates
incorrectly multiplied together.

For an even quadratic factor (u=0), the available no-additional-carry bounds
are consequently

| Prime exponent range | Usable carry levels | Conditional failure bound |
|---|---|---|
| `4/(J+1)<alpha<4/J`, `J>=6` | `2,...,J` | `2^(1-J)` |
| `2/3<alpha<4/5` | `2,4,5` | `1/8` |
| `4/5<alpha<1` | `2,4` | `1/4` |

For an odd quadratic factor (u nonzero), the last two bounds improve to
1/16 and 1/8, using levels 2,3,4,5 and 2,3,4 respectively. Every assertion
is made away from fixed-width exponent boundary strips, whose costs are then
sent to zero in the usual order of limits.

### 5.4 Numerical consequence for the first quartic

There are six even and two odd irreducible quadratic factors, each with
average root count one. Put

\[
T_6=\sum_{J\geq6}2^{1-J}\log\frac{J+1}{J}
 =0.008579551796005669\ldots.
\]

The total failure cost from primes strictly below the parameter height is

\[
\begin{aligned}
C_{<1}
&=8T_6+\left(\frac68+\frac2{16}\right)\log\frac65
       +\left(\frac64+\frac28\right)\log\frac54\\
&=0.6186689913626227\ldots.
\end{aligned}
\]

This improves the scope of the earlier elementary bound to all p<X^(1-eta).
It does not bound the failures from larger prime divisors.

### 5.5 The improved partial-splitting family

For the construction

\[
m=(y^2-1681)/552,\qquad F(y)=m(m+1)/210-3,
\]

the six irreducible quadratic factors consist of four even and two odd ones.
Their contribution below the parameter height is

\[
C_{\rm quad}=6T_6+\frac58\log\frac65+\frac54\log\frac54
 =0.44435772291501785\ldots.
\]

The four remaining linear factors can be handled with stable residue tables
for levels 2,3,4. On y=23971 modulo D=63987840 the construction agent's exact
joint bad counts are

\[
B_2=\frac53,\qquad B_{23}=\frac56,\qquad B_{234}=\frac{311}{1242}.
\]

All three counts were independently recomputed in this audit by a separate
C++ enumeration of all phi(D)=13,989,888 unit classes. The individual results
are as follows.

| Linear root s | Bad at level 2 | Bad at levels 2 and 3 | Bad at levels 2, 3 and 4 |
|---|---:|---:|---:|
| 41 | 49/99 | 49/198 | 94/759 |
| -41 | 1/3 | 1/6 | 0 |
| 97 | 1/3 | 1/6 | 0 |
| -97 | 50/99 | 25/99 | 1729/13662 |

For reproducibility, put z=23971-s and B=-2810. The three numerators whose
residues are enumerated are

\[
A_2=z^2(23971^2+2s\cdot23971+3s^2+B),\quad
A_3=z^3(23971+3s),\quad A_4=z^4.
\]

For each unit class v modulo D, use A_j v^j modulo D; inversion merely
permutes the unit classes, so this equals enumerating prime inverse powers.
The first omitted correction has sign s at levels 2 and 4, and is positive
at level 3. Thus a zero leading residue is bad for a positive correction and
good for a negative correction; a half-denominator residue has the opposite
endpoint rule. These signs were included in the independent enumeration.

The analytic use of those counts is

\[
\begin{aligned}
C_{\rm lin}={}&4\sum_{J\geq8}2^{1-J}\log\frac{J+1}{J}
 +B_2\sum_{J=6}^7 2^{2-J}\log\frac{J+1}{J}\\
&+\frac{B_{23}}4\log\frac98
 +\frac{B_{234}}2\log\frac{16}{15}
 +B_{234}\log\frac54\\
={}&0.1182099412734677\ldots.
\end{aligned}
\]

Thus the combined estimate below y is

\[
C_{\rm quad}+C_{\rm lin}=0.5625676641884855292\ldots.
\]

A union-bound proof of positive density would still require a rigorous bound
strictly less than 0.4374323358115144707 for all unresolved larger-prime failures.
No such bound has been established in this audit. In particular, an elementary
divisor-switching sieve on q(y)=ap sees progression length about X/a; replacing
its available logarithm log(X/a) by log(X^2/a) would be an unjustified gain.

### 5.6 Closing the endpoint ranges and stating the exact partial result

The fixed-margin estimates above imply a statement for every prime p<=y,
not merely for p<=Y^(1-eta). Here the parameter y belongs to
`A_Y={y in [Y,2Y]: y=23971 mod D}`, with D=63987840 fixed, and
`N_Y=|A_Y|=Y/D+O(1)`. Outside the fixed exceptional prime set, a root of any
one factor modulo p occurs among these parameters `N_Y/p+O(1)` times.
There are at most sixteen roots across the four linear and six quadratic
factors. Consequently the proportion of parameters having any such factor
divisible by a prime in the omitted upper strip is at most

\[
\begin{aligned}
\frac1{N_Y}\#\{y\in A_Y:\exists p\in(Y^{1-\eta},3Y]
\text{ dividing a factor}\}
&\leq 16\sum_{Y^{1-\eta}<p\leq3Y}\frac1p
       +O_D\!\left(\frac{\pi(3Y)}Y\right)\\
&=16\log\frac1{1-\eta}+o(1).
\end{aligned}
\]

This tends to zero as eta tends to zero after Y tends to infinity. The ceiling
3Y covers p<=y and every prime divisor of the four linear factors y+-41 and
y+-97 for all sufficiently large Y. Thus the endpoint argument also covers
linear-factor primes exceeding y. The other finitely many exponent boundary
strips at each fixed lower cutoff have vanishing harmonic prime mass by the
same argument.

The small-prime closures apply to this quartic as well. For a fixed prime the
carry-scarcity and valuation-tail argument of Section 4.2 has density zero.
Outside the fixed discriminant, resultant and denominator primes, every root
of F+i is simple, so the Hensel bijection and digit count of Section 2.4 apply
with this fixed progression modulus; the contribution from growing primes
p<=Y^delta tends to zero as delta tends to zero. Repeated nonfixed primes
have cost O(1/P)+o(1) by Section 4.2. At each fixed P and delta only finitely
many bands are used, so no interchange of an infinite sum with a nonuniform
error is required.

In particular, let E_Y consist of the parameters y in A_Y for which

\[
v_p\binom{2F(y)}{F(y)}
 \geq2\,v_p\!\left(\prod_{i=1}^4(F(y)+i)\right)
\]

holds for every prime p<=y and also for every prime divisor of one of the
four factors y-97, y-41, y+41 and y+97. The estimates prove

\[
\liminf_{Y\to\infty}\frac{|E_Y|}{N_Y}
 \geq1-C_{\rm quad}-C_{\rm lin}
 =0.43743233581151447079\ldots.
\]

This is a simultaneous partial valuation statement. Prime divisors exceeding
y in the quadratic factors are still outside its scope.

## 6. A new centered quartic with a nonzero linear term

The loss of level 3 for even quadratic factors suggests changing the quartic,
rather than assuming an unavailable prime-square equidistribution theorem.
Consider P(x)=x^4+Bx^2+Lx+E with L nonzero. A rational quadratic factorization

\[
P(x)+iD=(x^2+u_i x+v_i)(x^2-u_i x+w_i)
\]

forces `u_i(w_i-v_i)=L`, so every quadratic factor has a nonzero linear
coefficient. The exact factorization condition is

\[
v_i=\frac{B+u_i^2-L/u_i}{2},\quad
w_i=\frac{B+u_i^2+L/u_i}{2},\quad
v_iw_i=H(u_i):=\frac{(B+u_i^2)^2-L^2/u_i^2}{4}.
\]

Thus one may choose four distinct positive rational u_i and ask that the four
H(u_i) form an arithmetic progression. After subtracting the common B^2/4,
this is a linear system for B, L^2, the offset and the common difference.
The only additional condition is that the solved positive L^2 be a rational
square. This condition is not rare in a small exact search.

The accompanying `non_even_quartic_search.cpp` checks every ordered quadruple
of distinct positive integers at most 60, identifying a progression with its
reversal. It tests 5,851,620 orderings and finds 5,038 rational-square hits.
Many hits have a local obstruction to integer values, so this count is not a
count of usable integer families.

One admissible hit is u=(5,4,3,1), L=15, B=-1505/32. After setting x=y/8 and
clearing denominators it gives

\[
F(y)=\frac{y^4-3010y^2+7680y+3969}{483840}.
\]

With D=483840, direct multiplication gives

\[
\begin{aligned}
D(F+1)&=(y^2+40y-801)(y^2-40y-609),\\
D(F+2)&=(y-21)(y+53)(y^2-32y-873),\\
D(F+3)&=(y-27)(y+51)(y^2-24y-1057),\\
D(F+4)&=(y^2+8y-1953)(y^2-8y-993).
\end{aligned}
\]

The numerator has 128 roots modulo 512, 8 roots modulo 27, 4 modulo 5 and
4 modulo 7. Hence there are 16,384 admissible residue classes modulo D; for
example y=77 gives the integer F(y)=37. These observations establish
integrality and factorization, not the factorial divisibility condition.

The six irreducible quadratics have discriminant squareclasses
1201,1009,1129,1201,1969,1009; all have nonzero linear coefficients. Their
level-2 leading slopes `4u^2-4v+2B` are respectively
3584,2816,1568,512,2048,-1792, all nonzero. Consequently the averaged
equidistribution argument in Section 5 supplies levels 2,3,4 for all six
quadratics below the parameter height, and their total contribution is

\[
6T_6+\frac38\log\frac65+\frac34\log\frac54
 =0.2872055580594243\ldots.
\]

The four linear-factor residue tables and the genuinely large quadratic
prime divisors require separate analysis. This is another explicit
construction with a stronger partial carry bound, not a completed k=4
infinitude proof.

### 6.1 A second non-even family with eight linear factors

The hit u=(7,5,1,2), L=35, B=-405/8 gives, after setting x=y/4,

\[
F(y)=\frac{y^4-810y^2+2240y-41751}{40320}.
\]

For D=40320, independently expanded exact identities are

\[
\begin{aligned}
D(F+1)&=(y-1)(y-27)(y^2+28y-53),\\
D(F+2)&=(y-9)(y+29)(y^2-20y-149),\\
D(F+3)&=(y-13)(y+9)(y^2+4y-677),\\
D(F+4)&=(y-19)(y+27)(y^2-8y-233).
\end{aligned}
\]

The numerator has respectively 64,4,4,4 roots modulo 128,9,5,7, giving 4096
integral parameter classes modulo D. For instance, F(1)=-1 establishes the
admissible progression y=1 modulo D, on which F(y) tends to positive infinity.
This fact alone makes no assertion about factorial divisibility.

The four quadratic factors have discriminant squareclasses 249,249,681,249.
Their level-2 slopes are 1728,576,1152,-432, all nonzero. All have a nonzero
linear coefficient. Their total contribution below the parameter height is

\[
4T_6+\frac14\log\frac65+\frac12\log\frac54
 =S=0.1914703720396162\ldots.
\]

The construction agent's exact residue optimization nevertheless finds a
minimum very-large-prime bad coefficient of 1 for this family. The linear
cost at the best class with that coefficient is approximately 0.27935.
Consequently the additional splitting does not by itself give a sufficient
total density bound.

### 6.2 Wider exact searches

At u_i<=120 the exhaustive search tests 98,574,840 orderings and returns
25,196 rational-square hits. After removing common scalings of the u_i and
identical quartic coefficients, there are 17,852 distinct candidates. The
numbers with 0,1,2,3,4 split quadratic factors are respectively
17,642; 173; 32; 4; 1. None has five or more split quadratic factors.

The completed extension to u_i<=250 tests 1,906,593,000 orderings and returns
127,846 rational-square hits. There are 87,016 primitive distinct coefficient
tuples, with 0,1,2,3,4 split quadratic factors occurring respectively
86,735; 239; 36; 5; 1 times. Again there is no five-split hit. The additional
three-split representation is equivalent by a variable scaling and a different
partition of an already fully split fiber to an earlier candidate. All
candidates with at least three split quadratics, including their exact rational
coefficients and discriminants, are recorded in
`non-even-search-250-summary.json`.

The receipt can be reproduced without third-party Python packages:

```sh
g++ -O3 -std=c++17 research/k4/non_even_quartic_search.cpp -o /tmp/search-quartics
/tmp/search-quartics 250 > /tmp/quartic-hits-250.txt
python research/k4/summarize_non_even_search.py /tmp/quartic-hits-250.txt --bound 250 --output /tmp/quartic-summary-250.json
```

The more targeted `rational_t_split_search.cpp` tests two rational families:

\[
\begin{array}{c|c|c}
\text{ordered }u&L&B\\ \hline
(1,3,t,5)&15&(-t^6+336t^2+225)/(2t^2(t^2-17))\\
(7,5,1,t)&35&(-t^6-3024t^2+1225)/(2t^2(t^2+23)).
\end{array}
\]

These formulas make the four H(u_i) an arithmetic progression identically.
For all reduced t=p/q with 1<=p,q<=5000, the program tests 15,200,915
fractions in both families. The only hit with at least four split quadratic
factors is t=2 in the second family, giving the D=40320 construction above.
There is no five-split hit in this bounded search. The integer-square tests
are exact: floating-point square roots provide only initial guesses, then
integer multiplication and correction decide equality.

## 7. The supplementary downward 432 family and the remaining scope

The construction agent also found

\[
m=(y^2-14401)/432,\qquad F(y)=m(m+1)/210-4,
\qquad y=16273\pmod{90720}.
\]

Here m=174 modulo 210. With D=210*432^2=39191040, its exact factor pattern is

\[
\begin{aligned}
D(F+1)&=(y^2-216y+9143)(y^2+216y+9143),\\
D(F+2)&=(y^2-23041)(y-73)(y+73),\\
D(F+3)&=(y-143)(y+143)(y-89)(y+89),\\
D(F+4)&=(y^2-14401)(y^2-13969).
\end{aligned}
\]

Thus there are six linear factors, three even irreducible quadratics and two
odd irreducible quadratics. The quadratic part of the lower-prime estimate is

\[
5T_6+\frac12\log\frac65+\log\frac54
 =0.3572020886912154\ldots.
\]

The construction agent's exact six-linear residue computation gives

\[
B_2=242/81,\qquad B_{23}=148/81,\qquad B_{234}=2299/2916,
\]

and a linear contribution 0.30647868120367294. These data are recorded in
`432-linear-refinement-optimization.json.gz`; they have not been independently
re-enumerated in this audit. The level-3 table must retain its signed boundary
rule: at the linear root 73 the leading residue is zero and the correction is
positive, so that coordinate is always bad. Merely replacing B23 by B2/2
would be incorrect.

The combined lower-prime contribution is therefore
0.6636807698948883. The stabilized large-prime residue coefficient is 1/3,
but this coefficient is only an arithmetic classification of bad prime
classes. It is not an upper bound for the number of parameters with such a
prime divisor. This family consequently does not supply a completed bound
either.

### Final scope of the audit

The 552 family remains the main candidate in the accompanying review. Its
linear table and averaged quadratic-root argument have been checked, and the
0.5625676641884855292 lower-prime contribution is defensible with the stated
published inputs and order of limits. To finish a union-bound proof, the
unresolved larger-prime failures would need a rigorously established upper
density strictly smaller than 0.4374323358115144707. None of the searches,
residue tables, individual divisibility certificates, or theorems checked in
this audit supplies that bound. In particular, the O(1) error in a count of
roots modulo each p cannot be summed up to p of order X^2.

## References checked in this audit

1. O. Robert, *On van der Corput's k-th derivative test for exponential sums*,
   Indagationes Mathematicae 27 (2016), 559--589, Theorem 1.
   <https://doi.org/10.1016/j.indag.2015.11.009>
2. D. Keliher and E. S. Lee, *On the constants in Mertens' theorems for primes in
   arithmetic progressions*, arXiv:2306.09981v2 (2024), equation (6), p. 2.
   <https://arxiv.org/pdf/2306.09981>
3. A. van der Poorten, *Fermat's Four Squares Theorem*, arXiv:0712.3850 (2007),
   theorem p. 1 and descent p. 2.
   <https://arxiv.org/pdf/0712.3850>
4. J. Arias de Reyna, *Explicit van der Corput's d-th derivative estimate*,
   arXiv:2407.02094v1 (2024), equations (1)--(2), specialized to d=4.
   <https://arxiv.org/html/2407.02094v1>
5. A. Toth, *Cycle Integrals of Modular Forms* (2018 dissertation),
   Theorem 6.1.3, printed p. 96 / PDF p. 98, reproducing the quadratic-root
   theorem for prime arithmetic progressions.
   <https://real-d.mtak.hu/1120/7/dc_1553_18_doktori_mu.pdf>
6. W. Duke, J. B. Friedlander and H. Iwaniec, *Equidistribution of roots of a
   quadratic congruence to prime moduli*, Annals of Mathematics 141 (1995),
   423--441.
   <https://annals.math.princeton.edu/1995/141-2/p08>
7. T. D. Wooley, *Translation invariance, exponential sums, and Waring's
   problem*, Proceedings of the ICM 2014, Theorem 7.1 (classical Weyl bound).
   <https://www.math.purdue.edu/~twooley/publ/2014icm.pdf>
