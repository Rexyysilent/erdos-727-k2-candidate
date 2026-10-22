# Exact linear-factor analysis for the 552 family

22 October 2026. These results establish the linear-factor portion of a
possible k=4 argument. They do not control every prime divisor of the six
remaining irreducible quadratics, and hence do not establish infinitude.

## 1. Family and exact Taylor coefficients

Put

\[
D=210\cdot552^2=63987840,\qquad Y=23971,\qquad y\equiv Y\pmod D,
\]
\[
m=(y^2-1681)/552,\qquad n=F(y)=m(m+1)/210-3
 =\frac{y^4-2810y^2-190065671}{D}.
\]

The congruence ensures that m is an integer with m≡195 (mod210). The four
linear factors of the product of the four consecutive shifts are y−r for
r∈{−97,−41,41,97}. We have F(±97)=−2 and F(±41)=−3. If y=r+ap, then

\[
F(y)=\frac{a^4p^4+4ra^3p^3+E_2a^2p^2+E_1ap}{D}-i,
\tag{1}
\]

where i=2 for |r|=97 and i=3 for |r|=41, and

\[
E_2=6r^2-2810>0,\qquad E_1=4r^3-5620r.
\]

The magnitudes of (E₂,E₁) are (7276,45264) when |r|=41 and
(53644,3105552) when |r|=97. The sign of E₁ is the sign of r, and all these
coefficients are nonzero. This nonvanishing is essential to the Fourier
argument below.

Throughout this note p is coprime to D. A prime dividing a linear factor has
n≡−i (modp), so its first carry is present once p>2i. Outside a fixed finite
set of primes, different factors have no common prime divisor. Primes p≥X^δ
whose square divides a linear factor contribute o(X) exceptional y in an
interval [X,2X], for each fixed δ>0, by summing X/p²+O(1). Thus, in the
linear-factor analysis, it suffices to supply one further carry.

## 2. Stabilized carries at levels 2, 3 and 4

Let u be the inverse of p modulo D and define

\[
\begin{aligned}
K_2(r)&=(Y-r)^2(Y^2+2rY+3r^2-2810)\pmod D,\\
K_3(r)&=(Y-r)^3(Y+3r)\pmod D,\\
K_4(r)&=(Y-r)^4\pmod D.
\end{aligned}
\]

Since a≡(Y−r)u (modD), equation (1) gives the following limiting
fractional parts. For y∈[X,2X], the error estimates are uniform.

| Carry level j | Main fraction | Error size | Error sign for large X | Stable when |
|---|---|---|---|---|
| 2 | (K₂u² modD)/D | O(X/p²) | sign(r) | p>X^(1/2+η) |
| 3 | (K₃u³ modD)/D | O(X²/p³) | positive | p>X^(2/3+η) |
| 4 | (K₄u⁴ modD)/D | O(X³/p⁴) | sign(r) | p>X^(3/4+η) |

For example, the correction at level2 is E₁a/(Dp)−i/p². At level3 it is
E₂a²/(Dp)+E₁a/(Dp²)−i/p³. At level4 it is
4ra³/(Dp)+E₂a²/(Dp²)+E₁a/(Dp³)−i/p⁴. These formulas establish both the
error sizes and their signs. Since the leading fractions have fixed
denominator D, a nonboundary value is separated from both endpoints of
[0,1/2) by at least 1/D. The stated exponent margins therefore stabilize
the carry decisions.

The correct endpoint rule is the following. If s is the leading numerator
modulo D and σ is the sign of the vanishing correction, absence of a carry
means

\[
\operatorname{bad}(s,\sigma)=
\begin{cases}
1_{\sigma>0},&s=0,\\
1_{\sigma<0},&s=D/2,\\
1_{2s<D},&0<s<D,\ s\ne D/2.
\end{cases}
\]

The chosen Y has no endpoint at any of the three levels. Indeed none of the
listed K values is 0 or D/2, and multiplication by a unit cannot turn a
different residue into either of these two residues.

| r | K₂ | K₃ | K₄ |
|---:|---:|---:|---:|
| −97 | 6628416 | 7326080 | 44941696 |
| −41 | 63073728 | 1828224 | 34736256 |
| 41 | 4625760 | 8534960 | 54724240 |
| 97 | 19958112 | 36793008 | 51875856 |

## 3. Exact joint residue table

There are φ(D)=13989888 prime unit classes. The following proportions are
exact counts over these classes; they are not random-independence assumptions.

| r | Bad at 2 | Bad at 3 | Bad at 4 | Bad at both 2,3 | Bad at all 2,3,4 |
|---:|---:|---:|---:|---:|---:|
| −97 | 50/99 | 1/2 | 1139/2277 | 25/99 | 1729/13662 |
| −41 | 1/3 | 1/2 | 0 | 1/6 | 0 |
| 41 | 49/99 | 1/2 | 1139/2277 | 49/198 | 94/759 |
| 97 | 1/3 | 1/2 | 1/3 | 1/6 | 0 |
| Sum | **5/3** | 2 | 3037/2277 | **5/6** | **311/1242** |

Denote the three bold totals by B₂, B₂₃ and B₂₃₄.
The full enumeration is in `construction_linear_residue_table.cpp`, with
output `linear-residues-23971.json`. A second enumeration in
`construction_optimize_linear.cpp` uses a smaller exact space:

* D=2⁷·3³·5·7·23² has 218592 distinct unit-square classes, each with 64
  preimages in its unit group.
* The decisions at levels2 and4 depend only on u². The involution u↦−u
  leaves those two decisions unchanged and interchanges the level3 decision
  with its complement, because level3 has no endpoint.
* Consequently B₂₃=B₂/2 and B₂₃₄=B₂₄/2. Enumerating the 218592 square
  classes therefore verifies the full joint table, with a factor64 reduction
  in the enumeration size. The implementation explicitly handles hypothetical
  K₃=0 or D/2 cases as well.

The prime number theorem for fixed arithmetic progressions modulo D weights
the unit classes uniformly in a logarithmic prime sum. No assumption about
independence of the carry decisions is needed for this step.

## 4. Equidistribution of the other carry coordinates

We give details because a generic high-derivative bound is unnecessarily weak
here. Restrict to one root r and one prime p, and write

\[
a=a_0+Db,\qquad a_0\equiv(Y-r)p^{-1}\pmod D.
\]

The relevant b interval has length N≈X/(Dp). For fixed integer frequencies
hⱼ, apply Weyl's criterion to the polynomial phase

\[
\sum_jh_j\frac{F(r+ap)}{p^j}.
\]

Every term cₑaᵉp^(e−j)/D with e≥j has integer-valued nonconstant variation
as b varies: each nonconstant binomial term in (a₀+Db)ᵉ has a factor D.
Those terms may therefore be removed from the exponential sum, leaving a
constant factor of absolute value1. The highest remaining coefficient is as
follows.

| Nonzero frequency of greatest relevant type | Degree in b | Leading coefficient |
|---|---:|---|
| Some hⱼ≠0 with j≥5; let j₀ be the smallest such j | 4 | D³ Σ(j≥5)hⱼp^(4−j) ∼ D³hⱼ₀p^(4−j₀) |
| No such frequency, but h₄≠0 | 3 | 4rD²h₄/p |
| Only levels2,3 remain and h₃≠0 | 2 | E₂Dh₃/p |
| Only h₂≠0 remains | 1 | E₁h₂/p |

For degrees d≥2, the leading coefficient α tends to0 with a nonzero fixed
leading constant. Put q=round(1/|α|) and choose numerator sign(α). Then
gcd(sign(α),q)=1 and |α−sign(α)/q|≤q⁻² for all sufficiently large p.
Classical Weyl's inequality gives

\[
\left|\sum_{b\in I}e(P(b))\right|
 \ll N^{1+\varepsilon}
       (q^{-1}+N^{-1}+qN^{-d})^{1/2^{d-1}}.
\tag{2}
\]

An explicit reference for (2) is Trevor D. Wooley, *Translation invariance,
exponential sums, and Waring's problem*, Proceedings of the ICM2014,
Theorem7.1, [author's PDF](https://www.math.purdue.edu/~twooley/publ/2014icm.pdf),
also [arXiv:1404.3508](https://arxiv.org/abs/1404.3508).

For the quartic case q≈p^(j₀−4), hence q/N⁴≈p^j₀/X⁴. For the cubic and
quadratic cases the corresponding quantities are p⁴/X³ and p³/X².
The linear case is a geometric sum with relative bound O(p²/X). Therefore,
if p≥X^δ and fixed positive exponent margins are retained, there is a
uniform power saving in every nonzero fixed Fourier mode under these
respective upper bounds:

\[
p<X^{4/J-\eta}\quad(j\le J,\ j\ge5),\qquad
p<X^{3/4-\eta}\quad(j=4),\qquad
p<X^{2/3-\eta}\quad(j=3),\qquad
p<X^{1/2-\eta}\quad(j=2).
\]

Choose the ε in (2) after δ and η, sufficiently small that N^ε is absorbed
by the resulting fixed power saving. There is no loss at this step. With
fixed Fourier truncation the frequencies are fixed; a finite collection of
these uniform bounds proves joint equidistribution. Then let the truncation
tend to infinity. The coefficients at lower degrees are unrestricted in
Weyl's inequality, so variation of a₀ or of the interval's starting point
does not weaken the argument.

This proves the joint distribution needed in each band below. In particular,
it does not assert that stabilized coordinates themselves are equidistributed;
their exact dependence on the prime residue was already counted in Section3.

## 5. Resulting loss from the four linear factors

Write α=log(p)/log(X). Ignore fixed narrow strips around band endpoints,
then let their widths tend to zero. The bands and their valid joint bounds are

| Prime exponent α | Stabilized coordinates required to be bad | Equidistributed coordinates required to be bad | Summed weight over four roots |
|---|---|---|---:|
| 4/(J+1)<α<4/J, J≥8 | None | 2,…,J | 4·2^(1−J) |
| 4/(J+1)<α<4/J, J=6,7 | 2 | 3,…,J | B₂·2^(2−J) |
| 2/3<α<3/4 | 2,3 | 4,5 | B₂₃/4 |
| 3/4<α<4/5 | 2,3,4 | 5 | B₂₃₄/2 |
| 4/5<α<1 | 2,3,4 | None | B₂₃₄ |

For each prime, the linear congruence in y has X/(Dp)+O(1) points in an
interval of length X. Summing the O(1) errors for p≪X costs O(X/logX),
which is negligible compared with X/D, since D is fixed. Summing the uniform
Fourier errors is also negligible. The prime reciprocal sum on
X^α<p<X^β is asymptotic to log(β/α), with the specified residue weights.
The resulting upper bound for the union of the linear-prime failures is

\[
\begin{aligned}
C_{\rm lin}={}&4\sum_{J\ge8}2^{1-J}\log\frac{J+1}{J}\\
&+\frac53\sum_{J=6}^{7}2^{2-J}\log\frac{J+1}{J}
 +\frac{5}{24}\log\frac98\\
&+\frac{311}{2484}\log\frac{16}{15}
 +\frac{311}{1242}\log\frac54\\
={}&\boxed{0.1182099412734676820956881759066784995989\ldots}.
\end{aligned}
\]

The five contributions in order are 0.006703660166582722,
0.023012122514533293, 0.024538132428413220, 0.008080305987836004 and
0.055875720176102443. This is a union-bound contribution for the linear
factors. A complete infinitude argument still needs the fixed-small-prime
conditions and all quadratic-prime ranges.

## 6. Exhaustive optimization of the progression refinement

There are 64 residues y modulo115920 satisfying
y²≡1681+552·195 (mod115920). Among them, the smallest combined stabilized
bad coefficient for the new two quadratic factors of n+4 is 1/3. It is
attained in exactly eight coarse classes:

23971, 26491, 31469, 33989, 81931, 84451, 89429, 91949.

Each class has 552 refinements modulo D. The second exact enumerator checks
all 4416 such refinements, using the 218592 square classes. Exactly 2208
refinements have C_lin=0.11820994127346769…, and the other 2208 have
C_lin=0.19015166959828658…. Thus Y=23971 already minimizes this linear
bound among all refinements that minimize the new quadratic high-prime
coefficient. This is the exact finite optimization performed; it is not a
claim of global optimality over other quartic families.

The complete receipt is `linear-refinement-optimization.json.gz`. Reproduce with

```sh
g++ -O3 -std=c++17 -fopenmp research/k4/construction_optimize_linear.cpp -o /tmp/construction_optimize_linear
OMP_NUM_THREADS=8 /tmp/construction_optimize_linear | gzip -n > research/k4/linear-refinement-optimization.json.gz
```

## 7. Independent arithmetic verification

`construction_verify_linear.py` verifies all four Taylor identities
symbolically, the 12 signed endpoint cases by rational arithmetic, and 4608
stabilized carry predictions at 128 large coprime integer bases. The test
includes cofactor sizes adapted to each of the three distinct stabilization
thresholds. Its output is `linear-signature-verification.json`. These
identities and inequalities apply to every coprime integer base; primality is
needed only in the subsequent counting argument.

`construction_552.py` reconstructs every factorization and verifies factorial
valuations against explicit digit addition. For primes above2⁶⁴ it supplies a
recursive Lucas primality proof; the witness output optionally includes such
proofs for every prime. The finite scans give

| Progression | Interval | Successful parameters | First success |
|---|---|---:|---|
| y=23971+115920t | 0≤t<1000 | 332 | t=4, n=883769523993343 |
| y=23971+63987840t | 0≤t<100 | 31 | t=2, n=4195055316568514657639543 |
| y=23971+63987840t | 0≤t<1000 | 315 | t=2, n=4195055316568514657639543 |

The exact refined witness, with all Lucas prime certificates, is
`552-witness-refined-t2.json`. These counts are finite observations and imply
neither a limiting density nor infinitely many solutions.
