# Review and extensions: Erdős problem 727, k = 2

Prepared 17 September 2026. Basis: the supplied 12-page PDF dated 16 September
2026 and the older `Rexyysilent/erdos-727-k2-candidate` repository at
`3dec50a4d4d08d551aabb910368b849d2648a4d0`.

## Assessment and scope

I found no fatal gap in the supplied PDF's argument after checking its algebra,
fixed-prime reduction, simple-root digit count, exponential-sum estimates,
Fourier-to-box passage, prime-band coverage and order of limits against the
stated hypotheses. The proposed lower bound survives this review. This is an
AI-assisted mathematical assessment, not independent peer review, a formal proof
certificate, or a claim of priority. The original manuscript and README are
preserved. Finite calculations are reported separately from the analytic proof.

The result is specific to `F(t)=6t²+11t+3` and k=2. Nothing here proves the
assertion for every fixed k. No assumption about whether a problem website labels
the problem open or solved is used to assess the argument.

Three additional mathematical consequences are developed below:

1. A particular large-prime obstruction has exact limiting density `(log 3)/2`.
   Hence the upper limiting success proportion of this family is at most
   `1-(log 3)/2 = 0.4506938556…`.
2. The mean number of simple-prime obstructions tends to `1-rho`. The excess
   multiplicity of those obstructions exactly identifies the slack in the
   original union bound, up to a vanishing remainder.
3. The proposed dyadic lower bound implies at least
   `(rho/sqrt(6)-o(1))*sqrt(N)` solutions n≤N from this family.

These do **not** increase the numerical lower bound rho or establish that the
success density exists. They use the analytic inputs checked in the PDF.

A separate companion, `GENERALIZED_FAMILIES_2026-09-17.md`, proves proposed
stronger bounds in different quadratic families and a ten-family counting
corollary. It does not replace the original-family lower bound reviewed here.

## 1. Review of the supplied argument

### Sections 2 and 3: exact reduction and factorization

Cancelling `(n!)²` gives the necessary and sufficient criterion

`((n+2)!)² | (2n)!  <=>  ((n+1)(n+2))² | binom(2n,n)`.

For each denominator prime, the available valuation is

`C_p(n) = sum_{j>=1} 1_{ 2*(n mod p^j) >= p^j }`.

This equals the number of outgoing carries in the explicit base-p addition
`n+n`. It is not a probabilistic approximation. The four linear factors of
`D(F(t))` multiply exactly to `(F(t)+1)(F(t)+2)`. Their m-shifts differ by at most
five, so p≥7 cannot divide two different shifts.

In formula (3.4), the fractional parts of `a²/6` are respectively
`1/6, 0, 1/2, 2/3` for `c=0,1,-2,3`. The small correction is negative at c=1
and positive at c=-2. Keeping these signs is essential; they correctly give a
second carry for those two boundary cases.

### Section 4: fixed primes and repeated factors

The primitive-linear-factor bound on denominator valuations is valid also for
p=2 and p=3: a factor with nonunit leading coefficient has no root because its
constant coefficient is then a unit. The carry-count estimate is an upper bound,
which is all that is needed. For p=2,3, the polynomial is a permutation modulo p
with unit derivative and therefore lifts bijectively modulo every power.

For p≥5, after excluding `p^H | 12t+11`, write `12t+11=p^s u`, s<H. When
L>2H the valuation of the square fixes s; the unit square congruence has at most
two roots modulo `p^(L-2s)`, each with `p^s` lifts in the required modulus.
Thus the claimed `2p^H` fiber bound is conservative and sufficient.

For p>P, repeated factors contribute `O(X/P+sqrt(X))`. The error is summed over
prime-square divisibility incidences. Independence of the four linear factors
is neither asserted nor required. The later limit P→infinity
removes it after fixed primes have been handled.

### Sections 5 and 6: small-prime summation

At the four roots, `F'(t)` is `5,-5,-1,1` modulo p, so p≥7 gives simple roots.
The lifting map is bijective on each root class. Taking
`L=floor(log_p(X/12))` ensures that the interval count per residue is bounded by
an absolute multiple of `X/p^L`, including its endpoint error.

The weight `w(y)=y^-1 exp(-c0 log(X)/log(y))` need not be monotone. This is not
a gap: the proof uses the one-sided inequality
`-w'(y) <= y^-2 exp(-c0 log(X)/log(y))`, multiplied by the nonnegative prime
counting function. The resulting integral is bounded by
`delta*exp(-c0/delta)/c0`, with `c0=log(7/4)`.

The proof actually bounds the **sum of bad prime incidences**, a stronger fact
than the union bound stated in Lemma 6.1. This stronger form is used in Section 4
below and follows directly from the displayed per-prime estimate (6.2).

### Section 7: uniform distribution

For only h₂ nonzero, subtracting the integer-valued quadratic variation leaves
slope `h₂*d_c/p`, whose distance to the nearest integer is `|h₂*d_c|/p` once
p>2|h₂*d_c|. The normalized geometric-sum bound is `O_h(p²/X)`.

For a nonzero higher coordinate, let j₀≥3 be the first such coordinate. The
quadratic second derivative is

`12 sum_{j=3}^J h_j p^(2-j) = 12 h_{j0} p^(2-j0)(1+O_h(1/p))`.

No later term can cancel the leading fixed-frequency term for sufficiently large
p. The second-derivative test gives

`O_h(p^(-(j0-2)/2) + p^(j0/2)/X)`.

The stated prime bounds make both terms uniformly small. In particular, replacing
j₀ by J in the second term is valid because p≥1 and j₀≤J.

Uniform Fourier convergence alone should not be described as probabilistic digit
independence. Lemma 7.1 supplies the needed step: fixed upper and lower continuous
approximations to a box, followed by finite trigonometric approximation, give
uniform box counts. The finite frequency set is chosen before X tends to
infinity. The exact relation between consecutive digit levels has a frequency
that grows with p and does not contradict this argument.

Lemma 7.4 is the one-dimensional instance with second derivative `12h/p` and
normalized error `O_h(p^-1/2+p^(3/2)/X)`.

### Sections 8 through 10: coverage and limits

Above `X^(1/2+eta)`, only the c=0 factor can obstruct. The reciprocal-prime costs
above and below the `2/3` boundary add to `(log 3)/2`. Below the square-root
scale, the intervals `(2/(J+1),2/J)` cover all exponents except their endpoints.
The explicit finite boundary set also covers the clipped delta interval.

For fixed delta, eta and J, the conditional errors are uniform in p. Summing
`L_{p,c}` contributes the reciprocal-prime sum plus `O(pi(2X+3))/N_X=o(1)`.
The limit order is legitimate:

`X→infinity; then eta→0; then delta→0; then P→infinity`.

There is no hidden requirement for estimates uniform in unbounded J. The series
tail bound `2^(3-R)/(R+1)` is correct. Positivity is independent of decimal
computation; the Taylor lower bound used for exp(11/10) is
`36015101/12000000 > 3`.

### References checked

The original question is present on printed page 90 of Erdős, Graham, Ruzsa and
Straus (1975). Robert's Theorem 1 is the second-derivative estimate with a constant
depending only on the derivative-comparability factor, which is fixed here.
Goldmakher's page-2 proposition gives the reciprocal-prime estimate used in (5.1).
The stronger theorems of Ford and Konyagin are not needed for this proof.

## 2. New proposition: an exact large-prime obstruction density

Use the PDF's notation `M_X={m in [X,2X]:m=5 mod 6}`, `N_X=#M_X`,
`n=m(m+1)/6-2`. Define

`H_X={m in M_X: some prime p>sqrt(2X+3) divides m and C_p(n)=1}`.

**Proposition.** As X tends to infinity,

`#H_X/N_X → (log 3)/2`.

Consequently, writing g_X for the success proportion in M_X,

`limsup g_X <= 1-(log 3)/2`.

**Proof.** At most one prime exceeding sqrt(2X+3) can divide any m≤2X. Its
valuation in m is one and it divides no other shift; thus `v_p(D(n))=1`.
These prime incidences are disjoint, not merely bounded by a union estimate.

Write m=ap. Both a and p are units modulo six, so `a²/6` has fractional part
1/6. Since `a/p=m/p²<1`, for sufficiently large X,

`0 < {n/p²} = 1/6 + a/(6p) - 2/p² < 1/3 < 1/2`.

Also `n/p^4 < 1/2`, and therefore every level j≥4 has no carry. The initial
carry is forced. Consequently the exact obstruction criterion in this range is

`C_p(n)=1 <=> {n/p³}<1/2`.

Fix a small positive eta. In
`X^(1/2+eta) <= p <= X^(2/3-eta)`, Lemma 7.4 gives a conditional frequency
`1/2+o_eta(1)`, uniformly in p. Summation gives

`(1/2) log((2/3-eta)/(1/2+eta)) + o_eta(1)`.

For `p>X^(2/3+eta)`, the inequality p³>2n holds uniformly for sufficiently
large X. Every such incidence is an obstruction. Their density contribution is

`log(1/(2/3+eta)) + o_eta(1)`.

The omitted interval from sqrt(2X+3) to `X^(1/2+eta)` and the band around
`X^(2/3)` have total incidence density `O(eta)+o_eta(1)` by Mertens' estimate.
All relevant p are at most 2X. Let X tend to infinity first and eta decrease to
zero afterward. Disjointness turns these two contributions into the exact limit

`(1/2)log(4/3)+log(3/2)=(1/2)log 3`.

Every member of H_X violates the original divisibility criterion, proving the
upper bound. QED.

This result uses only the c=0 higher-carry distribution, the exact identities,
and reciprocal-prime summation. It does not depend on an independence heuristic
or on improving the lower-band estimates.

## 3. Two-sided consequence for the same quadratic family

Combining the proposed lower theorem with the proposition above gives

`rho <= liminf g_X <= limsup g_X <= 1-(log 3)/2`.

Numerically, these endpoints are

`0.2592234836263289458542163785…` and
`0.4506938556659451543023773815…`.

These are lower and upper limits, not a proof that a density exists. The upper
bound concerns this family only, not all solutions of the factorial question.

## 4. New proposition: the exact first moment and the union-bound slack

For m in M_X let

`W(m) = #{primes p>=7: p divides D(n), C_p(n)=1}`.

Each counted prime necessarily has denominator valuation one, because higher
valuation forces at least two initial carries.

**Proposition.** The average of W over M_X tends to

`tau = (log 3)/2 + 4 sum_{J>=4} 2^(1-J) log((J+1)/J) = 1-rho`.

**Proof.** In a lower prime band
`X^(max(delta,2/(J+1))+eta) < p < X^(2/J-eta)`, not only are the levels 2,...,J jointly distributed as stated in Lemma 7.2, but
there can be no carry at any level above J: uniformly,

`n/p^(J+1) = O(X^(-(J+1)*eta)) < 1/2`.

Thus the no-carry box in (7.2) counts `C_p(n)=1` exactly, rather than just
bounding it above. At fixed delta, each of the four shifts contributes, as eta
decreases to zero,

`2^(1-J) log((2/J)/max(delta,2/(J+1)))`, for `2/J>delta`.

The bottom band may be clipped: its contribution must not be replaced by the
full-band value in an equality. Each clipped term is bounded by
`2^(1-J) log((J+1)/J)` and tends to that value for every fixed J as delta decreases
to zero. The convergent dominating series justifies passage to the full sum.
Summation is linear and requires no independence between different primes.

Outside the vanishing exponent band about 1/2, Section 8 of the supplied
PDF eliminates the other three shifts. The large-prime contribution therefore
tends to `(log 3)/2` by Section 2 above; boundary strips can be included or omitted.
The deterministic exclusion of the other shifts is not being extended without
a margin all the way down to sqrt(2X+3). The per-prime proof of Lemma 6.1
bounds the remaining small-prime incidence sum by a quantity tending to zero
with delta. For fixed delta there are only finitely many lower bands. Take X
to infinity, then eta to zero, then delta to zero. This proves the limit. QED.

Now let r_X be the proportion of bad parameters with W(m)=0. Such a failure
arises at a prime p≤5 or at a prime with denominator valuation at least two.
For any fixed P, the former and other p≤P failures have density o_P(1); the
remaining repeated-prime set has upper density O(1/P). Letting P tend to
infinity proves r_X→0.

The elementary identity `1_{W>=1}=W-(W-1)_+` therefore gives

`g_X = rho + average_{M_X}(W-1)_+ + o(1)`.

This isolates the target for a stronger lower bound: prove a positive lower bound
for **excess obstruction multiplicity**. Finite overlap counts alone are not
such a proof. No conjectural overlap estimate is inserted into rho here.

## 5. Counting solutions up to N

A dyadic lower limit at least rho implies the same lower limit for prefix counts
of good parameters. Fix epsilon>0. Above a fixed threshold all dyadic intervals
have at least `(rho-epsilon)` times their length in good parameters. Decompose a
long prefix into dyadic intervals, discard its bounded initial piece, and absorb
the O(log T) endpoint overlaps. The same argument transfers the upper bound.

Since `F(t)=6t²+O(t)` is increasing, its inverse cutoff is

`floor((sqrt(24N+49)-11)/12) = sqrt(N/6)+O(1)`.

The number S_F(N) of successful values n≤N **from this family** consequently
satisfies

`liminf S_F(N)/sqrt(N) >= rho/sqrt(6) = 0.1058275440385359781503532597…`,

`limsup S_F(N)/sqrt(N) <= (1-(log 3)/2)/sqrt(6) = 0.1839949960981891318246561220…`.

In particular the lower bound also applies to the number of all successful n,
since the family gives distinct values. The upper bound does not transfer to all n.

## 6. Finite computations, separate from the proofs

All intervals below are half-open. The original GitHub verifier was reconstructed
byte-for-byte, and its Git blob hash was checked before execution. It reproduced
all three counts stated in the PDF. An independently implemented C++ scanner
agreed and extended the ranges.

| t interval | Parameters checked | Valid parameters | Fraction |
|---|---:|---:|---:|
| [0, 1,000) | 1,000 | 120 | 0.1200000 |
| [10,000, 20,000) | 10,000 | 1,879 | 0.1879000 |
| [100,000, 200,000) | 100,000 | 21,716 | 0.2171600 |
| [1,000,000, 2,000,000) | 1,000,000 | 234,157 | 0.2341570 |
| [10,000,000, 20,000,000) | 10,000,000 | 2,483,873 | 0.2483873 |

The final sample remains below rho. This does not contradict an asymptotic lower
limit with no effective threshold. It does show why these sample fractions must
not be advertised as a verification of rho.

Additional checks: literal factorial divisibility for every 0≤n≤500; agreement
between Legendre sums, fractional-part indicators and explicit digit carries;
2,000 independent high-parameter Python/C++ comparisons; small prime-power
permutation, fiber and lift checks; exact algebraic identities. The companion
suite has 38 passing tests. Those finite checks do not verify the analytic limits.

For the 10-million-parameter final interval, exact obstruction accounting gives:
8,987,644 simple-prime incidences; 1,495,440 excess incidences; and 23,923 bad
parameters without a simple-prime obstruction. Accordingly,

`10,000,000 - 8,987,644 + 1,495,440 - 23,923 = 2,483,873`.

The large-obstruction event occurred 6,467,619 times. Its finite proportion,
0.6467619, is not the limiting value `(log 3)/2`. None of the convergence rates
or implied constants needed for a finite error bar has been made effective here.

The constants file uses rational upper/lower bounds for logarithms, not merely
floating-point evaluation: `log x=2*sum z^(2r+1)/(2r+1)`, z=(x-1)/(x+1), with
an explicit geometric tail. The infinite J-series has its own certified tail.
Displayed decimal endpoints are rounded outward.

## 7. Repository update and next research branches

The remote has four files: the September 7 note, README, original verifier and
September 12 results. The supplied PDF is an expanded exposition, not a different
quadratic family or a numerically stronger lower bound. Its main improvements
are the explicit partial summation, the uniform box lemma, clipped prime bands,
complete error inequality and clearer limit order.

The PDF's Section 11 mentions `prism-uploads/verify_727_k2(1).py`. That path is
not in the GitHub tree; the repository's real entry point is `verify_727_k2.py`.
This is a packaging/reference correction, not a mathematical defect. The copied
PDF remains unchanged; this review supplies the cross-reference.

Suggested separate branches, each with a clear deliverable:

- **Proof exposition and review:** incorporate the September 16 exposition and
  this addendum with a lemma dependency map and independent reviewer comments.
- **Overlap estimates:** target a positive asymptotic lower bound on the mean of
  `(W-1)_+`; do not replace this by a product of presumed independent probabilities.
- **Quantitative discrepancy:** make Fourier truncation, derivative constants and
  boundary choices effective before claiming any finite threshold for rho.
- **Alternative families:** search factorizable consecutive quadratic values and
  analyze the residue/sign table first. A numerical family search is not a proof
  and gives no automatic extension to k≥3.
- **Formal arithmetic core:** formalize cancellation, the four polynomial
  identities and the carry identity separately from the analytic theorems.

> **Publication update, later on 17 September 2026.** Authenticated GitHub write
> actions subsequently became available and the additive continuation was published.
> The paragraph below records the environment at the time this review was prepared.

No automated email, public status claim, GitHub push or PR was made during the
review-preparation run. The connector exposed no write action at that time, and the
local environment had no authenticated GitHub CLI route. The additions were therefore
first assembled as a reviewable patch before later publication.

## Sources

- Supplied PDF: *A quadratic-family approach to the k=2 case of Erdős problem 727*,
  16 September 2026, Sections 2–10 and computational Section 11.
- Older source note and verifier, pinned tree:
  https://github.com/Rexyysilent/erdos-727-k2-candidate/tree/3dec50a4d4d08d551aabb910368b849d2648a4d0
- P. Erdős, R. L. Graham, I. Z. Ruzsa, E. G. Straus, *On the prime factors of
  binom(2n,n)*, Mathematics of Computation 29 (1975), 83–92, printed p. 90:
  https://renyi.hu/~p_erdos/1975-27.pdf
- O. Robert, *On van der Corput's k-th derivative test for exponential sums*,
  Indagationes Mathematicae 27 (2016), 559–589, Theorem 1.
  DOI: 10.1016/j.indag.2015.11.009. The publisher search extract supplied the
  theorem statement; direct article fetch was blocked.
- L. Goldmakher, *A quick proof of Mertens' theorem*, p. 2:
  https://web.williams.edu/Mathematics/lg5/mertens.pdf


## Companion updates made later in this continuation

The separate generalized-family note develops stronger lower bounds in other
families and a finite-union counting bound. A separate k=3 note supplies the
new six-factor argument for `210t²+391t+179`. Those proposed results were not
asserted in the supplied PDF and do not raise the original family's lower
constant. The combined arithmetic suite contains 91 passing tests (38 for the
original family, 35 for generalized k=2 families, and 18 for k=3).
