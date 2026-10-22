# A smoothness reduction for the k = 4 continuation

22 October 2026. This note proves a reduction, not the k = 4 infinitude assertion.

## 1. The relevant published input

J. Teräväinen, *On the Liouville function at polynomial arguments*, Definition
2.8, Proposition 2.10, and the proof of that proposition (printed pp. 5 and
27–28), records the following consequence of G. Harman's smooth-quadratic
theorem: an irreducible quadratic takes values whose prime factors are at most
`t^(4/5 + epsilon)` on a set of positive lower density. The proof explains the
extension from `t^2+1` to `t^2-D`, then to arbitrary irreducible quadratics, and
the insertion of a fixed congruence condition. Its exported Property S says
that the lower-density constant for `t`-smooth values can be chosen independently
of each fixed progression modulus and residue.

Primary text:

- https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-teravainen-FINAL.pdf
- G. Harman, *On values of n^2+1 free of large prime factors*, Archiv der
  Mathematik 90 (2008), 239–245 (Teräväinen's reference 20).

All congruence moduli are fixed before the parameter tends to infinity. This
input gives no uniform estimate for a modulus growing as `p^j` with the scale.
It also gives no simultaneous `t^(4/5+epsilon)` smoothness theorem for eight
distinct quadratic polynomials.

The argument below only needs the explicitly exported Property S, so it does
not require a stronger uniform version of the `4/5+epsilon` assertion.

## 2. One quadratic times finitely many linear factors has Property S

Let `Q(t)` be a fixed irreducible quadratic with positive leading coefficient,
and let `L_1(t),...,L_r(t)` be fixed nonconstant integral linear polynomials.
Assume their values are nonzero at all sufficiently large positive parameters.
Put

```
P(t) = Q(t) * product_i L_i(t).
```

Then `P` has Property S. In particular, there is `eta > 0` such that for every
fixed positive `q` and every fixed residue `b`,

```
liminf_{X -> infinity} (q/X) * # {1 <= t <= X :
    t = b (mod q), P^+(|P(t)|) <= t} >= eta.
```

Here `P^+(M)` denotes the largest prime factor of `M`, with `P^+(1)=1`.

**Proof using the exported theorem.** Let `eta_Q > 0` be the constant in
Property S for `Q`. For any fixed `gamma > 0`, remove `t < gamma X` at a
relative cost at most `gamma + o(1)`. On the remaining interval, any prime
factor of a linear value exceeding `t` belongs to `(gamma X, C X]`, for a
fixed `C`. For all sufficiently large `X`, these primes are coprime to the
fixed slope and to `q`; the simultaneous congruences `t = b (mod q)` and
`p | L_i(t)` form one class modulo `pq`. Therefore the number of these
incidences is at most

```
(X/q) * sum_{gamma X < p <= C X} 1/p + O(pi(C X)) = o(X/q).
```

Both the reciprocal-prime sum and `q*pi(CX)/X` tend to zero for fixed
`gamma,C,q`. Taking `gamma < eta_Q/2` proves the assertion, for example with
any `eta < eta_Q/2`. No assertion about independence between the factors is
used. This proof does not require the stronger `4/5+epsilon` smoothness level.

**Stronger optional version.** If the fixed-progression extension of Harman's
stronger statement is used, fix any `alpha` with `4/5 < alpha < 1` and let
`delta > 0` be its positive relative lower-density constant. One can then
choose `beta < 1`, as close to 1 as necessary and greater than `alpha`, such
that all factors of `P(t)` are simultaneously `t^beta`-smooth on a set of
positive relative lower density. After removing `t < gamma X`, the linear
bad-prime incidence cost is

```
r * log(1/beta) + o(1)
```

relative to the progression. Choose `gamma + r*log(1/beta) < delta/2`.
This stronger version uses the extension discussed in Teräväinen's proof;
the unqualified reduction below uses only Property S.

## 3. Application to the repository's quadratic

For

```
H(t) = 210t^2 + 391t + 179,
P4(t) = product_{i=1}^4 (H(t)+i),
```

the first three shifted values have the six linear factors already in the
repository. The fourth is

```
Q4(t) = 210t^2 + 391t + 183,
discriminant(Q4) = -839.
```

Thus `Q4` is irreducible over the rationals, and Section 2 applies to `P4`.
There is a positive lower-density set of parameters on which **every prime
dividing the k=4 denominator is at most the parameter t**. This is a proved
smoothness statement, with the published theorem as input. It is not a
factorial-divisibility statement.

Fix a progression `b mod q`. Work on an interval `[gamma X,X]`, where
`gamma > 0` is fixed and small enough that the smooth set still contributes
a positive multiple of `X/q`. The following further exclusions can be made:

1. A fixed finite set of primes fails the required carry inequalities on
   only `o(X/q)` parameters. The fixed-prime argument uses bounded-in-probability
   valuations and the bounded fibers of a nonconstant quadratic modulo powers
   of a fixed prime, just as in the repository's original argument.
2. For primes above a sufficiently large fixed `Y`, repeated factors cost
   `O(X/(qY)) + o(X/q)`. Outside the fixed discriminants and resultants,
   there are at most eight root classes modulo `p^2` across the six linear
   factors and `Q4`. If `p^2` divides one of these factor values, then
   `p <= C X`. Summing the progression count `X/(q p^2)+O(1)` gives the
   claimed bound, since `pi(CX)=o(X)` and `q` is fixed.
3. After removing those repeated factors, primes below `X^delta` that fail
   the second-carry test cost `epsilon(delta) X/q + o(X/q)`, with
   `epsilon(delta) -> 0` as `delta -> 0`. At every root of a denominator
   factor, `H'(t)` is nonzero modulo all but finitely many primes. The
   simple-root digit-lifting argument therefore applies, with at most eight
   roots instead of the four or six in the older proofs. Its weighted
   prime sum has the same exponentially decaying integral tail.

Choose `Y` large and `delta` small compared with the positive smoothness
constant. A positive multiple of `X/q` parameters remains after these three
exclusions. At all remaining primes dividing `P4(t)`, the valuation is one
and the first carry is forced by `H(t) = -i (mod p)` for `1 <= i <= 4`.

## 4. The precise remaining statement

On the surviving smooth set, every unresolved obstruction has the form

```
X^delta < p <= X,
p | P4(t),
{H(t)/p^j} < 1/2 for every j >= 2.
```

With the optional stronger smoothness version, the upper endpoint becomes
`X^beta` for some fixed `beta < 1`.

To finish the k=4 proof by this route, one must prove that the union of these
events does not cover the surviving smooth set. A lower bound for smooth
parameters and an unweighted upper bound for k=3 failures do not establish
this intersection. Positive-density subsets can be disjoint. Likewise,
fixed-progression smoothness does not justify inserting a growing collection
of carry restrictions modulo `p^j`.

A sufficient new input would be a carry-sensitive lower-bound sieve for
smooth `Q4(t)` values, retaining a positive proportion after the six linear
factor conditions and the extra quadratic's carry conditions. No such input
has been proved in this continuation.

## 5. Why other smoothness results do not finish this step

Balog and Wooley's 1998 paper, *On strings of consecutive integers with no
large prime factors*, gives arbitrarily smooth consecutive blocks by an
explicit cyclotomic construction:

- https://www.math.purdue.edu/~twooley/publ/1998%20consec.pdf

That removes large prime factors but does not force the requisite base-prime
carries or bound repeated valuations. Its thin parameter set cannot simply
be intersected with the repository's positive-density family.

Bober, Fretwell, Martin and Wooley's *Smooth values of polynomials* constructs
auxiliary substitutions making any one quadratic arbitrarily polysmooth:

- https://arxiv.org/abs/1710.01970

Again, polynomial factorization and smoothness alone are weaker than the
factorial-divisibility condition.

The general Dartyge–Martin–Tenenbaum positive-proportion theorem allows
smoothness level `X^(d - 1/K + epsilon)` when the maximal irreducible degree
is `d` and occurs `K` times. For eight quadratic factors its directly
available level is `X^(15/8 + epsilon)`, not `X`:

- https://personal.math.ubc.ca/~gerg/index.shtml?abstract=PVFLPF
- https://personal.math.ubc.ca/~gerg/papers/downloads/PVFLPF.pdf

The quartic construction therefore requires a separate argument for primes
above its parameter length; a theorem for one quadratic cannot be silently
applied to their product.

## 6. Further primary literature checked during the continuation

### A newer one-quadratic smoothness exponent

Glyn Harman, *Two problems on the greatest prime factor of n^2+1*, Acta
Arithmetica 213 (2024), 273–287, DOI `10.4064/aa230710-18-12`, improves the
positive-proportion smoothness exponent for `n^2+1` from `4/5+epsilon` to
`0.567`:

- https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/213/3/115467/two-problems-on-the-greatest-prime-factor-of-n-2-1

The publisher's abstract verifies this statement. The full paper was not
available through the accessible publisher page, so neither a general
quadratic extension nor a fixed-progression version of this improved
exponent has been used above. The improvement does not itself establish a
simultaneous or carry-weighted result.

### Equidistribution of quadratic roots over prime moduli

Duke, Friedlander and Iwaniec proved the negative-discriminant case in
*Equidistribution of roots of a quadratic congruence to prime moduli*,
Annals of Mathematics 141 (1995), 423–441:

- https://annals.math.princeton.edu/1995/141-2/p08

Tóth's *Roots of quadratic congruences*, IMRN 2000, 719–739, extends this to
nonsquare discriminants of either sign. The fixed-arithmetic-progression
form is stated explicitly as Theorem 6.1.3 on printed page 96 (PDF page 98)
of Tóth's own 2018 dissertation:

- https://real-d.mtak.hu/1120/7/dc_1553_18_doktori_mu.pdf

For a fixed quadratic `P` of nonsquare discriminant and a fixed arithmetic
progression of prime moduli, the normalized roots `nu/p` are equidistributed
in `[0,1]`, provided there are infinitely many such roots. This supports
fixed finite unions of intervals and harmonic weights `1/p` by partial
summation. It does not directly supply uniform estimates for intervals
shrinking with `p`, nor independence from higher Hensel digits. In the
quartic construction, finite-difference arguments can sometimes reduce
carry restrictions to this fixed-interval theorem; that reduction must be
proved separately for the actual polynomial and prime range.

### A relevant switching argument for prescribed prime classes

J. Brüdern and K. Soundararajan, *Common divisors of totients of polynomial
sequences*, Mathematische Zeitschrift 299 (2021), 527–542:

- https://link.springer.com/article/10.1007/s00209-020-02674-7

Their Theorem 4 proves that an irreducible quadratic without a fixed prime
factor has infinitely many values avoiding all primes `p = 1 (mod h)`,
for sufficiently large fixed `h`, while also avoiding sufficiently small
prime factors. Sections 5–7 use a switching argument and estimates for
quadratic roots over composite moduli to control primes near the parameter
length. Lemmas 7 and 8 give a useful precise framework for such estimates.
The stated upper bounds contain unspecified constants; the theorem does
not automatically apply to a particular union of residue classes of density
`1/6`, nor to the carry conditions here.

For orientation, elementary switching of `Q(y)=a p`, with `y` of size `Y`
and `p > Y^alpha`, gives `a < C Y^(2-alpha)`. Sieving the variable
`y = r + a ell` uses an interval of length `Y/a`, rather than the size
`Y^2/a` of the prime quotient. One cannot replace `log(Y/a)` in a basic
upper-bound sieve by `log(Y^2/a)`. The range `p` near `Y` is consequently
a substantive issue, not an endpoint error.
