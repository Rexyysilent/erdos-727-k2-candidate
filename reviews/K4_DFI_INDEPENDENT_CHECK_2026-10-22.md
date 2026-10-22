# Independent check of the averaged quadratic carry argument

Date: 22 October 2026.

Scope: Section 5 of `research/k4/analytic-audit.md`, including its
`0.5625676641884855292...` upper bound for the failure contribution of primes
below the parameter height in the 552 family. This is a separate mathematical
check by the research agent, not a proof-assistant certificate.

## Result

I found no mathematical error in the fixed-arc, finite-differencing argument
or in the resulting below-height bound. The statement is an average over
prime/root pairs; it does not assert equidistribution for each prime, nor
extend to primes above the parameter height.

## 1. Exact scope of Tóth's theorem

The primary source is Á. Tóth, *Cycle Integrals of Modular Forms* (2018),
Chapter 6, Theorem 6.1.3, printed p. 96 (PDF p. 98):

https://real-d.mtak.hu/1120/7/dc_1553_18_doktori_mu.pdf

The definitions on printed pp. 94–95 explicitly use normalized roots
`nu/p` in `[0,1]`, ordered by increasing prime denominator. The theorem
allows any fixed arithmetic progression of prime moduli and a fixed
quadratic of nonsquare discriminant. Both signs of the discriminant are
covered. A progression on which there are no roots contributes zero and
can simply be omitted.

For a fixed finite union of arcs E, let A_E(t) denote the root count up to
t in a fixed prime progression and R(t) the total root count there. The
theorem and the prime number theorem in progressions give

```
A_E(t) - |E| R(t) = o(t/log t).
```

Partial summation on `[X^alpha,X^beta]`, with fixed `0<alpha<beta`, gives
an error `o(1)` for the corresponding `1/p`-weighted counts: the integral
of `dt/(t log t)` over this interval is the fixed number
`log(beta/alpha)`. A rate uniform in shrinking arcs is not needed.

## 2. Root representatives and progression restrictions

Write `P(x)=x^4+B*x^2+L*x+C`, `F=P/D`, and suppose
`g=x^2+u*x+v` divides `P+iD`. Take a fixed parameter modulus M divisible
by D. For each prime outside the fixed exceptional set, imposing a root
modulo p and the parameter progression gives

```
x=R+M*p*b,  R=r+k*p,  0<=r<p,  g(r)=0 (mod p).
```

The representative k depends on p and r. This dependence is harmless:
the terms it changes in the relevant leading coefficients are integers.
The coefficients of the highest surviving powers of b, after subtracting
integer-valued terms, are exactly

```
level 2: (M/D)*((4u^2-4v+2B)*r+4uv+L)/p,
level 3: (M^2/D)*(-6u*r+B-6v)/p,
level 4: (4M^3/D)*r/p.
```

All effects of higher lifts occur in lower coefficients. The
finite-differencing lemma permits arbitrary lower coefficients and arbitrary
starting points of the b interval, so it introduces no hidden independence
assumption about Hensel digits.

## 3. Fixed differencing parameters are sufficient

For degree d and leading coefficient `A*r/p+B0/p`, with fixed nonzero
integer A, choose a fixed differencing cutoff H. After d-1 differencings,
the finitely many linear slopes are

```
d! * h1 * ... * h_(d-1) * (A*r/p+B0/p).
```

Exclude a fixed union of arcs where the distance of one of
`d! * h1 * ... * h_(d-1) * A*r/p` from an integer is below `2 epsilon`.
The union has measure `O_H(epsilon)`. Outside it, the linear sums have
normalized size `O_H(1/(N epsilon))` once p is large; `N` tends to infinity
uniformly for `p<=X^(1-eta)`. The differencing boundary errors are
`O_H(1/N)`, and the diagonal terms vanish after H tends to infinity.

The order of limits is `X -> infinity`, then `epsilon -> 0`, then
`H -> infinity`. Tóth bounds the exceptional arc mass, so the weighted
average of the absolute normalized Fourier sums tends to zero. Taking
the highest nonzero carry frequency establishes joint equidistribution,
not merely marginal equidistribution.

When u=0, the level-3 leading coefficient is a fixed constant divided by
p. The argument correctly falls back to the uniform second-derivative
estimate and its threshold `p<X^(2/3-eta)`; it does not claim the extra
level throughout the whole range.

For levels at least 5, the fourth-derivative estimate used in the audit
agrees with the primary theorem in J. Arias de Reyna (2024):

https://arxiv.org/html/2407.02094v1

With constant fourth derivative of size lambda, that theorem gives the
three quantities `N^(-1/8)`, `lambda^(1/14)`, and
`(lambda*N^4)^(-1/8)`, as used there.

## 4. Stable linear factors and exponent bands

For a linear root s of `P+iD`, set `y=s+a*p`. Direct expansion gives

```
D*F(y)=a^4*p^4+4s*a^3*p^3+(6s^2+B)*a^2*p^2
       +P'(s)*a*p-iD.
```

The errors after extracting the stable residues at carry levels 2, 3,
and 4 are respectively of sizes `O(X/p^2)`, `O(X^2/p^3)`, and
`O(X^3/p^4)`. Their thresholds are therefore `1/2`, `2/3`, and `3/4`,
exactly as in the audit.

For B=-2810 and s in `{41,-41,97,-97}`, the signs of the first omitted
terms are respectively `sign(s)`, positive, and `sign(s)`. Thus the
zero and half-denominator endpoint conventions in the residue table are
correct.

The integrations match the available coordinate sets:

- Below exponent 1/2, all coordinates 2 through J are used.
- From 1/2 to 2/3, stable level 2 is combined with equidistributed levels
  3 through J.
- From 2/3 to 3/4, stable levels 2 and 3 are combined with levels 4 and 5.
- From 3/4 to 4/5, stable levels 2,3,4 are combined with level 5.
- From 4/5 to 1, only the three stable levels are used.

At fixed lower cutoff there are finitely many bands. The omitted strips
have arbitrarily small harmonic prime mass. The O(1) errors in progression
counts sum to `O(pi(CX))=o(X)`, with the fixed parameter modulus absorbed
in the constant. This argument would fail for primes up to `CX^2`.

I independently recomputed the displayed sums and obtained

```
T6    = 0.008579551796005669...
Cquad = 0.44435772291501785...
Clin  = 0.1182099412734677...
total = 0.5625676641884855292...
```

The exact finite residue counts were already recomputed by two other
agents; this check verified their analytic use and the endpoint signs,
without adding a third enumeration of all 13,989,888 unit classes.

## 5. Separate caution about a possible large-prime input

Grimmelt–Merikoski (2025), Theorems 1.4 and 1.5,

https://arxiv.org/html/2505.00493v2

assume `a,h>=1` in `a*l^2+h`. Their discussion at the end of Section 1.2
explicitly restricts the paper to negative discriminants. These estimates
cannot be applied as stated to the newer positive-discriminant factors
`(y+c)^2-249` or `(y+c)^2-681`.

A relevant input that does allow either sign is de la Bretèche–Drappeau,
Section 8.4 of

https://arxiv.org/pdf/1703.03197

For fixed nonsquare Delta in `y^2-Delta`, prime-divisor scale P, parameter
height X, and sieve level D, their remainder is bounded by

```
X^(epsilon+O(eta)) *
  (D*P/X + sqrt(D*P) + sqrt(X)*P^theta*D^(1-2theta)).
```

Their allowed level is therefore

```
D <= X^(-K*eta) * min(X^2/P,
                     X^(1/(2-4theta))*P^(-theta/(1-2theta))).
```

At `theta=7/64`, `P=X^alpha`, the exponent is
`min(2-alpha, 0.64-0.14*alpha)`. In the range `1<=alpha<=1.5`, the second
term is smaller. Restrictions to the desired parameter progression,
prime residue classes, and carry weights still require a derivation;
the quoted remainder alone is not a completed upper bound for the
unresolved failures.
