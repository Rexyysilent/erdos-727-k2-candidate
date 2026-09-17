# Separate proposed k = 3 extension for Erdős problem 727

17 September 2026. **New proposed argument, seeking independent review.**

The supplied k=2 manuscript explicitly stops at k=2 because its four-factor
construction does not control `n+3`. This note records a separate six-factor
construction. It is not claimed to solve the problem for every fixed k, and no
claim for k>=4 is made.

## 1. Family and exact factorization

Take

```text
H(t) = 210t^2 + 391t + 179.
```

Equivalently let `m=210t+195`. Then

```text
n+3 = m(m+1)/210,
n+2 = (m-14)(m+15)/210,
n+1 = (m-20)(m+21)/210.
```

The three consecutive factors split completely:

```text
n+1 = (6t+5)(35t+36),
n+2 = (t+1)(210t+181),
n+3 = (14t+13)(15t+14).
```

Hence

```text
(n+1)(n+2)(n+3)
 = (6t+5)(35t+36)(t+1)(210t+181)(14t+13)(15t+14).
```

All six linear factors are primitive.

The factorial condition is exactly

```text
((n+3)!)^2 | (2n)!
<=> ((n+1)(n+2)(n+3))^2 | binom(2n,n).
```

## 2. Carry formulation

For every prime p,

```text
C_p(n) = v_p binom(2n,n)
       = sum_{j>=1} 1_{ {n/p^j} >= 1/2 }.
```

Thus divisibility holds exactly when

```text
C_p(n) >= 2 v_p((n+1)(n+2)(n+3))
```

for every prime dividing the product.

For large primes, each of the six shifts of `m` occurs in at most one factor,
and valuation one forces one initial carry. The remaining question is whether a
second carry appears.

## 3. Six simple roots

The relevant shifts of `m` are

```text
0, 1, -14, 15, -20, 21.
```

For primes larger than the shift separations and not dividing 210, the six roots
are distinct. Since

```text
H'(t) = 420t+391 = 2m+1,
```

at the root corresponding to a shift c the derivative is congruent to `1-2c`.
These residues are nonzero for the primes in the growing-prime analysis, giving
simple-root lifting exactly as in the k=2 construction.

At the fixed denominator primes 2,3,5,7, the polynomial has unit derivative and
is a permutation modulo the prime; lifting gives a permutation modulo all prime
powers. The same fixed-prime carry-count argument therefore removes every fixed
prime outside a density-zero exceptional set.

## 4. Large-prime residue table

For `m+c=ap`, substitution has the same form

```text
n/p^2 = a^2/210 + (1-2c)a/(210p) - i_c/p^2,
```

with `i_c` determined by whether p divides `n+1`, `n+2`, or `n+3`.

Averaging the signed second-carry obstruction over units modulo 210 gives the
exact rational value

```text
b_3 = 5/6.
```

The six unit-square classes modulo 210 are

```text
1, 79, 109, 121, 151, 169,
```

each with multiplicity 8. The corresponding numbers of bad shifts are

```text
3, 1, 0, 0, 0, 1.
```

Their weighted average is therefore `5/6`.

## 5. Proposed lower bound

Let

```text
S = 4 * sum_{J>=4} 2^(1-J) log((J+1)/J)
  = 0.191470372039616208448161003021...
```

The lower-prime analysis now has six root classes rather than four, so its cost
is `(3/2)S`. The intermediate and large-prime ranges contribute
`(5/12)log 3` from the residue average `5/6`.

The proposed theorem is

```text
liminf_{T->infinity} (1/T) * # {
  t in [T,2T] : ((H(t)+3)!)^2 | (2H(t))!
}
>= rho_3,
```

where

```text
rho_3 = 1 - (3/2)S - (5/12)log 3
      = 0.25503932166219664924...
      > 0.
```

Positivity does not depend on the decimal evaluation. Using `S<0.192` and
`log 3 < 1.1` already gives a positive rational margin.

The analytic proof follows the k=2 order of limits: fix all auxiliary cutoffs,
take the height to infinity, then remove boundary margins and finally the fixed
prime cutoff. No simultaneous growing-J or growing-modulus estimate is asserted.

## 6. Exact finite evidence

For `H(t)=210t²+391t+179`, exact counts are:

| t interval | valid parameters |
|---|---:|
| `[0,1000)` | 126 |
| `[10000,20000)` | 1,880 |
| `[100000,200000)` | 21,456 |

The first witness in this family is

```text
t = 61,
n = 805440,
((805443)!)^2 | 1610880!.
```

A prime-by-prime certificate verifies the required valuations for this n. A
second implementation agreed on all 1,000 selected parameters in
`[199000,200000)`, of which 240 pass. Literal factorial, central-binomial, and
carry formulations also agree for every `0<=n<=500`.

These computations do not prove the asymptotic theorem.

## 7. What must be independently checked

Independent review should focus on:

- the six-factor identities and fixed-prime permutation argument;
- simplicity and distinctness of all six growing-prime roots;
- the signed residue table giving `5/6`;
- the multiplication of the lower-band cost from four to six roots;
- the fixed-modulus prime-class summation;
- complete coverage of exponent bands and the stated limit order.

Even if this k=3 argument is accepted, extending to k>=4 requires new control of
additional consecutive factors. There is no automatic induction here.
