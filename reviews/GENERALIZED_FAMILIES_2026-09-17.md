# Generalized quadratic families for Erdős 727, k = 2

17 September 2026. **Proposed extension, seeking independent review.**

This note does not replace the original family `6t²+11t+3`. It records a
broader fixed-denominator construction and the strongest checked finite-family
consequence found in the current continuation. The argument is analytic and
remains unreviewed by an independent mathematician.

## 1. Family

Fix an integer `u >= 1`, let

```text
Q = u(u+1),
0 <= r < Q,
Q | r(r+1),
m = Qt+r,
```

and define

```text
F_{Q,r}(t) = Qt^2 + (2r+1)t + r(r+1)/Q - 2.
```

Then

```text
n+2 = m(m+1)/Q,
n+1 = (m-u)(m+u+1)/Q.
```

If `C={0,1,-u,u+1}` and `g_c=gcd(Q,r+c)`, the four integers

```text
L_c(t) = (Q/g_c)t + (r+c)/g_c
```

are primitive linear polynomials and

```text
(n+1)(n+2) = product_{c in C} L_c(t).
```

The identity `product_c g_c = Q^2` follows prime-power by prime-power from
`Q|r(r+1)` and `Q=u(u+1)`.

## 2. Residue coefficient

For a unit `v mod Q`, let `s_c(v)` be the least nonnegative residue of
`(r+c)^2 v^2 mod Q`, and set `d_c=1-2c`. Define the boundary-aware indicator

```text
chi(s,d) = 1  if 0 < s < Q/2
           1  if s = 0   and d > 0
           1  if 2s = Q and d < 0
           0  otherwise.
```

The exact rational coefficient is

```text
b(Q,r) = (1/phi(Q)) * sum_{v unit mod Q} sum_{c in C} chi(s_c(v),d_c).
```

This is finite residue algebra, not a fitted probability.

Write

```text
S = 4 * sum_{J>=4} 2^(1-J) log((J+1)/J)
  = 0.191470372039616208448161003021...
```

## 3. Proposed lower bound

For every fixed admissible pair `(Q,r)`, the continuation gives

```text
liminf_{T->infinity} (1/T) * # {
  t in [T,2T] : ((F_{Q,r}(t)+2)!)^2 | (2F_{Q,r}(t))!
}
>= 1 - S - (b(Q,r)/2) log 3.
```

The proof is an adaptation of the supplied k=2 argument:

1. cancellation reduces the factorial condition to central-binomial valuations;
2. fixed primes and repeated large factors have vanishing total density;
3. the four roots are simple for sufficiently large primes;
4. the same Fourier/box argument gives the lower-prime carry cost `S`;
5. for primes above the square-root parameter scale, the signed residue table
   determines which shifts can fail to supply the second carry;
6. fixed-modulus Mertens asymptotics in arithmetic progressions weight those
   unit residue classes.

All moduli are fixed before the height tends to infinity. No growing-modulus
uniformity or GRH-dependent estimate is used in the stated argument.

## 4. Strongest checked family in the bounded search

For

```text
u = 14,
Q = 210,
r = 209,
F(t) = 210t^2 + 419t + 207,
```

the exact residue calculation gives

```text
b(210,209) = 2/3.
```

Hence the proposed parameter lower bound is

```text
rho_210 = 1 - S - (1/3)log 3
        = 0.44232553173768056108...
```

The same coefficient occurs for `r=195`. This is a stronger lower **parameter
proportion in a different family**. It does not improve the original family’s
constant `rho=0.2592234836263289...`.

A bounded exact search through `1 <= u <= 150` examined 2,058 admissible
root classes. The minimum observed residue coefficient was `2/3`, attained at
`(u,r)=(14,195)` and `(14,209)`. This is a finite-search statement only.

## 5. Combining a fixed finite set of denominators

Suppose two distinct fixed denominators `Q,R` yield the same value:

```text
m(m+1)/Q - 2 = k(k+1)/R - 2.
```

With `y=2m+1` and `z=2k+1`, this becomes

```text
R y^2 - Q z^2 = R-Q.
```

After scaling, this is a fixed generalized Pell-type equation. Partitioning
solutions by residue class modulo the fixed nonzero right-hand side shows that
successive positive solutions in one class grow by a fixed factor. Therefore
any two fixed different denominators have only `O_{Q,R}(log N)` common values
below height `N`.

Consequently the contributions of any **fixed finite** collection of families
can be added up to `o(sqrt(N))` overlap loss.

Ten explicitly checked families give the proposed bound

```text
liminf A_2(N)/sqrt(N) >= 0.43890932340518345610...
```

where

```text
A_2(N) = # {1 <= n <= N : ((n+2)!)^2 | (2n)!}.
```

The fixed bounded enumeration through `u<=150` selects 296 positive-bound
families across 92 denominators and gives the larger proposed coefficient

```text
0.83733579699699106572...
```

for `A_2(N)/sqrt(N)`. This is **not** an 83.7% natural density. The finite
family set is fixed before `N->infinity`, and no global optimality is claimed.

## 6. Finite arithmetic checks

Selected exact counts, with upper endpoints excluded:

| Family | `[0,1000)` | `[10000,20000)` | `[100000,200000)` |
|---|---:|---:|---:|
| `6t²+7t` | 126 | 2,105 | 23,896 |
| `210t²+419t+207` | 223 | 2,849 | 31,711 |
| `210t²+391t+180` | 220 | 2,888 | 31,725 |

A second implementation checked a further 1,000 high-range parameters for each
of these three families. These computations validate arithmetic identities and
individual divisibility decisions; they do not establish the analytic liminf.

## 7. Review targets

The most valuable independent checks are:

- the fixed-prime argument for arbitrary fixed `Q`;
- the sign convention at residues `0` and `Q/2`;
- the fixed-modulus prime-class weighting;
- the cross-denominator `O(log N)` intersection argument;
- the rule that the finite family set is fixed before `N->infinity`.

A future improvement should classify or lower-bound `b(Q,r)` theoretically,
rather than infer an optimum from a bounded search.

## Sources

- Supplied 16 September 2026 expanded k=2 manuscript, reviewed locally.
- D. Keliher and E. S. Lee, *On the constants in Mertens' theorems for primes
  in arithmetic progressions*, arXiv:2306.09981v2, equation (6), fixed modulus.
- O. Robert, *On van der Corput's k-th derivative test for exponential sums*,
  Indagationes Mathematicae 27 (2016), Theorem 1.
- L. Goldmakher, *A quick proof of Mertens' theorem*.
