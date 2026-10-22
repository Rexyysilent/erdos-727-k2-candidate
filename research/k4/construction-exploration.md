# Construction branch for k = 4

22 October 2026. This is a mathematical exploration based on repository commit
`8c1e3e50040224b5d3394f91f8c491be8307cd47`. The identities and finite certificates
below are exact. They do **not** establish infinitely many solutions for k = 4.

Later developments are recorded in
[`construction-linear-analysis.md`](construction-linear-analysis.md), which
derives the exact552-family linear carry cost and complete residue refinement
optimization, and
[`construction-additional-families.md`](construction-additional-families.md),
which documents the downward432 family, signed all-six norm searches,
non-even quartics, and their finite certificates and limitations.

## 1. A rigorous barrier for every quadratic family

Let `F(t)=At²+Bt+C` with rational coefficients and `A != 0`. If all four
polynomials `F(t)+i`, `1 <= i <= 4`, split into rational linear factors, their
discriminants

`D_i = B²-4A(C+i)`

are rational squares. They form a nonconstant four-term arithmetic progression
with common difference `-4A`. Fermat's four-squares theorem prohibits such a
progression. Therefore **no nonconstant quadratic over Q has all four consecutive
shifts completely split over Q**. Changing the denominator, the residue class,
or the starting discriminant cannot evade this obstruction.

The fixed-denominator ansatz `n=m(m+1)/Q-4` would require
`1,1+4Q,1+8Q,1+12Q` to be squares, which is the same obstruction in a special
form. This blocks one particular proof architecture, not the original k=4
divisibility statement.

Primary reference: Alf van der Poorten, *Fermat's Four Squares Theorem*,
[arXiv:0712.3850](https://arxiv.org/abs/0712.3850), 2007.
An additional exposition is Keith Conrad, *Arithmetic Progressions of Four
Squares*, <https://kconrad.math.uconn.edu/blurbs/ugradnumthy/4squarearithprog.pdf>.

## 2. A quartic family controlling all four factors

Start from the existing k=3 form

`n=m(m+1)/210-3`.

Choose

`x=210t+14`, `m=x²+209`, `t>=0`.

Then `m == 195 mod 210`, so n is an integer, and

`n = (x⁴+419x²+43260)/210`.

All four consecutive integers split into quadratic factors:

| Shift | Exact factorization |
|---|---|
| `210(n+1)` | `(x²+189)(x²+230)` |
| `210(n+2)` | `(x²+195)(x²+224)` |
| `210(n+3)` | `(x²+209)(x²+210)` |
| `210(n+4)` | `(x²-x+210)(x²+x+210)` |

The first three identities follow by substituting `m=x²+209` into the existing
six-factor construction. For the fourth,

`m(m+1)+210 = x⁴+419x²+44100`

`= (x²+210)²-x²`.

Thus every prime dividing the four-shift product is bounded by `x²+230`, apart
from the harmless fixed denominator bookkeeping. Since `n~x⁴/210`, this gives
a simultaneous square-root-scale factorization. The factors are not linear in
the actual counting parameter t, and the existing analytic proof cannot simply
be transferred by restricting its m parameter to squares plus 209.

### General norm construction

For any rational quadratic `q(s)=A(s-h)²+D`, put

`g(z)=h+c z²+D/(Ac)-c/4`, with rational `c != 0`.

Then

`q(g(z))=Ac²(z²+z+D/(Ac²)+1/4)(z²-z+D/(Ac²)+1/4)`.

Applying this to `q(s)=210s²+391s+183=H(s)+4` with `c=1/210` gives
`g(z)=(z²+14)/210`, which is exactly the quartic above. This construction
explains why a quadratic substitution can split the formerly irreducible
fourth shift without pretending to obtain eight linear factors.

## 3. A proved deterministic carry lemma for the fourth shift

Let x be as above and define

`B(x)=(x³-x²+210x-420)/7`.

**Lemma.** If a prime p divides n+4 and `p²>B(x)`, then
`v_p(n+4)=1` and the doubling of n in base p has at least two carries.
Consequently this prime cannot obstruct the k=4 divisibility condition.

**Proof.** Put `P_σ=x²+σx+210` for `σ=+1,-1`. For x>=14, `B(x)` is larger
than both `P_+` and `P_-`. Thus `p²>B(x)` implies p>7, p² cannot divide either
factor, and p divides only one factor: a common prime divisor of the two
factors divides `2x` and 210, so is among 2,3,5,7. This proves valuation one.
Also `n == -4 mod p`, so the first carry is present.

Write the divisible factor as `P_σ=ap`. Since the other factor is
`P_{-σ}=ap-2σx`,

`n/p² = a²/210 - σ ax/(105p) - 4/p²`.                     (1)

When σ=+1, `P_+ == 0 mod 210` and p is a unit modulo 210, hence a is
divisible by 210. The right side of (1) is an integer minus

`δ_+ = (2xP_+ + 840)/(210p²)`.

The inequality `p²>B(x)` implies `0<δ_+<1/2`. Indeed

`B(x) - (2xP_+ +840)/105`

`= (13x³-17x²+2730x-7140)/105 > 0` for x>=14.

Its fractional part therefore exceeds 1/2, giving the second carry.

When σ=-1, `P_- == 182 mod 210`. The square of a unit modulo 210 lies
in `{1,79,109,121,151,169}`, and multiplication by `182²` sends these
classes into `{154,196}`. Therefore the fractional part of `a²/210` is
either `154/210` or `196/210`. Equation (1) adds

`δ_- = (2xP_- -840)/(210p²)`.

This is positive. The assumed bound is exactly

`2xP_- -840 < 14p²`,

so `0<δ_-<14/210`. The resulting fractional part remains between 1/2 and 1.
This supplies the second carry and completes the proof.

In particular all fourth-shift primes above approximately `x^(3/2)/sqrt(7)`
are automatically harmless. This is an exact sufficient bound, not a claim
that smaller primes necessarily fail.

## 4. Large-prime residue obstruction in the first three shifts

For a first-three-shift factor `x²+209+c=ap`, with

`c in {0,1,-14,15,-20,21}`,

the exact old identity remains

`n/p² = a²/210 + (1-2c)a/(210p) - i_c/p²`,

where `i_c=3,3,2,2,1,1` respectively. Since `m==195 mod210`, the residue
table in the existing K3 paper still applies for sufficiently large `p/x`.
For example `p>4x` is a uniform sufficient separation for x>=14.

Only the following factors and inverse-prime-square classes can then miss
the second carry:

| Quadratic factor | Bad classes of `p^(-2) mod210` |
|---|---|
| `x²+209` | `1,169` |
| `x²+195` | `1,79` |
| `x²+230` | `1` |

The other three factors give the second carry in this range. Merely imposing
the necessary quadratic-root condition does not improve the mean coefficient
`5/6`: the three bad factors have quadratic characters with conductors involving
the additional primes `{11,19}`, `{13}`, and `{23}`, respectively. Those
characters are not determined by the residue modulo 210. Over compatible
prime classes their mean root count remains one, just as in the old linear
incidence calculation. An actual sieve or divisor-switching theorem is needed
to exploit any further information.

If additionally `p³>2n`, absence of the second carry is an exact obstruction:
there are no higher carry levels available. This explains why the first three
quadratics' bad large primes cannot simply be ignored.

## 5. What remains unproved

The old proof counts `m+c==0 modp` in an interval of m of length X by
`X/p+O(1)`. Here m is restricted to `x²+209`. A quadratic factor vanishing
modulo p gives at most two root classes of x. Counting these in an x interval
of length T gives `2T/p+O(1)` per prime, but the primes can be as large as
`O(T²)`. Summing the `O(1)` errors is not negligible. The long range p>T
is the central new obstacle; it cannot be repaired just by reusing the old
six-root constant.

Moreover the existing k=3 lower proportion for all m in a residue class does
not imply any lower proportion on the sparse subset m=x²+209. Such a
restriction requires a new theorem. Finite scans below do not establish that
theorem.

Two concrete routes remain:

1. A sufficiently sharp upper-bound sieve or divisor-switching estimate for
   prime divisors of the three bad even quadratics, preserving the specified
   residue classes and carry restrictions.
2. A different higher-degree family whose large-prime residues are all forced
   into a good second-carry region, together with controlled estimates for the
   smaller primes.

## 6. Exact finite verification

`construction_quartic.py` verifies the factorization identities, reconstructs
the product's full prime factorization, independently factors n+1 through n+4,
and checks Legendre valuations against explicit base-p digit addition.

For t=8:

`x=1694`, `m=2869845`, `n=39219110444`,

and

`((39219110448)!)² | (78438220888)!`.

The complete prime certificate is `quartic-witness-t8.json`. The exact scan
`0<=t<1000` finds 157 successes, with the first at t=8; its full result is
`quartic-scan-1000.json`. This is a finite family scan, not a proof of positive
asymptotic density, infinitude, or smallest n globally.

Reproduce with Python and SymPy:

```sh
python research/k4/construction_quartic.py --t 8
python research/k4/construction_quartic.py --scan-stop 1000
```

## 7. An improved partial-splitting family

A further norm-parameter search yields

`m=(y²-1681)/552`, `n=m(m+1)/210-3`.

Choose `y == 23971 mod115920`; then `m == 195 mod210` and n is an integer.
The exact factorizations, with `D=210*552²=63987840`, are

| Shift | Exact factorization |
|---|---|
| `D(n+1)` | `(y²-12721)(y²+9911)` |
| `D(n+2)` | `(y-97)(y+97)(y²+6599)` |
| `D(n+3)` | `(y-41)(y+41)(y²-1129)` |
| `D(n+4)` | `(y²-138y+8117)(y²+138y+8117)` |

This splits **two of the three bad old quadratic factors** into linear factors:
`m` and `m-14`. Among the first three shifts, only `m+21` can therefore
produce a bad prime much larger than y. Its bad unit-square class is just
`p^(-2)==1 mod210`, occupying one sixth of the prime unit classes.

For the fourth shift put `Q_σ=y²+σ138y+8117`. If `Q_σ=ap`, then

`n/p² = a²/D - σ276ay/(Dp) -4/p²`.

Above `y^(3/2+eta)`, the correction tends to zero uniformly. On the chosen
y residue class, the leading fractions are as follows:

| Factor | Possible leading fractional parts | Bad fraction of prime unit classes |
|---|---|---|
| `Q_-` | `10/21,13/21,19/21` | `1/3` |
| `Q_+` | `21/40,29/40` | `0` |

Here the residue of `Q_σ² modD` is fixed by y modulo115920: if y is increased
by 115920, its square changes by a multiple of D because
`552 | 2Q_σ(2y+σ138)` on these classes. Reduction by the common divisor of
`Q_σ²` and D leaves the small moduli21 and40 in the table. The possible
fractions follow by enumerating unit squares in those moduli. None is a
boundary fraction0 or1/2, so no correction-sign ambiguity occurs.

The combined coefficient for genuinely large quadratic-prime obstructions is
therefore `1/6+1/3=1/2`, compared with `5/6` for the simpler family of Section2.
This is a residue-table coefficient only. It is **not** an established density
or a license to sum quadratic prime incidences as if they were linear ones.
The missing large-prime counting theorem remains missing.

### Search parametrization

For coprime positive integers P,Q define

`M=4PQ`, `N=839Q²-P²-2PQ`, `C=839Q²+P²`,

`m=(y²+N)/M`.

Then

`M²(m(m+1)+210)=(y²-2Py+C)(y²+2Py+C)`.

The choice `(P,Q)=(69,2)` gives the552 family. Complete splitting of
`m,m-14,m+21` in this parametrization would require

`-N`, `-N+14M`, `-N-21M`

all to be integer squares (rational squares of integers are integer squares).
The last condition forces `P/Q>91`; the first is equivalent to

`(P+Q-U)(P+Q+U)=840Q²`.

The exact bounded search `construction_three_bad_split_search.cpp` enumerates
the factor pairs on the right. For `1<=P<=10000000`, all possible Q are
covered by `Q<=109890`. It checked196650 admissible candidates and found no
all-three-splitting example. This is a bounded negative result, not an
impossibility theorem.

## 8. Completely split cubics: a bounded negative search

A rational cubic can be shifted and rescaled so that its roots, after clearing
denominators, are integer roots of `z³-Az+B`. Every split depressed cubic has
roots `r,s,-r-s`, and thus

`A=r²+rs+s²`, `B=rs(r+s)`.

Up to changing all root signs and permuting roots, one can take `r>=s>=0`.
For each A, the search records both signs of every available B and checks for
every nonconstant four-term arithmetic progression among these B values.
Such a progression would yield four consecutive completely split translates
after dividing the polynomial by its common difference.

`construction_cubic_split_search.cpp` exhaustively checked `A<=100000000`:
30237841 nonnegative root-pair records, 8738132 A values with multiple
representations, and zero four-term progression hits. The search covers
general four-term progressions, including but not restricted to symmetric
ones `{-3b,-b,b,3b}`. It does not prove that larger-coefficient examples do not
exist.
