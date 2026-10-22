# Additional quartic families and completed searches

22 October 2026. All identities and finite residue counts below are exact.
The numerical comparisons explicitly labeled heuristic are not bounds on
prime divisors of polynomial values. No k=4 infinitude theorem follows.

## 1. Both directions of the norm construction

Let P be a nonzero integer, Q a positive integer, and define

\[
M=4PQ,\quad N=-\Delta Q^2-P^2-2PQ,\quad C=P^2-\Delta Q^2,
\qquad m=(y^2+N)/M.
\]

With K=(1−Δ)/4, the identity is

\[
M^2\bigl(m(m+1)+K\bigr)
 =(y^2-2Py+C)(y^2+2Py+C).
\]

The choice Δ=−839 gives K=210 and n=m(m+1)/210−3. The choice Δ=2521
gives K=−630 and n=m(m+1)/210−4, extending the three old shifts downward.
Both signs of P, and hence of M, were searched. Replacing P by −P exchanges
m with −m−1 and leaves the quartic n itself unchanged; it permutes the old
factor labels. Both signs were retained explicitly so that signed endpoint
rules could not silently be carried over from positive m.

For example, at m≡195 (mod210), the positive-m factor m+15 has a leading
second-carry fraction 0 with a negative correction. For negative m that
correction becomes positive, so this factor becomes bad in every prime unit
class. Therefore using the old positive-m table unchanged would be incorrect.

For all current small candidates, the local search includes every one of the
16 classes m₀ (mod210) with m₀(m₀+1)≡0, every admissible y class in the
specified period, and all six old factors

\[
m+c,\qquad c\in\{0,1,-14,15,-20,21\}.
\]

A factor is declared split only when −N−cM is an actual nonnegative integer
square. Merely sharing a squareclass with another constant is insufficient.

For the new quadratic pair, putting

\[
A_\sigma=2Qm+P+Q+\sigma y
\]

gives y²+σ2Py+C=2PAσ, so its leading carry fraction reduces to

\[
\frac{A_\sigma^2p^{-2}}{840Q^2}\pmod1.
\]

This is why the new residue denominator is independent of P. A period
210|M|Q suffices for these squared residues and the old m₀ class. The program
`construction_norm_local_search.py` checks all admissible classes in that
period when it does not exceed the stated `--max-period`, default 3000000.
Its JSON records exactly which larger candidates were skipped. These
high-prime residue counts do not include a counting theorem for primes above
the parameter height.

## 2. The downward M=432 family

Choose P=108, Q=1 and Δ=2521. Then

\[
m=\frac{y^2-14401}{432},\qquad
n=\frac{m(m+1)}{210}-4
 =\frac{y^4-28370y^2+44403409}{39191040}.
\]

For y≡16273 (mod90720), we have m≡174 (mod210). Put D=39191040.
The four exact factorizations are

\[
\begin{aligned}
D(n+1)&=(y^2-216y+9143)(y^2+216y+9143),\\
D(n+2)&=(y^2-23041)(y-73)(y+73),\\
D(n+3)&=(y-143)(y+143)(y-89)(y+89),\\
D(n+4)&=(y^2-14401)(y^2-13969).
\end{aligned}
\]

Thus three old quadratic factors split: m−14, m+15 and m+21. There are six
linear factors and five irreducible quadratics, of which two have a nonzero
linear term. The remaining old quadratic discriminant squareclasses are
14401, 61·229 and 23041, with prime support outside210. The new pair has
discriminant squareclass 2521. Hence their root conditions do not remove the
bad classes in the following residue averages.

The stable high-prime coefficient is **1/3**. Among the three even quadratics
only m contributes, with fraction 1/6. For the two new factors, their possible
leading fractions on the displayed progression are

| Factor | Leading fractional parts | Bad proportion |
|---|---|---:|
| y²−216y+9143 | 6/35,19/35,24/35,26/35,31/35,34/35 | 1/6 |
| y²+216y+9143 | 5/6 | 0 |

The minimum 1/3 is attained in exactly four coarse classes modulo90720:
16273, 29087, 61633 and 74447. Every one has m≡174 (mod210).

### Exact linear cost

Refine y to a fixed class modulo D. For each of the four coarse classes,
there are 432 refinements. The exact enumeration checks all 1728 refinements,
using all 279936 unit-square classes modulo D. Every refinement has the same
three joint totals

\[
B_2=\frac{242}{81},\qquad B_{23}=\frac{148}{81},\qquad
B_{234}=\frac{2299}{2916}.
\]

At Y=16273 the individual rows are

| Root | Bad at2 | Bad at2,3 | Bad at2,3,4 |
|---:|---:|---:|---:|
| −143 | 1/2 | 1/4 | 0 |
| −89 | 1/3 | 1/6 | 1/12 |
| −73 | 40/81 | 20/81 | 10/81 |
| 73 | 2/3 | 2/3 | 1/3 |
| 89 | 1/2 | 1/4 | 91/729 |
| 143 | 40/81 | 20/81 | 10/81 |

The exceptional row 73 is important: its level 3 leading fraction is0 with a
positive correction, so that carry is always absent. A blanket factor 1/2
for level 3 would give an incorrect improvement.

The Weyl argument and bands in `construction-linear-analysis.md` apply with
six roots, with the actual signs of F′(r) and F″(r) at each root. The result is

\[
\begin{aligned}
C_{\rm lin,432}={}&6\sum_{J\ge8}2^{1-J}\log\frac{J+1}{J}
 +\frac{242}{81}\sum_{J=6}^7 2^{2-J}\log\frac{J+1}{J}\\
&+\frac{148}{324}\log\frac98
 +\frac{2299}{5832}\log\frac{16}{15}
 +\frac{2299}{2916}\log\frac54\\
={}&0.30647868120367294\ldots.
\end{aligned}
\]

The sources and complete integer counts are
`construction_432_optimize_linear.cpp` and
`432-linear-refinement-optimization.json.gz`.

The quadratic cost below the parameter height is 0.3572020886912154 under the
separate averaged quadratic-carry argument. The optimistic, unproved
high-prime calculation is 0.40112210801324655. Adding these to the exact
linear cost gives a heuristic total 1.0648028779081347. Thus the smaller high
residue coefficient does not improve the full union-bound calculation over
the developed552 family.

### Finite certificate

For y=16273+90720t, the exact scan 0≤t<1000 finds 241 successful parameters.
The first is t=1:

\[
y=106993,\quad m=26498814,\quad n=3343748428117.
\]

Consequently ((3343748428121)!)² divides (6687496856234)!. The certificate
`432-witness-coarse-t1.json` checks 20 relevant prime valuations and contains
31 recursive Lucas prime certificates, including the smaller primes needed
in those proofs. An independent readback verified every Lucas relation and
all 20 certified relevant primes. The scan is `432-scan-coarse-1000.json`.

Reproduce with Python and SymPy:

```sh
python research/k4/construction_norm_certificate.py --t 1 --include-prime-proofs
python research/k4/construction_norm_certificate.py --scan-stop 1000
g++ -O3 -std=c++17 -fopenmp research/k4/construction_432_optimize_linear.cpp -o /tmp/k4_432_linear
OMP_NUM_THREADS=8 /tmp/k4_432_linear
```

## 3. Completed search for three split old factors

The expanded discovery search checks every coprime pair with
1≤|P|≤10000 and 1≤Q≤10000, for each of Δ=2521 andΔ=−839. It tests all six
old constants, retaining any case in which at least three are squares.

The full outputs and captured stderr are

* `downward-all-six-search-10000.txt` and its `.stderr.txt` receipt;
* `upward-all-six-search-10000.txt` and its `.stderr.txt` receipt.

For Δ=2521 there are exactly six signed hits: P=±30,±52,±108, all with Q=1.
For Δ=−839 there are exactly four signed hits: P=±38,±172, again Q=1.
No hit has four or more old factors split. Up to P↦−P, these are five
distinct even quartics. This is a bounded classification only.

Reproduction:

```sh
g++ -O3 -std=c++17 -fopenmp research/k4/construction_downward_norm_search.cpp -o /tmp/k4_norm_search
OMP_NUM_THREADS=8 /tmp/k4_norm_search 10000 2521 allsix
OMP_NUM_THREADS=8 /tmp/k4_norm_search 10000 -839 allsix
```

All coarse integral y representatives were then checked for these five
families, with exact linear and quadratic residue counts. Those additional
comparisons are explicitly **coarse representative searches**, not complete
refinement searches. The best heuristic sums among those representatives are

| M magnitude | Best representative y | Exact linear cost | Old/new high coefficients | Heuristic total |
|---:|---:|---:|---:|---:|
| 120 | 5809 | 0.2927485596 | 1/3,1/6 | 1.1426237803 |
| 208 | 1591 | 0.1832568646 | 1/3,1/6 | 1.0331320853 |
| 432 | 16273 | 0.3064786812 | 1/6,1/6 | 1.0648028779 |
| 152 | 7843 | 0.2632622629 | 1/3,1/6 | 1.1131374836 |
| 688 | 35311 | 0.2704362127 | 1/2,1/6 | 1.2118624575 |

The receipts are `norm<M>-coarse-linear-and-high.json`. No favorable
heuristic sum below 1 was found in this particular finite comparison.

## 4. Non-even quartics and both directions of growth

The analytic branch found the two useful non-even constructions

\[
F(y)=\frac{y^4-3010y^2+7680y+3969}{483840},\qquad
G(y)=\frac{y^4-810y^2+2240y-41751}{40320}.
\]

F has four linear factors and six odd quadratics across its four shifts:

\[
\begin{aligned}
483840(F+1)&=(y^2+40y-801)(y^2-40y-609),\\
483840(F+2)&=(y-21)(y+53)(y^2-32y-873),\\
483840(F+3)&=(y-27)(y+51)(y^2-24y-1057),\\
483840(F+4)&=(y^2+8y-1953)(y^2-8y-993).
\end{aligned}
\]

G has eight linear factors and four odd quadratics. Its full formulas are in
the main review. All these identities were independently expanded exactly,
as recorded in `construction-additional-verification.json`.

The program `construction_noneven_residues.cpp` enumerates every integral
class modulo each denominator, and was run for both F(y) and F(−y), and both
G(y) and G(−y). Sign reversal changes boundary decisions and therefore cannot
be dismissed as a relabeling of positive-y residue classes.

| Family / direction | Integral classes | Minimum high quadratic coefficient | Linear cost at the best class attaining that minimum |
|---|---:|---:|---:|
| F(y), denominator 483840 | 16384 | 22/27 | 0.4865888162 |
| F(−y) | 16384 | 22/27 | 0.4865888162 |
| G(y), denominator 40320 | 4096 | 1 | 0.2793466409 |
| G(−y) | 4096 | 2/3 | 0.4011306039 |

There are better tradeoffs at larger high coefficients; the complete data
are retained in the four `noneven-*-residues.json` files. For instance,
F(−y) at y≡3411 has high coefficient 31/27 and linear cost 0.1205726873.
None of these comparisons closes the needed bound.

Independent full-unit enumerations checked the joint linear tables at three
representatives: G at 1177, G(−y) at 713, and F(−y) at 3411. These involve
9216,9216 and110592 unit classes respectively, and agree exactly with the
separate unit-square/pairing enumeration, including signed boundaries.

Reproduction examples:

```sh
g++ -O3 -std=c++17 -fopenmp research/k4/construction_noneven_residues.cpp -o /tmp/k4_family_residues
OMP_NUM_THREADS=8 /tmp/k4_family_residues 40320
OMP_NUM_THREADS=8 /tmp/k4_family_residues 40320 reverse
OMP_NUM_THREADS=8 /tmp/k4_family_residues 483840 reverse
OMP_NUM_THREADS=8 /tmp/k4_family_residues norm208
```
