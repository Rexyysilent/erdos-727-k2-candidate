# Higher-degree split translates: research notes

These notes record exact algebraic reductions and literature leads from the
2026-10-22 k=4 investigation. They do not establish the k=4 factorial
divisibility theorem.

## 1. The four-fiber target

One possible replacement for the quadratic k=3 family is a polynomial
`f in Q[t]` for which `f+1`, `f+2`, `f+3`, and `f+4` all split into linear
factors over `Q`. Then all denominator primes are at most a constant times
the parameter. The carry estimates would still need to be proved for that
family.

More generally, it is enough to find four completely split fibers whose
values form a nonconstant rational arithmetic progression: an affine change
of the output makes them four consecutive translates.

The root multisets of these fibers have all elementary symmetric functions
except the last in common. They are therefore four ideal
Prouhet–Tarry–Escott root sets, with the additional requirement that their
signed products form an arithmetic progression. Ordinary multigrade chains
do not automatically satisfy this extra requirement.

For quadratics, the four discriminants would be four distinct rational
squares in arithmetic progression, which is impossible. This obstruction
does not, by itself, extend to arbitrary degree.

An additional class is excluded: if, after an affine change of input,
`f(t)=g(t^2)`, then every completely split fiber has roots in opposite pairs.
Its constant term is a fixed sign times a rational square. Four split
fibers in arithmetic progression would again give four distinct rational
squares in arithmetic progression. This rules out the even quartic and
even Chebyshev constructions as candidates for this specific target.

## 2. A cubic reduction

After depressing a cubic, use `f(t)=t^3-A t`. A split fiber has roots
`r,s,-r-s`, with

```
A = r^2+r*s+s^2,
fiber value = -r*s*(r+s).
```

The split fibers occur in opposite pairs. Therefore, a sufficient special
case for four fibers in arithmetic progression is two triples satisfying

```
r^2+r*s+s^2 = u^2+u*v+v^2,
r*s*(r+s) = 3*u*v*(u+v),
u*v*(u+v) != 0.
```

This is a homogeneous intersection of a quadric and a cubic in projective
3-space. The smooth complete intersection has genus 4. A bounded failure
to find points in this special case is not a proof of nonexistence, nor
does it rule out a nonsymmetric four-term progression of split cubic fibers.

## 3. A genus-5 reduction for splitting the three bad old factors

The k=3 outer family can be written

```
n = m*(m+1)/210 - 3,
n+1 = (m-20)*(m+21)/210,
n+2 = (m-14)*(m+15)/210,
n+3 = m*(m+1)/210.
```

On the favorable residue `m=195 (mod 210)`, the old factors with potentially
bad stable second carries are `m`, `m-14`, and `m+21`. Consider making all
three split after a quadratic substitution

```
m = (y^2-a^2)/d,
m-14 = (y^2-b^2)/d,
m+21 = (y^2-c^2)/d.
```

These identities require

```
b^2 = a^2+14*d,
c^2 = a^2-21*d.
```

For the fourth shift, require the symmetric factorization

```
d^2*(m*(m+1)+210)
  = (y^2+h*y+v)*(y^2-h*y+v).
```

Comparing coefficients gives

```
v^2 = a^4-a^2*d+210*d^2,
h^2 = 2*v+2*a^2-d.
```

Assuming `h != 0`, eliminate `v` and set `z=d/h^2`. Then all three
expressions

```
A(z) = 1+2*z-839*z^2,
B(z) = 1+58*z-839*z^2,
C(z) = 1-82*z-839*z^2
```

must be rational squares, specifically `(2a/h)^2`, `(2b/h)^2`, and
`(2c/h)^2`. Conversely, a rational solution with `z != 0` yields the
coefficient identities above after choosing a nonzero rational `h` and
setting `d=z*h^2`.

The three quadratic polynomials have six distinct roots over the algebraic
closure; their square classes are independent. The resulting degree-8
cover of the `z`-line has six branch points, each with ramification index
2. Riemann–Hurwitz gives

```
2*g-2 = 8*(-2) + 6*4 = 8,
g = 5.
```

The point `z=0` is a degenerate solution (`d=0`) and does not define the
desired substitution. No nondegenerate rational point, nor a proof of its
nonexistence, was obtained in this continuation.

## 4. Primary literature leads

Arvind Suresh's 2022 seminar *Factoring translates of polynomials* explicitly
studies `S_f={a in Q : f(t)-a splits over Q}` and its Prouhet–Tarry–Escott
connection:

- https://colleges.claremont.edu/ccms/event/antc-talk-arvind-suresh-university-of-arizona-tucson/

His paper *Constructing curves of large rank using composite polynomials*,
sections 6–8, treats the Galois-cover and cyclic/dihedral constructions
underlying families of split fibers:

- https://arxiv.org/pdf/2102.02113

These results do not supply the required four-term arithmetic progression.

L. Hajdu, Á. Papp and R. Tijdeman, *The Prouhet–Tarry–Escott problem,
indecomposability of polynomials and Diophantine equations*, Ramanujan
Journal 58 (2022), 1075–1093, gives a precise relationship between ideal
PTE partitions and polynomial decompositions. Its main restrictions concern
root sets drawn from an almost complete arithmetic progression; they do
not establish the absence of the general four-fiber target here:

- https://link.springer.com/article/10.1007/s11139-022-00555-7

The numerical multigrade-chain catalog maintained by Chen Shuwen includes
the nonsymmetric quartic chain

```
{0,87,93,214}, {9,52,123,210}, {24,30,133,207}.
```

Each set has the same first three power sums, but the constant terms are
not the desired four-term arithmetic progression. This is a search lead,
not a construction used in a proof:

- https://eslpower.org/chains.htm

### Nonsymmetric quartics with four or five split fibers

The more substantial primary source is A. Choudhry and J. Wróblewski,
*Triads of integers with equal sums of squares and equal products and a
related multigrade chain*, Acta Arithmetica 178 (2017), 87–100,
DOI `10.4064/aa8574-9-2016`:

- https://www.impan.pl/shop/publication/transaction/download/product/92058

The authors give thirteen chains of four positive integer triads with
entries below 30000 (Table 2), and two chains of five triads with larger
entries (pp. 88–89). Within each chain, `x^2+y^2+z^2` and `xyz` are
constant. The four numbers

```
x+y+z, x-y-z, -x+y-z, -x-y+z
```

therefore have the same first three power sums for every triad in the
chain. They provide nonsymmetric quartics with four or five fully split
rational fibers. Thus the existence of four split fibers itself is not
the obstacle, even for an indecomposable quartic.

I checked the constant terms from all thirteen four-triad examples and
both five-triad examples printed in that paper; none contains a four-term
arithmetic progression. This is a finite check of the displayed examples,
not a theorem about the authors' undisplayed search data or the general
four-fiber problem. For the smallest displayed four-triad chain,

```
(4455,1596,902), (4428,1694,855),
(3895,2772,594), (3630,3116,567),
```

the quartic constant terms, in increasing order, are

```
-2588900337911, 40052136267529,
247817675921929, 263505254796169.
```

### Composition cannot create the missing outer fibers

If `f=g(h(t))` and `f(t)-c` splits completely over `Q`, then `g(u)-c`
also splits completely over `Q`: for every root `alpha` of `g(u)-c`,
every root of the nonconstant polynomial `h(t)-alpha` is a root of
`f(t)-c`, and hence is rational; evaluating `h` at one such root proves
`alpha` rational. Therefore four split AP fibers for a composite `f`
already give the same four split AP fibers for its outer polynomial `g`.

In particular, a quadratic outer polynomial is excluded by the
four-squares-in-AP obstruction. A quadratic inner polynomial makes the
composite even after a rational input translation, and is excluded by
the even-polynomial obstruction above. Consequently a degree-4 or
degree-6 example would have to be indecomposable over `Q`.

The focused literature search found no construction of four completely
split fibers in a nonconstant arithmetic progression for any degree at
least three, and no general impossibility theorem. This statement records
the search result only; it is not a claim that the question is known to
be open or that such an example cannot exist.
