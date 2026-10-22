#!/usr/bin/env python3
"""Independent exact checks of the 552-family linear carry identities."""
from fractions import Fraction
from math import gcd
from pathlib import Path
import json

from sympy import symbols, expand

D = 63987840
Y = 23971
ROOTS = (-97, -41, 41, 97)


def bad(s: int, sign: int) -> bool:
    if s == 0:
        return sign > 0
    if 2*s == D:
        return sign < 0
    return 2*s < D


def main() -> None:
    y,a,q = symbols("y a q")
    F = (y**4-2810*y**2-190065671)/D
    for r in ROOTS:
        i = 2 if abs(r)==97 else 3
        E2 = 6*r*r-2810
        E1 = 4*r**3-5620*r
        assert E2 > 0 and E1*r > 0
        assert expand(F.subs(y,r+a*q)-
            (a**4*q**4+4*r*a**3*q**3+E2*a*a*q*q+E1*a*q)/D+i) == 0

    # The endpoint convention itself is checked with exact rational arithmetic.
    eps = Fraction(1,100*D)
    for s in (0,1,D//2-1,D//2,D//2+1,D-1):
        for sign in (-1,1):
            x = Fraction(s,D)+sign*eps
            fraction = x-(x.numerator//x.denominator)
            assert bad(s,sign) == (fraction < Fraction(1,2))

    table = json.loads((Path(__file__).parent/"linear-residues-23971.json").read_text())
    phase_rows = {row["root"]:row for row in table["rows"]}
    checks = 0
    bases = []
    k = 0
    while len(bases)<128:
        base = 10**120 + 104729*k + 1
        k += 1
        if gcd(base,D)==1:
            bases.append(base)
    for base in bases:
        u = pow(base,-1,D)
        for r in ROOTS:
            ks = [(Y-r)**2*(Y*Y+2*r*Y+3*r*r-2810)%D,
                  (Y-r)**3*(Y+3*r)%D,(Y-r)**4%D]
            assert ks == [phase_rows[r][f"K{j}"] for j in (2,3,4)]
            # No boundary is met for this chosen progression; the generic
            # endpoint rule above is nevertheless independently tested.
            assert all(gcd(K,D) not in (D,D//2) for K in ks)
            a0 = (Y-r)*u%D
            for size, levels in ((0,(2,3,4)),(10**108,(2,)),
                                 (10**48,(2,3)),(10**30,(2,3,4))):
                cofactor = a0+D*(size//D)
                value = r+cofactor*base
                assert value%D == Y
                m = (value*value-1681)//552
                assert (value*value-1681)%552 == 0 and m%210 == 195
                n = m*(m+1)//210-3
                i = 2 if abs(r)==97 else 3
                assert n%base == base-i
                for j in levels:
                    sign = 1 if j==3 or r>0 else -1
                    residue = ks[j-2]*pow(u,j,D)%D
                    predicted = bad(residue,sign)
                    actual = 2*(n%base**j) < base**j
                    assert predicted == actual,(base,r,size,j,residue)
                    checks += 1
    print(json.dumps({"symbolic_taylor_identities":4,"endpoint_checks":12,
                      "large_coprime_bases":len(bases),
                      "exact_carry_predictions_checked":checks,
                      "result":"All checks passed; exact integer/rational arithmetic.",
                      "scope":"Verifies the algebra and stabilized signatures, not infinitude."},indent=2))


if __name__ == "__main__":
    main()
