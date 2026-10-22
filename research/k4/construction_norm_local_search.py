#!/usr/bin/env python3
"""Exact stabilized large-prime residue search for small signed norm families.

The six old factors use m+c, c=0,1,-14,15,-20,21. Actual rational-linear
factors are omitted from the quadratic high-prime budget. Both signs of m
and all integral residue classes are included. New-factor phases have the
reduced denominator 840*Q^2, independent of P.
"""
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import json

from sympy import factorint, kronecker_symbol

CS=(0,1,-14,15,-20,21)
U210=[u for u in range(210) if gcd(u,210)==1]


def bad(s,mod,sign):
    if s==0:
        return sign>0
    if 2*s==mod:
        return sign<0
    return 2*s<mod


def squareclasses(modulus):
    out=[0]
    for p,e in factorint(modulus).items():
        q=int(p**e)
        sq={u*u%q for u in range(1,q) if gcd(u,q)==1}
        basis=(modulus//q)*pow(modulus//q,-1,q)
        out=[(a+b*basis)%modulus for a in out for b in sq]
    return out


def character_on_210(constant):
    fs=factorint(abs(-constant))
    d=prod(int(p) for p,e in fs.items() if e%2)
    if -constant<0:
        d=-d
    disc=d if d%4==1 else 4*d
    if 210%abs(disc):
        return None
    return [int(kronecker_symbol(disc,u)) for u in U210]


def search(P,Q,N,M,C):
    period=210*abs(M)*Q
    R=840*Q*Q
    squares=squareclasses(R)
    constants=[N+c*M for c in CS]
    split=[z<=0 and isqrt(-z)**2==-z for z in constants]
    signs=[1 if (1-2*c)*M>0 else -1 for c in CS]
    characters=[None if s else character_on_210(z) for s,z in zip(split,constants)]
    old={}
    for m0 in range(210):
        if m0*(m0+1)%210:
            continue
        counts=[]
        for c,s,sign,chars in zip(CS,split,signs,characters):
            if s:
                counts.append(Fraction(0))
            else:
                count=sum(bad((m0+c)**2*u*u%210,210,sign)*(1 if chars is None else 1+chars[k]) for k,u in enumerate(U210))
                counts.append(Fraction(count,48))
        old[m0]=counts

    @lru_cache(maxsize=None)
    def newbad(K,sign):
        return Fraction(sum(bad(K*s%R,R,sign) for s in squares),len(squares))

    best=None
    bestrows=[]
    admissible=0
    for Y in range(period):
        if (Y*Y+N)%M:
            continue
        m=(Y*Y+N)//M
        m0=m%210
        if m0 not in old:
            continue
        admissible+=1
        oldcounts=old[m0]
        newcounts=[]
        for sigma in (-1,1):
            A=2*Q*m+P+Q+sigma*Y
            assert (Y*Y+sigma*2*P*Y+C)==2*P*A
            newcounts.append(newbad(A*A%R,-1 if sigma*P>0 else 1))
        total=sum(oldcounts)+sum(newcounts)
        if best is None or total<best:
            best=total;bestrows=[]
        if total==best and len(bestrows)<16:
            bestrows.append({"Y":Y,"m_mod210":m0,
                             "old_fractions":[str(x) for x in oldcounts],
                             "new_fractions":[str(x) for x in newcounts]})
    return {"P":P,"Q":Q,"M":M,"N":N,"C":C,"period":period,
            "new_phase_denominator":R,"integral_residues":admissible,
            "old_split_c":[c for c,s in zip(CS,split) if s],
            "minimum_high_coefficient":str(best),"best_residues_first16":bestrows}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("--max-period",type=int,default=3000000)
    args=p.parse_args()
    out=[];skipped=[]
    for line in args.input.read_text().splitlines():
        P,Q,M,N,C,*_=map(int,line.split())
        if 210*abs(M)*Q>args.max_period:
            skipped.append([P,Q])
            continue
        out.append(search(P,Q,N,M,C))
    out.sort(key=lambda row:Fraction(row["minimum_high_coefficient"]) if row["minimum_high_coefficient"]!="None" else 100)
    print(json.dumps({"source":str(args.input),"max_period":args.max_period,
                      "rows":out,"skipped_parameters":skipped,
                      "scope":"Exact stabilized high-prime residue coefficients only."},indent=2))


if __name__=="__main__":
    main()
