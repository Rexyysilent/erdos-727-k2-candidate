#!/usr/bin/env python3
"""Exact finite checks for the proposed k=3 extension, not proof certification."""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from math import comb,gcd
import json
from review_checks import natural,factor_trial,carry_count,factorial_v,product,dec
from family_search import unit_square_counts,misses_second_carry,constant_enclosures

SHIFTS=(0,1,-14,15,-20,21)
LINES=((6,5),(35,36),(1,1),(210,181),(14,13),(15,14))


def family(t: int) -> int:
    t=natural(t,'t')
    return 210*t*t+391*t+179


def certificate(t: int) -> dict:
    n=family(t);fs=Counter()
    for A,B in LINES:fs.update(factor_trial(A*t+B))
    assert product(p**e for p,e in fs.items())==(n+1)*(n+2)*(n+3)
    rows=[]
    for p,e in sorted(fs.items()):
        cp=carry_count(n,p)
        assert cp==factorial_v(2*n,p)-2*factorial_v(n,p)
        rows.append({'prime':p,'required':2*e,'available':cp,'passes':cp>=2*e})
    passes=all(x['passes'] for x in rows)
    if n<=10000:
        assert passes==(comb(2*n,n)%(((n+1)*(n+2)*(n+3))**2)==0)
    return {'k':3,'t':t,'m':210*t+195,'n':n,'passes':passes,'prime_checks':rows}


def residue_table() -> dict:
    squares=unit_square_counts(210);rows=[];total=0
    for s,w in sorted(squares.items()):
        fs=[((195+c)**2*s)%210 for c in SHIFTS]
        bad=[c for c,v in zip(SHIFTS,fs) if misses_second_carry(v,210,1-2*c)]
        total+=w*len(bad)
        rows.append({'inverse_prime_square_residue':s,'multiplicity':w,'fractional_numerators':fs,'bad_shifts':bad})
    return {'Q':210,'r':195,'shifts':SHIFTS,'phi':sum(squares.values()),
            'b_exact':str(Fraction(total,sum(squares.values()))),'rows':rows}


def constants() -> dict:
    a=constant_enclosures();S=a['lower_prime_cost_S'];rho=a['rho_original_D6']
    slo,shi=Fraction(S['lower']),Fraction(S['upper'])
    llo=1-shi-Fraction(rho['upper']);lhi=1-slo-Fraction(rho['lower'])
    lo=1-Fraction(3,2)*shi-Fraction(5,6)*lhi
    hi=1-Fraction(3,2)*slo-Fraction(5,6)*llo
    return {'formula':'rho3=1-(3/2)*S-(5/12)*log(3)',
        'rho3_interval':{'lower':dec(lo,45),'upper':dec(hi,45,upward=True)},
        'method':'outward rational arithmetic from previously certified logarithm/series intervals'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--t',type=int,default=0)
    p.add_argument('--residues',action='store_true');p.add_argument('--constants',action='store_true')
    args=p.parse_args()
    try:out=constants() if args.constants else residue_table() if args.residues else certificate(args.t)
    except ValueError as exc:p.error(str(exc))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
