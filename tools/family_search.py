#!/usr/bin/env python3
"""Exact finite residue search for the generalized quadratic-family argument.

This enumerates the coefficient in a proposed analytic lower bound, not measured
success rates. No claim of optimality beyond the explicitly searched u range.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from math import gcd
import json


def family_parameters(u: int, r: int) -> tuple[int, tuple[int, ...]]:
    if type(u) is not int or not 1 <= u <= 1000:
        raise ValueError('require integer 1 <= u <= 1000')
    D = u * (u + 1)
    if type(r) is not int or not 0 <= r < D or r*(r+1) % D:
        raise ValueError('r must be a root of r(r+1)=0 modulo D, in [0,D)')
    return D, (0, 1, -u, u+1)


def family_value(u: int, r: int, t: int) -> int:
    D, _ = family_parameters(u, r)
    if type(t) is not int or t < 0:
        raise ValueError('t must be a nonnegative integer')
    return D*t*t + (2*r+1)*t + r*(r+1)//D - 2


def primitive_factors(u: int, r: int) -> list[tuple[int, int]]:
    D, shifts = family_parameters(u, r)
    return [(D//gcd(D, r+c), (r+c)//gcd(D, r+c)) for c in shifts]


def misses_second_carry(residue: int, D: int, d: int) -> bool:
    if not 0 <= residue < D or not d:
        raise ValueError('invalid residue or zero correction coefficient')
    if residue == 0:
        return d > 0
    if 2*residue == D:
        return d < 0
    return 2*residue < D


def unit_square_counts(D: int) -> Counter[int]:
    return Counter(q*q % D for q in range(D) if gcd(q,D)==1)


def residue_analysis(u: int, r: int, squares: Counter[int] | None=None) -> dict:
    D, shifts = family_parameters(u,r)
    squares = unit_square_counts(D) if squares is None else squares
    phi = sum(squares.values())
    rows = []
    bad_sum = 0
    for s,multiplicity in sorted(squares.items()):
        residues = [((r+c)**2*s) % D for c in shifts]
        bad_shifts = [c for c,f in zip(shifts,residues)
                      if misses_second_carry(f,D,1-2*c)]
        bad_sum += multiplicity*len(bad_shifts)
        rows.append({'inverse_prime_square_residue':s, 'unit_multiplicity':multiplicity,
                     'fractional_numerators':residues,'bad_shifts':bad_shifts})
    b = Fraction(bad_sum,phi)
    return {'u':u,'D':D,'r':r,'polynomial_coefficients':[D,2*r+1,r*(r+1)//D-2],
            'primitive_linear_factors':primitive_factors(u,r),
            'shifts':list(shifts),'phi_D':phi,'b_numerator':b.numerator,
            'b_denominator':b.denominator,'b_exact':str(b),'rows':rows}


def search(max_u: int) -> dict:
    if type(max_u) is not int or not 1 <= max_u <= 300:
        raise ValueError('bounded search requires 1 <= max_u <= 300')
    compact = []
    best: Fraction | None = None
    winners = []
    roots_count = 0
    for u in range(1,max_u+1):
        D = u*(u+1)
        squares = unit_square_counts(D)
        for r in range(D):
            if r*(r+1) % D:
                continue
            roots_count += 1
            result = residue_analysis(u,r,squares)
            b = Fraction(result['b_numerator'],result['b_denominator'])
            if b <= Fraction(3,2):
                compact.append({k:result[k] for k in ('u','D','r','polynomial_coefficients',
                                                      'b_numerator','b_denominator','b_exact')})
            if best is None or b<best:
                best,winners=b,[result]
            elif b==best:
                winners.append(result)
    return {'u_min':1,'u_max_inclusive':max_u,'root_classes_checked':roots_count,
            'minimum_b_in_this_search':str(best),'minimizers':winners,
            'shortlist_b_at_most_3_over_2':compact,
            'scope':'Exact bounded residue enumeration only. No global optimality or finite success-density claim.'}


def certificate(u: int, r: int, t: int) -> dict:
    from review_checks import factor_trial, carry_count, factorial_v, product
    n=family_value(u,r,t)
    if n<0:
        raise ValueError('n must be nonnegative')
    fs=Counter()
    for A,B in primitive_factors(u,r):
        v=A*t+B
        if v<1:
            raise ValueError('linear factors must be positive')
        fs.update(factor_trial(v))
    assert product(p**e for p,e in fs.items())==(n+1)*(n+2)
    checks=[]
    for p,e in sorted(fs.items()):
        cp=carry_count(n,p)
        assert cp==factorial_v(2*n,p)-2*factorial_v(n,p)
        checks.append({'prime':p,'required':2*e,'available':cp,'passes':cp>=2*e})
    return {'u':u,'r':r,'t':t,'n':n,'passes':all(x['passes'] for x in checks),'prime_checks':checks}


def constant_enclosures() -> dict:
    from review_checks import ln_interval, dec
    from math import isqrt
    R=200; l3,h3=ln_interval(3,1)
    low=high=Fraction(0)
    for j in range(4,R+1):
        a,b=ln_interval(j+1,j);weight=Fraction(4,2**(j-1))
        low+=weight*a;high+=weight*b
    high+=Fraction(1,2**(R-3)*(R+1))
    def interval(a,b):return {'lower':dec(a),'upper':dec(b,upward=True)}
    rho=(1-h3/2-high,1-l3/2-low)
    rho210=(1-h3/3-high,1-l3/3-low)
    scale=10**70;s=isqrt(6*scale*scale)
    return {'method':'exact rational logarithm intervals; J<=200 plus rigorous series tail',
        'lower_prime_cost_S':interval(low,high),
        'rho_original_D6':interval(*rho),
        'rho_D210_b_2_over_3':interval(*rho210),
        'improvement_in_parameter_lower_bound':interval(l3/6,h3/6),
        'two_D6_families_sqrtN_lower_coefficient':interval(2*rho[0]/Fraction(s+1,scale),2*rho[1]/Fraction(s,scale)),
        'ten_families_sqrtN_lower_coefficient':combined_counting_interval(rho,rho210,interval,scale)}


def selected_families() -> list[tuple[int,int]]:
    return [(2,3),(2,5),(5,24),(5,29)]+[(14,r) for r in (20,35,174,189,195,209)]


def combined_counting_interval(rho, rho210, interval, scale):
    from math import isqrt
    low=high=Fraction(0)
    for Q,multiplicity,bound in [(6,2,rho),(30,2,rho),(210,4,rho),(210,2,rho210)]:
        s=isqrt(Q*scale*scale)
        low+=multiplicity*bound[0]/Fraction(s+1,scale)
        high+=multiplicity*bound[1]/Fraction(s,scale)
    return interval(low,high)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--u',type=int,default=14);p.add_argument('--r',type=int,default=209)
    p.add_argument('--search-max-u',type=int)
    p.add_argument('--constants',action='store_true')
    p.add_argument('--t',type=int)
    args=p.parse_args()
    try:
        result=(constant_enclosures() if args.constants else certificate(args.u,args.r,args.t)
                if args.t is not None else search(args.search_max_u)
                if args.search_max_u is not None else residue_analysis(args.u,args.r))
    except ValueError as exc:p.error(str(exc))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
