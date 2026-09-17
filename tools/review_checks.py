#!/usr/bin/env python3
"""Exact, independent finite checks for the 727 k=2 review.

No random claim, floating-point divisibility, network, or asymptotic certification.
The original repository verifier is left unchanged. This companion recomputes
prime factors by trial division and carries by explicit digit addition.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
from math import comb, factorial, isqrt


def natural(value: int, name: str = 'value') -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f'{name} must be a nonnegative integer')
    return value


def family(t: int) -> int:
    t = natural(t, 't')
    return 6*t*t + 11*t + 3


def factor_trial(n: int) -> Counter[int]:
    natural(n, 'n')
    if n == 0:
        raise ValueError('cannot factor zero')
    out: Counter[int] = Counter()
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] += 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] += 1
    return out


def carry_count(n: int, base: int) -> int:
    """Outgoing carries in n+n in any integer base >=2."""
    natural(n, 'n')
    if type(base) is not int or base < 2:
        raise ValueError('base must be an integer >=2')
    incoming = count = 0
    while n or incoming:
        n, digit = divmod(n, base)
        incoming = (2*digit + incoming) // base
        count += incoming
    return count


def factorial_v(n: int, p: int) -> int:
    """Legendre sum; callers supply a prime p from exact factorization."""
    natural(n, 'n')
    if type(p) is not int or p < 2:
        raise ValueError('p must be an integer >=2')
    total = 0
    while n:
        n //= p
        total += n
    return total


def certificate_for_t(t: int) -> dict:
    n = family(t)
    factors: Counter[int] = Counter()
    for v in (2*t+1, 3*t+4, t+1, 6*t+5):
        factors.update(factor_trial(v))
    assert product(p**e for p,e in factors.items()) == (n+1)*(n+2)
    rows = []
    for p,e in sorted(factors.items()):
        cp = carry_count(n,p)
        assert cp == factorial_v(2*n,p)-2*factorial_v(n,p)
        rows.append({'prime':p,'denominator_valuation':e,'required':2*e,'available':cp,'passes':cp>=2*e})
    good=all(r['passes'] for r in rows)
    if n<=10000:
        assert good == (comb(2*n,n) % (((n+1)*(n+2))**2) == 0)
    return {'t':t,'m':6*t+5,'n':n,'passes':good,'prime_checks':rows}


def product(values):
    result=1
    for v in values:
        result*=v
    return result


def large_obstruction(m: int, X: int) -> dict:
    """The addendum's exact finite event; X>=10, m in M_X."""
    natural(X,'X');natural(m,'m')
    if X<10 or not X<=m<=2*X or m%6!=5:
        raise ValueError('require X>=10 and m in [X,2X], m=5 mod 6')
    n=m*(m+1)//6-2
    ps=[p for p in factor_trial(m) if p>=7 and p*p>2*X+3]
    assert len(ps)<=1
    if not ps:
        return {'prime':None,'bad':False}
    p=ps[0]
    assert 2*(n%(p*p))<p*p
    assert 2*n<p**4
    bad=2*(n%(p**3))<p**3
    assert bad==(carry_count(n,p)==1)
    return {'prime':p,'bad':bad}


def ln_interval(num: int, den: int, *, terms: int=140, scale: int=10**70) -> tuple[Fraction,Fraction]:
    """Rational enclosure via log x = 2*atanh((x-1)/(x+1))."""
    if type(num) is not int or type(den) is not int or not 0<den<=num<=3*den:
        raise ValueError('require integers with 1 <= num/den <= 3')
    if type(terms) is not int or not 1<=terms<=1000 or type(scale) is not int or scale<1:
        raise ValueError('invalid series precision')
    a,b=num-den,num+den
    if a==0:return Fraction(0),Fraction(0)
    an,bn=a,b
    total=0
    for k in range(terms):
        total += (2*scale*an)//((2*k+1)*bn)
        an*=a*a;bn*=b*b
    tail=Fraction(2*an*b*b,(2*terms+1)*bn*(b*b-a*a))
    return Fraction(total,scale),Fraction(total+terms,scale)+tail


def dec(frac: Fraction, places: int=60, *, upward: bool=False) -> str:
    with localcontext() as ctx:
        ctx.prec=places+40
        ctx.rounding=ROUND_CEILING if upward else ROUND_FLOOR
        result=Decimal(frac.numerator)/Decimal(frac.denominator)
        return str(result.quantize(Decimal(1).scaleb(-places)))


def constants(R: int=200) -> dict:
    if type(R) is not int or not 4<=R<=400:
        raise ValueError('require 4 <= R <= 400')
    l3,u3=ln_interval(3,1)
    low=high=Fraction(0)
    for j in range(4,R+1):
        a,b=ln_interval(j+1,j)
        weight=Fraction(4,2**(j-1))
        low+=weight*a;high+=weight*b
    tail=Fraction(1,2**(R-3)*(R+1))
    rho_lo,rho_hi=1-u3/2-high-tail,1-l3/2-low
    u_lo,u_hi=1-u3/2,1-l3/2
    scale=10**70;s=isqrt(6*scale*scale)
    sqrt_lo,sqrt_hi=Fraction(s,scale),Fraction(s+1,scale)
    def interval(a,b):
        return {'lower':dec(a),'upper':dec(b,upward=True),'width_upper':dec(b-a,upward=True)}
    return {'method':'exact rational atanh log bounds plus explicit positive series tail; displayed decimals are outward-rounded',
        'series_last_J':R,'series_tail_upper':dec(tail,upward=True),
        'rho':interval(rho_lo,rho_hi),'large_obstruction_density':interval(l3/2,u3/2),
        'family_upper_density':interval(u_lo,u_hi),
        'simple_obstruction_mean':interval(1-rho_hi,1-rho_lo),
        'family_solution_sqrtN_lower_coefficient':interval(rho_lo/sqrt_hi,rho_hi/sqrt_lo),
        'family_solution_sqrtN_upper_coefficient':interval(u_lo/sqrt_hi,u_hi/sqrt_lo)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int);p.add_argument('--constants',action='store_true')
    args=p.parse_args()
    try:
        result=constants() if args.constants or args.t is None else certificate_for_t(args.t)
    except ValueError as exc:p.error(str(exc))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
