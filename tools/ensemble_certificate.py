#!/usr/bin/env python3
"""Revalidate a finite family ensemble and certify its sqrt(N) coefficient.

The analytic theorem and intersection lemma remain mathematical obligations.
This script certifies the bounded residue arithmetic and constant arithmetic,
not a formal proof of those analytic statements.
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from family_search import residue_analysis, family_parameters, constant_enclosures
from review_checks import dec


def direct_residue_average(u: int, r: int) -> tuple[Fraction,int]:
    D,C=family_parameters(u,r)
    total=units=0
    for q in range(D):
        if gcd(q,D)!=1:continue
        units+=1;p=10**6*D*D+q;inv=pow(p,-1,D)
        for c in C:
            a=((r+c)*inv)%D+D
            m=a*p-c
            assert (m-r)%D==0 and m*(m+1)%D==0
            n=m*(m+1)//D-2
            total+=2*(n%(p*p))<p*p
    return Fraction(total,units),units


def certify(path: Path) -> dict:
    if path.stat().st_size>8*1024*1024:
        raise ValueError('search input exceeds 8 MiB')
    source=json.loads(path.read_text(encoding='utf-8'))
    if source.get('u_min')!=1 or not 1<=source.get('u_max_inclusive',0)<=150:
        raise ValueError('require bounded source search starting at u=1, ending at <=150')
    constants=constant_enclosures();S=constants['lower_prime_cost_S']
    rho=constants['rho_original_D6']
    slo,shi=Fraction(S['lower']),Fraction(S['upper'])
    llo=1-shi-Fraction(rho['upper']);lhi=1-slo-Fraction(rho['lower'])
    assert 1-slo-Fraction(3,2)*llo<0
    selected=[];seen=set();unit_checks=0;scale=10**60
    total_lo=total_hi=0
    for row in source['shortlist_b_at_most_3_over_2']:
        u,r=row['u'],row['r']
        if (u,r) in seen:raise ValueError('duplicate family in supplied search')
        seen.add((u,r))
        if not 1<=u<=source['u_max_inclusive']:raise ValueError('out-of-scope u')
        b=Fraction(row['b_numerator'],row['b_denominator'])
        actual=residue_analysis(u,r)
        if b!=Fraction(actual['b_numerator'],actual['b_denominator']):
            raise ValueError('saved coefficient failed fresh square-class enumeration')
        lo,hi=1-shi-b*lhi,1-slo-b*llo
        if hi<=0:continue
        if lo<=0:raise ValueError('sign of bound is unresolved at the chosen precision')
        other,units=direct_residue_average(u,r);unit_checks+=units
        if other!=b:raise ValueError('direct finite fraction calculation disagrees')
        D=u*(u+1);s=isqrt(D*scale*scale)
        low=lo/Fraction(s+1,scale);high=hi/Fraction(s,scale)
        low_i=(low.numerator*scale)//low.denominator
        high_i=-((-high.numerator*scale)//high.denominator)
        total_lo+=low_i;total_hi+=high_i
        selected.append({'u':u,'D':D,'r':r,'b_exact':str(b),
            'parameter_lower_bound_interval':{'lower':dec(lo,45),'upper':dec(hi,45,upward=True)},
            'sqrtN_contribution_interval':{'lower':dec(Fraction(low_i,scale),45),
                                           'upper':dec(Fraction(high_i,scale),45,upward=True)},
            'direct_unit_classes_checked':units})
    return {'source_u_max_inclusive':source['u_max_inclusive'],
            'source_root_classes_checked':source['root_classes_checked'],
            'selected_families':len(selected),'distinct_denominators':len({x['D'] for x in selected}),
            'direct_unit_class_checks':unit_checks,'direct_fraction_checks':4*unit_checks,
            'coefficient_interval':{'lower':dec(Fraction(total_lo,scale),45),
                                    'upper':dec(Fraction(total_hi,scale),45,upward=True)},
            'families':selected,
            'scope':'Exact finite residue/constant certificate; analytic lower bounds and O(log N) intersections require the supplied proof. No claim of global optimality.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('search_json',type=Path)
    args=p.parse_args()
    try:result=certify(args.search_json)
    except (ValueError,OSError,KeyError,TypeError) as exc:p.error(str(exc))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
