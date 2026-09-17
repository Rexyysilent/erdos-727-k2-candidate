import math
import sys
from pathlib import Path
from fractions import Fraction
from collections import Counter
from decimal import Decimal
import pytest
TOOLS=Path(__file__).resolve().parents[1]/'tools'
sys.path.insert(0,str(TOOLS))
import family_search as f
import review_checks as c

@pytest.mark.parametrize('u,r',[(2,3),(2,5),(14,195),(14,209),(5,24),(5,29)])
def test_general_identities_and_primitive_factors(u,r):
    D,C=f.family_parameters(u,r)
    factors=f.primitive_factors(u,r)
    assert all(math.gcd(A,B)==1 for A,B in factors)
    for t in range(1,1000):
        n=f.family_value(u,r,t);m=D*t+r
        assert n==m*(m+1)//D-2
        assert (n+1)*D==(m-u)*(m+u+1)
        assert (n+1)*(n+2)==math.prod(A*t+B for A,B in factors)
        assert 4*D*n+8*D+1==(2*m+1)**2


def test_all_factor_gcd_products_through_u40():
    for u in range(1,41):
        D=u*(u+1)
        for r in range(D):
            if r*(r+1)%D:continue
            _,C=f.family_parameters(u,r)
            assert math.prod(math.gcd(D,r+c0) for c0 in C)==D*D

@pytest.mark.parametrize('u,r,b',[(2,3,Fraction(1)),(2,5,Fraction(1)),
                                (14,195,Fraction(2,3)),(14,209,Fraction(2,3))])
def test_exact_residue_coefficient(u,r,b):
    a=f.residue_analysis(u,r)
    assert Fraction(a['b_numerator'],a['b_denominator'])==b
    D,C=f.family_parameters(u,r)
    count=phi=0
    for q in range(D):
        if math.gcd(q,D)!=1:continue
        phi+=1
        for shift in C:
            s=(((r+shift)*q)**2)%D
            d=1-2*shift
            count+=f.misses_second_carry(s,D,d)
    assert Fraction(count,phi)==b


def test_D210_six_square_classes():
    assert f.unit_square_counts(210)==Counter({1:8,79:8,109:8,121:8,151:8,169:8})
    a=f.residue_analysis(14,209)
    assert [len(row['bad_shifts']) for row in a['rows']]==[2,1,0,0,0,1]

@pytest.mark.parametrize('u,r',[(2,3),(2,5),(14,195),(14,209)])
def test_signed_second_carry_rule_on_actual_integer_values(u,r):
    D,C=f.family_parameters(u,r)
    for q in range(D):
        if math.gcd(q,D)!=1:continue
        p=(1000000//D+1)*D+q
        while len(c.factor_trial(p))!=1 or c.factor_trial(p).get(p)!=1:p+=D
        inv=pow(p,-1,D)
        for shift in C:
            a=((r+shift)*inv)%D+D
            m=a*p-shift;t=(m-r)//D;n=f.family_value(u,r,t)
            assert m==D*t+r
            residue=a*a%D
            predicted=f.misses_second_carry(residue,D,1-2*shift)
            assert predicted==(2*(n%(p*p))<p*p)

@pytest.mark.parametrize('u,r,expected',[(2,3,126),(14,195,220),(14,209,223)])
def test_new_family_counts_1000_with_independent_carries(u,r,expected):
    assert sum(f.certificate(u,r,t)['passes'] for t in range(1000))==expected


def test_small_literal_factorial_checks():
    for u,r in [(2,3),(2,5),(14,195),(14,209)]:
        for t in range(8):
            n=f.family_value(u,r,t)
            if n>1000:continue
            assert f.certificate(u,r,t)['passes']==(math.factorial(2*n)%(math.factorial(n+2)**2)==0)


def test_distinct_D6_families():
    left={f.family_value(2,3,t) for t in range(1000)}
    right={f.family_value(2,5,t) for t in range(1000)}
    assert not left&right


def test_new_constants_are_certified():
    a=f.constant_enclosures()
    assert Decimal(a['rho_D210_b_2_over_3']['lower'])>Decimal('0.44232553173')
    assert Decimal(a['rho_D210_b_2_over_3']['upper'])<Decimal('0.44232553175')
    assert Decimal(a['two_D6_families_sqrtN_lower_coefficient']['lower'])>Decimal('0.21165508807')

@pytest.mark.parametrize('u,r',[(-1,0),(True,0),(2,1),(2,6),(2,-1),(1001,0)])
def test_invalid_family(u,r):
    with pytest.raises(ValueError):f.family_parameters(u,r)


def test_bounded_search20():
    a=f.search(20)
    assert a['minimum_b_in_this_search']=='2/3'
    assert {(x['u'],x['r']) for x in a['minimizers']}=={(14,195),(14,209)}
    with pytest.raises(ValueError):f.search(301)


def test_selected_ten_family_coefficients():
    pairs=f.selected_families()
    assert len(pairs)==len(set(pairs))==10
    for u,r in pairs:
        a=f.residue_analysis(u,r)
        expected=Fraction(2,3) if u==14 and r in (195,209) else Fraction(1)
        assert Fraction(a['b_numerator'],a['b_denominator'])==expected
    a=f.constant_enclosures()['ten_families_sqrtN_lower_coefficient']
    assert Decimal(a['lower'])>Decimal('0.4389093234051834561049915697270864665644')
    assert Decimal(a['upper'])<Decimal('0.4389093234051834561049915697270864665646')


def test_known_cross_denominator_intersections_need_deduplication():
    cases=[(143,2,5,4,14,174,0),(285,2,5,6,14,35,1),
           (13207,5,29,20,14,195,7),(20590,2,3,58,14,189,9)]
    for n,u,r,t,v,s,k in cases:
        Q=u*(u+1);R=v*(v+1)
        assert f.family_value(u,r,t)==f.family_value(v,s,k)==n
        y=2*(Q*t+r)+1;z=2*(R*k+s)+1
        assert R*y*y-Q*z*z==R-Q
        X=R*y;Y=z;d=Q*R;h=R*(R-Q)
        assert X*X-d*Y*Y==h


def test_pell_ratio_integral_unit_algebra():
    for d,h,left,right in [(2,1,(3,2),(17,12)),(2,-1,(1,1),(7,5)),
                           (2,4,(2,0),(6,4))]:
        x,y=left;X,Y=right
        assert x*x-d*y*y==X*X-d*Y*Y==h
        assert (X-x)%abs(h)==(Y-y)%abs(h)==0
        a_num=X*x-d*Y*y;b_num=Y*x-X*y
        assert a_num%h==b_num%h==0
        a,b=a_num//h,b_num//h
        assert a*a-d*b*b==1 and a>=2 and b>0


def test_small_ensemble_certificate(tmp_path):
    import json
    from ensemble_certificate import certify
    path=tmp_path/'search.json';path.write_text(json.dumps(f.search(2)))
    result=certify(path)
    assert result['selected_families']==2 and result['distinct_denominators']==1
    assert result['direct_fraction_checks']==16
    a=result['coefficient_interval']
    assert Decimal(a['lower'])>Decimal('0.2116550880770719563007065194241220574597')
    assert Decimal(a['upper'])<Decimal('0.2116550880770719563007065194241220574599')


def test_ensemble_rejects_modified_coefficient(tmp_path):
    import json
    from ensemble_certificate import certify
    source=f.search(2);source['shortlist_b_at_most_3_over_2'][0]['b_numerator']+=1
    path=tmp_path/'bad.json';path.write_text(json.dumps(source))
    with pytest.raises(ValueError):certify(path)


def test_ensemble_rejects_duplicate_family(tmp_path):
    import json
    from ensemble_certificate import certify
    source=f.search(2);source['shortlist_b_at_most_3_over_2']*=2
    path=tmp_path/'duplicate.json';path.write_text(json.dumps(source))
    with pytest.raises(ValueError):certify(path)
