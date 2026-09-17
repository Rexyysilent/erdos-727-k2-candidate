"""Finite arithmetic tests of the new k=3 companion, not its asymptotics."""
from pathlib import Path
import sys
from fractions import Fraction
from decimal import Decimal
from math import comb, factorial, gcd, prod
import pytest

TOOLS=Path(__file__).resolve().parents[1]/'tools'
sys.path.insert(0,str(TOOLS))
import verify_k3 as k
import review_checks as c


def test_six_linear_factorization_and_square_identity():
    assert all(gcd(A,B)==1 for A,B in k.LINES)
    for t in range(2000):
        n=k.family(t);m=210*t+195
        assert 210*(n+3)==m*(m+1)
        assert 210*(n+2)==(m-14)*(m+15)
        assert 210*(n+1)==(m-20)*(m+21)
        assert (6*t+5)*(35*t+36)==n+1
        assert (t+1)*(210*t+181)==n+2
        assert (14*t+13)*(15*t+14)==n+3
        assert 840*n+2521==(420*t+391)**2
        assert prod(A*t+B for A,B in k.LINES)==(n+1)*(n+2)*(n+3)

@pytest.mark.parametrize('p',[2,3,5,7])
def test_permutation_at_fixed_denominator_primes(p):
    q=p**3
    assert len({k.family(t)%q for t in range(q)})==q
    assert all((420*t+391)%p for t in range(p))

@pytest.mark.parametrize('p',[43,47,53])
def test_six_simple_roots_and_bijective_lifting(p):
    roots=[((-195-shift)*pow(210,-1,p))%p for shift in k.SHIFTS]
    assert len(set(roots))==6
    for shift,t0 in zip(k.SHIFTS,roots):
        i=3 if shift in (0,1) else 2 if shift in (-14,15) else 1
        assert k.family(t0)%p==(-i)%p
        assert (420*t0+391)%p==(1-2*shift)%p!=0
        values={k.family(t0+p*a)%(p**3) for a in range(p*p)}
        assert len(values)==p*p
        assert all(n%p==(-i)%p for n in values)


def test_residue_table_and_exact_average():
    out=k.residue_table()
    assert out['b_exact']=='5/6'
    assert out['phi']==48
    assert [row['inverse_prime_square_residue'] for row in out['rows']]==[1,79,109,121,151,169]
    assert [len(row['bad_shifts']) for row in out['rows']]==[3,1,0,0,0,1]
    assert all(row['multiplicity']==8 for row in out['rows'])
    count=0
    for v in range(210):
        if gcd(v,210)!=1:continue
        for shift in k.SHIFTS:
            s=((195+shift)*v)**2%210
            d=1-2*shift
            count+=0<s<105 or (s==0 and d>0) or (s==105 and d<0)
    assert Fraction(count,48)==Fraction(5,6)


def test_signed_fraction_rule_on_actual_integer_values():
    for v in range(210):
        if gcd(v,210)!=1:continue
        p=210*100000+v
        inv=pow(p,-1,210)
        for shift in k.SHIFTS:
            a=((195+shift)*inv)%210+210
            m=a*p-shift
            t=(m-195)//210
            assert m==210*t+195
            n=k.family(t)
            s=(a*a)%210;d=1-2*shift
            predicted=0<s<105 or (s==0 and d>0) or (s==105 and d<0)
            assert predicted==(2*(n%(p*p))<p*p)


def test_first_1000_and_first_witness():
    good=[k.certificate(t) for t in range(1000)]
    good=[row for row in good if row['passes']]
    assert len(good)==126
    assert (good[0]['t'],good[0]['n'])==(61,805440)
    assert [x['n'] for x in good[:8]]==[805440,1056549,1852832,1892493,2055337,2678211,3277623,3880214]
    assert all(x['available']>=x['required'] for x in good[0]['prime_checks'])


def test_literal_factorials_and_binomial_criterion_up_to_500():
    for n in range(501):
        E=(n+1)*(n+2)*(n+3)
        literal=factorial(2*n)%(factorial(n+3)**2)==0
        binomial=comb(2*n,n)%(E*E)==0
        fs=c.factor_trial(n+1)+c.factor_trial(n+2)+c.factor_trial(n+3)
        carries=all(c.carry_count(n,p)>=2*e for p,e in fs.items())
        assert literal==binomial==carries


def test_rho3_interval_is_positive_and_encloses_formula():
    result=k.constants()['rho3_interval']
    assert Decimal('0.25503932166219664924')<Decimal(result['lower'])
    assert Decimal(result['upper'])<Decimal('0.25503932166219664925')
    a=k.constant_enclosures()['lower_prime_cost_S']
    assert Fraction(a['upper'])<Fraction(24,125)
    assert c.ln_interval(3,1)[1]<Fraction(11,10)
    assert 1-Fraction(3,2)*Fraction(24,125)-Fraction(5,12)*Fraction(11,10)>Fraction(1,4)

@pytest.mark.parametrize('value',[-1,True,1.0,'61',None])
def test_invalid_parameter(value):
    with pytest.raises(ValueError):k.certificate(value)
