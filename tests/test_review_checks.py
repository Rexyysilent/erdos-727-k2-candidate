import importlib.util
import math
import sys
from pathlib import Path
from decimal import Decimal
from collections import Counter
import pytest

TOOL=Path(__file__).resolve().parents[1]/'tools/review_checks.py'
spec=importlib.util.spec_from_file_location('checks727',TOOL)
c=importlib.util.module_from_spec(spec);sys.modules[spec.name]=c;spec.loader.exec_module(c)

@pytest.mark.parametrize('t',[-1,True,1.5,'5'])
def test_invalid_parameters(t):
    with pytest.raises(ValueError):c.family(t)

@pytest.mark.parametrize('n,base',[(-1,2),(2,1),(2,True),(False,2)])
def test_invalid_carry_inputs(n,base):
    with pytest.raises(ValueError):c.carry_count(n,base)


def test_exact_polynomial_identities():
    for t in range(2000):
        n=c.family(t);m=6*t+5
        assert n==m*(m+1)//6-2
        assert n+1==(2*t+1)*(3*t+4)
        assert n+2==(t+1)*(6*t+5)
        assert 24*n+49==(12*t+11)**2

@pytest.mark.parametrize('p',[2,3,5,7,11,19,101])
def test_legendre_floor_and_explicit_digit_carries(p):
    for n in range(2000):
        v=c.factorial_v(2*n,p)-2*c.factorial_v(n,p)
        q=p;floor_sum=0
        while q<=2*n:
            floor_sum+=int(2*(n%q)>=q);q*=p
        assert c.carry_count(n,p)==v==floor_sum


def test_literal_factorials_all_n_through_500():
    for n in range(501):
        fs=c.factor_trial(n+1)+c.factor_trial(n+2)
        good=all(c.carry_count(n,p)>=2*e for p,e in fs.items())
        assert good==(math.factorial(2*n)%(math.factorial(n+2)**2)==0)
        assert good==(math.comb(2*n,n)%(((n+1)*(n+2))**2)==0)


def test_original_witness_and_first_eight():
    r=c.certificate_for_t(5)
    assert r['n']==208 and r['passes']
    assert [x['available'] for x in r['prime_checks']]==[3,3,2,2,2,2]
    good=[c.family(t) for t in range(1000) if c.certificate_for_t(t)['passes']]
    assert len(good)==120
    assert good[:8]==[208,3430,8175,17440,19435,34578,46374,47435]

@pytest.mark.parametrize('p',[2,3])
def test_fixed_prime_permutation_small_powers(p):
    for L in range(1,9):
        q=p**L
        assert len({c.family(t)%q for t in range(q)})==q

@pytest.mark.parametrize('p',[5,7])
def test_square_fiber_bound(p):
    H=1;L=3;q=p**L
    counts=Counter(c.family(t)%q for t in range(q) if (12*t+11)%p)
    assert max(counts.values())<=2*p**H

@pytest.mark.parametrize('p',[7,11,13])
def test_simple_root_lift_bijections(p):
    q=p**3
    roots=[t for t in range(p) if any(v%p==0 for v in (2*t+1,3*t+4,t+1,6*t+5))]
    assert len(roots)==4
    for t0 in roots:
        assert (12*t0+11)%p in {1,p-1,5,p-5}
        values={c.family(t)%q for t in range(t0,q,p)}
        assert len(values)==p*p
        assert len({v%p for v in values})==1

@pytest.mark.parametrize('p',[2,3,5,7])
def test_low_carry_digit_upper_bound(p):
    for L in range(1,5):
        for K in range(3):
            actual=sum(c.carry_count(n,p)<=K for n in range(p**L))
            bound=sum(math.comb(L,r)*(p//2)**r*((p+1)//2)**(L-r) for r in range(min(K,L)+1))
            assert actual<=bound

@pytest.mark.parametrize('X',[25,100,1000,10000])
def test_large_prime_event_exact_finite_equivalence(X):
    for m in range(X,2*X+1):
        if m%6!=5:continue
        r=c.large_obstruction(m,X)
        if r['bad']:
            assert not c.certificate_for_t((m-5)//6)['passes']


def test_rational_log_bounds_and_rho():
    a,b=c.ln_interval(3,1)
    assert a<b and b-a<c.Fraction(1,10**60)
    x=c.constants()
    assert Decimal(x['rho']['lower'])>Decimal('0.259223483626328945854216378517358130655998348891588306056952')
    assert Decimal(x['rho']['upper'])<Decimal('0.259223483626328945854216378517358130655998348891588306056955')
    a=c.constants(40);b=c.constants(100)
    assert Decimal(a['rho']['lower'])<=Decimal(b['rho']['lower'])
    assert Decimal(b['rho']['upper'])<=Decimal(a['rho']['upper'])

@pytest.mark.parametrize('args',[(0,1),(1,0),(4,1),(True,1)])
def test_log_interval_input_contract(args):
    with pytest.raises(ValueError):c.ln_interval(*args)
