#!/usr/bin/env python3
"""Finite exact certificates for signed norm-composition k=4 families.

Defaults to the improved downward M=432 family. The output proves only each
specified integer instance. It does not establish infinitely many instances.
"""
import argparse
from collections import Counter
import json
from math import isqrt, prod
import sys

from construction_552 import factor, factorial_valuation, doubling_carries, prove_prime, PROOFS


def certificate(t,P=108,Q=1,delta=2521,residue=16273,step=90720,proofs=False):
    assert t>=0 and Q>0 and P!=0 and delta in (2521,-839)
    M=4*P*Q
    N=-delta*Q*Q-P*P-2*P*Q
    C=P*P-delta*Q*Q
    D=210*M*M
    y=residue+step*t
    assert (y*y+N)%M==0
    m=(y*y+N)//M
    assert m*(m+1)%210==0
    offset=4 if delta==2521 else 3
    n=m*(m+1)//210-offset
    assert n>0
    groups={}
    for i,cs in ((offset,(0,1)),(offset-1,(-14,15)),(offset-2,(-20,21))):
        values=[]
        for c in cs:
            const=N+c*M
            if const<=0 and isqrt(-const)**2==-const:
                r=isqrt(-const)
                values.extend([y-r,y+r])
            else:
                values.append(y*y+const)
        groups[i]=values
    groups[1 if delta==2521 else 4]=[y*y-2*P*y+C,y*y+2*P*y+C]
    denominator=factor(D)
    factors=Counter()
    factorizations=[]
    for i,values in sorted(groups.items()):
        assert min(values)>0
        assert prod(values)==D*(n+i)
        local=Counter()
        for value in values:
            fs=factor(value)
            local.update(fs)
            factorizations.append({"shift":i,"value":value,
                                   "prime_factors":[[p,e] for p,e in sorted(fs.items())]})
        for p,e in denominator.items():
            local[p]-=e
        assert all(e>=0 for e in local.values())
        assert prod(p**e for p,e in local.items())==n+i
        factors.update({p:e for p,e in local.items() if e})
    assert prod(p**e for p,e in factors.items())==prod(n+i for i in range(1,5))
    rows=[]
    for p,e in sorted(factors.items()):
        available=factorial_valuation(2*n,p)-2*factorial_valuation(n,p)
        assert available==doubling_carries(n,p)
        margin=factorial_valuation(2*n,p)-2*factorial_valuation(n+4,p)
        assert margin==available-2*e
        rows.append({"p":p,"product_valuation":e,"required":2*e,
                     "available_carries":available,"factorial_margin":margin})
    out={"k":4,"P":P,"Q":Q,"delta":delta,"M":M,"N":N,"C":C,"D":D,
         "t":t,"step":step,"residue":residue,"y":y,"m":m,"n":n,
         "passes":all(row["factorial_margin"]>=0 for row in rows),
         "factorizations":factorizations,"prime_checks":rows,
         "scope":"Exact finite certificate; not an infinitude proof."}
    if proofs:
        for p in factors:
            assert prove_prime(p)
        needed=set(factors)
        todo=list(factors)
        while todo:
            p=todo.pop()
            for q,e in PROOFS[p].get("p_minus_one_factors",[]):
                if q not in needed:
                    needed.add(q);todo.append(q)
        out["lucas_prime_certificates"]=[PROOFS[p] for p in sorted(needed)]
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--t",type=int,default=0)
    p.add_argument("--P",type=int,default=108)
    p.add_argument("--Q",type=int,default=1)
    p.add_argument("--delta",type=int,default=2521)
    p.add_argument("--residue",type=int,default=16273)
    p.add_argument("--step",type=int,default=90720)
    p.add_argument("--scan-stop",type=int)
    p.add_argument("--include-prime-proofs",action="store_true")
    a=p.parse_args()
    if a.scan_stop is None:
        out=certificate(a.t,a.P,a.Q,a.delta,a.residue,a.step,a.include_prime_proofs)
    else:
        passed=[]
        for t in range(a.scan_stop):
            row=certificate(t,a.P,a.Q,a.delta,a.residue,a.step)
            if row["passes"]:
                passed.append({k:row[k] for k in ("t","y","m","n")})
            if (t+1)%100==0:
                print(f"Checked {t+1}; successes {len(passed)}",file=sys.stderr)
        out={"P":a.P,"Q":a.Q,"delta":a.delta,"residue":a.residue,"step":a.step,
             "interval":[0,a.scan_stop],"upper_endpoint_excluded":True,
             "count":len(passed),"witnesses":passed,
             "scope":"Finite exact scan only; no asymptotic density claim."}
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
