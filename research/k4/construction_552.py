#!/usr/bin/env python3
"""Exact finite verification for the partially split k=4 family with M=552.

This verifies finite integers, not infinitude. All factor products are checked.
SymPy's deterministic primality algorithm is used below 2^64. For larger
primes, a recursive Lucas certificate is produced and checked. The optional
--include-prime-proofs switch supplies Lucas certificates for every prime.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
import json
from math import prod
import sys

from sympy import factorint, isprime

D = 210 * 552**2
L = 210 * 552
Y0 = 23971
PROOFS: dict[int, dict] = {}


@lru_cache(maxsize=None)
def prove_prime(p: int) -> bool:
    """Converse of Fermat via a certified element of order p-1."""
    if p == 2:
        PROOFS[p] = {"p": 2, "base_case": True}
        return True
    fs = {int(q): int(e) for q, e in factorint(p - 1).items()}
    assert prod(q**e for q, e in fs.items()) == p - 1
    assert all(prove_prime(q) for q in fs)
    a = 2
    while a < p:
        if pow(a, p-1, p) == 1 and all(pow(a, (p-1)//q, p) != 1 for q in fs):
            PROOFS[p] = {"p": p, "lucas_witness": a,
                         "p_minus_one_factors": [[q, e] for q, e in sorted(fs.items())]}
            return True
        a += 1
    raise AssertionError(f"Composite alleged prime: {p}")


def factor(n: int) -> dict[int, int]:
    fs = {int(p): int(e) for p, e in factorint(n).items()}
    assert prod(p**e for p, e in fs.items()) == n
    assert all(e > 0 for e in fs.values())
    for p in fs:
        assert isprime(p)
        if p >= 2**64:
            assert prove_prime(p)
    return fs


def factorial_valuation(n: int, p: int) -> int:
    total = 0
    while n:
        n //= p
        total += n
    return total


def doubling_carries(n: int, p: int) -> int:
    total = incoming = 0
    while n or incoming:
        digit = n % p
        incoming = int(2 * digit + incoming >= p)
        total += incoming
        n //= p
    return total


def certificate(t: int, step: int = L, residue: int = Y0,
                include_prime_proofs: bool = False) -> dict:
    if t < 0 or step <= 0 or step % L:
        raise ValueError("t must be nonnegative and step a positive multiple of 115920")
    y = step*t + residue
    assert y > 113
    assert (y*y-1681) % 552 == 0
    m = (y*y-1681)//552
    assert m % 210 == 195
    assert m*(m+1) % 210 == 0
    n = m*(m+1)//210 - 3
    factors = [
        [y*y-12721, y*y+9911],
        [y-97, y+97, y*y+6599],
        [y-41, y+41, y*y-1129],
        [y*y-138*y+8117, y*y+138*y+8117],
    ]
    fac = Counter()
    factor_rows = []
    for i, values in enumerate(factors, 1):
        assert prod(values) == D*(n+i)
        shift = Counter()
        for value in values:
            fs = factor(value)
            shift.update(fs)
            factor_rows.append({"shift": i, "value": value,
                                "prime_factors": [[p,e] for p,e in sorted(fs.items())]})
        for p, e in {2:7, 3:3, 5:1, 7:1, 23:2}.items():
            shift[p] -= e
        assert all(e >= 0 for e in shift.values())
        assert prod(p**e for p,e in shift.items()) == n+i
        fac.update({p:e for p,e in shift.items() if e})
    assert prod(p**e for p,e in fac.items()) == prod(n+i for i in range(1,5))
    rows = []
    for p,e in sorted(fac.items()):
        available = factorial_valuation(2*n,p)-2*factorial_valuation(n,p)
        assert available == doubling_carries(n,p)
        margin = factorial_valuation(2*n,p)-2*factorial_valuation(n+4,p)
        assert margin == available-2*e
        rows.append({"p":p,"product_valuation":e,"required":2*e,
                     "available_carries":available,"factorial_margin":margin})
    out = {"k":4,"t":t,"step":step,"residue":residue,"y":y,"m":m,"n":n,
           "passes":all(row["factorial_margin"] >= 0 for row in rows),
           "factorizations":factor_rows,"prime_checks":rows,
           "scope":"Exact finite certificate; not an infinitude proof."}
    if include_prime_proofs:
        for p in fac:
            assert prove_prime(p)
        needed = set(fac)
        todo = list(fac)
        while todo:
            p = todo.pop()
            for q,e in PROOFS[p].get("p_minus_one_factors",[]):
                if q not in needed:
                    needed.add(q)
                    todo.append(q)
        out["lucas_prime_certificates"] = [PROOFS[p] for p in sorted(needed)]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=0)
    parser.add_argument("--step", type=int, default=L)
    parser.add_argument("--residue", type=int, default=Y0)
    parser.add_argument("--scan-stop", type=int)
    parser.add_argument("--include-prime-proofs", action="store_true")
    args = parser.parse_args()
    if args.scan_stop is None:
        out = certificate(args.t,args.step,args.residue,args.include_prime_proofs)
    else:
        if args.scan_stop < 0:
            parser.error("scan-stop must be nonnegative")
        passed = []
        for t in range(args.scan_stop):
            row = certificate(t,args.step,args.residue)
            if row["passes"]:
                passed.append({k:row[k] for k in ("t","y","m","n")})
            if (t+1)%100 == 0:
                print(f"Checked {t+1}; successes {len(passed)}",file=sys.stderr)
        out = {"interval":[0,args.scan_stop],"upper_endpoint_excluded":True,
               "step":args.step,"residue":args.residue,"count":len(passed),
               "witnesses":passed,
               "scope":"Finite exact scan only; no asymptotic density claim."}
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
