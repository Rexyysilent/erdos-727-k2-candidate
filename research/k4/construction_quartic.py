#!/usr/bin/env python3
"""Exact quartic-family arithmetic for Erdős 727, k=4.

This file certifies finite instances and factorization identities. It does not
prove infinitude. Requires SymPy for integer factorization.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from math import prod

from sympy import factorint


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


def certificate(t: int) -> dict:
    if t < 0:
        raise ValueError("t must be nonnegative")
    x = 210 * t + 14
    m = x * x + 209
    assert m % 210 == 195
    numerator = m * (m + 1)
    assert numerator % 210 == 0
    n = numerator // 210 - 3
    pairs = [
        (x*x + 189, x*x + 230),
        (x*x + 195, x*x + 224),
        (x*x + 209, x*x + 210),
        (x*x - x + 210, x*x + x + 210),
    ]
    fac = Counter()
    rows = []
    for i, pair in enumerate(pairs, 1):
        assert prod(pair) == 210 * (n+i)
        for value in pair:
            fac.update({int(p): int(e) for p, e in factorint(value).items()})
    for p in (2, 3, 5, 7):
        fac[p] -= 4
    fac = Counter({p: e for p, e in fac.items() if e})
    assert all(e > 0 for e in fac.values())
    assert prod(p**e for p, e in fac.items()) == prod(n+i for i in range(1, 5))
    # Independent factorization of the four original consecutive integers.
    original = Counter()
    for i in range(1, 5):
        original.update({int(p): int(e) for p, e in factorint(n+i).items()})
    assert original == fac
    for p, e in sorted(fac.items()):
        available = factorial_valuation(2*n, p) - 2*factorial_valuation(n, p)
        assert available == doubling_carries(n, p)
        factorial_margin = factorial_valuation(2*n, p) - 2*factorial_valuation(n+4, p)
        assert factorial_margin == available - 2*e
        rows.append({"p": p, "product_valuation": e, "required": 2*e,
                     "available_carries": available,
                     "factorial_margin": factorial_margin})
    return {"k": 4, "t": t, "x": x, "m": m, "n": n,
            "quadratic_pairs": pairs,
            "passes": all(row["factorial_margin"] >= 0 for row in rows),
            "prime_checks": rows,
            "scope": "Exact finite certificate; not an infinitude proof."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t", type=int, default=8)
    parser.add_argument("--scan-stop", type=int)
    args = parser.parse_args()
    if args.scan_stop is None:
        result = certificate(args.t)
    else:
        if args.scan_stop < 0:
            parser.error("scan-stop must be nonnegative")
        passed = []
        for t in range(args.scan_stop):
            record = certificate(t)
            if record["passes"]:
                passed.append({key: record[key] for key in ("t", "x", "n")})
        result = {"interval": [0, args.scan_stop], "upper_endpoint_excluded": True,
                  "count": len(passed), "witnesses": passed,
                  "scope": "Finite exact scan only; no asymptotic density claim."}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
