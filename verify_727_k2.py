#!/usr/bin/env python3
"""Exact numerical checks for a candidate k=2 argument for Erdős problem 727.

This is a numerical verifier, NOT a formal verification of the asymptotic proof.
It uses only the Python standard library and never constructs factorials.

Examples:
    python verify_727_k2.py --n 208
    python verify_727_k2.py --start 0 --stop 1000
    python verify_727_k2.py --start 10000 --stop 20000

The parameterized family is n = 6*t*t + 11*t + 3, t >= 0.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
import json
from math import comb, isqrt
from typing import Optional


def factorial_valuation(n: int, p: int) -> int:
    """Return v_p(n!) for a prime p, by Legendre's formula."""
    if n < 0 or p < 2:
        raise ValueError("Require n >= 0 and p >= 2.")
    total = 0
    while n:
        n //= p
        total += n
    return total


def make_spf(limit: int) -> array:
    """A zero entry for an integer >1 means that integer is prime."""
    if limit < 1:
        raise ValueError("The sieve limit must be positive.")
    spf = array('I', [0]) * (limit + 1)
    for p in range(2, isqrt(limit) + 1):
        if not spf[p]:
            for multiple in range(p * p, limit + 1, p):
                if not spf[multiple]:
                    spf[multiple] = p
    return spf


def factor(n: int, spf: Optional[array] = None) -> Counter[int]:
    if n < 1:
        raise ValueError("Factorization requires a positive integer.")
    result: Counter[int] = Counter()
    if spf is not None:
        if n >= len(spf):
            raise ValueError("Sieve too short for this input.")
        while n > 1:
            p = int(spf[n]) or n
            while n % p == 0:
                result[p] += 1
                n //= p
        return result
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] += 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        result[n] += 1
    return result


def certificate(n: int, factors: Optional[Counter[int]] = None) -> dict:
    """Check ((n+2)!)^2 | (2n)! by checking just the denominator primes.

    Equivalently, ((n+1)*(n+2))^2 divides binomial(2*n, n).
    All other prime valuations are nonnegative because that binomial is integer.
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")
    if factors is None:
        factors = factor(n + 1) + factor(n + 2)
    product = 1
    for p, exponent in factors.items():
        product *= p ** exponent
    if product != (n + 1) * (n + 2):
        raise ValueError("The supplied factorization is not the required product.")
    rows = []
    for p, exponent in sorted(factors.items()):
        available = factorial_valuation(2 * n, p) - 2 * factorial_valuation(n, p)
        rows.append({"prime": p, "required": 2 * exponent,
                     "available": available, "passes": available >= 2 * exponent})
    passes = all(row["passes"] for row in rows)
    # An independent literal-binomial cross-check for modest inputs.
    if n <= 10000:
        literal = comb(2 * n, n) % (((n + 1) * (n + 2)) ** 2) == 0
        if literal != passes:
            raise AssertionError("Valuation and literal-binomial checks disagree.")
    return {"n": n, "passes": passes, "prime_checks": rows}


def scan(start: int, stop: int) -> dict:
    if not 0 <= start < stop:
        raise ValueError("Require 0 <= start < stop (stop is exclusive).")
    spf = make_spf(max(5, 6 * (stop - 1) + 5))
    count = 0
    witnesses = []
    for t in range(start, stop):
        n = 6 * t * t + 11 * t + 3
        linear_values = (2 * t + 1, 3 * t + 4, t + 1, 6 * t + 5)
        factors: Counter[int] = Counter()
        for value in linear_values:
            factors.update(factor(value, spf))
        result = certificate(n, factors)
        if result["passes"]:
            count += 1
            if len(witnesses) < 20:
                witnesses.append({"t": t, "m": 6 * t + 5, "n": n})
    return {"family": "n = 6*t^2 + 11*t + 3", "start_inclusive": start,
            "stop_exclusive": stop, "parameters_checked": stop - start,
            "valid_parameters": count, "observed_fraction": count / (stop - start),
            "first_witnesses": witnesses,
            "caution": "Finite computation does not prove infinitude or the density bound."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, help="Check a particular n instead of scanning the family.")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=1000)
    args = parser.parse_args()
    try:
        result = certificate(args.n) if args.n is not None else scan(args.start, args.stop)
    except (ValueError, MemoryError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
