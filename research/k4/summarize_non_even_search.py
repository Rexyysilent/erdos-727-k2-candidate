#!/usr/bin/env python3
"""Summarize exact rational quartic hits from non_even_quartic_search.cpp.

The search output supplies ordered u values and a rational linear coefficient L.
This script recovers all coefficients with fractions.Fraction, verifies the
four factorization constants form an arithmetic progression, removes common
scalings and duplicate quartics, and counts rationally split quadratic factors.
It does not assert integer admissibility or factorial divisibility.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import comb, gcd, isqrt
from pathlib import Path
import json


def second_difference(values: list[Fraction], index: int) -> Fraction:
    return values[index] - 2 * values[index + 1] + values[index + 2]


def rational_square(value: Fraction) -> bool:
    return (
        value > 0
        and isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


def summarize(path: Path) -> dict:
    counts: Counter[tuple[int, int]] = Counter()
    seen: set[tuple[Fraction, ...]] = set()
    candidates = []
    raw_hits = 0
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        raw_hits += 1
        shifts = tuple(map(int, fields[:4]))
        if gcd(*shifts) > 1:
            continue
        linear = Fraction(fields[4].removeprefix("C="))
        squares = [Fraction(u * u) for u in shifts]
        index = 0 if second_difference(squares, 0) else 1
        quadratic = (
            linear**2 * second_difference([1 / a for a in squares], index)
            - second_difference([a * a for a in squares], index)
        ) / (2 * second_difference(squares, index))
        constants = [
            ((a + quadratic) ** 2 - linear**2 / a) / 4 for a in squares
        ]
        difference = constants[1] - constants[0]
        assert difference != 0
        assert all(
            constants[i] == constants[0] + i * difference for i in range(4)
        )
        if difference < 0:
            shifts = shifts[::-1]
            constants = constants[::-1]
            difference = -difference
        constant = constants[0] - difference
        key = quadratic, linear, difference, constant
        if key in seen:
            continue
        seen.add(key)
        discriminants = [
            -u * u - 2 * quadratic + sign * 2 * linear / u
            for u in shifts
            for sign in (-1, 1)
        ]
        split = sum(rational_square(value) for value in discriminants)
        zero = sum(value == 0 for value in discriminants)
        counts[split, zero] += 1
        if split >= 3:
            candidates.append(
                {
                    "u": shifts,
                    "B": str(quadratic),
                    "L": str(linear),
                    "D": str(difference),
                    "E": str(constant),
                    "split_quadratics": split,
                    "zero_discriminants": zero,
                    "discriminants": [str(value) for value in discriminants],
                }
            )
    candidates.sort(
        key=lambda row: (-row["split_quadratics"], max(row["u"]), sum(row["u"]))
    )
    return {
        "scope": "Bounded exact rational-factorization search, not an infinitude proof.",
        "raw_hits": raw_hits,
        "primitive_distinct_quartics": len(seen),
        "split_count_distribution": [
            {"split_quadratics": split, "zero_discriminants": zero, "count": count}
            for (split, zero), count in sorted(counts.items())
        ],
        "candidates_with_at_least_three_split_quadratics": candidates,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bound", type=int, help="Generator's exhaustive u bound")
    args = parser.parse_args()
    result = summarize(args.input)
    if args.bound is not None:
        if not 4 <= args.bound <= 300:
            parser.error("--bound must be between 4 and 300")
        result["exhaustive_u_bound"] = args.bound
        result["orderings_tested"] = 12 * comb(args.bound, 4)
    report = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(report)
    else:
        print(report, end="")
