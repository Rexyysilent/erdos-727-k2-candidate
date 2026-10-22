"""Independent finite checks for the k=4 research continuation.

These tests certify identities, saved witnesses, and the finite scanner.
They do not certify the analytic estimates or infinitude.
"""
from collections import Counter
import json
from math import isqrt, prod
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "research" / "k4"


def load(name):
    return json.loads((DATA / name).read_text())


def valuation_factorial(n, p):
    # Independent power-sum form of Legendre's formula.
    power, total = p, 0
    while power <= n:
        total += n // power
        power *= p
    return total


def carries_from_fractional_parts(n, p):
    power, total = p, 0
    while power <= 2 * n:
        total += 2 * (n % power) >= power
        power *= p
    return total


def prime_by_trial(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def check_valuation_rows(record):
    n = record["n"]
    rows = record["prime_checks"]
    assert prod(row["p"] ** row["product_valuation"] for row in rows) == prod(
        n + i for i in range(1, 5)
    )
    for row in rows:
        p, e = row["p"], row["product_valuation"]
        assert e > 0
        available = carries_from_fractional_parts(n, p)
        margin = valuation_factorial(2 * n, p) - 2 * valuation_factorial(n + 4, p)
        assert available == row["available_carries"]
        assert row["required"] == 2 * e
        assert margin == available - 2 * e == row["factorial_margin"]
    assert record["passes"] == all(row["factorial_margin"] >= 0 for row in rows)


def test_simple_quartic_witness_independently():
    record = load("quartic-witness-t8.json")
    assert (record["t"], record["x"], record["n"]) == (8, 1694, 39219110444)
    for i, pair in enumerate(record["quadratic_pairs"], 1):
        assert prod(pair) == 210 * (record["n"] + i)
    assert all(prime_by_trial(row["p"]) for row in record["prime_checks"])
    check_valuation_rows(record)
    assert record["passes"]


def test_552_witness_and_recursive_primality_proofs():
    record = load("552-witness-refined-t2.json")
    assert record["n"] == 4195055316568514657639543
    assert record["step"] == 63987840 and record["t"] == 2
    certified = set()
    for proof in record["lucas_prime_certificates"]:
        p = proof["p"]
        if p == 2:
            assert proof["base_case"]
        else:
            factors = proof["p_minus_one_factors"]
            assert all(q in certified and e >= 1 for q, e in factors)
            assert prod(q ** e for q, e in factors) == p - 1
            a = proof["lucas_witness"]
            assert pow(a, p - 1, p) == 1
            assert all(pow(a, (p - 1) // q, p) != 1 for q, _ in factors)
        certified.add(p)
    assert all(row["p"] in certified for row in record["prime_checks"])
    shifts = Counter()
    for row in record["factorizations"]:
        assert prod(p ** e for p, e in row["prime_factors"]) == row["value"]
        assert all(p in certified or p in (2, 3, 5, 7, 23)
                   for p, _ in row["prime_factors"])
        shifts[row["shift"]] = shifts.get(row["shift"], 1) * row["value"]
    for i in range(1, 5):
        assert shifts[i] == 63987840 * (record["n"] + i)
    check_valuation_rows(record)
    assert record["passes"]


def test_quartic_factorization_identities():
    # Degree-four identities are determined by their values at five points.
    for x in range(-3, 4):
        simple = x ** 4 + 419 * x * x + 43260
        pairs = [(x*x+189, x*x+230), (x*x+195, x*x+224),
                 (x*x+209, x*x+210), (x*x-x+210, x*x+x+210)]
        assert all(a * b == simple + 210 * i
                   for i, (a, b) in enumerate(pairs, 1))
        d = 63987840
        improved = (x*x-1681) * (x*x-1129) - 3 * d
        factorizations = [
            (x*x-12721) * (x*x+9911),
            (x-97) * (x+97) * (x*x+6599),
            (x-41) * (x+41) * (x*x-1129),
            (x*x-138*x+8117) * (x*x+138*x+8117),
        ]
        assert all(value == improved + d * i
                   for i, value in enumerate(factorizations, 1))
        eight = x**4 - 810*x*x + 2240*x - 41751
        factorizations = [
            (x-1)*(x-27)*(x*x+28*x-53),
            (x-9)*(x+29)*(x*x-20*x-149),
            (x-13)*(x+9)*(x*x+4*x-677),
            (x-19)*(x+27)*(x*x-8*x-233),
        ]
        assert all(value == eight + 40320 * i
                   for i, value in enumerate(factorizations, 1))


def test_deterministic_fourth_shift_large_prime_lemma():
    tested = 0
    for t in range(40):
        x = 210 * t + 14
        n = (x**4 + 419*x*x + 43260) // 210
        bound = (x**3 - x*x + 210*x - 420) // 7
        for original in (x*x-x+210, x*x+x+210):
            assert original < bound
            value, p, primes = original, 2, set()
            while p*p <= value:
                if value % p == 0:
                    primes.add(p)
                    while value % p == 0:
                        value //= p
                p += 1
            if value > 1:
                primes.add(value)
            for p in primes:
                if p*p > bound:
                    tested += 1
                    assert (n + 4) % p == 0 and (n + 4) % (p*p) != 0
                    assert 2 * (n % (p*p)) >= p*p
                    assert carries_from_fractional_parts(n, p) >= 2
    assert tested > 0


@pytest.fixture(scope="module")
def scanner(tmp_path_factory):
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("GCC with unsigned 128-bit support is not installed")
    binary = tmp_path_factory.mktemp("k4_scan") / "scan"
    subprocess.run([compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
                    str(ROOT / "tools" / "scan_k4_quartic.cpp"), "-o", str(binary)],
                   check=True, capture_output=True, text=True)
    return binary


def test_cpp_scanner_matches_independent_python_scan(scanner):
    actual = json.loads(subprocess.check_output([str(scanner), "0", "1000"], text=True))
    saved_cpp = load("quartic-cpp-scan-0-1000.json")
    saved_python = load("quartic-scan-1000.json")
    assert actual == saved_cpp
    assert actual["valid"] == saved_python["count"] == 157
    assert [r["t"] for r in actual["first_witnesses"]] == [
        r["t"] for r in saved_python["witnesses"][:10]
    ]


def test_cpp_scanner_smallest_interval_has_complete_trial_basis(scanner):
    # At t=0, x=14 and x^2+195=17*23. A sieve stopping at x+2 is too short.
    actual = json.loads(subprocess.check_output([str(scanner), "0", "1"], text=True))
    assert actual["checked"] == 1 and actual["valid"] == 0
    assert actual["has_failing_prime_above_x"] == 1
