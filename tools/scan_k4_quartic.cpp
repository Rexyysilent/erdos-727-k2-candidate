// Exact finite scanner for the k=4 quartic family. No asymptotic claim.
// g++ -std=c++17 -O3 tools/scan_k4_quartic.cpp -o /tmp/scan_k4_quartic
// /tmp/scan_k4_quartic START STOP
// GCC/Clang unsigned 128-bit arithmetic is required. No external libraries.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using U64 = std::uint64_t;
using U128 = __uint128_t;

std::string decimal(U128 n) {
    if (!n) return "0";
    std::string s;
    while (n) {
        s.push_back(static_cast<char>('0' + n % 10));
        n /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

U64 number(const char* arg) {
    std::string s(arg);
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument("arguments must be nonnegative integers");
    return std::stoull(s);
}

U128 factorial_valuation(U128 n, U64 p) {
    U128 v = 0;
    while (n) {
        n /= p;
        v += n;
    }
    return v;
}

unsigned doubling_carries(U128 n, U64 p) {
    unsigned total = 0;
    U64 incoming = 0;
    while (n || incoming) {
        U128 digit = n % p;
        incoming = (2 * digit + incoming >= p);
        total += incoming;
        n /= p;
    }
    return total;
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::invalid_argument("use START STOP");
        const U64 start = number(argv[1]), stop = number(argv[2]);
        if (start >= stop || stop > 100000)
            throw std::invalid_argument("require 0 <= START < STOP <= 100000");

        const U64 max_x = 210 * (stop - 1) + 14;
        const U64 max_factor = std::max(max_x * max_x + 230,
                                       max_x * max_x + max_x + 210);
        U64 limit = static_cast<U64>(std::sqrt(static_cast<long double>(max_factor)));
        while (limit + 1 <= max_factor / (limit + 1)) ++limit;
        while (limit > max_factor / limit) --limit;
        std::vector<bool> composite(limit + 1);
        std::vector<U64> primes;
        for (U64 p = 2; p <= limit; ++p) {
            if (composite[p]) continue;
            primes.push_back(p);
            if (p <= limit / p)
                for (U64 j = p * p; j <= limit; j += p) composite[j] = true;
        }
        auto add_factors = [&](U64 value, std::map<U64, int>& factors) {
            for (U64 p : primes) {
                if (p > value / p) break;
                while (value % p == 0) {
                    ++factors[p];
                    value /= p;
                }
            }
            if (value > 1) ++factors[value];
        };

        U64 good = 0, small_failure = 0, large_failure = 0, both_failure = 0;
        std::vector<std::pair<U64, U128>> first;
        for (U64 t = start; t < stop; ++t) {
            const U64 x = 210 * t + 14, x2 = x * x, m = x2 + 209;
            const U128 n = U128(m) * (m + 1) / 210 - 3;
            if (m % 210 != 195 || U128(m) * (m + 1) % 210)
                throw std::logic_error("integrality failure");
            const std::array<std::pair<U64, U64>, 4> pairs = {{
                {x2 + 189, x2 + 230}, {x2 + 195, x2 + 224},
                {x2 + 209, x2 + 210}, {x2 - x + 210, x2 + x + 210}
            }};
            std::map<U64, int> factors;
            for (unsigned i = 0; i < 4; ++i) {
                const auto [a, b] = pairs[i];
                if (U128(a) * b != 210 * (n + i + 1))
                    throw std::logic_error("factorization mismatch");
                add_factors(a, factors);
                add_factors(b, factors);
            }
            for (U64 p : {2, 3, 5, 7}) factors[p] -= 4;

            bool bad_small = false, bad_large = false;
            for (const auto [p, e] : factors) {
                if (e < 0) throw std::logic_error("negative exponent");
                if (!e) continue;
                const U128 available = factorial_valuation(2 * n, p)
                                       - 2 * factorial_valuation(n, p);
                if (available != doubling_carries(n, p))
                    throw std::logic_error("carry mismatch");
                const bool pass = available >= U128(2 * e);
                if (pass != (factorial_valuation(2 * n, p)
                             >= 2 * factorial_valuation(n + 4, p)))
                    throw std::logic_error("factorial mismatch");
                if (!pass) {
                    if (p <= x) bad_small = true;
                    else bad_large = true;
                }
            }
            if (!bad_small && !bad_large) {
                ++good;
                if (first.size() < 10) first.push_back({t, n});
            }
            small_failure += bad_small;
            large_failure += bad_large;
            both_failure += bad_small && bad_large;
        }
        std::cout << "{\"family\":\"x=210t+14; n=(x^4+419x^2+43260)/210\","
                  << "\"k\":4,\"start_inclusive\":" << start
                  << ",\"stop_exclusive\":" << stop << ",\"checked\":" << stop - start
                  << ",\"valid\":" << good
                  << ",\"has_failing_prime_at_most_x\":" << small_failure
                  << ",\"has_failing_prime_above_x\":" << large_failure
                  << ",\"has_both\":" << both_failure << ",\"first_witnesses\":[";
        for (unsigned i = 0; i < first.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << "{\"t\":" << first[i].first
                      << ",\"n\":" << decimal(first[i].second) << '}';
        }
        std::cout << "],\"scope\":\"Finite exact arithmetic only; "
                     "not an infinitude proof.\"}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
