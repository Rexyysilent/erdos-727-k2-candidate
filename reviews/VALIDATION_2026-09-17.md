# Validation receipt and limits

## Actual source state

The update was based on `Rexyysilent/erdos-727-k2-candidate` at
`3dec50a4d4d08d551aabb910368b849d2648a4d0`. The existing README, September 7
manuscript, original verifier and September 12 results were preserved unchanged.

The original `verify_727_k2.py` was reconstructed and checked byte-for-byte
against Git blob `44584200d7a2f7042b852f3ddcf74fe3bc7fa07e` before execution.
The supplied 16 September PDF was reviewed locally without edits. Its SHA-256 is
`52db8eecadf4a041c8c7e93c922b92593818b0dac89bb2293bbc4cffec3fbd61`.
The binary PDF is not stored in this text-first repository update.

## Executed checks

- **91 passing tests** on Python 3.13.5 / pytest 9.0.2.
- 38 tests cover the original k=2 family.
- 35 tests cover generalized k=2 family algebra, residue coefficients,
  constants, overlap algebra and finite-family calculations.
- 18 tests cover the separate k=3 construction.
- Literal factorial divisibility for every integer `0<=n<=500` agrees with
  independent carry and central-binomial formulations.
- The original verifier reproduces the supplied-PDF counts exactly:
  `120/1000`, `1879/10000`, `21716/100000`.
- A separate C++ implementation gives `234157/1000000` on
  `[10^6,2*10^6)` and `2483873/10000000` on `[10^7,2*10^7)`.
- 2,000 selected high-parameter original-family cases agree between independent
  C++ SPF/Legendre and Python trial-division/digit-carry paths.
- Three additional k=2 families were cross-checked on 3,000 selected high-range
  parameters in total.
- The bounded search through `1<=u<=150` examines 2,058 admissible root classes;
  the minimum observed `b` is `2/3` at `(14,195)` and `(14,209)`.
- The fixed positive-bound ensemble contains 296 families across 92 denominators.
  Its residue algebra was separately checked over 650,180 unit classes and
  2,600,720 direct fractional-part calculations.
- Logarithms, square roots and infinite-series constants were enclosed with exact
  rational/integer arithmetic and explicit tails. Displayed decimals are rounded
  summaries, not substitutes for the enclosures.

## Separate k=3 checks

For `210t²+391t+179`, exact counts are `126/1000`, `1880/10000`, and
`21456/100000` on the corresponding half-open intervals. A Python implementation
agrees with C++ on all 1,000 selected cases in `[199000,200000)`.

The first witness found in that family is `t=61,n=805440`; its prime-by-prime
valuation certificate verifies `((805443)!)² | 1610880!`. It is not asserted to
be the smallest k=3 solution outside this family.

## Analytic review boundary

The supplied argument's cancellation, four-factor algebra, fixed-prime reduction,
root lifting, small-prime summation, exponential sums, Fourier-to-box passage,
prime-band coverage and sequential limits were checked. No fatal gap was identified
under its stated analytic inputs. This is an AI-assisted mathematical review, not
independent peer review or a formal proof certificate.

The generalized-family and k=3 notes are new proposed arguments. Their most
important independent-review targets are stated in those files.

## Not established

This repository update does not establish:

- an exact success density for the original family;
- an effective convergence threshold for its asymptotic lower bound;
- a global optimum over all quadratic families;
- a result for every fixed `k`;
- priority or the current resolution status of Erdős problem 727.

Finite computations verify arithmetic statements and regression behavior. They do
not replace the analytic arguments.
