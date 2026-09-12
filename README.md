# Erdős problem 727, k = 2: candidate argument seeking review

This repository shares a proposed quadratic-family argument for the **k = 2 case** of [Erdős problem 727](https://www.erdosproblems.com/727), together with an exact Python verifier and reproducible numerical results.

**Status: unreviewed candidate.** The argument was developed with ChatGPT. It has not been independently reviewed or formally verified. Publication here is a request for mathematical assessment, not a claim of an accepted solution. The general question for every fixed k ≥ 2 is not resolved by this note.

## Proposed result

For the family

$$n=F(t)=6t^2+11t+3,$$

the note proposes that a positive asymptotic proportion of integer parameters satisfy

$$((n+2)!)^2\mid(2n)!.$$

More specifically, it proposes

$$\liminf_{T\to\infty}\frac{1}{T}\#\{t\in[T,2T]\cap\mathbb Z:((F(t)+2)!)^2\mid(2F(t))!\}\geq\rho,$$

where the claimed lower bound is approximately 0.25922348. This is an asymptotic claim, not a bound on the success rate in every finite sample.

## Read and reproduce

- [Complete candidate argument and references](erdos_727_k2_candidate.md)
- [Exact numerical verifier](verify_727_k2.py), using only the Python standard library
- [Saved results](verification_results.json), rerun on 12 September 2026
- [Earlier public Gist](https://gist.github.com/Rexyysilent/55c6f58e3a8b349ef720604ce5042ff2)

With Python 3.9 or later:

```sh
python verify_727_k2.py --n 208
python verify_727_k2.py --start 0 --stop 1000
```

The scan finds **120 valid parameters out of 1,000** for 0 ≤ t < 1000. At t = 5, n = 208 satisfies (210!)² | 416!. The verifier uses prime valuations and also cross-checks against a literal binomial coefficient for modest inputs.

These computations verify individual examples. They do not verify the infinitude argument, the proposed density bound, or the analytic estimates.

## Request for mathematical review

The most useful feedback would be a precise assessment of these steps in the note:

1. **Section 7, Lemma 3, equations (19)–(24):** the uniform exponential-sum estimates and their passage to equidistribution over the stated prime ranges.
2. **Section 6, equations (14)–(15):** the simple-root digit count and summation over small growing primes.
3. **Sections 8–10:** the prime-range coverage, accumulation of uniform errors, and order of limits.

Please [open an issue](https://github.com/Rexyysilent/erdos-727-k2-candidate/issues/new) with the relevant lemma or equation, the objection or justification, and any counterexample or missing hypothesis. Partial reviews are welcome; identify which steps you have checked. Numerical evidence alone should not be treated as verification of the proof.

The supplied note and verifier are preserved unchanged in this repository. The note is dated 7 September 2026; this repository and the saved verification results were prepared on 12 September 2026.
