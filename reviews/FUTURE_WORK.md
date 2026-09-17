# Reviewable next branches

| Branch | Question | Concrete acceptance condition |
|---|---|---|
| `independent-analytic-review` | Do the generalized root, Fourier and residue-weighted prime sums justify their stated uniformity? | A human review identifies exact equations and dependencies; objections become revisions, not suppressed caveats |
| `residue-coefficient-optimization` | Can `b(Q,r)` be classified or bounded beyond a finite search? | A proof, or an explicitly bounded new enumeration with independent checks; no global optimum inferred from a search |
| `obstruction-overlap` | Can the original family's positive excess obstruction count be bounded asymptotically? | A proved positive lower bound for `E[(W-1)+]`, with all joint-prime errors controlled; finite overlaps do not suffice |
| `effective-discrepancy` | Can the nested limits become an explicit finite threshold? | Explicit constants, finite Fourier approximation, all error terms and a checked parameter choice |
| `arithmetic-formalization` | Can the exact identities be machine-checked separately from the analytic proof? | A formal development of cancellation, valuations, factorization and boundary signs; no implication of full formalization until analytic lemmas are handled |

The cross-denominator overlap lemma in the generalized-family note permits a
fixed finite family collection to be combined without an independence claim.
A growing-modulus or growing-family argument still requires new uniform estimates.

For repository usefulness, keep a stable review index, exact test commands, small
regression witnesses, machine-readable certificates where useful, and a record of
which lemmas external reviewers actually checked. Numerical scans should remain
clearly separated from proof verification.

## k=3 before k>=4

`K3_EXTENSION_2026-09-17.md` gives a proposed six-linear-factor k=3 argument.
Its signed residue table, simple-root hypotheses and fixed-modulus prime sum should
receive independent scrutiny before being treated as established. Extending to
k>=4 requires new control of further consecutive factors; the k=2 and k=3
constructions do not provide an automatic induction.
