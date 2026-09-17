# Checked continuation: Erdős 727, k=2 and a separate k=3 extension

17 September 2026. **Proposed arguments seeking independent mathematical review.**
This index is additive. The existing README, September 7 manuscript, original
verifier and September 12 results are preserved unchanged.

**Publication note.** The dated review PDFs were rendered before authenticated GitHub
write actions became available in this session. Any embedded “not pushed” wording
records that preparation-time environment; this index and Git history are the current
publication record. `README.md` and the pre-existing proof/verifier/results remain
unchanged.

## Read in this order

1. `reviews/REVIEW_AND_EXTENSIONS_2026-09-17.md`: review of the supplied
   16 September expanded proof, an original-family upper density bound, an
   obstruction-count identity, and extended finite checks.
2. `reviews/GENERALIZED_FAMILIES_2026-09-17.md`: continuation for other
   quadratic families, including a proposed parameter lower bound
   `0.44232553173768…` for `210t²+419t+207` and fixed finite-family counting
   bounds.
3. `reviews/K3_EXTENSION_2026-09-17.md`: a separate proposed `k=3` proof for
   `210t²+391t+179`, with lower parameter proportion
   `0.2550393216621966…`.
4. `reviews/VALIDATION_2026-09-17.md`: exact test and evidence boundaries.
5. `reviews/FUTURE_WORK.md`: concrete next proof and verification branches.

The supplied 12-page expanded PDF and locally rendered review PDFs were checked
during preparation but are not stored in this text-first repository commit. The
supplied-PDF SHA-256 and the validation boundary are recorded in
`reviews/VALIDATION_2026-09-17.md`.

## Principal distinctions

The original family's lower constant remains `rho=0.2592234836263289…`.
The review derives the complementary upper limiting proportion
`1-(log 3)/2=0.450693855665945…` for that same family. It does not establish
existence of its success density.

The generalized-family theorem uses `Q=u(u+1)` and `m=Qt+r`, with
`Q | r(r+1)`. An exact residue average `b(Q,r)` gives a candidate lower bound
`1-S-b(Q,r)*(log 3)/2`, where
`S=4*sum(J>=4, 2^(1-J)*log((J+1)/J))`.
For `Q=210,r=195 or 209`, `b=2/3`. This improves the parameter bound for
**different families**, not the original polynomial.

A proved-in-the-addendum intersection lemma shows that distinct fixed
denominators have only `O(log N)` common values below height `N`. It permits
any fixed finite collection of these family lower bounds to be added. Ten
explicitly listed families yield a candidate coefficient `0.43890932340518…`
for counting all successful `n<=N`. A bounded enumeration of 296 families
through `u=150` yields `0.83733579699699…`. These are coefficients of `sqrt(N)`,
not natural densities among integers. The family set never grows with `N` in
the stated proof. Neither coefficient is claimed globally optimal.

## Separate k=3 continuation

The identities `210=14*15` and `420=20*21` allow all of `n+1,n+2,n+3`
to factor into six primitive linear terms when `n=210t²+391t+179`. The
six-class residue average is `5/6`. The proposed lower proportion is
`rho3=1-(3/2)*S-(5/12)*log(3)=0.2550393216621966…`. The proof writes
out the six-root lifting, changed prime costs, residue signs and limit order.
This is a new proposed k=3 result, not a claim present in the old PDF. No
assertion for k>=4 or every fixed k is made.

## Reproduce

The companion tools require Python 3.11+; the Python arithmetic itself uses only
the standard library. Install pytest for tests. The original verifier's previous
requirements are unchanged.

```sh
python -m pytest tests -q
python tools/review_checks.py --constants
python tools/family_search.py --u 14 --r 209
python tools/verify_k3.py --constants
python tools/verify_k3.py --t 61
```

Exact finite scanning with GCC/Clang (unsigned 128-bit extension required):

```sh
g++ -std=c++17 -O3 tools/scan_family.cpp -o scan_family
./scan_family 10000000 20000000
g++ -std=c++17 -O3 tools/scan_general_family.cpp -o scan_general_family
./scan_general_family 14 209 100000 200000
g++ -std=c++17 -O3 tools/scan_k3.cpp -o scan_k3
./scan_k3 100000 200000
```

The large original-family scan uses about 480 MB; the generalized scanner caps
its sieve at 50 million entries, about 200 MB. Counts are exact integer
arithmetic. Displayed fractions are rounded decimal summaries.

The checked tools and regression tests are included below `tools/` and `tests/`.
The larger machine-readable evidence bundle is retained separately from this compact
text-first repository update. The 91 tests and finite scans do not formally verify
the analytic arguments. No claim of priority, independent peer approval, a current
problem-status change, or a result for all fixed `k>=2` is made.
