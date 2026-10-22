# Receipts for the construction searches

These notes distinguish captured output from statistics retained in the
research session. All negative results below are bounded searches.

## Completely split cubics, A≤100000000

Source: `construction_cubic_split_search.cpp`.
The numeric search argument was `100000000`; hits were redirected to
`cubic-split-search-100000000.txt`. The program writes progress and summary
statistics to stderr, so the empty result file records only the absence of
hits. The following stderr statistics were observed during the completed run
and retained in the session; this receipt was reconstructed afterward:

```text
pairs 30237841; max A 100000000
multiple-representation groups 8738132; AP hits 0
```

The search enumerates every pair r≥s≥0 with r≥1 and r²+rs+s²≤A_max,
groups them by A=r²+rs+s², inserts both signs of B=rs(r+s), and tests every
four-term arithmetic progression among the distinct B values in each group.
The displayed statistics were not inferred from the empty output file.

Reproduction command, with both streams preserved on a future run:

```sh
g++ -O3 -std=c++17 research/k4/construction_cubic_split_search.cpp -o /tmp/construction_cubic_split_search
/tmp/construction_cubic_split_search 100000000 > research/k4/cubic-split-search-100000000.txt 2> research/k4/cubic-split-search-100000000.stderr.txt
```

The earlier `cubic-split-search-20000000.txt` is also empty. Detailed stderr
counts for that earlier run were not retained here, and are not invented.

## Three bad factors simultaneously split, P≤10000000

Source: `construction_three_bad_split_search.cpp`.
The completed run used P_max=10000000 and reported 196650 admissible
candidates, zero passes at the second square condition, and zero all-three
hits. The file `three-bad-split-search-10000000.txt` contains hits only and
is empty. These statistics are retained session observations, recorded here
afterward, rather than a captured copy of stderr.

The condition that the third required square be nonnegative implies P/Q>91.
Thus Q≤109890 covers the stated P bound. The first square condition is
enumerated through factor pairs of840Q², as described in the source and
`construction-exploration.md`.

## Linear residue refinement optimization

Source: `construction_optimize_linear.cpp`.
This run used an exact enumeration, not random samples: 4416 admissible
refinements and218592 unit-square classes modulo63987840.
The captured JSON is `linear-refinement-optimization.json.gz`. The completed
run printed

```text
Exact candidates 4416; square classes 218592; best Y 23971; best linear cost 0.11821
```

The JSON retains the full precision of the cost calculation and the exact
integer class counts for every refinement. There are2208 instances of each
of the two costs0.11820994127346769 and0.19015166959828658.

## Expanded signed all-six norm discovery

These completed runs preserve both stdout and stderr in the repository:

```sh
OMP_NUM_THREADS=8 /tmp/construction_downward_norm_search 10000 2521 allsix > research/k4/downward-all-six-search-10000.txt 2> research/k4/downward-all-six-search-10000.stderr.txt
OMP_NUM_THREADS=8 /tmp/construction_downward_norm_search 10000 -839 allsix > research/k4/upward-all-six-search-10000.txt 2> research/k4/upward-all-six-search-10000.stderr.txt
```

The captured summaries report six signed downward hits and four signed upward
hits. Their root-column order is c=0,−14,21,1,15,−20. Each negative entry−1
marks a nonsquare. Every retained candidate has exactly three squares; none
has four. All hits have Q=1. Details and the resulting families are in
`construction-additional-families.md`.

## Downward432 refinement optimization

The completed exact run of `construction_432_optimize_linear.cpp` checked1728
refinements, with279936 unit-square classes per refinement, and reported

```text
Exact candidates 1728; square classes 279936; best Y 16273; best linear cost 0.306479
```

Every refinement has the same cost; the complete integer counts and full
numeric precision are in `432-linear-refinement-optimization.json.gz`.
