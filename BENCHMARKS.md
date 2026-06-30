# BitMako Benchmark Results

All measurements on a single workstation, Windows 11, 32 GB RAM.
Rust release build: `lto = "fat"`, `codegen-units = 1`, `opt-level = 3`, `target-cpu = native`.
Corpus: Enamine REAL, **1,364,304,490 compounds**.

---

## Linear scan baseline

Sequential Tanimoto scan over the flat fingerprint store (`compounds.fp`, 163 GB),
no index, no pruning — equivalent to a brute-force RDKit similarity screen.

| Metric | Value |
|---|---|
| Throughput | 1,576,536 compounds / sec |
| Full-corpus time (extrapolated from 10M sample) | **865 s (14.4 min)** |

Rate is dominated by the POPCNT inner loop (`target-cpu=native`, AVX2 autovectorized).
Any similarity search that avoids a full scan is competing against this baseline.

---

## WAND similarity search

Block-Max WAND with streaming skip-indexed cursors.
Each query returns the exact top-10 most Tanimoto-similar compounds above threshold —
no approximation.

| Compound | Query size | Threshold | WAND time | Docs evaluated | Eval fraction | vs. linear |
|---|---|---|---|---|---|---|
| Ethanol | 5 bits | 0.15 | 2.7 s | 92 | 6.7 per billion | **320×** |
| Benzoic acid | 17 bits | 0.20 | 7.7 s | 27 | 2.0 per billion | **112×** |
| Paracetamol | 22 bits | 0.25 | 8.4 s | 50 | 3.7 per billion | **103×** |
| Aspirin | 24 bits | 0.25 | 17 s | 26 | 1.9 per billion | **51×** |
| Ibuprofen | 28 bits | 0.30 | 14 s | 29 | 2.1 per billion | **62×** |

**Docs evaluated** is the number of compounds for which exact Tanimoto was computed.
All other compounds were pruned by WAND's pivot rule, block-level upper bounds, or
the popcount upper bound — without ever fetching their fingerprint.

The dominant cost is posting-list traversal (advancing streaming cursors through
the 70 GB inverted index), not fingerprint evaluation. WAND spends most of its time
ruling out candidates, not computing similarity scores.

---

## Property-filtered search (prop-store vs. Lance post-filter)

When `--prop-store compounds.prop` is supplied, property filters (MW, LogP, etc.)
are applied inside the pivot loop as an O(1) mmap read before the fingerprint fetch.
Compounds failing a filter never reach exact Tanimoto evaluation, so the dynamic
threshold rises faster on genuinely qualifying results.

| Query | Filters | Mode | Time |
|---|---|---|---|
| Benzoic acid @ 0.20, top-10 | mw ≤ 350, logP ≤ 3 | `--prop-store` (in-loop) | **8 s** |
| Benzoic acid @ 0.20, top-10 | mw ≤ 350, logP ≤ 3 | `--lance` post-filter (20× over-fetch) | 30 s |

**3.7× faster** than the Lance post-filter path for equivalent results.

---

## Pruning rate summary

WAND evaluates between **26 and 92 documents** out of 1.36 billion to find top-10
results — a pruning rate of **>99.999999%**. Put differently, the engine computes
exact Tanimoto for roughly 1 compound per 10 million in the corpus.

This is the core property of Block-Max WAND: the min-shared-bits pivot rule and
per-block maximum-popcount bounds together eliminate virtually all of the corpus
before any fingerprint is read, and the dynamic threshold tightens pruning further
as the top-k heap fills.

---

## Methodology notes

- **Linear baseline** is measured on a 10M-compound sample and extrapolated to the
  full 1.36B corpus (rate is independent of corpus size for a sequential scan).
- **WAND times** are wall-clock, cold page cache, single query, single thread.
- Speedup = linear extrapolated time / WAND time.
- All WAND results are **exact** — the same top-k as an exhaustive scan.
- Queries with 0 results still exercise the full engine; the early-exit mechanism
  (50M stale pivot iterations) caps worst-case search time at ~35 s.
