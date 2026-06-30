# BitMako

Ultra-high-throughput chemoinformatics ETL and **Block-Max WAND** similarity search
engine, in pure Rust. Built to ingest and search the full
[Enamine REAL](https://enamine.net/compound-collections/real-compounds/real-database-subsets)
library — **1.36 billion compounds** — on a single workstation.

## What it does

1. **Ingests** a bzip2-compressed CXSMILES dump, computing an ECFP4-style 1024-bit
   Morgan fingerprint and molecular properties (MW, LogP, rotatable bonds, heavy
   atoms, ring count) for every compound, and writes them to a columnar
   [Lance](https://lancedb.github.io/lance/) dataset.
2. **Builds an inverted index** (1024 posting lists, one per fingerprint bit) with
   per-block max-popcount metadata for Block-Max WAND pruning.
3. **Builds a flat fingerprint store** — a memory-mappable file of every
   fingerprint, addressed directly by `doc_id`.
4. **Builds a flat property store** — a memory-mappable file of MW, LogP, rotatable
   bonds, heavy atoms, and ring count per compound, also addressed by `doc_id`.
5. **Searches** for the top-k most Tanimoto-similar compounds to a query SMILES,
   using BMW dynamic pruning to avoid scanning the full corpus. Property filters
   (`--mw-max`, `--logp-max`) are evaluated inside the pivot loop via O(1) mmap
   reads — before the fingerprint fetch — so the engine returns exactly top-k
   filtered results with no over-fetch.

## Scale (this build)

| Artifact | Value |
|---|---|
| Compounds ingested | 1,364,304,490 (0 parse failures) |
| Lance dataset | ~333 GB |
| Inverted index | ~70 GB |
| Flat fingerprint store | ~163 GB (128 bytes/compound) |
| Flat property store | ~27 GB (20 bytes/compound) |
| Ingest throughput | ~270k compounds/sec (single bzip2 stream) |

## Architecture

```
 .bz2 CXSMILES
      │  stream-decompress (crossbeam backpressure)
      ▼
 raw line chunks ──rayon──▶ SMILES parse + Morgan FP + properties
      │
      ▼
 Arrow RecordBatches ──▶ Lance dataset  (compounds.lance)
      │
      ├── build-index ──────▶ inverted index    (compounds.bitmako)
      │        └── build-skip ──▶ skip index    (compounds.skip)
      ├── build-fp-store ───▶ flat fingerprints (compounds.fp)
      └── build-prop-store ─▶ flat properties   (compounds.prop)
                                   │                    │
 query SMILES ──▶ BMW engine ◀─────┘                    │
                     │  1. WAND pivot: min-shared-bits   │
                     │  2. Property pre-screen ◀─────────┘ O(1) mmap
                     │  3. Fingerprint fetch + exact Tanimoto
                     ▼
              top-k (doc_id, Tanimoto) — already filtered
```

### Fingerprints

ECFP4 (radius 2), 1024 bits stored as `[u64; 16]`. A pure-Rust SMILES parser builds
the molecular graph; atom environments are hashed with CRC32 over two iterations.

### Tanimoto similarity

`|A ∩ B| / |A ∪ B|` = `popcount(A & B) / popcount(A | B)`, computed over the 16
words with POPCNT. Compiled with `target-cpu=native` so LLVM autovectorizes the
inner loop.

### WAND search with streaming posting lists

Posting lists are delta + LEB128-varint encoded in 128-doc blocks. Search runs
document-at-a-time over the active-bit lists with two key properties:

1. **Min-shared-bits pivoting.** A doc sharing `c` of the query's bits has
   Tanimoto `c/(P+K−c) ≤ c/P`, so reaching threshold `t` requires `c ≥ ⌈t·P⌉`
   shared bits. The engine jumps directly to docs appearing in `≥ θ` posting
   lists (the WAND pivot), and reads the exact intersection size `c` for free
   from how many cursors align at the pivot — a tight upper bound before any
   fingerprint fetch. `θ` rises as the top-k heap fills, tightening pruning.
2. **Streaming, skip-indexed cursors.** Each cursor decodes one 128-doc block at
   a time straight from the mmap'd index. `advance_to` consults the **skip index**
   (`compounds.skip`) to binary-search the block containing the target and
   decodes only that block — so jumping over a list of hundreds of millions of
   postings costs O(log blocks + one block), and resident memory per query is a
   few blocks rather than the whole lists. This is what keeps common-fragment
   queries from OOM-ing (the naive "decode every active list" approach needs
   tens of GB and can exceed RAM).

For each surviving pivot, the engine optionally checks molecular properties from the
mmap'd flat property store (O(1), 20 bytes) before fetching the fingerprint (128
bytes). Compounds failing a property filter never reach the exact-Tanimoto step and
never enter the heap — so `top_k` results come back already filtered and no
over-fetch heuristic is needed.

Surviving pivots then fetch the fingerprint from the mmap'd flat fingerprint store,
compute exact Tanimoto, and update the top-k min-heap.

The streaming cursor is verified against an exhaustive linear scan in the test
suite (identical results across many queries × thresholds).

## Usage

Requires Rust (with the MSVC toolchain on Windows) and
[`protoc`](https://github.com/protocolbuffers/protobuf) (Lance uses Protocol
Buffers). Set `PROTOC` to the compiler path if it isn't on `PATH`.

```bash
cargo build --release

# 1. Ingest the bz2 dump → Lance dataset
bitmako ingest        --input REAL.cxsmiles.bz2 --output compounds.lance

# 2. Build the inverted index (multi-pass, memory-bounded)
bitmako build-index   --lance compounds.lance   --output compounds.bitmako --bits-per-pass 64

# 3. Build the skip index from the index (single pass, no rebuild)
bitmako build-skip    --index compounds.bitmako --output compounds.skip

# 4. Build the flat fingerprint store
bitmako build-fp-store   --lance compounds.lance  --output compounds.fp

# 5. Build the flat property store (~27 GB, one Lance scan)
bitmako build-prop-store --lance compounds.lance  --output compounds.prop

# 6. Search (doc_id + Tanimoto only)
bitmako search --index compounds.bitmako --skip compounds.skip --fp-store compounds.fp \
    --query "OC(=O)c1ccccc1" --threshold 0.3 --top-k 10

# 6b. Search with SMILES / property lookup (add --lance)
bitmako search --index compounds.bitmako --skip compounds.skip --fp-store compounds.fp \
    --lance compounds.lance \
    --query "OC(=O)c1ccccc1" --threshold 0.3 --top-k 10

# 6c. Search with property filters — fast path (in-loop, no over-fetch)
bitmako search --index compounds.bitmako --skip compounds.skip --fp-store compounds.fp \
    --prop-store compounds.prop \
    --query "OC(=O)c1ccccc1" --threshold 0.1 --top-k 10 \
    --mw-max 350 --logp-max 3

# 6d. Search with property filters + SMILES output (combine --prop-store and --lance)
bitmako search --index compounds.bitmako --skip compounds.skip --fp-store compounds.fp \
    --prop-store compounds.prop --lance compounds.lance \
    --query "OC(=O)c1ccccc1" --threshold 0.1 --top-k 10 \
    --mw-max 350 --logp-max 3
```

> **Tanimoto is size-sensitive.** A small query (few set bits) can't reach a high
> Tanimoto against the large building-block compounds in REAL — even full
> containment of a ~12-bit fragment in a ~40-bit compound scores only ~0.3. Pick a
> threshold appropriate to your query's size, or query molecules of similar size to
> the corpus.

### `build-index` is memory-bounded

The builder runs in multiple passes, holding only `--bits-per-pass` posting lists in
RAM at once (default 32, ~33 scans). Higher values are faster but use more RAM:
64 bits/pass ≈ 17 scans on ~10 GB peak. The index uses 64-bit posting offsets, so the
posting section can exceed 4 GiB (the full build is ~70 GB).

## On-disk formats

**Index (`*.bitmako`, v2):**

```
[8B magic "BITMAKO1"][4B version=2][4B num_compounds][4B num_bits=1024]
[8B × 1024 posting-list offsets]
[1B × num_compounds compound popcounts][pad to 4B]
[posting-list data: per bit, 128-doc blocks of delta+varint doc_ids]
```

**Skip index (`*.skip`):** sidecar for streaming traversal. Per block, per bit, it
stores the decoder base and byte offset so a cursor can jump to the block holding a
target doc_id without decoding the list. Built in one pass over the index (~9 min
for the 70 GB index; ~6 GB output) — no rebuild required.

```
[8B magic "BMSKIP01"][4B version][4B num_bits]
[num_bits × { num_blocks: u32, data_off: u64 }]          // directory
[per bit: num_blocks × { base: u32, byte_offset: u64 }]
```

**Fingerprint store (`*.fp`):** flat array of `num_compounds × [u64; 16]`, little
-endian, indexed directly by `doc_id`. Memory-mapped at search time.

**Property store (`*.prop`):** flat array of `num_compounds × 5 × f32` (mw, logp,
rot_bonds, heavy_atoms, ring_count), 20 bytes per compound, little-endian, indexed
directly by `doc_id`. Memory-mapped at search time. Integer-valued fields (rot_bonds
etc.) are stored as f32 for uniform layout; all values round-trip exactly. At 1.36B
compounds the file is ~27 GB.

## Search performance

Measured on the full 1.36B-compound build, single workstation (32 GB RAM, cold
cache), before and after the streaming skip-cursor:

| Query | Before (full decode) | After (streaming + skip) |
|---|---|---|
| Aspirin @ 0.8 | 1097 s, ~26 GB RAM | 210 s, ~1 GB RAM |
| Medium drug-like @ 0.7 | **OOM crash** (8 GB alloc) | 341 s, ~1.3 GB RAM |
| Benzoic acid @ 0.2 (10 hits) | — | 105 s, ~2.8 GB RAM |

The decode/OOM bottleneck is gone — peak memory is now bounded by the query, not
the corpus. Remaining latency for low-threshold or common-fragment queries is
dominated by cold random fingerprint fetches from the 163 GB flat store (one disk
seek per surviving candidate); these benefit directly from page-cache warmth, an
SSD, or more RAM.

### Property-filtered search

When `--prop-store` is supplied, molecular properties are screened in O(1) inside
the WAND pivot loop — before the fingerprint fetch — so rejected compounds never
reach the exact-Tanimoto step and the dynamic threshold rises faster on genuinely
qualifying results:

| Query | Mode | Time |
|---|---|---|
| Benzoic acid @ 0.2, mw≤350, logp≤3, top-10 | `--prop-store` (in-loop) | **8 s** |
| Benzoic acid @ 0.2, mw≤350, logp≤3, top-10 | `--lance` post-filter (20× over-fetch) | 30 s |

3.7× faster than the Lance post-filter path for the same results.

## Known limitations

- **Cold fingerprint fetches dominate large candidate sets.** Each pivot that
  survives the count-based upper bound triggers a random 128-byte read from the 163 GB
  flat store. For low thresholds or common-fragment queries this is the main cost;
  it benefits directly from page-cache warmth, an SSD, or more RAM.
- **Property filters without `--prop-store` over-fetch 20×.** Passing `--mw-max` /
  `--logp-max` with `--lance` but no `--prop-store` falls back to a post-filter pass:
  WAND fetches `top_k × 20` candidates and drops non-matches. With very restrictive
  filters you may receive fewer than `top_k` results; increase `--top-k` to
  compensate, or build `compounds.prop` to use the fast in-loop path.

## Benchmarks

Full benchmark results — pruning rates, speedup vs. linear scan, property-filter
comparison — are in [BENCHMARKS.md](BENCHMARKS.md).

## Testing

```bash
cargo test --release
```

## License

MIT
