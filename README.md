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
4. **Searches** for the top-k most Tanimoto-similar compounds to a query SMILES,
   using BMW dynamic pruning to avoid scanning the full corpus.

## Scale (this build)

| Artifact | Value |
|---|---|
| Compounds ingested | 1,364,304,490 (0 parse failures) |
| Lance dataset | ~333 GB |
| Inverted index | ~70 GB |
| Flat fingerprint store | ~163 GB (128 bytes/compound) |
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
      ├── build-index ─────▶ inverted index   (compounds.bitmako)
      └── build-fp-store ──▶ flat fingerprints (compounds.fp)
                                   │
 query SMILES ──▶ BMW engine ◀─────┘
                     │  (decodes only the query's active-bit posting lists,
                     │   block-max pruning, popcount upper bounds)
                     ▼
              top-k (doc_id, Tanimoto)
```

### Fingerprints

ECFP4 (radius 2), 1024 bits stored as `[u64; 16]`. A pure-Rust SMILES parser builds
the molecular graph; atom environments are hashed with CRC32 over two iterations.

### Tanimoto similarity

`|A ∩ B| / |A ∪ B|` = `popcount(A & B) / popcount(A | B)`, computed over the 16
words with POPCNT. Compiled with `target-cpu=native` so LLVM autovectorizes the
inner loop.

### Block-Max WAND

Posting lists are delta + LEB128-varint encoded in 128-doc blocks, each carrying the
max compound popcount in that block. Search:

1. Decode posting lists **only for the query's active bits** (memory stays
   proportional to the query, not the corpus).
2. Walk candidates in `doc_id` order; use the block max-popcount to compute a
   Tanimoto upper bound and **skip whole blocks** that can't beat the current
   threshold.
3. Apply a per-doc popcount upper bound before fetching the full fingerprint.
4. Fetch the fingerprint from the mmap'd flat store, compute exact Tanimoto, and
   update the top-k min-heap — raising the threshold as it fills.

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

# 3. Build the flat fingerprint store
bitmako build-fp-store --lance compounds.lance  --output compounds.fp

# 4. Search
bitmako search --index compounds.bitmako --fp-store compounds.fp \
    --query "CC(=O)Oc1ccccc1C(=O)O" --threshold 0.8 --top-k 10 \
    [--mw-max 500 --logp-max 5]
```

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

**Fingerprint store (`*.fp`):** flat array of `num_compounds × [u64; 16]`, little
-endian, indexed directly by `doc_id`. Memory-mapped at search time.

## Known limitations / next steps

- **Common-fragment queries are slow.** The engine fully decodes each active-bit
  posting list into a `Vec<u32>` before walking it. For a query whose fragments are
  common across the corpus (e.g. aspirin's benzene + carbonyl bits), those lists run
  to hundreds of millions of entries, so decoding dominates — an aspirin query over
  the full 1.36B set took ~19 min (cold cache) and ~25 GB RAM. The fix is to iterate
  the varint stream directly and use block-max metadata to skip blocks without
  decoding them. Selective (drug-like, rarer-fragment) queries are far cheaper.
- **Results are `doc_id`s, not SMILES.** Mapping back to compound IDs requires a join
  against the Lance dataset; not yet wired into the `search` command.
- **Property filters** (`--mw-max`, `--logp-max`) are parsed into the query but the
  property columns aren't yet loaded during search.

## Testing

```bash
cargo test --release
```

## License

MIT
