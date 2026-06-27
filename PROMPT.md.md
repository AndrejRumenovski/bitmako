# Role & Objective
You are a Staff Systems and Data Platform Engineer specializing in high-performance, zero-copy, columnar data architectures in Rust. Your goal is to build a production-grade, ultra-high-throughput chemoinformatics ETL and query platform capable of processing a 1.4-billion compound Enamine REAL subset (~40GB .bz2 compressed file). 

The final system must achieve zero-copy parsing, compute-ready compressed columnar storage, and a custom Block-Max WAND (BMW) similarity and property search engine compiled on the STABLE Rust toolchain.

---

## Technical Specifications & Architecture

### Phase 1: High-Performance ETL (.bz2 to Lance/Vortex-Array)
1. **Streaming Ingestion:** Stream the `.bz2` file on-the-fly using `bzip2` or `flate2` into a buffered reader. Do NOT decompress to disk or read the entire file into RAM.
2. **Zero-Copy Parsing:** Parse the SMILES strings and properties column-by-column. Leverage string slices (`&str`) and lifetime boundaries to ensure no arbitrary `String` allocations happen during the parsing phase.
3. **Parallel Ingestion:** Use `rayon` to implement a chunk-based work-stealing parallel processing loop. Accumulate blocks of rows (e.g., 100,000 compounds per chunk), parse structural data, compute basic properties (MW, LogP, Rotatable Bonds), and generate 1024-bit Morgan fingerprints.
4. **Columnar Array Layout (Vortex & Lance):** Integrate the `vortex-array` crate. Convert the chunk data directly into Apache Arrow `RecordBatch` formats, mapping the fingerprints to a `FixedSizeList(UInt64, 16)`. Use Vortex’s compressed, compute-ready memory representations (like FastLanes or ALP) to compress columns while keeping them accessible to SIMD instructions without full decompression. Write the final materialized fragments out to a `lance` dataset layout.
5. **General Property Querying:** The storage architecture must support filter push-down capabilities so users can query structures based on scalar chemical properties (e.g., `LogP <= 5.0 AND MW < 500`) alongside structural similarity.

### Phase 2: Inverted Index & Block-Max WAND Engine
1. **Stable Toolchain Refactor:** Clone or pull down the `rise-rs` inverted index search engine crate/repository. Audit the codebase, strip out any nightly-only features (such as experimental SIMD intrinsics or unstable features), and refactor them to use standard, high-performance stable Rust structures (e.g., `core::simd` where stable, or explicit compiler-vectorized loops).
2. **Columnar Execution Engine:** Build the search execution loop to operate over column arrays rather than row iterators. Postings lists for fingerprint bits must be tightly packed in memory (utilizing bit-packing or frame-of-reference compression matching the Vortex/Lance ecosystem).
3. **Block-Max WAND (BMW) Implementation:** - Segment the postings lists into fixed-size blocks (e.g., 64 or 128 elements).
   - Pre-calculate and store the maximum possible Tanimoto score contribution per block.
   - Implement the BMW dynamic pruning loop. The engine must track the running top-$k$ threshold score and aggressively skip non-viable postings blocks using direct byte-offset jumps (`.advance_to(target_id)`) without reading the underlying document payloads from disk/RAM.
4. **Streaming Execution:** Ensure queries process lazily via streaming iterators over Lance dataset fragments, maintaining a tight, constant memory overhead regardless of database size.

---

## Quality & Performance Requirements (No Slop)
* **Zero Unwraps:** All errors must use explicit error propagation. Implement clean custom errors using `thiserror`.
* **Hardware Optimization:** Guide LLVM by using vectorizable loops or native instructions (`.count_ones()`). Ensure code compiles natively with `RUSTFLAGS="-C target-cpu=native"`.
* **Microbenchmarking:** Provide a robust benchmarking harness using `criterion` or `divan` inside the `benches/` directory. Benchmark search speed and block-skipping efficiency against varying similarity thresholds (e.g., Tanimoto >= 0.7 vs >= 0.9) and complex property intersection filters.
* **Diagnostics:** Use the `tracing` crate for asynchronous, structured logging instead of raw `println!` macros.

Generate the complete project layout, including `Cargo.toml`, the Phase 1 streaming ETL module, the Phase 2 stable-compatible BMW execution engine, and the benchmarking architecture.