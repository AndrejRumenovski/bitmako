//! ETL pipeline benchmarks.
//!
//! Run with: cargo bench --bench etl_bench -- --output-format bencher
//!
//! Benchmarks cover:
//!   - ECFP4 fingerprint computation throughput
//!   - Property estimation throughput
//!   - Parallel chunk parsing (simulated)
//!   - Arrow RecordBatch serialization

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use bitmako::etl::fingerprint::compute_morgan_fp;
use bitmako::etl::parser::{parse_chunk_parallel, ParsedCompound};
use bitmako::etl::properties::compute_properties;
use bitmako::etl::reader::{RawLine, ReaderConfig};
use bitmako::etl::writer::build_record_batch;

/// Representative SMILES from Enamine REAL-like compounds
const BENCH_SMILES: &[&str] = &[
    "CCO",
    "c1ccccc1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
    "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O",
    "O=C(O)[C@@H](N)Cc1ccc(O)cc1",
    "CC1=CC(=O)c2ccccc2C1=O",
    "c1ccc2c(c1)ccc3ccccc23",
    "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21",
    "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "O=C(O)c1ccc(cc1)N",
    "CC(=O)c1ccc(cc1)O",
    "Fc1ccc(cc1)C(=O)c1ccc(F)cc1",
    "CC1(C)OCC(O1)C(=O)O",
    "O=C(O)c1ccc(Cl)cc1",
    "NC(=O)c1ccc(cc1)Br",
    "O=C(O)C1CCCCC1",
    "CN1CCN(CC1)c1ncnc2c1ccc(n2)Nc1ccc(cc1)C(F)(F)F",
    "CC(C)(C)c1ccc(cc1)C(=O)Cl",
];

fn bench_fingerprint_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("fingerprint");
    group.throughput(Throughput::Elements(BENCH_SMILES.len() as u64));

    group.bench_function("ecfp4_batch", |b| {
        b.iter(|| {
            for smiles in BENCH_SMILES {
                black_box(compute_morgan_fp(black_box(smiles)));
            }
        });
    });

    // Benchmark individual SMILES by complexity
    for smiles in BENCH_SMILES.iter().take(5) {
        group.bench_with_input(
            BenchmarkId::new("ecfp4_single", smiles.len()),
            smiles,
            |b, s| b.iter(|| compute_morgan_fp(black_box(s))),
        );
    }
    group.finish();
}

fn bench_property_throughput(c: &mut Criterion) {
    let mut group = c.benchmark_group("properties");
    group.throughput(Throughput::Elements(BENCH_SMILES.len() as u64));

    group.bench_function("compute_properties_batch", |b| {
        b.iter(|| {
            for smiles in BENCH_SMILES {
                black_box(compute_properties(black_box(smiles)));
            }
        });
    });
    group.finish();
}

fn bench_parallel_parsing(c: &mut Criterion) {
    let config = ReaderConfig::default();

    let mut group = c.benchmark_group("parallel_parsing");

    for chunk_size in [1_000usize, 10_000, 50_000] {
        let chunk: Vec<RawLine> = BENCH_SMILES
            .iter()
            .cycle()
            .take(chunk_size)
            .enumerate()
            .map(|(i, s)| RawLine {
                line_num: i + 1,
                raw: format!("{}\tZ{:09}", s, i),
            })
            .collect();

        group.throughput(Throughput::Elements(chunk_size as u64));
        group.bench_with_input(
            BenchmarkId::new("rayon_chunk", chunk_size),
            &chunk,
            |b, chunk| {
                b.iter(|| {
                    parse_chunk_parallel(chunk.clone(), &config)
                });
            },
        );
    }
    group.finish();
}

fn bench_record_batch_build(c: &mut Criterion) {
    let config = ReaderConfig::default();
    let chunk: Vec<RawLine> = BENCH_SMILES
        .iter()
        .cycle()
        .take(10_000)
        .enumerate()
        .map(|(i, s)| RawLine {
            line_num: i + 1,
            raw: format!("{}\tZ{:09}", s, i),
        })
        .collect();

    let (compounds, _) = parse_chunk_parallel(chunk, &config);

    c.bench_function("build_record_batch_10k", |b| {
        b.iter(|| build_record_batch(black_box(&compounds)));
    });
}

criterion_group!(
    benches,
    bench_fingerprint_throughput,
    bench_property_throughput,
    bench_parallel_parsing,
    bench_record_batch_build,
);
criterion_main!(benches);
