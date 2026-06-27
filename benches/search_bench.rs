//! BMW search engine benchmarks.
//!
//! Measures:
//!   - Tanimoto similarity computation throughput
//!   - Block-Max WAND search at varying thresholds (0.5, 0.7, 0.9)
//!   - Posting list serialization/deserialization
//!   - Block-skipping efficiency as a function of threshold

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use bitmako::etl::fingerprint::{compute_morgan_fp, fp_popcount, Fingerprint};
use bitmako::index::builder::IndexBuilder;
use bitmako::index::posting_list::PostingList;
use bitmako::index::IndexReader;
use bitmako::search::query::SimilarityQuery;
use bitmako::search::tanimoto::{tanimoto, tanimoto_upper_bound};
use bitmako::search::wand::BmwEngine;

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
];

/// Build a test index with synthetic compounds (from repeated BENCH_SMILES)
fn build_bench_index(n_compounds: usize) -> (Vec<Fingerprint>, tempfile::NamedTempFile) {
    let fps: Vec<Fingerprint> = BENCH_SMILES
        .iter()
        .cycle()
        .take(n_compounds)
        .map(|s| compute_morgan_fp(s))
        .collect();

    let mut builder = IndexBuilder::new();
    for (i, fp) in fps.iter().enumerate() {
        builder.add_compound(i as u32, fp);
    }

    let tmp = tempfile::NamedTempFile::new().unwrap();
    builder.write_index(tmp.path()).unwrap();
    (fps, tmp)
}

fn bench_tanimoto(c: &mut Criterion) {
    let fps: Vec<Fingerprint> = BENCH_SMILES.iter().map(|s| compute_morgan_fp(s)).collect();
    let n = fps.len();

    let mut group = c.benchmark_group("tanimoto");
    group.throughput(Throughput::Elements((n * n) as u64));

    group.bench_function("all_pairs_10x10", |b| {
        b.iter(|| {
            let mut sum = 0.0f32;
            for a in &fps {
                for b_fp in &fps {
                    sum += tanimoto(black_box(a), black_box(b_fp));
                }
            }
            black_box(sum)
        });
    });

    group.bench_function("upper_bound_10x10", |b| {
        let pops: Vec<u32> = fps.iter().map(fp_popcount).collect();
        b.iter(|| {
            let mut sum = 0.0f32;
            for &pa in &pops {
                for &pb in &pops {
                    sum += tanimoto_upper_bound(black_box(pa), black_box(pb));
                }
            }
            black_box(sum)
        });
    });
    group.finish();
}

fn bench_posting_list_roundtrip(c: &mut Criterion) {
    let doc_ids: Vec<u32> = (0u32..10_000).collect();
    let pops = vec![30u8; 79]; // ceil(10000/128) blocks
    let pl = PostingList { doc_ids: doc_ids.clone(), block_max_pop: pops };

    let mut group = c.benchmark_group("posting_list");
    group.throughput(Throughput::Elements(10_000));

    group.bench_function("serialize_10k", |b| {
        b.iter(|| pl.serialize().unwrap())
    });

    let bytes = pl.serialize().unwrap();
    group.bench_function("deserialize_10k", |b| {
        b.iter(|| PostingList::deserialize(black_box(&bytes)).unwrap())
    });

    group.bench_function("advance_to_sequential", |b| {
        b.iter(|| {
            let mut pos = 0;
            for target in (0u32..10_000).step_by(10) {
                pos = pl.advance_to(pos, target);
            }
            black_box(pos)
        });
    });
    group.finish();
}

fn bench_bmw_search(c: &mut Criterion) {
    let n_compounds = 10_000;
    let (fps, tmp) = build_bench_index(n_compounds);
    let index = IndexReader::open(tmp.path()).unwrap();
    let engine = BmwEngine::new(&index);

    let query_fp = compute_morgan_fp("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
    let fps_ref = &fps;

    let mut group = c.benchmark_group("bmw_search");
    group.throughput(Throughput::Elements(n_compounds as u64));

    for threshold in [0.5f32, 0.7, 0.9] {
        group.bench_with_input(
            BenchmarkId::new("threshold", format!("{:.1}", threshold)),
            &threshold,
            |b, &thresh| {
                let query = SimilarityQuery::new(query_fp, thresh, 50);
                b.iter(|| {
                    engine.search(black_box(&query), |doc_id| {
                        fps_ref.get(doc_id as usize).copied()
                    })
                    .unwrap()
                });
            },
        );
    }
    group.finish();
}

fn bench_bmw_scale(c: &mut Criterion) {
    let mut group = c.benchmark_group("bmw_scale");
    let query_fp = compute_morgan_fp("c1ccccc1");

    for n_compounds in [1_000usize, 10_000, 100_000] {
        let (fps, tmp) = build_bench_index(n_compounds);
        let index = IndexReader::open(tmp.path()).unwrap();
        let engine = BmwEngine::new(&index);
        let fps_ref = fps.clone();
        let query = SimilarityQuery::new(query_fp, 0.7, 50);

        group.throughput(Throughput::Elements(n_compounds as u64));
        group.bench_with_input(
            BenchmarkId::new("n_compounds", n_compounds),
            &n_compounds,
            |b, _| {
                b.iter(|| {
                    engine.search(black_box(&query), |doc_id| {
                        fps_ref.get(doc_id as usize).copied()
                    })
                    .unwrap()
                });
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_tanimoto,
    bench_posting_list_roundtrip,
    bench_bmw_search,
    bench_bmw_scale,
);
criterion_main!(benches);
