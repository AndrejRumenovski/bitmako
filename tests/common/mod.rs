//! Shared black-box test helpers: build a small in-memory-backed `Searcher`
//! from a list of SMILES, without ever needing the real 1.4B-compound data
//! files on disk. Used by every file in `tests/` that needs a working
//! `Searcher` over a synthetic corpus.
//!
//! `cargo test` compiles this `mod common;` fresh into every integration-test
//! binary that includes it, and each binary only calls the subset of helpers
//! it needs — so dead_code warnings here are expected noise, not a signal
//! about the crate under test.
#![allow(dead_code)]

use std::io::Write;

use bitmako::etl::fingerprint::{compute_morgan_fp, Fingerprint};
use bitmako::etl::properties::{compute_properties, MolecularProperties};
use bitmako::index::builder::IndexBuilder;
use bitmako::index::skip::SkipIndex;
use bitmako::index::IndexReader;
use bitmako::search::fp_store::{FpStore, FP_BYTES};
use bitmako::search::prop_store::{PropStore, PROP_BYTES};
use bitmako::search::Searcher;
use tempfile::NamedTempFile;

/// A synthetic corpus plus a ready-to-query `Searcher` over it.
///
/// The temp files are kept alive for the lifetime of the struct: `Searcher`
/// holds them memory-mapped, and on some platforms a mapped file can't be
/// deleted out from under the mapping.
pub struct Corpus {
    pub fps: Vec<Fingerprint>,
    pub props: Vec<MolecularProperties>,
    pub searcher: Searcher,
    _index_tmp: NamedTempFile,
    _fp_tmp: NamedTempFile,
    _prop_tmp: NamedTempFile,
}

/// Build a corpus (fingerprints, properties, full index/skip/fp/prop stores,
/// and an attached `Searcher`) from a list of SMILES, indexed in the order given.
pub fn build_corpus(smiles: &[&str]) -> Corpus {
    let fps: Vec<Fingerprint> = smiles.iter().map(|s| compute_morgan_fp(s)).collect();
    let props: Vec<MolecularProperties> = smiles.iter().map(|s| compute_properties(s)).collect();

    let mut builder = IndexBuilder::new();
    for (i, fp) in fps.iter().enumerate() {
        builder.add_compound(i as u32, fp);
    }
    let index_tmp = NamedTempFile::new().expect("create temp index file");
    builder.write_index(index_tmp.path()).expect("write index");
    let index = IndexReader::open(index_tmp.path()).expect("open index");
    let skip = SkipIndex::build_in_memory(&index).expect("build skip index");

    let fp_tmp = write_fp_store(&fps);
    let fp_store = FpStore::open(fp_tmp.path()).expect("open fp store");

    let prop_tmp = write_prop_store(&props);
    let prop_store = PropStore::open(prop_tmp.path()).expect("open prop store");

    let searcher = Searcher::open_from_index(index, skip, fp_store).with_prop_store(prop_store);

    Corpus {
        fps,
        props,
        searcher,
        _index_tmp: index_tmp,
        _fp_tmp: fp_tmp,
        _prop_tmp: prop_tmp,
    }
}

fn write_fp_store(fps: &[Fingerprint]) -> NamedTempFile {
    let mut tmp = NamedTempFile::new().expect("create temp fp store");
    for fp in fps {
        let mut buf = [0u8; FP_BYTES];
        for (i, &w) in fp.iter().enumerate() {
            buf[i * 8..i * 8 + 8].copy_from_slice(&w.to_le_bytes());
        }
        tmp.write_all(&buf).expect("write fp record");
    }
    tmp.flush().expect("flush fp store");
    tmp
}

fn write_prop_store(props: &[MolecularProperties]) -> NamedTempFile {
    let mut tmp = NamedTempFile::new().expect("create temp prop store");
    for p in props {
        let mut buf = [0u8; PROP_BYTES];
        buf[0..4].copy_from_slice(&p.mw.to_le_bytes());
        buf[4..8].copy_from_slice(&p.logp.to_le_bytes());
        buf[8..12].copy_from_slice(&(p.rot_bonds as f32).to_le_bytes());
        buf[12..16].copy_from_slice(&(p.heavy_atoms as f32).to_le_bytes());
        buf[16..20].copy_from_slice(&(p.ring_count as f32).to_le_bytes());
        tmp.write_all(&buf).expect("write prop record");
    }
    tmp.flush().expect("flush prop store");
    tmp
}

/// Build a corpus directly from pre-computed fingerprints, skipping SMILES
/// entirely (used for synthetic/perf tests where only the bit patterns
/// matter). Properties are left at their default (all-zero) values.
pub fn build_corpus_from_fingerprints(fps: Vec<Fingerprint>) -> Corpus {
    let props = vec![MolecularProperties::default(); fps.len()];

    let mut builder = IndexBuilder::new();
    for (i, fp) in fps.iter().enumerate() {
        builder.add_compound(i as u32, fp);
    }
    let index_tmp = NamedTempFile::new().expect("create temp index file");
    builder.write_index(index_tmp.path()).expect("write index");
    let index = IndexReader::open(index_tmp.path()).expect("open index");
    let skip = SkipIndex::build_in_memory(&index).expect("build skip index");

    let fp_tmp = write_fp_store(&fps);
    let fp_store = FpStore::open(fp_tmp.path()).expect("open fp store");

    let prop_tmp = write_prop_store(&props);
    let prop_store = PropStore::open(prop_tmp.path()).expect("open prop store");

    let searcher = Searcher::open_from_index(index, skip, fp_store).with_prop_store(prop_store);

    Corpus {
        fps,
        props,
        searcher,
        _index_tmp: index_tmp,
        _fp_tmp: fp_tmp,
        _prop_tmp: prop_tmp,
    }
}

/// Deterministic xorshift64 PRNG — fixed seed, no external `rand` dependency,
/// reproducible across runs/machines for synthetic corpus generation.
pub struct Xorshift64(u64);

impl Xorshift64 {
    pub fn new(seed: u64) -> Self {
        Xorshift64(seed | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

/// Generate a synthetic sparse fingerprint with roughly `num_bits` bits set,
/// drawn uniformly from the 1024-bit space using `rng`.
pub fn random_fingerprint(rng: &mut Xorshift64, num_bits: usize) -> Fingerprint {
    let mut fp = [0u64; 16];
    for _ in 0..num_bits {
        let bit = (rng.next_u64() % 1024) as usize;
        fp[bit / 64] |= 1u64 << (bit % 64);
    }
    fp
}

/// Exhaustive Tanimoto scan over the whole corpus — ground truth to check
/// WAND's pruning against. Mirrors WAND semantics: only docs with positive
/// overlap (score > 0) are candidates, since a doc sharing zero query bits
/// appears in none of the query's posting lists.
pub fn brute_force(query_fp: &Fingerprint, fps: &[Fingerprint], threshold: f32) -> Vec<(u32, f32)> {
    use bitmako::search::tanimoto::tanimoto;
    let mut out: Vec<(u32, f32)> = fps
        .iter()
        .enumerate()
        .map(|(i, fp)| (i as u32, tanimoto(query_fp, fp)))
        .filter(|(_, s)| *s >= threshold && *s > 0.0)
        .collect();
    out.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    out
}

/// A chemically diverse SMILES set spanning small/large, aromatic/aliphatic,
/// and ring/chain structures, so query popcounts range from small to large
/// and posting lists range from rare to common bits.
pub const DIVERSE_SMILES: &[&str] = &[
    "CCO", "CCCO", "CCCCO", "CCCCCO", "c1ccccc1", "Cc1ccccc1",
    "CN", "CC(=O)O", "CNC(=O)c1ccccc1", "c1ccncc1", "OCC(O)CO",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(=O)Oc1ccccc1C(=O)O",
    "C1CCCCC1", "c1ccc2ccccc2c1", "NCCO", "CCN(CC)CC",
    "CC(N)C(=O)O", "OC(=O)c1ccccc1", "Clc1ccccc1", "CC(C)O",
    "c1ccc(cc1)c1ccccc1", "O=C(O)c1ccccc1O",
];
