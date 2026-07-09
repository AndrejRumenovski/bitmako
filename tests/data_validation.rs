//! Data-validation tests: malformed/corrupt/truncated/misaligned on-disk
//! files must return `Err`, never panic, and query construction must reject
//! invalid parameters cleanly. These guard the boundary between "real files
//! on disk" and the in-process data structures that mmap them.

use bitmako::etl::fingerprint::compute_morgan_fp;
use bitmako::etl::properties::MolecularProperties;
use bitmako::index::skip::SkipIndex;
use bitmako::index::IndexReader;
use bitmako::search::fp_store::FpStore;
use bitmako::search::prop_store::PropStore;
use bitmako::search::query::{validate_query, PropertyField, PropertyFilter, SimilarityQuery};
use tempfile::NamedTempFile;

fn write_bytes(bytes: &[u8]) -> NamedTempFile {
    let tmp = NamedTempFile::new().unwrap();
    std::fs::write(tmp.path(), bytes).unwrap();
    tmp
}

// ---- IndexReader ----

#[test]
fn index_reader_rejects_bad_magic() {
    let tmp = write_bytes(b"NOTBITMAKO00000000000");
    assert!(IndexReader::open(tmp.path()).is_err());
}

#[test]
fn index_reader_rejects_empty_file() {
    let tmp = write_bytes(b"");
    assert!(IndexReader::open(tmp.path()).is_err());
}

#[test]
fn index_reader_rejects_header_shorter_than_minimum() {
    let tmp = write_bytes(b"BITMAKO1"); // magic only, no version/counts
    assert!(IndexReader::open(tmp.path()).is_err());
}

#[test]
fn index_reader_rejects_unsupported_version() {
    let mut data = Vec::new();
    data.extend_from_slice(b"BITMAKO1");
    data.extend_from_slice(&42u32.to_le_bytes());
    data.extend_from_slice(&0u32.to_le_bytes());
    data.extend_from_slice(&0u32.to_le_bytes());
    let tmp = write_bytes(&data);
    assert!(IndexReader::open(tmp.path()).is_err());
}

#[test]
fn index_reader_rejects_offset_table_beyond_file_length() {
    let mut data = Vec::new();
    data.extend_from_slice(b"BITMAKO1");
    data.extend_from_slice(&2u32.to_le_bytes());
    data.extend_from_slice(&1u32.to_le_bytes()); // num_compounds
    data.extend_from_slice(&1024u32.to_le_bytes()); // num_bits — offsets table alone is 8KB
    // File ends right after the header; no offset table, pops, or postings.
    let tmp = write_bytes(&data);
    assert!(IndexReader::open(tmp.path()).is_err());
}

// ---- SkipIndex ----

#[test]
fn skip_index_rejects_bad_magic() {
    let tmp = write_bytes(b"NOTASKIP12345678");
    assert!(SkipIndex::open(tmp.path()).is_err());
}

#[test]
fn skip_index_rejects_truncated_directory() {
    let mut data = Vec::new();
    data.extend_from_slice(b"BMSKIP01");
    data.extend_from_slice(&1u32.to_le_bytes());
    data.extend_from_slice(&1000u32.to_le_bytes()); // claims 1000 bits, no directory follows
    let tmp = write_bytes(&data);
    assert!(SkipIndex::open(tmp.path()).is_err());
}

#[test]
fn skip_index_rejects_empty_file() {
    let tmp = write_bytes(b"");
    assert!(SkipIndex::open(tmp.path()).is_err());
}

// ---- FpStore ----

#[test]
fn fp_store_rejects_size_not_a_multiple_of_fp_bytes() {
    let tmp = write_bytes(&[0u8; 100]); // 100 is not a multiple of 128
    assert!(FpStore::open(tmp.path()).is_err());
}

#[test]
fn fp_store_accepts_empty_file() {
    // Zero fingerprints is a degenerate but valid (size % FP_BYTES == 0) store.
    let tmp = write_bytes(&[]);
    let store = FpStore::open(tmp.path()).unwrap();
    assert_eq!(store.len(), 0);
    assert!(store.is_empty());
    assert!(store.get(0).is_none());
}

// ---- PropStore ----

#[test]
fn prop_store_rejects_size_not_a_multiple_of_prop_bytes() {
    let tmp = write_bytes(&[0u8; 17]); // PROP_BYTES = 20
    assert!(PropStore::open(tmp.path()).is_err());
}

#[test]
fn prop_store_accepts_empty_file() {
    let tmp = write_bytes(&[]);
    let store = PropStore::open(tmp.path()).unwrap();
    assert_eq!(store.len(), 0);
    assert!(store.get(0).is_none());
}

// ---- validate_query ----

#[test]
fn validate_query_rejects_threshold_above_one() {
    let fp = compute_morgan_fp("CCO");
    let query = SimilarityQuery::new(fp, 1.5, 10);
    assert!(validate_query(&query).is_err());
}

#[test]
fn validate_query_rejects_negative_threshold() {
    let fp = compute_morgan_fp("CCO");
    let query = SimilarityQuery::new(fp, -0.1, 10);
    assert!(validate_query(&query).is_err());
}

#[test]
fn validate_query_accepts_threshold_boundaries() {
    let fp = compute_morgan_fp("CCO");
    assert!(validate_query(&SimilarityQuery::new(fp, 0.0, 10)).is_ok());
    assert!(validate_query(&SimilarityQuery::new(fp, 1.0, 10)).is_ok());
}

#[test]
fn validate_query_rejects_zero_top_k() {
    let fp = compute_morgan_fp("CCO");
    let query = SimilarityQuery::new(fp, 0.5, 0);
    assert!(validate_query(&query).is_err());
}

#[test]
fn validate_query_rejects_all_zero_fingerprint() {
    let query = SimilarityQuery::new([0u64; 16], 0.5, 10);
    assert!(validate_query(&query).is_err());
}

// ---- PropertyFilter boundary behavior ----

#[test]
fn property_filter_max_is_inclusive() {
    let filter = PropertyFilter { field: PropertyField::MolWeight, min: None, max: Some(100.0) };
    let at_max = MolecularProperties { mw: 100.0, logp: 0.0, rot_bonds: 0, heavy_atoms: 0, ring_count: 0 };
    let just_over = MolecularProperties { mw: 100.001, ..at_max };
    assert!(filter.passes(&at_max), "value exactly at max should pass (inclusive)");
    assert!(!filter.passes(&just_over));
}

#[test]
fn property_filter_min_is_inclusive() {
    let filter = PropertyFilter { field: PropertyField::LogP, min: Some(1.0), max: None };
    let at_min = MolecularProperties { mw: 0.0, logp: 1.0, rot_bonds: 0, heavy_atoms: 0, ring_count: 0 };
    let just_under = MolecularProperties { logp: 0.999, ..at_min };
    assert!(filter.passes(&at_min), "value exactly at min should pass (inclusive)");
    assert!(!filter.passes(&just_under));
}

#[test]
fn property_filter_with_no_bounds_always_passes() {
    let filter = PropertyFilter { field: PropertyField::RingCount, min: None, max: None };
    let props = MolecularProperties { mw: 1e9, logp: -1e9, rot_bonds: u32::MAX, heavy_atoms: 0, ring_count: 0 };
    assert!(filter.passes(&props));
}

#[test]
fn property_filter_range_excludes_outside_both_bounds() {
    let filter = PropertyFilter { field: PropertyField::HeavyAtomCount, min: Some(5.0), max: Some(10.0) };
    let below = MolecularProperties { mw: 0.0, logp: 0.0, rot_bonds: 0, heavy_atoms: 4, ring_count: 0 };
    let inside = MolecularProperties { heavy_atoms: 7, ..below };
    let above = MolecularProperties { heavy_atoms: 11, ..below };
    assert!(!filter.passes(&below));
    assert!(filter.passes(&inside));
    assert!(!filter.passes(&above));
}
