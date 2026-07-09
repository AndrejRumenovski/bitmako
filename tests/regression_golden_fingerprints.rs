//! Golden/regression tests: hardcoded known-answer values for canonical
//! molecules, computed once and locked in. If `compute_morgan_fp`,
//! `compute_properties`, or `tanimoto` ever change behavior — intentionally
//! or by accident — these tests catch it immediately, since the numbers here
//! are not recomputed at test time from the same code path being tested.
//!
//! Values were captured by running the current implementation once (see git
//! history for the one-off script) — they are a snapshot of *current*
//! behavior, not independently-verified ground truth. If a future change
//! deliberately improves fingerprinting/property-estimation accuracy, these
//! goldens must be regenerated and reviewed as part of that change, not
//! silently patched to make tests pass.

use bitmako::etl::fingerprint::{compute_morgan_fp, fp_popcount};
use bitmako::etl::properties::compute_properties;
use bitmako::search::tanimoto::tanimoto;

const ETHANOL: &str = "CCO";
const BENZENE: &str = "c1ccccc1";
const ASPIRIN: &str = "CC(=O)Oc1ccccc1C(=O)O";
const CAFFEINE: &str = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
const IBUPROFEN: &str = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";

#[test]
fn golden_fingerprint_popcounts() {
    assert_eq!(fp_popcount(&compute_morgan_fp(ETHANOL)), 9);
    assert_eq!(fp_popcount(&compute_morgan_fp(BENZENE)), 3);
    assert_eq!(fp_popcount(&compute_morgan_fp(ASPIRIN)), 27);
    assert_eq!(fp_popcount(&compute_morgan_fp(CAFFEINE)), 27);
    assert_eq!(fp_popcount(&compute_morgan_fp(IBUPROFEN)), 28);
}

#[test]
fn golden_molecular_properties() {
    let e = compute_properties(ETHANOL);
    assert!((e.mw - 50.1000).abs() < 1e-3);
    assert!((e.logp - 0.0820).abs() < 1e-3);
    assert_eq!(e.rot_bonds, 0);
    assert_eq!(e.heavy_atoms, 3);
    assert_eq!(e.ring_count, 0);

    let b = compute_properties(BENZENE);
    assert!((b.mw - 96.2556).abs() < 1e-3);
    assert!((b.logp - 0.9486).abs() < 1e-3);
    assert_eq!(b.rot_bonds, 1);
    assert_eq!(b.heavy_atoms, 6);
    assert_eq!(b.ring_count, 1);

    let a = compute_properties(ASPIRIN);
    assert!((a.mw - 216.4426).abs() < 1e-3);
    assert!((a.logp - 0.5561).abs() < 1e-3);
    assert_eq!(a.rot_bonds, 4);
    assert_eq!(a.heavy_atoms, 13);
    assert_eq!(a.ring_count, 1);

    let c = compute_properties(CAFFEINE);
    assert!((c.mw - 232.4932).abs() < 1e-3);
    assert!((c.logp - (-1.1120)).abs() < 1e-3);
    assert_eq!(c.rot_bonds, 2);
    assert_eq!(c.heavy_atoms, 14);
    assert_eq!(c.ring_count, 2);

    let i = compute_properties(IBUPROFEN);
    assert!((i.mw - 244.5834).abs() < 1e-3);
    assert!((i.logp - 1.5449).abs() < 1e-3);
    assert_eq!(i.rot_bonds, 5);
    assert_eq!(i.heavy_atoms, 15);
    assert_eq!(i.ring_count, 1);
}

#[test]
fn golden_pairwise_tanimoto_scores() {
    let ethanol = compute_morgan_fp(ETHANOL);
    let benzene = compute_morgan_fp(BENZENE);
    let aspirin = compute_morgan_fp(ASPIRIN);
    let caffeine = compute_morgan_fp(CAFFEINE);
    let ibuprofen = compute_morgan_fp(IBUPROFEN);

    let cases: [(&str, &bitmako::etl::fingerprint::Fingerprint, &bitmako::etl::fingerprint::Fingerprint, f32); 10] = [
        ("ethanol/benzene", &ethanol, &benzene, 0.000000),
        ("ethanol/aspirin", &ethanol, &aspirin, 0.058824),
        ("ethanol/caffeine", &ethanol, &caffeine, 0.058824),
        ("ethanol/ibuprofen", &ethanol, &ibuprofen, 0.121212),
        ("benzene/aspirin", &benzene, &aspirin, 0.071429),
        ("benzene/caffeine", &benzene, &caffeine, 0.071429),
        ("benzene/ibuprofen", &benzene, &ibuprofen, 0.033333),
        ("aspirin/caffeine", &aspirin, &caffeine, 0.148936),
        ("aspirin/ibuprofen", &aspirin, &ibuprofen, 0.195652),
        ("caffeine/ibuprofen", &caffeine, &ibuprofen, 0.100000),
    ];

    for (label, a, b, expected) in cases {
        let got = tanimoto(a, b);
        assert!(
            (got - expected).abs() < 1e-5,
            "tanimoto({label}) drifted: expected {expected}, got {got}"
        );
    }
}

#[test]
fn golden_self_tanimoto_is_always_one() {
    for smiles in [ETHANOL, BENZENE, ASPIRIN, CAFFEINE, IBUPROFEN] {
        let fp = compute_morgan_fp(smiles);
        assert!((tanimoto(&fp, &fp) - 1.0).abs() < 1e-6, "{smiles} self-tanimoto must be 1.0");
    }
}
