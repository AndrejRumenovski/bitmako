//! Morgan/ECFP fingerprint generation from SMILES strings.
//!
//! Implements ECFP4 (radius=2) producing 1024-bit fingerprints stored as [u64; 16].
//! Uses CRC32-based hashing for atom environment identifiers.

use crate::error::{BitMakoError, Result};
use crc32fast::Hasher as Crc32Hasher;

pub const FP_WORDS: usize = 16;
pub const FP_BITS: usize = FP_WORDS * 64; // 1024

/// Human-readable identifier for the fingerprint scheme every compound in the
/// corpus was encoded with. Single source of truth for display purposes (e.g.
/// the HTTP API's `/health` endpoint) — if BitMako ever supports more than one
/// fingerprint type, this is the one place that label needs to change.
pub const FINGERPRINT_KIND: &str = "ECFP4 (1024-bit Morgan)";

pub type Fingerprint = [u64; FP_WORDS];

/// Atomic numbers for organic subset atoms.
#[inline]
fn atomic_number(symbol: &str) -> u8 {
    match symbol {
        "H" => 1,
        "B" => 5,
        "C" | "c" => 6,
        "N" | "n" => 7,
        "O" | "o" => 8,
        "F" => 9,
        "P" | "p" => 15,
        "S" | "s" => 16,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        _ => 0,
    }
}

#[inline]
fn is_aromatic_symbol(c: char) -> bool {
    matches!(c, 'c' | 'n' | 'o' | 's' | 'p')
}

/// Bond order encoding
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Bond {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl Bond {
    fn order(self) -> u8 {
        match self {
            Bond::Single => 1,
            Bond::Double => 2,
            Bond::Triple => 3,
            Bond::Aromatic => 4,
        }
    }
}

/// Lightweight atom representation for fingerprint computation
#[derive(Clone, Debug)]
struct Atom {
    atomic_num: u8,
    is_aromatic: bool,
    charge: i8,
    h_count: u8,
}

/// Molecular graph: atoms + adjacency list
struct Molecule {
    atoms: Vec<Atom>,
    neighbors: Vec<Vec<(usize, Bond)>>,
}

impl Molecule {
    fn new() -> Self {
        Molecule {
            atoms: Vec::new(),
            neighbors: Vec::new(),
        }
    }

    fn add_atom(&mut self, atom: Atom) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.neighbors.push(Vec::new());
        idx
    }

    fn add_bond(&mut self, a: usize, b: usize, bond: Bond) {
        self.neighbors[a].push((b, bond));
        self.neighbors[b].push((a, bond));
    }
}

/// Parse SMILES into a molecular graph.
/// Handles: organic subset atoms, bracket atoms, branches, ring closures, explicit bonds.
fn parse_smiles(smiles: &str) -> Result<Molecule> {
    let mut mol = Molecule::new();
    let bytes = smiles.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    // Stack for branch points: stores (atom_index, current_bond)
    let mut branch_stack: Vec<(Option<usize>, Bond)> = Vec::new();
    // Ring closure map: ring_num -> (atom_idx, bond)
    let mut ring_closures: std::collections::HashMap<u32, (usize, Bond)> = std::collections::HashMap::new();

    let mut prev_atom: Option<usize> = None;
    let mut next_bond = Bond::Single;
    let mut explicit_bond = false;

    while i < len {
        let ch = bytes[i] as char;
        match ch {
            // Branch open
            '(' => {
                branch_stack.push((prev_atom, next_bond));
                next_bond = Bond::Single;
                explicit_bond = false;
                i += 1;
            }
            // Branch close
            ')' => {
                if let Some((pa, pb)) = branch_stack.pop() {
                    prev_atom = pa;
                    next_bond = pb;
                    explicit_bond = false;
                }
                i += 1;
            }
            // Explicit bonds
            '-' => {
                next_bond = Bond::Single;
                explicit_bond = true;
                i += 1;
            }
            '=' => {
                next_bond = Bond::Double;
                explicit_bond = true;
                i += 1;
            }
            '#' => {
                next_bond = Bond::Triple;
                explicit_bond = true;
                i += 1;
            }
            ':' => {
                next_bond = Bond::Aromatic;
                explicit_bond = true;
                i += 1;
            }
            '/' | '\\' => {
                // Stereo bonds treated as single
                next_bond = Bond::Single;
                explicit_bond = true;
                i += 1;
            }
            // Bracket atom: [C@@H], [NH3+], etc.
            '[' => {
                let start = i + 1;
                let end = smiles[start..]
                    .find(']')
                    .ok_or_else(|| BitMakoError::InvalidSmiles {
                        line: 0,
                        smiles: smiles.to_string(),
                    })?
                    + start;
                let inner = &smiles[start..end];
                let atom = parse_bracket_atom(inner)?;
                let idx = mol.add_atom(atom);
                if let Some(prev) = prev_atom {
                    let bond = infer_bond(&mol.atoms[prev], &mol.atoms[idx], next_bond, explicit_bond);
                    mol.add_bond(prev, idx, bond);
                }
                prev_atom = Some(idx);
                next_bond = Bond::Single;
                explicit_bond = false;
                i = end + 1;
            }
            // Ring closure (single digit)
            '0'..='9' => {
                let ring_num = (ch as u8 - b'0') as u32;
                handle_ring_closure(&mut mol, &mut ring_closures, ring_num, prev_atom, next_bond, explicit_bond)?;
                next_bond = Bond::Single;
                explicit_bond = false;
                i += 1;
            }
            // Two-digit ring closure: %12
            '%' if i + 2 < len => {
                let ring_num = (bytes[i + 1] - b'0') as u32 * 10 + (bytes[i + 2] - b'0') as u32;
                handle_ring_closure(&mut mol, &mut ring_closures, ring_num, prev_atom, next_bond, explicit_bond)?;
                next_bond = Bond::Single;
                explicit_bond = false;
                i += 3;
            }
            // Organic atoms
            'B' | 'C' | 'N' | 'O' | 'F' | 'P' | 'S' | 'I'
            | 'c' | 'n' | 'o' | 's' | 'p' => {
                let (sym, advance) = if ch == 'B' && i + 1 < len && bytes[i + 1] == b'r' {
                    ("Br", 2)
                } else if ch == 'C' && i + 1 < len && bytes[i + 1] == b'l' {
                    ("Cl", 2)
                } else {
                    (&smiles[i..i + 1], 1)
                };
                let aromatic = is_aromatic_symbol(ch);
                let atom = Atom {
                    atomic_num: atomic_number(sym),
                    is_aromatic: aromatic,
                    charge: 0,
                    h_count: default_h_count(atomic_number(sym), aromatic),
                };
                let idx = mol.add_atom(atom);
                if let Some(prev) = prev_atom {
                    let bond = infer_bond(&mol.atoms[prev], &mol.atoms[idx], next_bond, explicit_bond);
                    mol.add_bond(prev, idx, bond);
                }
                prev_atom = Some(idx);
                next_bond = Bond::Single;
                explicit_bond = false;
                i += advance;
            }
            // Skip stereo/unknown chars
            _ => { i += 1; }
        }
    }

    Ok(mol)
}

fn infer_bond(prev: &Atom, curr: &Atom, explicit: Bond, is_explicit: bool) -> Bond {
    if is_explicit {
        return explicit;
    }
    if prev.is_aromatic && curr.is_aromatic {
        Bond::Aromatic
    } else {
        Bond::Single
    }
}

/// Implicit hydrogen count used for ECFP atom invariants.
///
/// Deliberately separate from `properties::default_h_count`: this module and
/// `etl::properties` each run their own lightweight SMILES scan tuned for their
/// own purpose (fingerprint invariants vs. MW/LogP estimation), so their default-H
/// tables take different inputs and are not meant to be unified.
fn default_h_count(atomic_num: u8, aromatic: bool) -> u8 {
    match atomic_num {
        6 => if aromatic { 0 } else { 4 }, // C: 4 valence, c: aromatic handled later
        7 => 3,
        8 => 2,
        16 => 2,
        15 => 3,
        5 => 3,
        _ => 0,
    }
}

fn parse_bracket_atom(inner: &str) -> Result<Atom> {
    let mut chars = inner.chars().peekable();
    let mut atomic_num = 0u8;
    let mut is_aromatic = false;
    let mut charge = 0i8;
    let mut h_count = 0u8;

    // Skip isotope number
    while chars.peek().map(|c| c.is_ascii_digit()).unwrap_or(false) {
        chars.next();
    }

    // Element symbol (1 or 2 chars)
    if let Some(c) = chars.peek().copied() {
        if c.is_alphabetic() {
            chars.next();
            let mut sym = c.to_string();
            if let Some(next_c) = chars.peek().copied() {
                if next_c.is_lowercase() && next_c.is_alphabetic() && next_c != 'h' {
                    chars.next();
                    sym.push(next_c);
                }
            }
            is_aromatic = sym.chars().next().map(|c| c.is_lowercase()).unwrap_or(false);
            atomic_num = atomic_number(&sym);
        }
    }

    // Parse remaining: stereo (@), H count, charge
    for ch in chars {
        match ch {
            '@' => {}
            'H' => { h_count = 1; }
            '+' => { charge += 1; }
            '-' => { charge -= 1; }
            '0'..='9' => {}
            _ => {}
        }
    }

    Ok(Atom { atomic_num, is_aromatic, charge, h_count })
}

fn handle_ring_closure(
    mol: &mut Molecule,
    ring_closures: &mut std::collections::HashMap<u32, (usize, Bond)>,
    ring_num: u32,
    prev_atom: Option<usize>,
    bond: Bond,
    explicit_bond: bool,
) -> Result<()> {
    let Some(prev) = prev_atom else { return Ok(()); };
    match ring_closures.remove(&ring_num) {
        Some((open_idx, open_bond)) => {
            let b = if explicit_bond { bond } else { open_bond };
            let resolved = infer_bond(&mol.atoms[open_idx], &mol.atoms[prev], b, explicit_bond);
            mol.add_bond(open_idx, prev, resolved);
        }
        None => {
            ring_closures.insert(ring_num, (prev, bond));
        }
    }
    Ok(())
}

/// ECFP initial atom invariant hash
#[inline]
fn atom_invariant_hash(atom: &Atom, degree: usize) -> u32 {
    let mut h = Crc32Hasher::new();
    h.update(&[atom.atomic_num, degree as u8, atom.h_count, (atom.charge + 10) as u8, atom.is_aromatic as u8]);
    h.finalize()
}

/// Single ECFP iteration: hash each atom's environment
fn ecfp_iterate(
    mol: &Molecule,
    identifiers: &[u32],
) -> Vec<u32> {
    let mut new_ids = identifiers.to_vec();
    for (idx, atom) in mol.atoms.iter().enumerate() {
        let mut neighbors: Vec<(u8, u32)> = mol.neighbors[idx]
            .iter()
            .map(|(nidx, bond)| (bond.order(), identifiers[*nidx]))
            .collect();
        neighbors.sort_unstable();

        let mut h = Crc32Hasher::new();
        h.update(&identifiers[idx].to_le_bytes());
        h.update(&(atom.atomic_num as u32).to_le_bytes());
        for (order, nbr_id) in &neighbors {
            h.update(&[*order]);
            h.update(&nbr_id.to_le_bytes());
        }
        new_ids[idx] = h.finalize();
    }
    new_ids
}

/// Set a bit in a 1024-bit fingerprint stored as [u64; 16]
#[inline]
fn set_fp_bit(fp: &mut Fingerprint, feature_hash: u32) {
    let bit_pos = (feature_hash as usize) % FP_BITS;
    let word = bit_pos / 64;
    let bit = bit_pos % 64;
    fp[word] |= 1u64 << bit;
}

/// Compute ECFP4 (radius=2) 1024-bit Morgan fingerprint.
/// Returns zeroed fingerprint on unparseable SMILES rather than failing.
pub fn compute_morgan_fp(smiles: &str) -> Fingerprint {
    let mut fp = [0u64; FP_WORDS];

    let mol = match parse_smiles(smiles) {
        Ok(m) if !m.atoms.is_empty() => m,
        _ => return fp,
    };

    // Initial identifiers from atom invariants
    let mut identifiers: Vec<u32> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(idx, atom)| atom_invariant_hash(atom, mol.neighbors[idx].len()))
        .collect();

    // Fold radius-0 features (individual atoms)
    for &id in &identifiers {
        set_fp_bit(&mut fp, id);
    }

    // Iterate up to radius 2 (ECFP4)
    for _radius in 1..=2 {
        identifiers = ecfp_iterate(&mol, &identifiers);
        for &id in &identifiers {
            set_fp_bit(&mut fp, id);
        }
    }

    fp
}

/// Popcount across the full 1024-bit fingerprint
#[inline]
pub fn fp_popcount(fp: &Fingerprint) -> u32 {
    fp.iter().map(|w| w.count_ones()).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_smiles() {
        let fp = compute_morgan_fp("");
        assert_eq!(fp, [0u64; 16]);
    }

    #[test]
    fn test_ethanol_not_empty() {
        let fp = compute_morgan_fp("CCO");
        assert!(fp_popcount(&fp) > 0);
    }

    #[test]
    fn test_identical_smiles_same_fp() {
        let fp1 = compute_morgan_fp("c1ccccc1");
        let fp2 = compute_morgan_fp("c1ccccc1");
        assert_eq!(fp1, fp2);
    }

    #[test]
    fn test_different_molecules_differ() {
        let fp1 = compute_morgan_fp("CCO");
        let fp2 = compute_morgan_fp("CNC");
        assert_ne!(fp1, fp2);
    }
}
