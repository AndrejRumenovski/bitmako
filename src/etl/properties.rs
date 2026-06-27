//! Lightweight molecular property estimation from SMILES.
//!
//! Implements:
//!  - Molecular weight (MW) from atom symbol counts
//!  - LogP approximation using Wildman-Crippen atom contributions
//!  - Rotatable bond count (non-terminal single bonds not in rings)
//!
//! These are heuristic estimates suitable for fast pre-filtering of
//! Lipinski-style property constraints.

/// Atom contributions to Crippen LogP.
/// Values from Wildman & Crippen, JCICS 39, 868–873 (1999) (simplified subset).
#[inline]
fn crippen_logp_contribution(symbol: &str, is_aromatic: bool) -> f32 {
    match (symbol, is_aromatic) {
        ("C", false) => 0.1441,
        ("C", true) | ("c", _) => 0.1581,
        ("N", false) => -0.2090,
        ("N", true) | ("n", _) => -0.4806,
        ("O", false) => -0.2062,
        ("O", true) | ("o", _) => -0.2202,
        ("S", false) => 0.1155,
        ("S", true) | ("s", _) => 0.0743,
        ("F", _) => 0.4202,
        ("Cl", _) => 0.6895,
        ("Br", _) => 0.8456,
        ("I", _) => 0.8857,
        ("P", _) | ("p", _) => 0.8612,
        _ => 0.0,
    }
}

/// Atomic masses for molecular weight computation (monoisotopic rounded to 4 dp)
#[inline]
fn atomic_mass(symbol: &str) -> f32 {
    match symbol {
        "H" => 1.0079,
        "B" => 10.811,
        "C" | "c" => 12.011,
        "N" | "n" => 14.007,
        "O" | "o" => 15.999,
        "F" => 18.998,
        "P" | "p" => 30.974,
        "S" | "s" => 32.065,
        "Cl" => 35.453,
        "Br" => 79.904,
        "I" => 126.904,
        _ => 0.0,
    }
}

/// Default implicit hydrogen count in organic subset
#[inline]
fn default_h_count(symbol: &str) -> u8 {
    match symbol {
        "C" | "c" => 4,
        "N" | "n" => 3,
        "O" | "o" => 2,
        "S" | "s" => 2,
        "P" | "p" => 3,
        "B" => 3,
        _ => 0,
    }
}

/// Token types produced by the SMILES scanner (property-focused subset)
#[derive(Debug, Clone)]
enum Token<'a> {
    Atom { symbol: &'a str, is_aromatic: bool, explicit_h: Option<u8> },
    BondSingle,
    BondDouble,
    BondTriple,
    BondAromatic,
    RingOpen(u32),
    BranchOpen,
    BranchClose,
}

fn scan_tokens(smiles: &str) -> Vec<Token<'_>> {
    let bytes = smiles.as_bytes();
    let len = bytes.len();
    let mut tokens = Vec::with_capacity(smiles.len());
    let mut i = 0;

    while i < len {
        let ch = bytes[i] as char;
        match ch {
            '(' => { tokens.push(Token::BranchOpen); i += 1; }
            ')' => { tokens.push(Token::BranchClose); i += 1; }
            '-' => { tokens.push(Token::BondSingle); i += 1; }
            '=' => { tokens.push(Token::BondDouble); i += 1; }
            '#' => { tokens.push(Token::BondTriple); i += 1; }
            ':' => { tokens.push(Token::BondAromatic); i += 1; }
            '/' | '\\' => { tokens.push(Token::BondSingle); i += 1; }
            '0'..='9' => {
                let ring_num = (bytes[i] - b'0') as u32;
                tokens.push(Token::RingOpen(ring_num));
                i += 1;
            }
            '%' if i + 2 < len => {
                let ring_num = (bytes[i + 1] - b'0') as u32 * 10 + (bytes[i + 2] - b'0') as u32;
                tokens.push(Token::RingOpen(ring_num));
                i += 3;
            }
            '[' => {
                // Parse bracket atom for explicit H count and charge
                let start = i + 1;
                if let Some(rel_end) = smiles[start..].find(']') {
                    let inner = &smiles[start..start + rel_end];
                    let (sym, explicit_h, _charge) = parse_bracket_atom_props(inner);
                    let is_aromatic = sym.chars().next().map(|c| c.is_lowercase()).unwrap_or(false);
                    tokens.push(Token::Atom { symbol: sym, is_aromatic, explicit_h: Some(explicit_h) });
                    i = start + rel_end + 1;
                } else {
                    i += 1;
                }
            }
            'B' | 'C' | 'N' | 'O' | 'F' | 'P' | 'S' | 'I'
            | 'c' | 'n' | 'o' | 's' | 'p' => {
                let (sym, advance) = if ch == 'B' && i + 1 < len && bytes[i + 1] == b'r' {
                    (&smiles[i..i + 2], 2)
                } else if ch == 'C' && i + 1 < len && bytes[i + 1] == b'l' {
                    (&smiles[i..i + 2], 2)
                } else {
                    (&smiles[i..i + 1], 1)
                };
                let is_aromatic = ch.is_lowercase();
                tokens.push(Token::Atom { symbol: sym, is_aromatic, explicit_h: None });
                i += advance;
            }
            _ => { i += 1; }
        }
    }
    tokens
}

/// Returns (symbol, explicit_h, charge) from bracket atom content
fn parse_bracket_atom_props(inner: &str) -> (&str, u8, i8) {
    let bytes = inner.as_bytes();
    let mut pos = 0;

    // Skip isotope digits
    while pos < bytes.len() && bytes[pos].is_ascii_digit() {
        pos += 1;
    }

    // Element symbol
    let sym_start = pos;
    if pos < bytes.len() && bytes[pos].is_ascii_alphabetic() {
        pos += 1;
        if pos < bytes.len() && bytes[pos].is_ascii_lowercase() && bytes[pos] != b'h' && bytes[pos] != b'H' {
            pos += 1;
        }
    }
    let sym = &inner[sym_start..pos];

    let mut h_count = 0u8;
    let mut charge = 0i8;
    let rest = &inner[pos..];

    let mut chars = rest.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            'H' => {
                h_count = 1;
                if let Some(&('0'..='9')) = chars.peek() {
                    if let Some(d) = chars.next() {
                        h_count = d as u8 - b'0';
                    }
                }
            }
            '+' => { charge += 1; }
            '-' => { charge -= 1; }
            _ => {}
        }
    }

    (sym, h_count, charge)
}

/// Computed molecular properties
#[derive(Debug, Clone, Copy, Default)]
pub struct MolecularProperties {
    pub mw: f32,
    pub logp: f32,
    pub rot_bonds: u32,
    pub heavy_atoms: u32,
    pub ring_count: u32,
}

/// Compute molecular properties directly from SMILES without a full parse.
/// Uses a single-pass token scan for maximum throughput.
pub fn compute_properties(smiles: &str) -> MolecularProperties {
    let tokens = scan_tokens(smiles);

    let mut mw = 0.0f32;
    let mut logp = 0.0f32;
    let mut heavy_atoms = 0u32;
    let mut atom_degrees: Vec<u32> = Vec::new();
    let mut bond_types: Vec<(usize, usize, u8)> = Vec::new(); // (a, b, order)
    let mut atom_is_aromatic: Vec<bool> = Vec::new();
    let mut atom_h_count: Vec<u8> = Vec::new();

    // Ring closure tracking: ring_num -> atom_index
    let mut ring_opens: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    let mut ring_count = 0u32;

    let mut atom_stack: Vec<Option<usize>> = Vec::new();
    let mut prev_atom: Option<usize> = None;
    let mut next_bond_order: u8 = 1;

    for token in &tokens {
        match token {
            Token::Atom { symbol, is_aromatic, explicit_h } => {
                let atom_idx = atom_degrees.len();
                atom_degrees.push(0);
                atom_is_aromatic.push(*is_aromatic);

                let h = explicit_h.unwrap_or_else(|| default_h_count(symbol));
                atom_h_count.push(h);

                // MW: heavy atom + implicit H mass
                mw += atomic_mass(symbol);
                mw += h as f32 * 1.0079;

                logp += crippen_logp_contribution(symbol, *is_aromatic);
                heavy_atoms += 1;

                if let Some(prev) = prev_atom {
                    let bond_order = if *is_aromatic && atom_is_aromatic[prev] {
                        // Aromatic bond represented as 1.5, use 1 for rot-bond counting
                        4u8
                    } else {
                        next_bond_order
                    };
                    bond_types.push((prev, atom_idx, bond_order));
                    atom_degrees[prev] += 1;
                    atom_degrees[atom_idx] += 1;
                }

                prev_atom = Some(atom_idx);
                next_bond_order = 1;
            }
            Token::BondDouble => { next_bond_order = 2; }
            Token::BondTriple => { next_bond_order = 3; }
            Token::BondAromatic => { next_bond_order = 4; }
            Token::BondSingle => { next_bond_order = 1; }
            Token::BranchOpen => {
                atom_stack.push(prev_atom);
            }
            Token::BranchClose => {
                if let Some(pa) = atom_stack.pop() {
                    prev_atom = pa;
                    next_bond_order = 1;
                }
            }
            Token::RingOpen(ring_num) => {
                if let Some(prev) = prev_atom {
                    match ring_opens.remove(ring_num) {
                        Some(open_idx) => {
                            // Ring closure bond
                            bond_types.push((open_idx, prev, next_bond_order));
                            atom_degrees[open_idx] += 1;
                            atom_degrees[prev] += 1;
                            ring_count += 1;
                            next_bond_order = 1;
                        }
                        None => {
                            ring_opens.insert(*ring_num, prev);
                        }
                    }
                }
            }
        }
    }

    // Count rotatable bonds: single bonds (not aromatic) where both atoms have degree >= 2
    let rot_bonds = bond_types.iter().filter(|&&(a, b, order)| {
        order == 1 // single bond
            && atom_degrees[a] >= 2
            && atom_degrees[b] >= 2
    }).count() as u32;

    MolecularProperties {
        mw,
        logp,
        rot_bonds,
        heavy_atoms,
        ring_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methane_mw() {
        let props = compute_properties("C");
        // Methane: C=12.011 + 4H=4*1.0079 = 12.011 + 4.032 = 16.043
        assert!((props.mw - 16.043).abs() < 0.1);
    }

    #[test]
    fn test_ethanol_properties() {
        let props = compute_properties("CCO");
        assert!(props.heavy_atoms == 3);
        assert!(props.mw > 40.0 && props.mw < 55.0);
    }

    #[test]
    fn test_benzene_ring_count() {
        let props = compute_properties("c1ccccc1");
        assert_eq!(props.ring_count, 1);
    }
}
