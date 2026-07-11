//! Bemis-Murcko scaffold extraction: reduces a molecule to its ring systems
//! plus the linkers connecting them, discarding terminal substituents — the
//! chemical "core" two molecules share even when their side chains differ.
//!
//! Like `search::analysis`, this is a pure post-processing step over the
//! (small) top-k result set a search already returned, never the corpus, and
//! it needs the candidate's SMILES text (only available via `--lance`), not
//! just its fingerprint. It reuses `etl::fingerprint`'s SMILES→graph parser
//! rather than adding a third independent SMILES scanner (the repo already
//! has two, in `fingerprint.rs` and `properties.rs`, each tuned for its own
//! narrow purpose) — scaffold extraction needs the *exact* graph the
//! fingerprinter sees, not a lighter approximation.
//!
//! ## Algorithm
//! 1. **2-core**: iteratively strip every atom with degree ≤ 1. What survives
//!    is exactly the ring systems plus the atoms bridging them — this *is*
//!    the Bemis-Murcko framework, no separate ring-perception pass needed.
//! 2. **Bridge detection** (Tarjan low-link) over the framework subgraph
//!    classifies each framework edge as a ring bond or a linker bond, which
//!    gives ring/linker atom counts and the number of distinct ring systems.
//! 3. **Canonical labeling** via iterative color refinement (1-Weisfeiler-Leman),
//!    the same "hash each atom's neighborhood, repeat" idea `fingerprint.rs`'s
//!    `ecfp_iterate` already uses. Because each round's color assignment is
//!    derived purely from sorted *intrinsic* atom/neighbor features — never
//!    from the atom's original index — two differently-labeled SMILES for an
//!    isomorphic scaffold converge to the same color histogram. `scaffold_key`
//!    (used for grouping results by shared scaffold) is hashed from that
//!    histogram plus the edge-color multiset, so it's stable across
//!    relabelings even where the *display* SMILES text might not be.
//! 4. **SMILES serialization**: a DFS rooted at the lowest-ranked atom, tied
//!    off with ring-closure digits for back edges. This produces a valid,
//!    single-computation-deterministic SMILES string — not a full RDKit-grade
//!    canonicalizer (it doesn't search over automorphism choices to find a
//!    minimal string), which is why `scaffold_key`, not the string itself, is
//!    the grouping key.
//! 5. **R-group decomposition** (`decompose`): every atom the 2-core pass
//!    strips is, by construction, part of an acyclic fragment touching the
//!    scaffold at exactly one bond (a second attachment point would put that
//!    fragment's atoms on a cycle, so they'd have survived the 2-core
//!    instead) — so "one dead connected component" and "one substituent"
//!    coincide exactly, with no extra bookkeeping to tell them apart. Each
//!    substituent's attachment position is tagged with its scaffold atom's
//!    color from step 3, which is what lets `r_group_tables` align
//!    substituents from *different* molecules sharing a scaffold into
//!    consistent R1/R2/… columns — an SAR table.

use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::etl::fingerprint::{parse_smiles, Atom, Bond, Molecule};

/// Bemis-Murcko scaffold summary for a single molecule.
#[derive(Debug, Clone, PartialEq)]
pub struct ScaffoldAnalysis {
    /// Scaffold rendered as SMILES; empty if the molecule is acyclic (no ring
    /// scaffold — e.g. ethanol) or unparseable.
    pub scaffold_smiles: String,
    /// Number of distinct fused-ring clusters (e.g. biphenyl's two separate
    /// rings joined by a linker count as 2; naphthalene's fused pair counts as 1).
    pub ring_systems: u32,
    /// Total independent (SSSR-style) ring count across all ring systems.
    pub ring_count: u32,
    /// Atoms retained in the scaffold (ring atoms + linker atoms).
    pub scaffold_atoms: u32,
    /// Scaffold atoms that are not themselves part of any ring (bridges
    /// connecting separate ring systems).
    pub linker_atoms: u32,
    /// Total heavy atoms in the original (unstripped) molecule.
    pub heavy_atoms: u32,
    /// `scaffold_atoms / heavy_atoms` — how much of the molecule is core vs.
    /// side chain. 0.0 for acyclic/unparseable input.
    pub framework_fraction: f32,
    /// Stable hash of the scaffold's canonical graph structure (colors +
    /// edge multiset, not the display string) — identical for two molecules
    /// sharing the same scaffold regardless of how their SMILES was written.
    /// Molecules with no ring scaffold all share `scaffold_key == 0`.
    pub scaffold_key: u64,
}

/// One shared-scaffold bucket in a grouped result set — "N results share this core."
#[derive(Debug, Clone, PartialEq)]
pub struct ScaffoldGroup {
    pub scaffold_smiles: String,
    pub scaffold_key: u64,
    pub count: u32,
}

/// One substituent removed from a molecule to reach its Bemis-Murcko
/// scaffold, plus where it attaches. Since any ring survives the 2-core
/// reduction, every stripped connected component is acyclic and touches the
/// scaffold at exactly one bond (see module doc comment) — so "one dead
/// component" and "one R-group" coincide exactly, no ambiguity to resolve.
#[derive(Debug, Clone, PartialEq)]
pub struct RGroup {
    /// `color_refine` color of the scaffold atom this substituent attaches
    /// to. Stable across any two molecules sharing a scaffold — the same
    /// invariant `scaffold_key` relies on — so it's the key `r_group_tables`
    /// aligns substituents by column on.
    pub attach_color: u32,
    /// Element symbol of the scaffold attachment atom (lowercase if
    /// aromatic), e.g. "c", "N" — enough context to read one substituent on
    /// its own, without its scaffold alongside it.
    pub attach_symbol: String,
    /// The substituent as SMILES, rooted at a `*` dummy marking the
    /// attachment bond — e.g. `"*F"`, `"*OC(C)=O"`.
    pub smiles: String,
}

/// A molecule split into its Bemis-Murcko scaffold plus the substituents
/// stripped to reach it. `r_groups` is empty and `scaffold_key == 0` for
/// acyclic or unparseable input, mirroring `analyze`'s degrade-gracefully
/// contract.
#[derive(Debug, Clone, PartialEq)]
pub struct RGroupDecomposition {
    pub scaffold_key: u64,
    pub scaffold_smiles: String,
    pub r_groups: Vec<RGroup>,
}

/// One aligned column of an `RGroupTable`: a scaffold attachment position
/// shared by every member of the group, labeled `R1`, `R2`, … in ascending
/// `attach_color` order (an arbitrary but deterministic numbering — the
/// colors themselves aren't meaningful outside this module).
#[derive(Debug, Clone, PartialEq)]
pub struct RGroupColumn {
    pub label: String,
    pub attach_symbol: String,
}

/// One molecule's row in an `RGroupTable`: `cells[i]` holds the substituent
/// SMILES attached at `columns[i]`'s position — empty if unsubstituted there
/// (i.e. that position is a bare hydrogen), more than one entry if the
/// molecule has two distinct substituents at symmetry-equivalent positions
/// that collapsed into the same column (see `r_group_tables`).
#[derive(Debug, Clone, PartialEq)]
pub struct RGroupRow {
    /// Index into the `decomps` slice `r_group_tables` was called with —
    /// callers map this back to whichever result/rank ordering they used.
    pub member_index: usize,
    pub cells: Vec<Vec<String>>,
}

/// A scaffold-group's substituents aligned into a table: one column per
/// distinct attachment position across the group, one row per member.
#[derive(Debug, Clone, PartialEq)]
pub struct RGroupTable {
    pub scaffold_key: u64,
    pub scaffold_smiles: String,
    pub columns: Vec<RGroupColumn>,
    pub rows: Vec<RGroupRow>,
}

/// Compact framework-only graph (ring + linker atoms only, reindexed 0..k),
/// built by inducing the original molecule's 2-core. All of this module's
/// graph algorithms after stripping operate on this smaller structure.
struct SubGraph {
    atomic_num: Vec<u8>,
    is_aromatic: Vec<bool>,
    charge: Vec<i8>,
    adj: Vec<Vec<(usize, Bond)>>,
}

impl SubGraph {
    fn len(&self) -> usize {
        self.atomic_num.len()
    }
}

/// Extract the Bemis-Murcko scaffold for a single SMILES string.
///
/// Never panics: unparseable or empty SMILES, and molecules with no ring
/// system, both degrade to a neutral all-zero/empty result — the same
/// degrade-gracefully contract `compute_morgan_fp` uses.
pub fn analyze(smiles: &str) -> ScaffoldAnalysis {
    let mol = match parse_smiles(smiles) {
        Ok(m) if !m.atoms.is_empty() => m,
        _ => return empty_analysis(0),
    };
    let heavy_atoms = mol.atoms.len() as u32;

    let alive = two_core(&mol);
    if !alive.iter().any(|&a| a) {
        return empty_analysis(heavy_atoms);
    }

    let (sub, _map) = induce_subgraph(&mol, &alive);
    let n = sub.len();

    let bridges = find_bridges(&sub);
    let ring_atom: Vec<bool> = (0..n)
        .map(|i| sub.adj[i].iter().any(|&(j, _)| !bridges.contains(&edge_key(i, j))))
        .collect();
    let ring_systems = count_ring_systems(&sub, &bridges, &ring_atom);
    let linker_atoms = ring_atom.iter().filter(|&&r| !r).count() as u32;

    let edge_count: usize = sub.adj.iter().map(|a| a.len()).sum::<usize>() / 2;
    let components = count_components(&sub);
    let ring_count = (edge_count + components).saturating_sub(n) as u32;

    let colors = color_refine(&sub);
    let scaffold_key = compute_scaffold_key(&sub, &colors);
    let scaffold_smiles = write_smiles(&sub, &colors);

    ScaffoldAnalysis {
        scaffold_smiles,
        ring_systems,
        ring_count,
        scaffold_atoms: n as u32,
        linker_atoms,
        heavy_atoms,
        framework_fraction: if heavy_atoms > 0 { n as f32 / heavy_atoms as f32 } else { 0.0 },
        scaffold_key,
    }
}

/// Group a batch of already-computed analyses by shared scaffold, sorted by
/// descending frequency (ties broken alphabetically for determinism).
pub fn group(analyses: &[ScaffoldAnalysis]) -> Vec<ScaffoldGroup> {
    let mut buckets: HashMap<u64, (String, u32)> = HashMap::new();
    for a in analyses {
        let entry = buckets.entry(a.scaffold_key).or_insert_with(|| (a.scaffold_smiles.clone(), 0));
        entry.1 += 1;
    }
    let mut groups: Vec<ScaffoldGroup> = buckets
        .into_iter()
        .map(|(scaffold_key, (scaffold_smiles, count))| ScaffoldGroup { scaffold_smiles, scaffold_key, count })
        .collect();
    groups.sort_by(|a, b| b.count.cmp(&a.count).then_with(|| a.scaffold_smiles.cmp(&b.scaffold_smiles)));
    groups
}

/// Given scaffold keys for an already score-ranked result list (best first),
/// returns the indices to keep so that at most one result survives per
/// distinct scaffold — the highest-scoring one, since it's first in rank
/// order. Acyclic hits (`scaffold_key == 0`, i.e. no ring scaffold at all)
/// are never collapsed into each other: they're structurally unrelated
/// molecules that only share the *absence* of a scaffold, so deduplicating
/// them would discard real diversity rather than add it.
pub fn diverse_indices(scaffold_keys: &[u64]) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut keep = Vec::new();
    for (i, &key) in scaffold_keys.iter().enumerate() {
        if key == 0 || seen.insert(key) {
            keep.push(i);
        }
    }
    keep
}

/// Split `smiles` into its Bemis-Murcko scaffold plus the substituents
/// removed to reach it. Never panics — degrades the same way `analyze` does.
pub fn decompose(smiles: &str) -> RGroupDecomposition {
    let mol = match parse_smiles(smiles) {
        Ok(m) if !m.atoms.is_empty() => m,
        _ => return empty_decomposition(),
    };
    let alive = two_core(&mol);
    if !alive.iter().any(|&a| a) {
        return empty_decomposition();
    }

    let (sub, map) = induce_subgraph(&mol, &alive);
    let colors = color_refine(&sub);
    let scaffold_key = compute_scaffold_key(&sub, &colors);
    let scaffold_smiles = write_smiles(&sub, &colors);

    let n = mol.atoms.len();
    let mut visited_dead = vec![false; n];
    let mut r_groups = Vec::new();
    for s in 0..n {
        if !alive[s] {
            continue;
        }
        for &(t, bond) in &mol.neighbors[s] {
            if alive[t] || visited_dead[t] {
                continue;
            }
            let frag_atoms = flood_fill_fragment(&mol, t, &alive, &mut visited_dead);
            let smiles = fragment_smiles(&mol, &frag_atoms, bond);
            r_groups.push(RGroup {
                attach_color: colors[map[s]],
                attach_symbol: element_symbol(&mol.atoms[s]),
                smiles,
            });
        }
    }
    // Sorted (not deduplicated by value) so genuinely distinct substituents
    // at symmetry-equivalent positions — e.g. resorcinol's two -OH groups,
    // which land in the same `attach_color` column — both survive; collapsing
    // same-(color, SMILES) entries would silently understate substitution.
    r_groups.sort_by(|a, b| a.attach_color.cmp(&b.attach_color).then_with(|| a.smiles.cmp(&b.smiles)));
    RGroupDecomposition { scaffold_key, scaffold_smiles, r_groups }
}

fn empty_decomposition() -> RGroupDecomposition {
    RGroupDecomposition { scaffold_key: 0, scaffold_smiles: String::new(), r_groups: Vec::new() }
}

/// Collects one maximal connected component of atoms stripped by the 2-core
/// pass, starting from `start` — exactly one R-group's atoms. Marks every
/// atom it visits in `visited` so the same fragment is never built twice even
/// though it may be reachable from more than one scaffold-atom neighbor scan
/// in the caller (relevant only for the degenerate single-atom-fragment case,
/// where the fragment has no other alive neighbors to rediscover it from).
fn flood_fill_fragment(mol: &Molecule, start: usize, alive: &[bool], visited: &mut [bool]) -> Vec<usize> {
    visited[start] = true;
    let mut atoms = vec![start];
    let mut stack = vec![start];
    while let Some(u) = stack.pop() {
        for &(v, _) in &mol.neighbors[u] {
            if !alive[v] && !visited[v] {
                visited[v] = true;
                atoms.push(v);
                stack.push(v);
            }
        }
    }
    atoms
}

/// Render one substituent fragment as SMILES: a standalone tiny `SubGraph`
/// with a dummy `*` atom (atomic_num 0) at index 0, bonded via the original
/// crossing bond to the fragment's copy of `frag_atoms[0]` (the atom directly
/// bonded to the scaffold), followed by the rest of the fragment with its
/// internal bonds. Canonically colored and rooted at the dummy so the
/// attachment point is always written first — reuses the same
/// `color_refine`/`write_smiles_rooted` the scaffold writer itself uses.
fn fragment_smiles(mol: &Molecule, frag_atoms: &[usize], attach_bond: Bond) -> String {
    let mut map = HashMap::new();
    let mut atomic_num = vec![0u8]; // index 0: the dummy attachment atom.
    let mut is_aromatic = vec![false];
    let mut charge = vec![0i8];
    for &orig in frag_atoms {
        map.insert(orig, atomic_num.len());
        atomic_num.push(mol.atoms[orig].atomic_num);
        is_aromatic.push(mol.atoms[orig].is_aromatic);
        charge.push(mol.atoms[orig].charge);
    }

    let mut adj = vec![Vec::new(); atomic_num.len()];
    let root_idx = map[&frag_atoms[0]];
    adj[0].push((root_idx, attach_bond));
    adj[root_idx].push((0, attach_bond));
    for &orig in frag_atoms {
        let u = map[&orig];
        for &(v_orig, bond) in &mol.neighbors[orig] {
            if let Some(&v) = map.get(&v_orig) {
                if u < v {
                    adj[u].push((v, bond));
                    adj[v].push((u, bond));
                }
            }
        }
    }

    let frag_sub = SubGraph { atomic_num, is_aromatic, charge, adj };
    let frag_colors = color_refine(&frag_sub);
    write_smiles_rooted(&frag_sub, &frag_colors, 0)
}

fn element_symbol(atom: &Atom) -> String {
    let base = symbol_for(atom.atomic_num);
    if atom.is_aromatic { base.to_lowercase() } else { base.to_string() }
}

/// Align a batch of already-computed decompositions into SAR tables — one
/// per scaffold shared by 2+ members (scaffold-free entries and singleton
/// groups have nothing to align, so they're skipped). Tables are ordered by
/// descending membership, then scaffold SMILES, matching `group`.
pub fn r_group_tables(decomps: &[RGroupDecomposition]) -> Vec<RGroupTable> {
    let mut buckets: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, d) in decomps.iter().enumerate() {
        if d.scaffold_key != 0 {
            buckets.entry(d.scaffold_key).or_default().push(i);
        }
    }

    let mut tables: Vec<RGroupTable> = buckets
        .into_iter()
        .filter(|(_, members)| members.len() >= 2)
        .map(|(scaffold_key, members)| build_table(decomps, scaffold_key, members))
        .collect();
    tables.sort_by(|a, b| b.rows.len().cmp(&a.rows.len()).then_with(|| a.scaffold_smiles.cmp(&b.scaffold_smiles)));
    tables
}

/// Build one scaffold group's SAR table: columns are the group's distinct
/// `attach_color` values in ascending order (an arbitrary but deterministic
/// R1..Rn numbering), rows are the group's members with each substituent
/// slotted into its color's column.
fn build_table(decomps: &[RGroupDecomposition], scaffold_key: u64, members: Vec<usize>) -> RGroupTable {
    let scaffold_smiles = decomps[members[0]].scaffold_smiles.clone();

    let mut colors: Vec<u32> = members
        .iter()
        .flat_map(|&i| decomps[i].r_groups.iter().map(|r| r.attach_color))
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    colors.sort_unstable();

    // `attach_symbol` is a function of the scaffold atom's element, so it's
    // identical for every occurrence of a given color in this group — take
    // it from wherever the color first turns up.
    let columns: Vec<RGroupColumn> = colors
        .iter()
        .enumerate()
        .map(|(idx, &color)| {
            let attach_symbol = members
                .iter()
                .find_map(|&i| {
                    decomps[i].r_groups.iter().find(|r| r.attach_color == color).map(|r| r.attach_symbol.clone())
                })
                .unwrap_or_default();
            RGroupColumn { label: format!("R{}", idx + 1), attach_symbol }
        })
        .collect();

    let rows: Vec<RGroupRow> = members
        .iter()
        .map(|&member_index| {
            let mut cells = vec![Vec::new(); colors.len()];
            for r in &decomps[member_index].r_groups {
                if let Ok(col) = colors.binary_search(&r.attach_color) {
                    cells[col].push(r.smiles.clone());
                }
            }
            RGroupRow { member_index, cells }
        })
        .collect();

    RGroupTable { scaffold_key, scaffold_smiles, columns, rows }
}

fn empty_analysis(heavy_atoms: u32) -> ScaffoldAnalysis {
    ScaffoldAnalysis {
        scaffold_smiles: String::new(),
        ring_systems: 0,
        ring_count: 0,
        scaffold_atoms: 0,
        linker_atoms: 0,
        heavy_atoms,
        framework_fraction: 0.0,
        scaffold_key: 0,
    }
}

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    if a < b { (a, b) } else { (b, a) }
}

/// Iteratively strip every degree-≤1 atom. What remains is the 2-core: every
/// ring system plus the atoms on the shortest paths linking them — exactly
/// the Bemis-Murcko framework, with no separate ring-perception step needed.
fn two_core(mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut alive = vec![true; n];
    let mut degree: Vec<usize> = (0..n).map(|i| mol.neighbors[i].len()).collect();
    let mut queue: std::collections::VecDeque<usize> = (0..n).filter(|&i| degree[i] <= 1).collect();

    while let Some(u) = queue.pop_front() {
        if !alive[u] || degree[u] > 1 {
            continue;
        }
        alive[u] = false;
        for &(v, _) in &mol.neighbors[u] {
            if alive[v] {
                degree[v] -= 1;
                if degree[v] <= 1 {
                    queue.push_back(v);
                }
            }
        }
    }
    alive
}

/// Build the framework-only subgraph: reindex surviving atoms 0..k and keep
/// only edges between two surviving atoms (this is what actually drops the
/// stripped side chains from the eventual SMILES output). Also returns the
/// original-molecule-index → subgraph-index map (`usize::MAX` for stripped
/// atoms) — `decompose` needs it to look up a scaffold atom's canonical color
/// starting from an original atom index.
fn induce_subgraph(mol: &Molecule, alive: &[bool]) -> (SubGraph, Vec<usize>) {
    let n = mol.atoms.len();
    let mut map = vec![usize::MAX; n];
    let mut atomic_num = Vec::new();
    let mut is_aromatic = Vec::new();
    let mut charge = Vec::new();
    for i in 0..n {
        if alive[i] {
            map[i] = atomic_num.len();
            atomic_num.push(mol.atoms[i].atomic_num);
            is_aromatic.push(mol.atoms[i].is_aromatic);
            charge.push(mol.atoms[i].charge);
        }
    }

    let mut adj = vec![Vec::new(); atomic_num.len()];
    for i in 0..n {
        if !alive[i] {
            continue;
        }
        for &(j, bond) in &mol.neighbors[i] {
            // Each undirected edge appears once in i's list and once in j's;
            // only add it when encountered from the lower original index, so
            // it lands in the new adjacency exactly once per side.
            if alive[j] && i < j {
                let (a, b) = (map[i], map[j]);
                adj[a].push((b, bond));
                adj[b].push((a, bond));
            }
        }
    }
    (SubGraph { atomic_num, is_aromatic, charge, adj }, map)
}

/// Tarjan bridge-finding over the framework subgraph. A bridge is a linker
/// edge (its removal would disconnect the graph); every non-bridge edge is
/// part of some ring.
fn find_bridges(sub: &SubGraph) -> HashSet<(usize, usize)> {
    let n = sub.len();
    let mut disc: Vec<Option<usize>> = vec![None; n];
    let mut low: Vec<usize> = vec![0; n];
    let mut timer = 0usize;
    let mut bridges = HashSet::new();

    for start in 0..n {
        if disc[start].is_none() {
            bridge_dfs(sub, start, None, &mut disc, &mut low, &mut timer, &mut bridges);
        }
    }
    bridges
}

fn bridge_dfs(
    sub: &SubGraph,
    u: usize,
    parent: Option<usize>,
    disc: &mut [Option<usize>],
    low: &mut [usize],
    timer: &mut usize,
    bridges: &mut HashSet<(usize, usize)>,
) {
    disc[u] = Some(*timer);
    low[u] = *timer;
    *timer += 1;

    let mut parent_edge_used = false;
    for &(v, _bond) in &sub.adj[u] {
        if Some(v) == parent && !parent_edge_used {
            parent_edge_used = true;
            continue;
        }
        match disc[v] {
            None => {
                bridge_dfs(sub, v, Some(u), disc, low, timer, bridges);
                low[u] = low[u].min(low[v]);
                if low[v] > disc[u].unwrap() {
                    bridges.insert(edge_key(u, v));
                }
            }
            Some(d) => {
                low[u] = low[u].min(d);
            }
        }
    }
}

/// Connected components counted over ring atoms + ring (non-bridge) edges
/// only — this is what distinguishes "two separate rings joined by a linker"
/// (2 systems) from "one fused polycyclic system" (1 system), even though
/// both are a single connected component of the *whole* framework.
fn count_ring_systems(sub: &SubGraph, bridges: &HashSet<(usize, usize)>, ring_atom: &[bool]) -> u32 {
    let n = sub.len();
    let mut visited = vec![false; n];
    let mut count = 0u32;
    for i in 0..n {
        if ring_atom[i] && !visited[i] {
            count += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(u) = stack.pop() {
                for &(v, _) in &sub.adj[u] {
                    if ring_atom[v] && !visited[v] && !bridges.contains(&edge_key(u, v)) {
                        visited[v] = true;
                        stack.push(v);
                    }
                }
            }
        }
    }
    count
}

/// Connected components of the whole framework subgraph (rings + linkers),
/// used by the `edges - atoms + components` circuit-rank formula for
/// `ring_count` — that formula holds regardless of how many bridge/tree
/// edges are mixed in, so it needs the *whole* framework's component count,
/// not just the ring-only one `count_ring_systems` computes.
fn count_components(sub: &SubGraph) -> usize {
    let n = sub.len();
    let mut visited = vec![false; n];
    let mut count = 0;
    for i in 0..n {
        if !visited[i] {
            count += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(u) = stack.pop() {
                for &(v, _) in &sub.adj[u] {
                    if !visited[v] {
                        visited[v] = true;
                        stack.push(v);
                    }
                }
            }
        }
    }
    count
}

/// Iterative color refinement (1-Weisfeiler-Leman): starts from an intrinsic
/// per-atom invariant (element, aromaticity, charge, framework degree — never
/// the atom's original index) and repeatedly refines each atom's color from
/// its sorted `(bond order, neighbor color)` multiset, exactly the same
/// "hash your neighborhood, repeat" idea `fingerprint.rs::ecfp_iterate` uses.
/// Because every round's key is a function of structure alone, this converges
/// to the same color histogram for any two differently-labeled instances of
/// an isomorphic graph — the property `compute_scaffold_key` relies on.
fn color_refine(sub: &SubGraph) -> Vec<u32> {
    let n = sub.len();
    if n == 0 {
        return Vec::new();
    }

    let initial_keys: Vec<(u8, bool, i8, usize)> =
        (0..n).map(|i| (sub.atomic_num[i], sub.is_aromatic[i], sub.charge[i], sub.adj[i].len())).collect();
    let mut color = dense_rank(&initial_keys);
    let mut num_colors = distinct_count(&color);

    for _ in 0..n {
        let new_keys: Vec<(u32, Vec<(u8, u32)>)> = (0..n)
            .map(|i| {
                let mut nbrs: Vec<(u8, u32)> =
                    sub.adj[i].iter().map(|&(j, bond)| (bond.order(), color[j])).collect();
                nbrs.sort_unstable();
                (color[i], nbrs)
            })
            .collect();
        let new_color = dense_rank(&new_keys);
        let new_num_colors = distinct_count(&new_color);
        color = new_color;
        if new_num_colors == num_colors || new_num_colors == n {
            break;
        }
        num_colors = new_num_colors;
    }
    color
}

fn distinct_count(colors: &[u32]) -> usize {
    colors.iter().collect::<HashSet<_>>().len()
}

/// Dense-rank a slice of sortable keys into `0..distinct_count`, ties getting
/// the same rank. The rank assigned to a given key value depends only on
/// where it falls in sorted order among the keys present — never on the
/// index it came from — which is what keeps this canonical across relabelings.
fn dense_rank<T: Ord>(keys: &[T]) -> Vec<u32> {
    let n = keys.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut rank = vec![0u32; n];
    let mut next = 0u32;
    for w in 0..n {
        if w > 0 && keys[order[w]] != keys[order[w - 1]] {
            next += 1;
        }
        rank[order[w]] = next;
    }
    rank
}

/// Hash the scaffold's canonical structure — the sorted color histogram plus
/// the sorted (color, color, bond order) edge multiset — rather than the
/// display SMILES string, so it's stable across different original atom
/// orderings of an isomorphic scaffold (see module doc comment).
fn compute_scaffold_key(sub: &SubGraph, colors: &[u32]) -> u64 {
    let mut color_multiset = colors.to_vec();
    color_multiset.sort_unstable();

    let mut edges: Vec<(u32, u32, u8)> = Vec::new();
    for i in 0..sub.len() {
        for &(j, bond) in &sub.adj[i] {
            if i < j {
                let (a, b) = if colors[i] <= colors[j] { (colors[i], colors[j]) } else { (colors[j], colors[i]) };
                edges.push((a, b, bond.order()));
            }
        }
    }
    edges.sort_unstable();

    let mut hasher = DefaultHasher::new();
    color_multiset.hash(&mut hasher);
    edges.hash(&mut hasher);
    hasher.finish()
}

/// Serialize the framework subgraph to SMILES: a DFS rooted at the
/// lowest-`(color, index)` atom, visiting neighbors in the same order,
/// emitting ring-closure digits for back edges. Deterministic for a single
/// computation (same input graph in, same string out) but not a full
/// canonicalizer — see module doc comment.
fn write_smiles(sub: &SubGraph, colors: &[u32]) -> String {
    let n = sub.len();
    if n == 0 {
        return String::new();
    }
    let order_key = |i: usize| (colors[i], i);
    let root = (0..n).min_by_key(|&i| order_key(i)).unwrap();
    write_smiles_rooted(sub, colors, root)
}

/// Same DFS-to-SMILES serialization as `write_smiles`, but rooted at a caller-
/// chosen atom instead of the lowest-`(color, index)` one — used by
/// `decompose` to root a substituent fragment's SMILES at its dummy `*` atom
/// so the attachment point is always written first.
fn write_smiles_rooted(sub: &SubGraph, colors: &[u32], root: usize) -> String {
    let n = sub.len();
    if n == 0 {
        return String::new();
    }

    let mut visited = vec![false; n];
    let mut children: Vec<Vec<(usize, Bond)>> = vec![Vec::new(); n];
    // (ring-closure digit, bond type, partner atom index) per atom.
    let mut ring_digits: Vec<Vec<(u32, Bond, usize)>> = vec![Vec::new(); n];
    let mut assigned: HashMap<(usize, usize), u32> = HashMap::new();
    let mut next_digit = 1u32;

    build_tree(sub, root, None, colors, &mut visited, &mut children, &mut ring_digits, &mut assigned, &mut next_digit);

    let mut out = String::new();
    emit(sub, root, None, &children, &ring_digits, &mut out);
    out
}

/// Recursive DFS building the spanning tree (`children`) and, for every back
/// edge encountered, a ring-closure digit recorded at *both* endpoints.
///
/// A given undirected edge can be "discovered" as a back edge from either
/// endpoint's neighbor scan depending on traversal order (whichever side's
/// scan reaches it after the other end is already visited) — `assigned`
/// dedupes so each physical edge gets exactly one digit, pushed once to each
/// endpoint's list, regardless of which side notices it first.
#[allow(clippy::too_many_arguments)]
fn build_tree(
    sub: &SubGraph,
    u: usize,
    parent: Option<usize>,
    colors: &[u32],
    visited: &mut [bool],
    children: &mut [Vec<(usize, Bond)>],
    ring_digits: &mut [Vec<(u32, Bond, usize)>],
    assigned: &mut HashMap<(usize, usize), u32>,
    next_digit: &mut u32,
) {
    visited[u] = true;
    let mut nbrs = sub.adj[u].clone();
    nbrs.sort_by_key(|&(v, _)| (colors[v], v));

    let mut parent_edge_used = false;
    for (v, bond) in nbrs {
        if Some(v) == parent && !parent_edge_used {
            parent_edge_used = true;
            continue;
        }
        if !visited[v] {
            children[u].push((v, bond));
            build_tree(sub, v, Some(u), colors, visited, children, ring_digits, assigned, next_digit);
        } else {
            let key = edge_key(u, v);
            if let std::collections::hash_map::Entry::Vacant(e) = assigned.entry(key) {
                let d = *next_digit;
                *next_digit += 1;
                e.insert(d);
                ring_digits[u].push((d, bond, v));
                ring_digits[v].push((d, bond, u));
            }
        }
    }
}

/// Emit the SMILES text for the subtree rooted at `u`, given the spanning
/// tree and ring-closure digits `build_tree` already computed.
fn emit(
    sub: &SubGraph,
    u: usize,
    incoming: Option<(Bond, bool)>,
    children: &[Vec<(usize, Bond)>],
    ring_digits: &[Vec<(u32, Bond, usize)>],
    out: &mut String,
) {
    if let Some((bond, parent_aromatic)) = incoming {
        out.push_str(bond_symbol(parent_aromatic, sub.is_aromatic[u], bond));
    }
    out.push_str(&atom_token(sub, u));

    for &(digit, bond, partner) in &ring_digits[u] {
        out.push_str(bond_symbol(sub.is_aromatic[u], sub.is_aromatic[partner], bond));
        push_ring_digit(out, digit);
    }

    let kids = &children[u];
    for (idx, &(v, bond)) in kids.iter().enumerate() {
        let is_last = idx + 1 == kids.len();
        if !is_last {
            out.push('(');
        }
        emit(sub, v, Some((bond, sub.is_aromatic[u])), children, ring_digits, out);
        if !is_last {
            out.push(')');
        }
    }
}

fn push_ring_digit(out: &mut String, digit: u32) {
    if digit < 10 {
        out.push((b'0' + digit as u8) as char);
    } else {
        out.push('%');
        out.push_str(&format!("{digit:02}"));
    }
}

/// Bond symbol to emit between two atoms, mirroring `fingerprint.rs::infer_bond`'s
/// parsing convention in reverse: a bond is left implicit exactly when a
/// parser would infer it by default (single between non-aromatic-pair atoms,
/// aromatic between an aromatic-pair) — anything else needs an explicit glyph,
/// including a single bond written between two aromatic atoms (e.g. the
/// biphenyl linkage), which would otherwise be misread as aromatic.
fn bond_symbol(a_aromatic: bool, b_aromatic: bool, bond: Bond) -> &'static str {
    match bond {
        Bond::Single => if a_aromatic && b_aromatic { "-" } else { "" },
        Bond::Aromatic => if a_aromatic && b_aromatic { "" } else { ":" },
        Bond::Double => "=",
        Bond::Triple => "#",
    }
}

/// Render one atom's SMILES token. Charged atoms use bracket notation; the
/// implicit-H count inside those brackets isn't recomputed for the
/// post-stripping valence (a documented, rare-case simplification — charged
/// ring atoms are uncommon, and getting this exactly right would need
/// per-atom valence bookkeeping this module otherwise avoids).
fn atom_token(sub: &SubGraph, i: usize) -> String {
    let base = symbol_for(sub.atomic_num[i]);
    let sym = if sub.is_aromatic[i] { base.to_lowercase() } else { base.to_string() };
    if sub.charge[i] != 0 {
        let sign = if sub.charge[i] > 0 { '+' } else { '-' };
        let mag = sub.charge[i].unsigned_abs();
        if mag == 1 {
            format!("[{sym}{sign}]")
        } else {
            format!("[{sym}{sign}{mag}]")
        }
    } else {
        sym
    }
}

/// `0` is never a real element — it's `decompose`'s dummy attachment atom
/// (or, in principle, a literal SMILES `*` wildcard), so it renders as `*`.
fn symbol_for(atomic_num: u8) -> &'static str {
    match atomic_num {
        0 => "*",
        1 => "H",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        // Unreachable in practice: atomic_num is always sourced from
        // fingerprint::atomic_number, which only ever produces the values
        // above (or 0 for unrecognized symbols, filtered out upstream by
        // the 2-core step never keeping a lone unparseable atom in a ring).
        // Fall back to carbon rather than panicking on any future symbol.
        _ => "C",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benzene_is_its_own_scaffold() {
        let a = analyze("c1ccccc1");
        assert_eq!(a.scaffold_smiles, "c1ccccc1");
        assert_eq!(a.ring_systems, 1);
        assert_eq!(a.ring_count, 1);
        assert_eq!(a.scaffold_atoms, 6);
        assert_eq!(a.linker_atoms, 0);
        assert_eq!(a.heavy_atoms, 6);
        assert!((a.framework_fraction - 1.0).abs() < 1e-6);
    }

    #[test]
    fn acyclic_molecule_has_no_scaffold() {
        let a = analyze("CCO");
        assert_eq!(a.scaffold_smiles, "");
        assert_eq!(a.ring_systems, 0);
        assert_eq!(a.ring_count, 0);
        assert_eq!(a.scaffold_atoms, 0);
        assert_eq!(a.heavy_atoms, 3);
        assert_eq!(a.scaffold_key, 0);
    }

    #[test]
    fn empty_and_unparseable_smiles_do_not_panic() {
        for s in ["", "[", "C(", "***"] {
            let a = analyze(s);
            assert_eq!(a.scaffold_smiles, "");
            assert_eq!(a.scaffold_key, 0);
        }
    }

    #[test]
    fn naphthalene_is_one_fused_ring_system() {
        // Two fused aromatic rings sharing an edge: 10 atoms, 2 independent rings,
        // but one connected ring *system* (not two separate ones).
        let a = analyze("c1ccc2ccccc2c1");
        assert_eq!(a.scaffold_atoms, 10);
        assert_eq!(a.ring_systems, 1);
        assert_eq!(a.ring_count, 2);
        assert_eq!(a.linker_atoms, 0);
    }

    #[test]
    fn aspirin_scaffold_is_the_bare_benzene_ring() {
        // Every substituent (acetyl ester, carboxylic acid) is an acyclic chain
        // hanging off the ring — all of it should strip away in the 2-core pass.
        let a = analyze("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(a.scaffold_smiles, "c1ccccc1");
        assert_eq!(a.ring_systems, 1);
        assert_eq!(a.scaffold_atoms, 6);
        assert_eq!(a.heavy_atoms, 13);
        assert!((a.framework_fraction - 6.0 / 13.0).abs() < 1e-5);
    }

    #[test]
    fn two_rings_joined_by_a_linker_are_two_systems() {
        // p-tolyl-CH2-phenyl: a methyl side chain (stripped), two aromatic
        // rings, and a one-atom CH2 linker bridging them (kept, not a ring atom).
        let a = analyze("Cc1ccc(Cc2ccccc2)cc1");
        assert_eq!(a.scaffold_atoms, 13); // 6 + 1 (linker) + 6
        assert_eq!(a.ring_systems, 2);
        assert_eq!(a.ring_count, 2);
        assert_eq!(a.linker_atoms, 1);
        assert_eq!(a.heavy_atoms, 14);
    }

    #[test]
    fn scaffold_of_a_scaffold_is_itself() {
        // The writer's output must be valid, re-parseable SMILES whose own
        // scaffold is unchanged — otherwise the displayed string would be lying
        // about what it represents.
        let first = analyze("CC(=O)Oc1ccccc1C(=O)O");
        let second = analyze(&first.scaffold_smiles);
        assert_eq!(first.scaffold_smiles, second.scaffold_smiles);
        assert_eq!(first.scaffold_key, second.scaffold_key);
        assert_eq!(second.scaffold_atoms, second.heavy_atoms);
    }

    #[test]
    fn unrelated_molecules_sharing_a_ring_share_a_scaffold_key() {
        // The whole point of grouping: two chemically different molecules
        // (aspirin, toluene) whose *only* ring is a bare, unsubstituted
        // benzene must land in the same scaffold bucket as plain benzene
        // itself — scaffold_key is a graph invariant, not a text hash, so
        // this holds regardless of how differently each SMILES was written.
        let benzene = analyze("c1ccccc1");
        let aspirin = analyze("CC(=O)Oc1ccccc1C(=O)O");
        let toluene = analyze("Cc1ccccc1");
        assert_eq!(benzene.scaffold_key, aspirin.scaffold_key);
        assert_eq!(benzene.scaffold_key, toluene.scaffold_key);
        assert_ne!(benzene.scaffold_key, 0);
    }

    #[test]
    fn different_scaffolds_get_different_keys() {
        let benzene = analyze("c1ccccc1");
        let naphthalene = analyze("c1ccc2ccccc2c1");
        assert_ne!(benzene.scaffold_key, naphthalene.scaffold_key);
    }

    #[test]
    fn group_buckets_by_shared_scaffold_and_sorts_by_frequency() {
        let analyses: Vec<ScaffoldAnalysis> = [
            "c1ccccc1",              // benzene ring alone
            "Cc1ccccc1",             // toluene -> same ring
            "CC(=O)Oc1ccccc1C(=O)O", // aspirin -> same ring
            "c1ccc2ccccc2c1",        // naphthalene -> different scaffold
            "CCO",                   // acyclic -> no scaffold
        ]
        .iter()
        .map(|s| analyze(s))
        .collect();

        let groups = group(&analyses);
        assert_eq!(groups.len(), 3); // benzene-ring bucket, naphthalene, acyclic
        assert_eq!(groups[0].scaffold_smiles, "c1ccccc1");
        assert_eq!(groups[0].count, 3);
    }

    #[test]
    fn charged_ring_atom_does_not_panic() {
        // Pyridinium-style charged aromatic nitrogen in a ring — exercises the
        // bracket-atom write path without asserting exact valence fidelity
        // (documented simplification, see `atom_token`).
        let a = analyze("Cc1cc[n+](C)cc1");
        assert!(!a.scaffold_smiles.is_empty());
        assert_eq!(a.ring_systems, 1);
        // Must itself re-parse without panicking.
        let _ = analyze(&a.scaffold_smiles);
    }

    #[test]
    fn diverse_indices_keeps_first_hit_per_scaffold() {
        // Rank order is score order: benzene-ring hit at 0 outranks the
        // later toluene/aspirin hits sharing the same scaffold_key.
        let keys = [10, 20, 10, 30, 20, 10];
        assert_eq!(diverse_indices(&keys), vec![0, 1, 3]);
    }

    #[test]
    fn diverse_indices_never_collapses_scaffold_free_hits() {
        // Every zero-key (acyclic) hit is its own scaffold, so all survive,
        // even though a real (non-zero) repeated key in the same list still
        // gets deduplicated down to its first occurrence.
        let keys = [0, 5, 0, 0, 5];
        assert_eq!(diverse_indices(&keys), vec![0, 1, 2, 3]);
    }

    #[test]
    fn diverse_indices_on_real_molecules_preserves_rank_and_dedups() {
        let analyses: Vec<ScaffoldAnalysis> = [
            "CC(=O)Oc1ccccc1C(=O)O", // aspirin -> benzene scaffold (best-ranked)
            "Cc1ccccc1",             // toluene -> same benzene scaffold, ranked lower
            "c1ccc2ccccc2c1",        // naphthalene -> distinct scaffold
            "CCO",                   // acyclic -> always kept
        ]
        .iter()
        .map(|s| analyze(s))
        .collect();
        let keys: Vec<u64> = analyses.iter().map(|a| a.scaffold_key).collect();
        assert_eq!(diverse_indices(&keys), vec![0, 2, 3]);
    }

    #[test]
    fn toluene_decomposes_to_one_methyl_r_group() {
        let d = decompose("Cc1ccccc1");
        assert_eq!(d.scaffold_smiles, "c1ccccc1");
        assert_eq!(d.r_groups.len(), 1);
        assert_eq!(d.r_groups[0].attach_symbol, "c");
        assert_eq!(d.r_groups[0].smiles, "*C");
    }

    #[test]
    fn aspirin_decomposes_to_two_r_groups_on_the_ring() {
        // Acetyl ester + carboxylic acid, both stripped down to the bare
        // benzene scaffold — matches `aspirin_scaffold_is_the_bare_benzene_ring`.
        let d = decompose("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(d.scaffold_smiles, "c1ccccc1");
        assert_eq!(d.r_groups.len(), 2);
        for r in &d.r_groups {
            assert_eq!(r.attach_symbol, "c");
            assert!(r.smiles.starts_with('*'));
        }
        let mut smiles: Vec<&str> = d.r_groups.iter().map(|r| r.smiles.as_str()).collect();
        smiles.sort_unstable();
        assert_eq!(smiles, vec!["*C(O)=O", "*OC(C)=O"]);
    }

    #[test]
    fn acyclic_and_unparseable_smiles_decompose_to_nothing() {
        for s in ["CCO", "", "[", "***"] {
            let d = decompose(s);
            assert_eq!(d.scaffold_smiles, "");
            assert_eq!(d.scaffold_key, 0);
            assert!(d.r_groups.is_empty());
        }
    }

    #[test]
    fn decompose_is_deterministic() {
        let a = decompose("CC(=O)Oc1ccccc1C(=O)O");
        let b = decompose("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(a, b);
    }

    #[test]
    fn single_substituent_molecules_agree_with_scaffold_analysis() {
        // Every r_group's attach_color must be a real color assigned to some
        // scaffold atom -- cross-checked against `analyze`'s independently
        // computed scaffold_key for the same molecule.
        let smiles = "CC(=O)Oc1ccccc1C(=O)O";
        let a = analyze(smiles);
        let d = decompose(smiles);
        assert_eq!(a.scaffold_key, d.scaffold_key);
        assert_eq!(a.scaffold_smiles, d.scaffold_smiles);
    }

    #[test]
    fn r_group_tables_aligns_positional_isomers_into_separate_columns() {
        // 2- and 3-methylpyridine share a bare pyridine scaffold, but the
        // ring position each methyl attaches to has a different
        // `color_refine` color -- pyridine's N breaks the ring's full
        // symmetry -- so the table gets two columns, each molecule filling
        // one and leaving the other as H (an empty cell).
        let decomps: Vec<RGroupDecomposition> = ["Cc1ccccn1", "Cc1cccnc1"].iter().map(|s| decompose(s)).collect();
        let tables = r_group_tables(&decomps);
        assert_eq!(tables.len(), 1);
        let t = &tables[0];
        assert_eq!(t.columns.len(), 2);
        assert_eq!(t.rows.len(), 2);
        for row in &t.rows {
            let filled: Vec<usize> = (0..row.cells.len()).filter(|&c| !row.cells[c].is_empty()).collect();
            assert_eq!(filled.len(), 1, "expected exactly one filled column per row");
        }
        let col_of = |row: &RGroupRow| (0..row.cells.len()).find(|&c| !row.cells[c].is_empty()).unwrap();
        assert_ne!(col_of(&t.rows[0]), col_of(&t.rows[1]));
    }

    #[test]
    fn r_group_tables_collapses_symmetric_positions_into_one_column() {
        // Benzene's ring positions are all graph-automorphic, so toluene's
        // methyl and benzoic acid's carboxyl -- attached to different
        // specific ring atoms in each molecule -- land on the same color,
        // hence the same single column.
        let decomps: Vec<RGroupDecomposition> =
            ["Cc1ccccc1", "c1ccccc1C(=O)O"].iter().map(|s| decompose(s)).collect();
        let tables = r_group_tables(&decomps);
        assert_eq!(tables.len(), 1);
        let t = &tables[0];
        assert_eq!(t.columns.len(), 1);
        assert_eq!(t.columns[0].label, "R1");
        assert_eq!(t.rows.len(), 2);
        assert_eq!(t.rows[0].cells[0], vec!["*C".to_string()]);
        assert_eq!(t.rows[1].cells[0], vec!["*C(O)=O".to_string()]);
    }

    #[test]
    fn r_group_tables_skips_scaffold_free_and_singleton_groups() {
        // "CCO" has no scaffold at all; naphthalene is the only member of
        // its own scaffold group -- neither has anything to align against.
        let decomps: Vec<RGroupDecomposition> =
            ["CCO", "Cc1ccccc1", "c1ccc2ccccc2c1"].iter().map(|s| decompose(s)).collect();
        assert!(r_group_tables(&decomps).is_empty());
    }

    #[test]
    fn r_group_tables_member_index_maps_back_to_input_order() {
        let decomps: Vec<RGroupDecomposition> =
            ["Cc1ccccc1", "c1ccccc1C(=O)O"].iter().map(|s| decompose(s)).collect();
        let tables = r_group_tables(&decomps);
        let mut indices: Vec<usize> = tables[0].rows.iter().map(|r| r.member_index).collect();
        indices.sort_unstable();
        assert_eq!(indices, vec![0, 1]);
    }
}
