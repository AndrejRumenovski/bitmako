//! HTTP API — wraps a `Searcher` (and optionally a Lance dataset for SMILES/property
//! resolution) in an Axum server so similarity search is network-queryable instead
//! of requiring a CLI process launch per query.
//!
//! The `Searcher` and Lance `Dataset` are loaded once at startup and shared across
//! all requests behind `Arc` — every field involved is a read-only mmap or an
//! async-safe handle, so concurrent requests need no locking.

use std::sync::Arc;

use axum::extract::State;
use axum::http::StatusCode;
use axum::response::Html;
use axum::routing::{get, post};
use axum::{Json, Router};
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::error::{BitMakoError, LanceResultExt, Result};
use crate::etl::fingerprint::compute_morgan_fp;
use crate::search::query::SimilarityQuery;
use crate::search::scaffold;
use crate::search::Searcher;

struct AppState {
    searcher: Searcher,
    lance: Option<lance::dataset::Dataset>,
}

/// The single-page search UI, embedded at compile time so the server ships as
/// one binary with no separate static-file deployment step.
const INDEX_HTML: &str = include_str!("../static/index.html");

async fn handle_index() -> Html<&'static str> {
    Html(INDEX_HTML)
}

#[derive(Deserialize)]
struct SearchRequest {
    smiles: String,
    #[serde(default = "default_threshold")]
    threshold: f32,
    #[serde(default = "default_top_k")]
    top_k: usize,
    mw_max: Option<f32>,
    logp_max: Option<f32>,
    /// Collapse same-scaffold hits down to their best-scoring representative
    /// before returning — requires `lance` to be attached (scaffolds need
    /// SMILES). See `search::scaffold::diverse_indices`.
    #[serde(default)]
    diverse: bool,
}

fn default_threshold() -> f32 {
    0.7
}
fn default_top_k() -> usize {
    50
}

#[derive(Serialize, Default)]
struct SearchResultItem {
    doc_id: u32,
    score: f32,
    #[serde(skip_serializing_if = "Option::is_none")]
    compound_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    smiles: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mw: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    logp: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    rot_bonds: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    heavy_atoms: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ring_count: Option<u32>,
    /// Similarity Analysis panel fields — always present, independent of
    /// `--lance`, since they come from the fingerprints alone. See
    /// `search::analysis`.
    shared_bits: u32,
    query_unique_bits: u32,
    candidate_unique_bits: u32,
    explanation: String,
    /// Bemis-Murcko scaffold fields — only present when `--lance` is
    /// attached, since extraction needs the candidate's SMILES text, not
    /// just its fingerprint. See `search::scaffold`.
    #[serde(skip_serializing_if = "Option::is_none")]
    scaffold_smiles: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ring_systems: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    scaffold_atoms: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    framework_fraction: Option<f32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    scaffold_key: Option<u64>,
    /// R-group decomposition — same gating as the scaffold fields above
    /// (only present with `--lance`). See `search::scaffold::decompose`.
    #[serde(skip_serializing_if = "Option::is_none")]
    r_groups: Option<Vec<RGroupItem>>,
}

#[derive(Serialize)]
struct ScaffoldGroupItem {
    scaffold_smiles: String,
    scaffold_key: u64,
    count: u32,
}

#[derive(Serialize)]
struct RGroupItem {
    attach_symbol: String,
    smiles: String,
}

#[derive(Serialize)]
struct RGroupColumnItem {
    label: String,
    attach_symbol: String,
}

#[derive(Serialize)]
struct RGroupRowItem {
    member_index: usize,
    cells: Vec<Vec<String>>,
}

#[derive(Serialize)]
struct RGroupTableItem {
    scaffold_key: u64,
    scaffold_smiles: String,
    columns: Vec<RGroupColumnItem>,
    rows: Vec<RGroupRowItem>,
}

#[derive(Serialize)]
struct SearchResponse {
    query_smiles: String,
    query_pop: u32,
    results: Vec<SearchResultItem>,
    docs_evaluated: u64,
    eval_fraction_pct: f64,
    /// Wall-clock time spent in `Searcher::search_with_stats` — the WAND
    /// search itself, not JSON serialization or (when `--lance` is attached)
    /// the SMILES/property lookup that happens afterward. Search Statistics
    /// panel field; see `handle_search`.
    search_time_ms: f64,
    /// "N results span M distinct scaffolds" grouping — empty without
    /// `--lance` (same gating as the per-result scaffold fields above).
    scaffold_groups: Vec<ScaffoldGroupItem>,
    /// SAR tables — one per scaffold shared by 2+ returned results, with
    /// substituents aligned into R1/R2/… columns by attachment position.
    /// Built from the results actually being returned (post-diversity-picking,
    /// if `diverse` was set) — empty without `--lance`.
    rgroup_tables: Vec<RGroupTableItem>,
}

#[derive(Serialize)]
struct HealthResponse {
    status: &'static str,
    compounds: u32,
    lance_attached: bool,
    prop_store_attached: bool,
    fingerprint_type: &'static str,
}

#[derive(Serialize)]
struct ErrorResponse {
    error: String,
}

type ApiError = (StatusCode, Json<ErrorResponse>);

fn bad_request(msg: impl Into<String>) -> ApiError {
    (StatusCode::BAD_REQUEST, Json(ErrorResponse { error: msg.into() }))
}

fn internal_error(e: impl std::fmt::Display) -> ApiError {
    (StatusCode::INTERNAL_SERVER_ERROR, Json(ErrorResponse { error: e.to_string() }))
}

async fn handle_health(State(state): State<Arc<AppState>>) -> Json<HealthResponse> {
    Json(HealthResponse {
        status: "ok",
        compounds: state.searcher.num_compounds(),
        lance_attached: state.lance.is_some(),
        prop_store_attached: state.searcher.has_prop_store(),
        fingerprint_type: crate::etl::fingerprint::FINGERPRINT_KIND,
    })
}

async fn handle_search(
    State(state): State<Arc<AppState>>,
    Json(req): Json<SearchRequest>,
) -> std::result::Result<Json<SearchResponse>, ApiError> {
    if !(0.0..=1.0).contains(&req.threshold) {
        return Err(bad_request("threshold must be in [0.0, 1.0]"));
    }
    if req.top_k == 0 {
        return Err(bad_request("top_k must be >= 1"));
    }

    let has_filters = req.mw_max.is_some() || req.logp_max.is_some();
    if has_filters && !state.searcher.has_prop_store() {
        return Err(bad_request(
            "mw_max/logp_max require the server to be started with --prop-store",
        ));
    }
    if req.diverse && state.lance.is_none() {
        return Err(bad_request("diverse requires the server to be started with --lance"));
    }

    let query_fp = compute_morgan_fp(&req.smiles);
    let query = SimilarityQuery::new(query_fp, req.threshold, req.top_k)
        .with_mw_logp_max(req.mw_max, req.logp_max);
    let query_pop = query.query_pop;

    // Timed strictly around the search itself (Search Statistics panel's
    // "execution time" field) — excludes JSON serialization and, when
    // `--lance` is attached, the SMILES/property lookup that happens after.
    let search_started = std::time::Instant::now();
    let (results, stats) = state
        .searcher
        .search_with_stats(&query)
        .map_err(|e| match e {
            BitMakoError::Query(msg) => bad_request(msg),
            other => internal_error(other),
        })?;
    let search_time_ms = search_started.elapsed().as_secs_f64() * 1000.0;

    let mut items = if let Some(dataset) = &state.lance {
        resolve_via_lance(dataset, &results).await.map_err(internal_error)?
    } else {
        results
            .iter()
            .map(|(doc_id, score)| SearchResultItem { doc_id: *doc_id, score: *score, ..Default::default() })
            .collect()
    };

    // Similarity Analysis: a cheap post-processing pass over the top-k results
    // already found above, not a second search — see `search::analysis`.
    for (item, analysis) in items.iter_mut().zip(state.searcher.analyze_results(&query, &results)) {
        item.shared_bits = analysis.shared_bits;
        item.query_unique_bits = analysis.query_unique_bits;
        item.candidate_unique_bits = analysis.candidate_unique_bits;
        item.explanation = analysis.explanation;
    }

    // Scaffold Analysis + R-group decomposition: only possible when `--lance`
    // is attached (both need the candidate's SMILES, not just its
    // fingerprint — see `search::scaffold`). Cheap post-processing passes
    // over the same top-k set, not a second search.
    let mut scaffold_groups = Vec::new();
    let mut rgroups: Vec<scaffold::RGroupDecomposition> = Vec::new();
    if state.lance.is_some() {
        let smiles_list: Vec<String> = items.iter().map(|it| it.smiles.clone().unwrap_or_default()).collect();
        let scaffolds = state.searcher.scaffold_results(&smiles_list);
        for (item, sc) in items.iter_mut().zip(scaffolds.iter()) {
            item.scaffold_smiles = Some(sc.scaffold_smiles.clone());
            item.ring_systems = Some(sc.ring_systems);
            item.scaffold_atoms = Some(sc.scaffold_atoms);
            item.framework_fraction = Some(sc.framework_fraction);
            item.scaffold_key = Some(sc.scaffold_key);
        }
        scaffold_groups = scaffold::group(&scaffolds)
            .into_iter()
            .map(|g| ScaffoldGroupItem { scaffold_smiles: g.scaffold_smiles, scaffold_key: g.scaffold_key, count: g.count })
            .collect();

        rgroups = state.searcher.rgroup_results(&smiles_list);
        for (item, d) in items.iter_mut().zip(rgroups.iter()) {
            item.r_groups = Some(
                d.r_groups
                    .iter()
                    .map(|r| RGroupItem { attach_symbol: r.attach_symbol.clone(), smiles: r.smiles.clone() })
                    .collect(),
            );
        }
    }

    // Diversity picking: keep only the best-scoring (first, since `items` is
    // still score-ranked) hit per distinct scaffold. Post-processes the same
    // top-k WAND already returned rather than over-fetching for a fuller
    // MaxMin-style pick, so a diverse response can come back shorter than
    // `top_k` — same tradeoff as the `search --diverse` CLI flag. `rgroups`
    // is filtered in lockstep so the SAR tables built below only ever
    // reflect results actually being returned.
    if req.diverse {
        let keys: Vec<u64> = items.iter().map(|it| it.scaffold_key.unwrap_or(0)).collect();
        let keep = scaffold::diverse_indices(&keys);
        let mut new_items = Vec::with_capacity(keep.len());
        let mut new_rgroups = Vec::with_capacity(keep.len());
        for &i in &keep {
            new_items.push(std::mem::take(&mut items[i]));
            if let Some(d) = rgroups.get(i) {
                new_rgroups.push(d.clone());
            }
        }
        items = new_items;
        rgroups = new_rgroups;
    }

    let rgroup_tables = scaffold::r_group_tables(&rgroups)
        .into_iter()
        .map(|t| RGroupTableItem {
            scaffold_key: t.scaffold_key,
            scaffold_smiles: t.scaffold_smiles,
            columns: t
                .columns
                .into_iter()
                .map(|c| RGroupColumnItem { label: c.label, attach_symbol: c.attach_symbol })
                .collect(),
            rows: t.rows.into_iter().map(|r| RGroupRowItem { member_index: r.member_index, cells: r.cells }).collect(),
        })
        .collect();

    Ok(Json(SearchResponse {
        query_smiles: req.smiles,
        query_pop,
        results: items,
        docs_evaluated: stats.docs_evaluated,
        eval_fraction_pct: stats.eval_fraction() * 100.0,
        search_time_ms,
        scaffold_groups,
        rgroup_tables,
    }))
}

async fn resolve_via_lance(
    dataset: &lance::dataset::Dataset,
    results: &[(u32, f32)],
) -> Result<Vec<SearchResultItem>> {
    let doc_ids: Vec<u32> = results.iter().map(|(doc_id, _)| *doc_id).collect();
    let resolved = crate::search::lance_lookup::resolve_compounds(dataset, &doc_ids).await?;

    Ok(results
        .iter()
        .zip(resolved)
        .map(|(&(doc_id, score), r)| SearchResultItem {
            doc_id,
            score,
            compound_id: Some(r.compound_id),
            smiles: Some(r.smiles),
            mw: Some(r.properties.mw),
            logp: Some(r.properties.logp),
            rot_bonds: Some(r.properties.rot_bonds),
            heavy_atoms: Some(r.properties.heavy_atoms),
            ring_count: Some(r.properties.ring_count),
            // Similarity Analysis fields are filled in by the caller (`handle_search`),
            // which has the query fingerprint needed to compute them.
            ..Default::default()
        })
        .collect())
}

/// Start the HTTP API on `bind:port`, serving until the process is killed.
///
/// `searcher` and the optional Lance dataset (opened from `lance_path` if given)
/// are loaded once and shared across all requests via `Arc<AppState>`.
pub async fn run_server(searcher: Searcher, lance_path: Option<String>, bind: &str, port: u16) -> Result<()> {
    let lance = match lance_path {
        Some(p) => {
            let dataset = lance::dataset::Dataset::open(&p).await.lance_err()?;
            info!("Lance dataset attached for SMILES/property resolution: {}", p);
            Some(dataset)
        }
        None => None,
    };

    let state = Arc::new(AppState { searcher, lance });

    let app = Router::new()
        .route("/", get(handle_index))
        .route("/health", get(handle_health))
        .route("/search", post(handle_search))
        .with_state(state);

    let addr = format!("{}:{}", bind, port);
    let listener = tokio::net::TcpListener::bind(&addr)
        .await
        .map_err(BitMakoError::Io)?;
    info!("BitMako HTTP API listening on http://{}", addr);
    info!("Search UI: http://{}/", addr);
    axum::serve(listener, app).await.map_err(BitMakoError::Io)?;
    Ok(())
}
