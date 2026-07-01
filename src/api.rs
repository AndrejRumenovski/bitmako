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

use crate::error::{BitMakoError, Result};
use crate::etl::fingerprint::compute_morgan_fp;
use crate::search::query::{PropertyField, PropertyFilter, SimilarityQuery};
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
}

#[derive(Serialize)]
struct SearchResponse {
    query_smiles: String,
    query_pop: u32,
    results: Vec<SearchResultItem>,
    docs_evaluated: u64,
    eval_fraction_pct: f64,
}

#[derive(Serialize)]
struct HealthResponse {
    status: &'static str,
    compounds: u32,
    lance_attached: bool,
    prop_store_attached: bool,
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

    let query_fp = compute_morgan_fp(&req.smiles);
    let mut query = SimilarityQuery::new(query_fp, req.threshold, req.top_k);
    if let Some(max) = req.mw_max {
        query = query.with_filter(PropertyFilter { field: PropertyField::MolWeight, min: None, max: Some(max) });
    }
    if let Some(max) = req.logp_max {
        query = query.with_filter(PropertyFilter { field: PropertyField::LogP, min: None, max: Some(max) });
    }
    let query_pop = query.query_pop;

    let (results, stats) = state
        .searcher
        .search_with_stats(&query)
        .map_err(|e| match e {
            BitMakoError::Query(msg) => bad_request(msg),
            other => internal_error(other),
        })?;

    let items = if let Some(dataset) = &state.lance {
        resolve_via_lance(dataset, &results).await.map_err(internal_error)?
    } else {
        results
            .iter()
            .map(|(doc_id, score)| SearchResultItem { doc_id: *doc_id, score: *score, ..Default::default() })
            .collect()
    };

    Ok(Json(SearchResponse {
        query_smiles: req.smiles,
        query_pop,
        results: items,
        docs_evaluated: stats.docs_evaluated,
        eval_fraction_pct: stats.eval_fraction() * 100.0,
    }))
}

async fn resolve_via_lance(
    dataset: &lance::dataset::Dataset,
    results: &[(u32, f32)],
) -> Result<Vec<SearchResultItem>> {
    use arrow_array::cast::AsArray;
    use arrow_array::types::{Float32Type, UInt32Type};

    if results.is_empty() {
        return Ok(Vec::new());
    }

    let row_indices: Vec<u64> = results.iter().map(|(d, _)| *d as u64).collect();
    let projection = dataset
        .schema()
        .project(&["compound_id", "smiles", "mw", "logp", "rot_bonds", "heavy_atoms", "ring_count"])
        .map_err(|e| BitMakoError::Lance(e.to_string()))?;
    let batch = dataset
        .take(&row_indices, projection)
        .await
        .map_err(|e| BitMakoError::Lance(e.to_string()))?;

    let cid_col = batch.column_by_name("compound_id").unwrap().as_string::<i32>();
    let smi_col = batch.column_by_name("smiles").unwrap().as_string::<i32>();
    let mw_col = batch.column_by_name("mw").unwrap().as_primitive::<Float32Type>();
    let logp_col = batch.column_by_name("logp").unwrap().as_primitive::<Float32Type>();
    let rot_col = batch.column_by_name("rot_bonds").unwrap().as_primitive::<UInt32Type>();
    let heavy_col = batch.column_by_name("heavy_atoms").unwrap().as_primitive::<UInt32Type>();
    let ring_col = batch.column_by_name("ring_count").unwrap().as_primitive::<UInt32Type>();

    Ok(results
        .iter()
        .enumerate()
        .map(|(i, (doc_id, score))| SearchResultItem {
            doc_id: *doc_id,
            score: *score,
            compound_id: Some(cid_col.value(i).to_string()),
            smiles: Some(smi_col.value(i).to_string()),
            mw: Some(mw_col.value(i)),
            logp: Some(logp_col.value(i)),
            rot_bonds: Some(rot_col.value(i)),
            heavy_atoms: Some(heavy_col.value(i)),
            ring_count: Some(ring_col.value(i)),
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
            let dataset = lance::dataset::Dataset::open(&p)
                .await
                .map_err(|e| BitMakoError::Lance(e.to_string()))?;
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
