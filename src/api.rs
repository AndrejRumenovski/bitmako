//! HTTP API — wraps a `Searcher` (and optionally a Lance dataset for SMILES/property
//! resolution) in an Axum server so similarity search is network-queryable instead
//! of requiring a CLI process launch per query.
//!
//! The `Searcher` and Lance `Dataset` are loaded once at startup and shared across
//! all requests behind `Arc` — every field involved is a read-only mmap or an
//! async-safe handle, so concurrent requests need no locking.

use std::collections::HashMap;
use std::net::{IpAddr, SocketAddr};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use axum::extract::{ConnectInfo, State};
use axum::http::StatusCode;
use axum::middleware::{self, Next};
use axum::response::{Html, IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use serde::{Deserialize, Serialize};
use tracing::info;

use crate::error::{BitMakoError, LanceResultExt, Result};
use crate::etl::fingerprint::compute_morgan_fp;
use crate::search::query::SimilarityQuery;
use crate::search::Searcher;

/// Hard ceiling on `top_k`, independent of whatever the UI's `<input max>`
/// suggests — a public-facing deployment can't trust the client to respect
/// that hint, since the request body is just JSON anyone can POST directly.
const MAX_TOP_K: usize = 500;

/// Requests allowed per IP per rolling window on `/search`. Generous enough
/// for a person clicking around the demo UI, tight enough that a script
/// can't hammer a public instance into a compute bill.
const RATE_LIMIT_MAX_REQUESTS: u32 = 30;
const RATE_LIMIT_WINDOW: Duration = Duration::from_secs(60);

struct AppState {
    searcher: Searcher,
    lance: Option<lance::dataset::Dataset>,
    /// Optional banner text shown in the UI (e.g. "demo running a 50M-compound
    /// subset of the full 1.36B corpus") — `None` for a normal/local deployment.
    demo_notice: Option<String>,
    limiter: RateLimiter,
}

/// Minimal fixed-window, per-IP request limiter. Deliberately not a crate
/// dependency: a `Mutex<HashMap>` is all a single-process demo server needs,
/// and it keeps `Cargo.toml` unchanged for what is otherwise a one-file concern.
struct RateLimiter {
    windows: Mutex<HashMap<IpAddr, (Instant, u32)>>,
}

impl RateLimiter {
    fn new() -> Self {
        RateLimiter { windows: Mutex::new(HashMap::new()) }
    }

    /// Returns `true` if `ip` is still within its rate limit for the current
    /// window (and records this call), `false` if it must be rejected.
    fn check(&self, ip: IpAddr) -> bool {
        let mut windows = self.windows.lock().unwrap();
        let now = Instant::now();
        let entry = windows.entry(ip).or_insert((now, 0));
        if now.duration_since(entry.0) > RATE_LIMIT_WINDOW {
            *entry = (now, 0);
        }
        entry.1 += 1;
        entry.1 <= RATE_LIMIT_MAX_REQUESTS
    }
}

async fn rate_limit_middleware(
    ConnectInfo(addr): ConnectInfo<SocketAddr>,
    State(state): State<Arc<AppState>>,
    request: axum::extract::Request,
    next: Next,
) -> Response {
    if state.limiter.check(addr.ip()) {
        next.run(request).await
    } else {
        (
            StatusCode::TOO_MANY_REQUESTS,
            Json(ErrorResponse {
                error: format!(
                    "rate limit exceeded: max {} requests per {}s per IP",
                    RATE_LIMIT_MAX_REQUESTS,
                    RATE_LIMIT_WINDOW.as_secs()
                ),
            }),
        )
            .into_response()
    }
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
    /// Similarity Analysis panel fields — always present, independent of
    /// `--lance`, since they come from the fingerprints alone. See
    /// `search::analysis`.
    shared_bits: u32,
    query_unique_bits: u32,
    candidate_unique_bits: u32,
    explanation: String,
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
}

#[derive(Serialize)]
struct HealthResponse {
    status: &'static str,
    compounds: u32,
    lance_attached: bool,
    prop_store_attached: bool,
    fingerprint_type: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    demo_notice: Option<String>,
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
        demo_notice: state.demo_notice.clone(),
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
    if req.top_k > MAX_TOP_K {
        return Err(bad_request(format!("top_k must be <= {MAX_TOP_K}")));
    }

    let has_filters = req.mw_max.is_some() || req.logp_max.is_some();
    if has_filters && !state.searcher.has_prop_store() {
        return Err(bad_request(
            "mw_max/logp_max require the server to be started with --prop-store",
        ));
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

    Ok(Json(SearchResponse {
        query_smiles: req.smiles,
        query_pop,
        results: items,
        docs_evaluated: stats.docs_evaluated,
        eval_fraction_pct: stats.eval_fraction() * 100.0,
        search_time_ms,
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
/// are loaded once and shared across all requests via `Arc<AppState>`. `demo_notice`,
/// when set, is surfaced via `/health` and shown as a banner in the UI — intended
/// for a public demo instance running a subset of the full corpus (see
/// `bitmako extract-subset`), so visitors aren't misled about scale.
pub async fn run_server(
    searcher: Searcher,
    lance_path: Option<String>,
    bind: &str,
    port: u16,
    demo_notice: Option<String>,
) -> Result<()> {
    let lance = match lance_path {
        Some(p) => {
            let dataset = lance::dataset::Dataset::open(&p).await.lance_err()?;
            info!("Lance dataset attached for SMILES/property resolution: {}", p);
            Some(dataset)
        }
        None => None,
    };

    let state = Arc::new(AppState { searcher, lance, demo_notice, limiter: RateLimiter::new() });

    // Rate limiting only applies to /search (the expensive route) — /health and
    // the static UI page are cheap and harmless to hit freely.
    let search_route = axum::Router::new()
        .route("/search", post(handle_search))
        .route_layer(middleware::from_fn_with_state(state.clone(), rate_limit_middleware));

    let app = Router::new()
        .route("/", get(handle_index))
        .route("/health", get(handle_health))
        .merge(search_route)
        .with_state(state);

    let addr = format!("{}:{}", bind, port);
    let listener = tokio::net::TcpListener::bind(&addr)
        .await
        .map_err(BitMakoError::Io)?;
    info!("BitMako HTTP API listening on http://{}", addr);
    info!("Search UI: http://{}/", addr);
    axum::serve(
        listener,
        app.into_make_service_with_connect_info::<SocketAddr>(),
    )
    .await
    .map_err(BitMakoError::Io)?;
    Ok(())
}
