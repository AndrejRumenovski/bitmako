//! BitMako: Ultra-high-throughput chemoinformatics ETL and Block-Max WAND
//! similarity search engine for billion-scale compound libraries.
//!
//! # Architecture
//!
//! ```text
//! .bz2 file ──▶ StreamingReader ──▶ [rayon chunk workers] ──▶ LanceWriter
//!                                            │
//!                                     ECFP4 FP + props
//!                                            │
//!                                      IndexBuilder
//!                                            │
//!                                    [1024 PostingLists]
//!                                            │
//!                                       BMW Engine
//!                                     (top-k Tanimoto)
//! ```

pub mod error;
pub mod etl;
pub mod index;
pub mod search;

pub use error::{BitMakoError, Result};
