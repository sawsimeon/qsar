#![warn(missing_docs)]
//! qsar — a small, focused Rust library for QSAR-related utilities.
//!
//! This crate provides:
//! - `descriptors` — a lightweight exact molecular-weight calculator from SMILES/InChI (small, pure-Rust parser).
//! - `data_io` — small CSV helpers to load descriptor matrices and target vectors.
//! - `models` — helpers to convert descriptor vectors into `ndarray` and a minimal Linfa example.
//!
//! The initial release intentionally avoids heavy native FFI dependencies (Open Babel / RDKit) so
//! the crate stays portable and easy to publish. Full parser/toolkit integration can be added later
//! behind optional Cargo features.
pub mod descriptors;
pub mod data_io;
pub mod models;

// Re-exports for convenience (stable API surface)
pub use descriptors::molecular_weight;
