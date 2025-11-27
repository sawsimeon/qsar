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
//!
//! # Quick examples
//!
//! Compute exact molecular weight (SMILES):
//!
//! ```no_run
//! use qsar::molecular_weight;
//!
//! let mw = qsar::descriptors::molecular_weight("CCO").expect("compute MW");
//! println!("Ethanol MW: {}", mw);
//! ```
//!
//! Convert vectors to `ndarray` and train a simple linear model (see `models`):
//!
//! ```no_run
//! use qsar::models::train_and_predict_example;
//! let pred = train_and_predict_example().expect("train and predict");
//! println!("Prediction: {:?}", pred);
//! ```
//!
//! Load descriptors from CSV (see `data_io`):
//!
//! ```no_run
//! use qsar::data_io::read_csv_descriptors;
//! let (x, y) = read_csv_descriptors("data/dataset.csv", &["mol_wt", "logp"], "pIC50")?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
pub mod descriptors;
pub mod data_io;
pub mod models;

// Re-exports for convenience (stable API surface)
pub use descriptors::molecular_weight;
pub use data_io::read_csv_descriptors;
pub use models::{to_ndarrays, train_and_predict_example};
