#![warn(missing_docs)]
//! qsar — a small, focused, zero-dependency Rust library for QSAR modeling.
//!
//! This crate provides fast, pure-Rust implementations of common QSAR building blocks:
//!
//! - **descriptors** — molecular descriptors from SMILES (no RDKit/OpenBabel needed)
//!   - Exact molecular weight
//!   - Classic physicochemical properties (LogP, TPSA, H-bond donors/acceptors, etc.)
//!   - Constitutional, topological, and fingerprint descriptors (growing!)
//! - **data_io** — lightweight CSV/JSON loaders for descriptor matrices
//! - **models** — helpers for `ndarray` ↔ `linfa` and simple training examples
//!
//! Everything is designed to be portable, fast to compile, and easy to embed in
//! larger chem-informatics pipelines. Heavy native dependencies can be added later
//! behind optional Cargo features.
//!
//! # Quick examples
//!
//! ### Physicochemical descriptors (most common in QSAR)
//! ```
//! use qsar::physchem_descriptors;
//!
//! let props = physchem_descriptors("CCO").unwrap();  // ethanol
//! println!("MolWt: {:.3}, LogP: {:.2}, TPSA: {:.1}, HBD: {}, HBA: {}",
//!          props.mol_wt, props.mol_log_p, props.tpsa,
//!          props.h_bond_donors, props.h_bond_acceptors);
//! ```
//!
//! ### Classic molecular weight (still available)
//! ```
//! use qsar::molecular_weight;
//! let mw = molecular_weight("CC(=O)Oc1ccccc1C(=O)O").unwrap();  // aspirin
//! println!("Aspirin MW: {:.3}", mw);
//! ```
//!
//! ### Load a dataset and train a model
//! ```no_run
//! use qsar::{data_io::read_csv_descriptors, models::to_ndarrays};
//! use linfa::prelude::*;
//! use linfa_linear::LinearRegression;
//!
//! let (desc, targets) = read_csv_descriptors("data.csv", &["mol_wt", "logp", "tpsa"], "pIC50")?;
//! let (x, y) = to_ndarrays(&desc, &targets);
//! let dataset = Dataset::new(x, y);
//! let model = LinearRegression::new().fit(&dataset)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod data_io;
pub mod models;

// ─────────────────────────────────────────────────────────────────────────────
// Descriptor modules (nested under a single `descriptors` module)
// ─────────────────────────────────────────────────────────────────────────────
pub mod descriptors {
    //! All molecular descriptor implementations live here.
    //!
    //! Currently available:
    //! - `physicochemical` — MolWt, LogP, TPSA, H-bond counts
    //! - `molecular_weight` — legacy exact mass function (still public at crate root)
    //! - `constitutional`, `topological`, `fingerprint` — coming soon / under development

    pub mod physicochemical;
    pub mod constitutional;
    pub mod topological;
    pub mod fingerprint;

    // Keep the old standalone function available for backward compatibility
    pub use physicochemical::molecular_weight;
}

// ─────────────────────────────────────────────────────────────────────────────
// Convenience re-exports (the nice, short names users love)
// ─────────────────────────────────────────────────────────────────────────────
pub use descriptors::physicochemical::{
    physchem_descriptors,
    PhysChemDescriptors,
};
pub use descriptors::molecular_weight; // still works exactly as before

pub use data_io::read_csv_descriptors;
pub use models::{to_ndarrays, train_and_predict_example};