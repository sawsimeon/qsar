#![warn(missing_docs)]
//! qsar — a fast, zero-dependency, pure-Rust QSAR toolbox.
//!
//! Ultra-lightweight molecular descriptor calculation and modeling utilities —
//! no RDKit, no OpenBabel, no Python. Just Rust. Ideal for embedded use, WASM,
//! serverless functions, or when you want full control and blazing compile times.
//!
//! ### Currently available descriptors
//!
//! | Category            | Functions / Types                                 | Status     |
//! |---------------------|----------------------------------------------------|------------|
//! | Physicochemical     | `physchem_descriptors`, `PhysChemDescriptors`      | Complete   |
//! | Constitutional      | `constitutional_descriptors`                       | Complete   |
//! | Fingerprints        | `maccs`, `ecfp4`, `atom_pairs`, `topological_torsion` | Complete |
//! | Molecular weight    | `molecular_weight` (legacy)                        | Complete   |
//! | Topological         | (Wiener, Zagreb, Balaban J — coming soon)          | In progress |
//!
//! ### Quick examples
//!
//! ```
//! use qsar::{
//!     physchem_descriptors,
//!     constitutional_descriptors,
//!     ecfp4,
//!     maccs,
//! };
//!
//! let smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"; // aspirin
//!
//! // Classic physicochemical + Lipinski
//! let phys = physchem_descriptors(smiles).unwrap();
//! println!("MW: {:.2}, LogP: {:.2}, TPSA: {:.1}", phys.mol_wt, phys.mol_log_p, phys.tpsa);
//! assert!(phys.lipinski_ro5());
//!
//! // Simple counts
//! let cons = constitutional_descriptors(smiles).unwrap();
//! println!("Heavy atoms: {}, Rotatable bonds: {}", cons.heavy_atom_count, cons.num_rotatable_bonds);
//!
//! // State-of-the-art fingerprint
//! let fp = ecfp4(smiles);
//! println!("ECFP4 active bits: {}", fp.len());
//!
//! // Industry-standard substructure keys
//! let maccs_fp = maccs(smiles);
//! println!("MACCS fingerprint: {} bytes", maccs_fp.len());
//! ```

pub mod data_io;
pub mod models;

// ─────────────────────────────────────────────────────────────────────────────
// All descriptor modules — grouped under `descriptors`
// ─────────────────────────────────────────────────────────────────────────────
pub mod descriptors {
    //! Core molecular descriptor implementations.
    //!
    //! All functions take a SMILES string and return pure Rust values.
    //! No native dependencies. No allocation-heavy toolkits.

    pub mod physicochemical;
    pub mod constitutional;
    pub mod topological;    // placeholder — will be filled soon
    pub mod fingerprint;

    // Legacy molecular_weight function — kept public at crate root for old users
    pub use physicochemical::molecular_weight;
}

// ─────────────────────────────────────────────────────────────────────────────
// Beautiful, ergonomic top-level re-exports
// ─────────────────────────────────────────────────────────────────────────────

// Physicochemical descriptors (most used in QSAR)
pub use descriptors::physicochemical::{
    physchem_descriptors,
    PhysChemDescriptors,
};

// Constitutional descriptors
pub use descriptors::constitutional::{
    constitutional_descriptors,
    ConstitutionalDescriptors,
};

// Fingerprints — the stars of modern ML-based QSAR
pub use descriptors::fingerprint::{
    maccs,
    ecfp4,
    atom_pairs,
    topological_torsion,
};

// Legacy molecular weight — still available at top level (no breaking change)
pub use descriptors::molecular_weight;

// Data I/O and modeling helpers
pub use data_io::read_csv_descriptors;
pub use models::{to_ndarrays, train_and_predict_example};
