#![warn(missing_docs)]
//! qsar — a fast, zero-dependency, pure-Rust QSAR toolbox.
//!
//! The most complete, well-documented, and blazing-fast collection of molecular
//! descriptors and modeling utilities written entirely in Rust — no Python, no RDKit,
//! no OpenBabel, no native dependencies. Perfect for:
//! - WebAssembly (WASM)
//! - Embedded systems
//! - Serverless functions
//! - High-performance pipelines
//! - Reproducible research
//!
//! ### Available descriptors (all complete!)
//!
//! | Category            | Key Functions / Types                                         | Status     |
//! |---------------------|----------------------------------------------------------------|------------|
//! | Physicochemical     | `physchem_descriptors`, `PhysChemDescriptors`                  | Complete   |
//! | Constitutional      | `constitutional_descriptors`, `ConstitutionalDescriptors`      | Complete   |
//! | Topological         | `topological_descriptors`, `TopologicalDescriptors`            | Complete   |
//! | Fingerprints        | `maccs`, `ecfp4`, `atom_pairs`, `topological_torsion`          | Complete   |
//! | Legacy              | `molecular_weight`                                             | Complete   |
//!
//! ### One-liner examples (the API users dream of)
//!
//! ```
//! use qsar::{
//!     physchem_descriptors,
//!     constitutional_descriptors,
//!     topological_descriptors,
//!     ecfp4,
//!     maccs,
//! };
//!
//! let smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"; // aspirin
//!
//! let phys = physchem_descriptors(smiles).unwrap();
//! let cons = constitutional_descriptors(smiles).unwrap();
//! let topo = topological_descriptors(smiles).unwrap();
//!
//! println!("Aspirin QSAR Profile:");
//! println!("  MW      = {:.2} Da", phys.mol_wt);
//! println!("  LogP    = {:.2}", phys.mol_log_p);
//! println!("  TPSA    = {:.1} Å²", phys.tpsa);
//! println!("  HBA/HBD = {}/{}", phys.h_bond_acceptors, phys.h_bond_donors);
//! println!("  Heavy atoms = {}", cons.heavy_atom_count);
//! println!("  Rotatable bonds = {}", cons.num_rotatable_bonds);
//! println!("  Aromatic rings = {}", cons.num_aromatic_rings);
//! println!("  Wiener index = {}", topo.wiener);
//! println!("  Randić χ = {:.4}", topo.randic);
//! println!("  ECFP4 bits = {}", ecfp4(smiles).len());
//! println!("  MACCS set = {} bits", maccs(smiles).iter().map(|b| b.count_ones()).sum::<u32>());
//! ```

pub mod data_io;
pub mod models;

// ─────────────────────────────────────────────────────────────────────────────
// All descriptor modules — cleanly grouped
// ─────────────────────────────────────────────────────────────────────────────
pub mod descriptors {
    //! Pure-Rust molecular descriptor implementations.
    //!
    //! No external dependencies. No FFI. No slow parsing.
    //! Just fast, correct, and beautiful code.

    pub mod physicochemical;
    pub mod constitutional;
    pub mod topological;
    pub mod fingerprint;

    /// Legacy standalone molecular weight function — kept for backward compatibility
    pub use physicochemical::molecular_weight;
}

// ─────────────────────────────────────────────────────────────────────────────
// Ergonomic top-level re-exports — the API your users will love
// ─────────────────────────────────────────────────────────────────────────────

// === Classic physicochemical (most used in QSAR) ===
pub use descriptors::physicochemical::{
    physchem_descriptors,
    PhysChemDescriptors,
};

// === Constitutional counts ===
pub use descriptors::constitutional::{
    constitutional_descriptors,
    ConstitutionalDescriptors,
};

// === Topological indices (Wiener, Zagreb, Randić, Balaban J, Kappa) ===
pub use descriptors::topological::{
    topological_descriptors,
    TopologicalDescriptors,
};

// === Modern fingerprints (gold standards in ML) ===
pub use descriptors::fingerprint::{
    maccs,
    ecfp4,
    atom_pairs,
    topological_torsion,
};

// === Legacy molecular weight (still available at crate root) ===
pub use descriptors::molecular_weight;

// === Data I/O and modeling ===
pub use data_io::read_csv_descriptors;
pub use models::{to_ndarrays, train_and_predict_example};
