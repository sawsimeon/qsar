# qsar

[![Crates.io](https://img.shields.io/crates/v/qsar.svg)](https://crates.io/crates/qsar)
[![docs.rs](https://docs.rs/qsar/badge.svg)](https://docs.rs/qsar)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE)

qsar is a lightweight Rust library for computing common molecular descriptors and integrating descriptor data with Linfa for basic QSAR modeling. The initial release focuses on an exact molecular weight calculator implemented in pure Rust (no native Open Babel / RDKit dependencies), convenient ndarray conversion utilities, and a minimal Linfa example.

Goals
- Provide a small, easy-to-publish crate for common QSAR descriptor tasks.
- Keep the core pure Rust and dependency-light so it is easy to use across platforms.
- Provide clear extension points for optional, feature-gated integrations with native chemoinformatics toolkits.

Features
- Exact molecular weight calculation from simple SMILES or InChI strings.
- Conversion helpers to prepare descriptor matrices for Linfa.
- Minimal CSV loader helper for descriptor datasets.
- Example demonstrating training and predicting with linfa_linear::LinearRegression.

Quickstart

Add to your `Cargo.toml`:
```toml
[dependencies]
qsar = "0.1.0"
```

Compute molecular weight:
```rust
use qsar::descriptors::molecular_weight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mw = molecular_weight("CCO")?; // ethanol
    println!("Ethanol exact molecular weight: {}", mw);
    Ok(())
}
```

Convert descriptor vectors and train a linear model:
```rust
use qsar::models::{to_ndarrays, train_and_predict_example};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Example: run the built-in Linfa example
    let pred = train_and_predict_example()?;
    println!("Prediction for sample [5.0, 6.0]: {}", pred[0]);

    // Example: convert your own data
    let descriptors = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
    let targets = vec![3.0, 7.0];
    let (x, y) = to_ndarrays(descriptors, targets)?;
    println!("Features shape: {:?}", x.dim());

    Ok(())
}
```

Loading descriptors from CSV
```rust
use qsar::data_io::read_csv_descriptors;

let (descriptors, targets) = read_csv_descriptors("data/my_dataset.csv", &["mol_wt", "logp"], "pIC50")?;
```

Limitations and roadmap
- The bundled SMILES parser supports a useful subset (linear molecules, `=`, `#`, bracketed atoms with H counts). It does not yet fully support rings, branches, or stereochemistry. This is deliberate to keep the first release small and portable.
- Future work:
  - Add an optional Cargo feature to enable Open Babel / RDKit FFI bindings for full parsing.
  - Add more descriptors (logP, TPSA, rotatable bonds, fingerprinting).
  - Add cross-validation and model selection helpers (Linfa pipelines).
  - Improve error messages and parsing robustness.

Contributing
Contributions are welcome. Please open issues and pull requests on the repository. The project follows standard Rust contribution practices: format with `cargo fmt`, lint with `cargo clippy`, and run tests with `cargo test`.

Recommended local validation
- cargo fmt
- cargo clippy --all-targets -- -D warnings
- cargo test
- cargo doc --no-deps --open

License
This project is dual-licensed under MIT OR Apache-2.0. See the LICENSE file for details.
