//! CSV/JSON data I/O helpers for QSAR descriptor datasets.
//!
//! This module contains small, lightweight utilities to load descriptor matrices
//! and target vectors from CSV files. The routines are intentionally minimal and
//! return plain `Vec` types so you can convert them to `ndarray` with the helpers
//! in `src/models.rs`.
use std::error::Error;
use std::path::Path;

/// Read a CSV file and extract descriptor columns and a target column.
///
/// - `path` is the filesystem path to the CSV file.
/// - `feature_cols` is a slice of header names (in order) to use as features.
/// - `target_col` is the header name for the target variable.
///
/// Returns `(descriptors, targets)` where:
/// - `descriptors` is `Vec<Vec<f64>>` with shape (n_samples x n_features)
/// - `targets` is `Vec<f64>` with length n_samples
///
/// Errors are returned if the CSV cannot be read, required columns are missing,
/// or values fail to parse as `f64`.
///
/// Example usage:
/// ```no_run
/// use qsar::data_io::read_csv_descriptors;
/// let (x, y) = read_csv_descriptors("data/dataset.csv", &["mol_wt", "logp"], "pIC50")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_csv_descriptors<P: AsRef<Path>>(
    path: P,
    feature_cols: &[&str],
    target_col: &str,
) -> Result<(Vec<Vec<f64>>, Vec<f64>), Box<dyn Error>> {
    let mut rdr = csv::Reader::from_path(&path)?;
    let headers = rdr
        .headers()?
        .clone();

    // Map requested header names to column indices
    let mut feature_idxs: Vec<usize> = Vec::with_capacity(feature_cols.len());
    for &col in feature_cols {
        let pos = headers
            .iter()
            .position(|h| h == col)
            .ok_or_else(|| format!("feature column '{}' not found in CSV headers", col))?;
        feature_idxs.push(pos);
    }
    let target_idx = headers
        .iter()
        .position(|h| h == target_col)
        .ok_or_else(|| format!("target column '{}' not found in CSV headers", target_col))?;

    let mut descriptors: Vec<Vec<f64>> = Vec::new();
    let mut targets: Vec<f64> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        // build feature row
        let mut row: Vec<f64> = Vec::with_capacity(feature_idxs.len());
        for &idx in &feature_idxs {
            let v = record
                .get(idx)
                .ok_or_else(|| format!("missing field at index {}", idx))?;
            let parsed: f64 = v.trim().parse().map_err(|e| {
                format!(
                    "failed to parse value '{}' in column index {}: {}",
                    v, idx, e
                )
            })?;
            row.push(parsed);
        }
        // parse target
        let t_str = record
            .get(target_idx)
            .ok_or_else(|| format!("missing target at index {}", target_idx))?;
        let t: f64 = t_str.trim().parse().map_err(|e| {
            format!(
                "failed to parse target '{}' in column index {}: {}",
                t_str, target_idx, e
            )
        })?;

        descriptors.push(row);
        targets.push(t);
    }

    Ok((descriptors, targets))
}

/// Convenience: load CSV from a reader (useful for tests and in-memory data).
pub fn read_csv_descriptors_from_reader(
    reader: impl std::io::Read,
    feature_cols: &[&str],
    target_col: &str,
) -> Result<(Vec<Vec<f64>>, Vec<f64>), Box<dyn Error>> {
    let mut rdr = csv::Reader::from_reader(reader);
    let headers = rdr
        .headers()?
        .clone();

    let mut feature_idxs: Vec<usize> = Vec::with_capacity(feature_cols.len());
    for &col in feature_cols {
        let pos = headers
            .iter()
            .position(|h| h == col)
            .ok_or_else(|| format!("feature column '{}' not found in CSV headers", col))?;
        feature_idxs.push(pos);
    }
    let target_idx = headers
        .iter()
        .position(|h| h == target_col)
        .ok_or_else(|| format!("target column '{}' not found in CSV headers", target_col))?;

    let mut descriptors: Vec<Vec<f64>> = Vec::new();
    let mut targets: Vec<f64> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        let mut row: Vec<f64> = Vec::with_capacity(feature_idxs.len());
        for &idx in &feature_idxs {
            let v = record
                .get(idx)
                .ok_or_else(|| format!("missing field at index {}", idx))?;
            let parsed: f64 = v.trim().parse().map_err(|e| {
                format!(
                    "failed to parse value '{}' in column index {}: {}",
                    v, idx, e
                )
            })?;
            row.push(parsed);
        }
        let t_str = record
            .get(target_idx)
            .ok_or_else(|| format!("missing target at index {}", target_idx))?;
        let t: f64 = t_str.trim().parse().map_err(|e| {
            format!(
                "failed to parse target '{}' in column index {}: {}",
                t_str, target_idx, e
            )
        })?;

        descriptors.push(row);
        targets.push(t);
    }

    Ok((descriptors, targets))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_csv_from_reader_example() {
        let data = "mol_wt,logp,pIC50\n18.0156, -0.67, 6.0\n46.0684, -0.18, 5.0\n";
        let (x, y) = read_csv_descriptors_from_reader(data.as_bytes(), &["mol_wt", "logp"], "pIC50")
            .expect("read CSV");
        assert_eq!(x.len(), 2);
        assert_eq!(x[0].len(), 2);
        assert_eq!(y.len(), 2);
        // basic sanity checks
        assert!((x[0][0] - 18.0156).abs() < 1e-6);
        assert!((y[0] - 6.0).abs() < 1e-6);
    }
}
