//! Utilities to convert descriptor vectors into ndarray and a minimal Linfa example.
//!
//! This module contains:
//! - `to_ndarrays` — convert `Vec<Vec<f64>>` and `Vec<f64>` into `ndarray::Array2<f64>` and `ndarray::Array1<f64>`.
//! - `train_and_predict_example` — minimal example that builds a `Dataset`, trains a `linfa_linear::LinearRegression`,
//!   and predicts on a new sample.

use std::error::Error;
use ndarray::{Array1, Array2};
use linfa::prelude::*;
use linfa_linear::LinearRegression;

/// Convert descriptor and target vectors into ndarray arrays suitable for Linfa.
///
/// - `descriptors` is a Vec of samples, each sample is a Vec of features (n_samples x n_features).
/// - `targets` is a Vec of target values (length n_samples).
///
/// Returns `(features: Array2<f64>, targets: Array1<f64>)` on success or an `Err` describing the problem.
pub fn to_ndarrays(
    descriptors: Vec<Vec<f64>>,
    targets: Vec<f64>,
) -> Result<(Array2<f64>, Array1<f64>), Box<dyn Error>> {
    let n_samples = descriptors.len();
    if n_samples == 0 {
        return Err("descriptors is empty".into());
    }
    let n_features = descriptors[0].len();

    // Ensure all rows have the same length and flatten into a single Vec
    let mut flat: Vec<f64> = Vec::with_capacity(n_samples * n_features);
    for row in &descriptors {
        if row.len() != n_features {
            return Err("inconsistent feature lengths in descriptors".into());
        }
        flat.extend_from_slice(&row[..]);
    }

    // Build Array2 in row-major order: shape = (n_samples, n_features)
    let x = Array2::from_shape_vec((n_samples, n_features), flat)
        .map_err(|e| format!("failed to construct Array2: {}", e))?;

    let y = Array1::from_vec(targets);
    if y.len() != n_samples {
        return Err("targets length does not match number of descriptor rows".into());
    }

    Ok((x, y))
}

/// Minimal, complete example:
/// - converts small in-memory data to ndarrays,
/// - creates a `Dataset`,
/// - fits `LinearRegression` (default),
/// - predicts on a new sample and returns the prediction.
///
/// This function is intended as a documented example and smoke-test for integration with Linfa.
pub fn train_and_predict_example() -> Result<Array1<f64>, Box<dyn Error>> {
    // Simple synthetic dataset (y ≈ x0 + x1)
    let descriptors = vec![
        vec![1.0_f64, 2.0_f64],
        vec![2.0_f64, 3.0_f64],
        vec![3.0_f64, 4.0_f64],
        vec![4.0_f64, 5.0_f64],
    ];
    let targets = vec![3.0_f64, 5.0_f64, 7.0_f64, 9.0_f64];

    // Convert to ndarray
    let (x, y) = to_ndarrays(descriptors, targets)?;

    // Build a Linfa Dataset
    let dataset = Dataset::from((x.clone(), y.clone()));

    // Fit linear regression (ordinary least squares)
    // `fit` returns a Result<FittedLinearRegression, LinearError>, so propagate errors with `?`
    let model = LinearRegression::default().fit(&dataset)?;

    // Prepare a new sample for prediction (shape: 1 x n_features)
    let new_sample = Array2::from_shape_vec((1, x.ncols()), vec![5.0_f64, 6.0_f64])?;

    // Predict (now `model` is the fitted model)
    let prediction = model.predict(&new_sample);

    // `predict` returns an Array1 with one element (for single sample)
    Ok(prediction)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_ndarrays_and_train_example_runs() {
        let pred = train_and_predict_example().expect("train_and_predict_example failed");
        // For our synthetic data y = x0 + x1, sample [5,6] => expected 11
        let value = pred[0];
        assert!((value - 11.0).abs() < 1e-6, "unexpected prediction: {}", value);
    }

    #[test]
    fn conversion_checks_shapes() {
        let desc = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let tgt = vec![3.0, 7.0];
        let (x, y) = to_ndarrays(desc, tgt).unwrap();
        assert_eq!(x.shape(), &[2, 2]);
        assert_eq!(y.len(), 2);
    }
}
