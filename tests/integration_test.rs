use std::error::Error;

#[test]
fn integration_smiles_and_model() -> Result<(), Box<dyn Error>> {
    // Test molecular weight for water (SMILES "O")
    let mw = qsar::descriptors::molecular_weight("O")?;
    let expected_mw = 2.0 * 1.00782503223 + 15.99491461956;
    assert!(
        (mw - expected_mw).abs() < 1e-8,
        "water MW mismatch: got {}, expected {}",
        mw,
        expected_mw
    );

    // Run the Linfa train & predict example and check the prediction
    let pred = qsar::models::train_and_predict_example()?;
    assert!(
        (pred[0] - 11.0).abs() < 1e-6,
        "unexpected prediction: {} (expected ~11.0)",
        pred[0]
    );

    // Test CSV reader helper (from in-memory reader)
    let csv_data = "mol_wt,logp,pIC50\n18.0156,-0.67,6.0\n46.0684,-0.18,5.0\n";
    let (x, y) = qsar::data_io::read_csv_descriptors_from_reader(
        csv_data.as_bytes(),
        &["mol_wt", "logp"],
        "pIC50",
    )?;
    assert_eq!(x.len(), 2);
    assert_eq!(x[0].len(), 2);
    assert_eq!(y.len(), 2);
    assert!((x[0][0] - 18.0156).abs() < 1e-6);
    assert!((y[0] - 6.0).abs() < 1e-6);

    Ok(())
}
