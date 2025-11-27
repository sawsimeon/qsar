//! Molecular descriptor calculations.
//!
//! This module provides lightweight, dependency-free routines for computing
//! molecular descriptors required by qsar. The initial implementation focuses
//! on computing the exact molecular weight (monoisotopic-style atomic masses)
//! from simple SMILES or InChI inputs. The SMILES support covers common linear
//! SMILES (implicit single bonds between consecutive atoms, simple double/triple
//! bonds, bracketed atoms with explicit hydrogen counts). This is intentionally
//! small and robust for common QSAR descriptor use-cases; you can later add a
//! feature gated integration with Open Babel or other toolkits for full parsing.
use std::collections::HashMap;
use thiserror::Error;

/// Errors returned by descriptor functions.
///
/// This enum describes failures returned by descriptor calculation routines:
/// - `ParseError`: input SMILES/InChI is unsupported or malformed.
/// - `UnknownElement`: an element symbol was encountered which is not present
///   in the built-in atomic mass table.
///
/// Use these variants to surface parsing errors or missing element data to callers.
#[derive(Debug, Error)]
pub enum DescriptorError {
    /// The provided input could not be parsed by the simple parser.
    #[error("unsupported or invalid molecule string: {0}")]
    ParseError(String),

    /// An element was found for which we don't have an atomic weight.
    #[error("unknown element: {0}")]
    UnknownElement(String),
}

/// Compute the exact molecular weight (sum of atomic masses including implicit
/// hydrogens) for a molecule given as a SMILES or InChI string.
///
/// Supported input forms:
/// - SMILES (basic subset): e.g., "O", "CCO", "C=O", "C#N", bracketed atoms like "[NH4+]" with explicit H counts.
/// - InChI strings starting with "InChI=": in that case we attempt to extract
///   a formula from the InChI and compute the mass from the formula. This is
///   a light-weight fallback and does not fully parse all InChI layers.
///
/// This function is conservative: it's designed for simple, common molecules
/// used in descriptor examples. For full cheminformatics parsing consider adding
/// an optional dependency on Open Babel or RDKit (via FFI) and gating that code
/// behind a Cargo feature.
///
/// # Examples
///
/// ```
/// use qsar::descriptors::molecular_weight;
///
/// // water (SMILES "O") -> H2O
/// let mw = molecular_weight("O").expect("calculate MW");
/// assert!((mw - (2.0 * 1.00782503223 + 15.99491461956)).abs() < 1e-6);
/// ```
///
/// # Errors
///
/// Returns `DescriptorError::ParseError` for unsupported or malformed inputs,
/// or `DescriptorError::UnknownElement` if the molecule contains an element not
/// present in the built-in atomic weights table.
pub fn molecular_weight(mol: &str) -> Result<f64, DescriptorError> {
    let weights = atomic_weights();

    // If InChI, try to extract formula layer "InChI=1S/..." where after version
    // the first slash-separated layer is the formula in many InChIs.
    if mol.starts_with("InChI=") {
        // naive extraction: split on '/' and take the first element after "InChI=..."
        if let Some(rest) = mol.splitn(2, '/').nth(1) {
            // the formula is the first slash part (until next slash or end)
            let formula = rest.split('/').next().unwrap_or(rest);
            return molecular_weight_from_formula(formula, &weights);
        } else {
            return Err(DescriptorError::ParseError(mol.to_string()));
        }
    }

    // Otherwise treat as SMILES (simple parser)
    // Tokenize atoms in a linear pass; handle bracketed atoms and aromatic lower-case
    // element symbols as their uppercase equivalents (c -> C, n -> N, etc).
    #[derive(Debug)]
    struct Atom {
        element: String,
        explicit_bonds: usize,
        explicit_h: Option<usize>, // Some(n) if [NH2] provided
    }

    let mut atoms: Vec<Atom> = Vec::new();
    let chars: Vec<char> = mol.chars().collect();
    let mut i = 0usize;

    // Helper: element parsing handled inline in the main loop (supports 1- and 2-letter symbols)

    // For implicit single bonds between consecutive atoms we will add 1 bond to each neighbor.
    // We will also respect explicit bond symbols: = (double), # (triple), : (aromatic treat as single)
    // For this small parser, double counts as two bonds, triple as three.
    let mut pending_bond: Option<usize> = None; // number of bonds for next connection

    while i < chars.len() {
        let c = chars[i];

        match c {
            // bracketed atom e.g. [NH4+], [Cl-], [NH3], [C@H], [13CH4]
            '[' => {
                let mut j = i + 1;
                let mut content = String::new();
                while j < chars.len() && chars[j] != ']' {
                    content.push(chars[j]);
                    j += 1;
                }
                if j >= chars.len() {
                    return Err(DescriptorError::ParseError(mol.to_string()));
                }
                // move i to after ]
                i = j + 1;

                // From content try to extract element symbol and explicit H count if present
                // Content formats vary; we look for an element at start (uppercase + optional lowercase)
                let content_chars: Vec<char> = content.chars().collect();
                if content_chars.is_empty() {
                    return Err(DescriptorError::ParseError(mol.to_string()));
                }
                // element
                let mut el = content_chars[0].to_string();
                if content_chars.len() >= 2 && content_chars[1].is_lowercase() {
                    el.push(content_chars[1]);
                }
                // find Hn pattern like H, H2, etc
                let mut explicit_h: Option<usize> = None;
                if let Some(pos) = content.find('H') {
                    // parse digits after H
                    let digits: String = content[pos + 1..]
                        .chars()
                        .take_while(|ch| ch.is_ascii_digit())
                        .collect();
                    if digits.is_empty() {
                        explicit_h = Some(1);
                    } else if let Ok(n) = digits.parse() {
                        explicit_h = Some(n);
                    }
                }

                let mut explicit_bonds = 0usize;
                if let Some(b) = pending_bond {
                    explicit_bonds += b;
                    pending_bond = None;
                }

                atoms.push(Atom {
                    element: el,
                    explicit_bonds,
                    explicit_h,
                });
                continue;
            }
            '=' => {
                pending_bond = Some(2);
                i += 1;
                continue;
            }
            '#' => {
                pending_bond = Some(3);
                i += 1;
                continue;
            }
            ':' => {
                // aromatic, treat as single bond for our simple purposes
                pending_bond = Some(1);
                i += 1;
                continue;
            }
            '-' => {
                pending_bond = Some(1);
                i += 1;
                continue;
            }
            // branch start or ring or digits indicates complexity we don't support in full
            '(' | ')' | '%' => {
                return Err(DescriptorError::ParseError(
                    "SMILES with branches/rings not supported in simple parser".to_string(),
                ));
            }
            ch if ch.is_whitespace() => {
                i += 1;
                continue;
            }
            _ => {
                // parse element symbol (uppercase or aromatic lowercase)
                let mut el = String::new();
                if c.is_uppercase() {
                    el.push(c);
                    if i + 1 < chars.len() && chars[i + 1].is_lowercase() {
                        el.push(chars[i + 1]);
                        i += 1;
                    }
                } else if c.is_lowercase() {
                    // aromatic lowercase: c, n, o, s, p -> map to uppercase equivalents
                    let mapped = match c {
                        'c' => "C",
                        'n' => "N",
                        'o' => "O",
                        's' => "S",
                        'p' => "P",
                        'b' => "B",
                        'r' => "R", // not a real element, but keep as-is (will error)
                        _ => {
                            return Err(DescriptorError::ParseError(format!(
                                "unknown aromatic element: {}",
                                c
                            )))
                        }
                    };
                    el.push_str(mapped);
                } else {
                    return Err(DescriptorError::ParseError(format!(
                        "unexpected character in SMILES: {}",
                        c
                    )));
                }

                let mut explicit_bonds = 0usize;
                if let Some(b) = pending_bond {
                    explicit_bonds += b;
                    pending_bond = None;
                }

                atoms.push(Atom {
                    element: el,
                    explicit_bonds,
                    explicit_h: None,
                });
                i += 1;
                continue;
            }
        }
    }

    if atoms.is_empty() {
        return Err(DescriptorError::ParseError(mol.to_string()));
    }

    // Add implicit bonds between consecutive atoms (single bond unless already counted)
    for idx in 0..atoms.len() - 1 {
        // each consecutive pair implies at least a single bond
        atoms[idx].explicit_bonds += 1;
        atoms[idx + 1].explicit_bonds += 1;
    }

    // Standard valences for common elements (used to compute implicit hydrogens).
    let valences: HashMap<&str, usize> = [
        ("H", 1),
        ("C", 4),
        ("N", 3),
        ("O", 2),
        ("S", 2), // many possible valences, we use common / minimal
        ("P", 3),
        ("F", 1),
        ("Cl", 1),
        ("Br", 1),
        ("I", 1),
        ("B", 3),
        ("Si", 4),
    ]
    .iter()
    .cloned()
    .collect();

    // For each atom compute implicit hydrogens and sum masses.
    let mut total_mass = 0f64;
    for atom in atoms.iter() {
        let el = atom.element.as_str();
        let mass = weights
            .get(el)
            .ok_or_else(|| DescriptorError::UnknownElement(el.to_string()))?;

        let explicit_bonds = atom.explicit_bonds;

        let implicit_h = if let Some(explicit_h) = atom.explicit_h {
            explicit_h
        } else {
            // try to compute: hydrogens = valence - explicit_bonds
            if let Some(&val) = valences.get(el) {
                if val >= explicit_bonds {
                    val - explicit_bonds
                } else {
                    // more bonds than common valence: fallback to 0 implicit H
                    0usize
                }
            } else {
                // unknown valence: assume 0 implicit hydrogens
                0usize
            }
        };

        // add mass of atom itself
        total_mass += *mass;
        // add hydrogens mass (if H exists in weights)
        if implicit_h > 0 {
            let h_mass = *weights
                .get("H")
                .ok_or_else(|| DescriptorError::UnknownElement("H".to_string()))?;
            total_mass += (implicit_h as f64) * h_mass;
        }
    }

    Ok(total_mass)
}

/// Compute mass from a plain molecular formula like "C2H6O" or "H2O".
fn molecular_weight_from_formula(
    formula: &str,
    weights: &HashMap<&'static str, f64>,
) -> Result<f64, DescriptorError> {
    // very small parser for element symbols followed by optional count
    let mut i = 0usize;
    let chars: Vec<char> = formula.chars().collect();
    let mut total = 0f64;

    while i < chars.len() {
        let c = chars[i];
        if !c.is_uppercase() {
            return Err(DescriptorError::ParseError(formula.to_string()));
        }
        let mut el = c.to_string();
        if i + 1 < chars.len() && chars[i + 1].is_lowercase() {
            el.push(chars[i + 1]);
            i += 1;
        }
        i += 1;
        // parse number
        let mut digits = String::new();
        while i < chars.len() && chars[i].is_ascii_digit() {
            digits.push(chars[i]);
            i += 1;
        }
        let count: usize = if digits.is_empty() { 1 } else { digits.parse().unwrap_or(1) };
        let mass = weights
            .get(el.as_str())
            .ok_or_else(|| DescriptorError::UnknownElement(el.clone()))?;
        total += (*mass) * (count as f64);
    }

    Ok(total)
}

/// Atomic weights used by the simple descriptor functions.
///
/// These are monoisotopic-like exact masses for common elements used in QSAR.
/// They are listed with reasonable precision for descriptor calculations.
fn atomic_weights() -> HashMap<&'static str, f64> {
    let mut m = HashMap::new();
    m.insert("H", 1.00782503223);
    m.insert("C", 12.0);
    m.insert("N", 14.00307400443);
    m.insert("O", 15.99491461956);
    m.insert("S", 31.9720711744);
    m.insert("P", 30.97376199842);
    m.insert("F", 18.998403163);
    m.insert("Cl", 34.968852682);
    m.insert("Br", 78.9183376);
    m.insert("I", 126.90447);
    m.insert("B", 11.00930536);
    m.insert("Si", 27.976926532);
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    #[test]
    fn water_molecular_weight_smiles_o() {
        // SMILES "O" should be interpreted as water (implicit H2)
        let mw = molecular_weight("O").expect("calculate MW");
        // expected: 2*H + O
        let expected = 2.0 * 1.00782503223 + 15.99491461956;
        assert!(
            (mw - expected).abs() < EPS,
            "water MW mismatch: got {}, expected {}",
            mw,
            expected
        );
    }

    #[test]
    fn ethanol_molecular_weight_smiles_cco() {
        // SMILES "CCO" should be ethanol C2H6O
        let mw = molecular_weight("CCO").expect("calculate MW");
        let expected = 2.0 * 12.0 + 6.0 * 1.00782503223 + 15.99491461956;
        assert!(
            (mw - expected).abs() < EPS,
            "ethanol MW mismatch: got {}, expected {}",
            mw,
            expected
        );
    }

    #[test]
    fn formula_parsing_h2o() {
        let mw = molecular_weight_from_formula("H2O", &atomic_weights()).unwrap();
        let expected = 2.0 * 1.00782503223 + 15.99491461956;
        assert!((mw - expected).abs() < EPS);
    }
}
