use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use crate::{resolve as resolve_rs, normalize as normalize_rs};

/// Resolve a chemical name to SMILES. Returns the SMILES string.
#[pyfunction]
fn resolve_to_smiles(name: &str) -> PyResult<String> {
    resolve_rs(name)
        .map(|r| r.smiles)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Resolve a chemical name and return a dict with smiles, canonical_name, source,
/// molecular_formula and molecular_weight fields.
#[pyfunction]
fn resolve_full(py: Python<'_>, name: &str) -> PyResult<PyObject> {
    let r = resolve_rs(name).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("smiles", &r.smiles)?;
    dict.set_item("canonical_name", &r.canonical_name)?;
    dict.set_item("source", format!("{:?}", r.source))?;
    dict.set_item("molecular_formula", r.molecular_formula.as_deref().unwrap_or(""))?;
    match r.molecular_weight {
        Some(mw) => dict.set_item("molecular_weight", mw)?,
        None => dict.set_item("molecular_weight", py.None())?,
    }
    Ok(dict.into())
}

/// Normalize a chemical name (CJK-safe).
#[pyfunction]
fn normalize_name(name: &str) -> String {
    normalize_rs(name).into_owned()
}

/// Resolve a list of names to SMILES strings.
/// Returns a list where each element is the SMILES string or None if resolution failed.
#[pyfunction]
fn resolve_batch(py: Python<'_>, names: Vec<String>) -> PyObject {
    let results: Vec<PyObject> = names.iter()
        .map(|n| {
            match resolve_rs(n) {
                Ok(r) => r.smiles.into_py(py),
                Err(_) => py.None(),
            }
        })
        .collect();
    results.into_py(py)
}

/// chem_name_resolver Python module.
#[pymodule]
fn chem_name_resolver(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(resolve_to_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(resolve_full, m)?)?;
    m.add_function(wrap_pyfunction!(normalize_name, m)?)?;
    m.add_function(wrap_pyfunction!(resolve_batch, m)?)?;
    Ok(())
}
