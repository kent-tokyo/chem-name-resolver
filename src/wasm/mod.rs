#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

/// Resolve a chemical name to SMILES. Returns JS null on failure.
#[cfg(feature = "wasm")]
#[wasm_bindgen]
pub fn resolve_to_smiles(name: &str) -> Option<String> {
    crate::resolve(name).ok().map(|r| r.smiles)
}

/// Resolve a chemical name and return full result as a JSON string.
/// Returns JS null on failure.
/// Fields: smiles, canonical_name, source, molecular_formula (nullable), molecular_weight (nullable).
#[cfg(feature = "wasm")]
#[wasm_bindgen]
pub fn resolve_full(name: &str) -> Option<String> {
    let result = crate::resolve(name).ok()?;
    serde_json::to_string(&result).ok()
}

/// Normalize a chemical name (CJK-safe). Never panics.
#[cfg(feature = "wasm")]
#[wasm_bindgen]
pub fn normalize_name(name: &str) -> String {
    crate::normalizer::normalize_lowercase(name)
}
