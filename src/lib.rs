pub mod dict;
pub mod normalizer;
pub mod parser;
mod error;

#[cfg(feature = "wasm")]
mod wasm;

#[cfg(feature = "python")]
mod py;

pub use error::ResolveError;
pub use parser::MolGraph;

use std::borrow::Cow;

#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub enum ResolveSource {
    Dictionary,
    Parser,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,
    /// Molecular formula in Hill notation (e.g. "C2H6O").
    /// None when resolved via DirectSmiles (no MolGraph available).
    pub molecular_formula: Option<String>,
    /// Molecular weight in g/mol.
    /// None when resolved via DirectSmiles (no MolGraph available).
    pub molecular_weight: Option<f64>,
}

/// Resolve any chemical name to SMILES.
///
/// Pipeline: normalize → dictionary → IUPAC parser
pub fn resolve(name: &str) -> Result<ResolveResult, ResolveError> {
    let normalized = normalizer::normalize_lowercase(name);

    if let Some(entry) = dict::lookup_synonym(&normalized) {
        match entry {
            dict::DictEntry::CanonicalName(canonical) => {
                let graph = parser::parse_iupac(canonical)?;
                let smiles = parser::smiles::to_smiles(&graph);
                let molecular_formula = Some(parser::formula::molecular_formula(&graph));
                let molecular_weight = Some(parser::formula::molecular_weight(&graph));
                return Ok(ResolveResult {
                    smiles,
                    canonical_name: canonical.to_string(),
                    source: ResolveSource::Dictionary,
                    molecular_formula,
                    molecular_weight,
                });
            }
            dict::DictEntry::DirectSmiles(smiles) => {
                return Ok(ResolveResult {
                    smiles: smiles.to_string(),
                    canonical_name: normalized.clone(),
                    source: ResolveSource::Dictionary,
                    molecular_formula: None,
                    molecular_weight: None,
                });
            }
        }
    }

    let graph = parser::parse_iupac(&normalized)?;
    let smiles = parser::smiles::to_smiles(&graph);
    let molecular_formula = Some(parser::formula::molecular_formula(&graph));
    let molecular_weight = Some(parser::formula::molecular_weight(&graph));
    Ok(ResolveResult {
        smiles,
        canonical_name: normalized,
        source: ResolveSource::Parser,
        molecular_formula,
        molecular_weight,
    })
}

/// Resolve a batch of chemical names, returning one result per input.
pub fn resolve_batch(names: &[&str]) -> Vec<Result<ResolveResult, ResolveError>> {
    names.iter().map(|&n| resolve(n)).collect()
}

/// Normalize a chemical name (CJK-safe, zero-copy when already normalized).
pub fn normalize(name: &str) -> Cow<'_, str> {
    normalizer::normalize(name)
}

/// Dictionary lookup only (no parsing).
pub fn lookup(name: &str) -> Option<dict::DictEntry<'static>> {
    dict::lookup_synonym(name)
}
