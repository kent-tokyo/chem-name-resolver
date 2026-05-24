//! # chem-name-resolver
//!
//! Pure-Rust library for resolving chemical names → SMILES, IUPAC names, InChI identifiers,
//! and molecular properties. Supports CJK (Japanese Katakana, Chinese Hanzi) synonyms,
//! trivial names, and systematic IUPAC nomenclature. WASM-compatible (no C/C++ deps).
//!
//! ## Quick start
//!
//! ```rust
//! use chem_name_resolver::resolve;
//!
//! let r = resolve("ethanol").unwrap();
//! assert_eq!(r.smiles, "CCO");
//! assert_eq!(r.canonical_name, "ethanol");
//! assert!(r.confidence >= 0.9);
//! ```
//!
//! ## Resolution pipeline
//!
//! 1. **Normalize** — CJK full-width, katakana prolonged mark, Greek letters, whitespace
//! 2. **Dictionary lookup** — Perfect-hash map of ~200 common synonyms (incl. CJK)
//! 3. **IUPAC parser** — Recursive-descent parser for systematic names (C1-C20, basic func. groups)

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

/// Source that produced a [`ResolveResult`].
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub enum ResolveSource {
    /// Name was found in the built-in synonym dictionary.
    Dictionary,
    /// Name was parsed by the IUPAC systematic-name parser.
    Parser,
}

/// The result of resolving a chemical name.
///
/// All fields that require a molecular graph (formula, weight, InChI) are `None`
/// when the name was resolved via a `DirectSmiles` dictionary entry, because no
/// graph is built in that path.
#[derive(Debug, Clone, serde::Serialize)]
pub struct ResolveResult {
    /// Canonical SMILES string.
    pub smiles: String,
    /// Normalized or dictionary-canonical name.
    pub canonical_name: String,
    /// Which resolution step produced this result.
    pub source: ResolveSource,
    /// Molecular formula in Hill notation (e.g. `"C2H6O"`).
    /// `None` for DirectSmiles dictionary entries.
    pub molecular_formula: Option<String>,
    /// Molecular weight in g/mol (IUPAC 2021 standard atomic weights).
    /// `None` for DirectSmiles dictionary entries.
    pub molecular_weight: Option<f64>,
    /// Confidence score in `0.0..=1.0`.
    ///
    /// | Scenario | Score |
    /// |----------|-------|
    /// | Exact dict match, DirectSmiles | 1.00 |
    /// | Exact dict match, CanonicalName | 0.95 |
    /// | Dict match after normalization  | 0.90 |
    /// | IUPAC parser                   | 0.85 |
    pub confidence: f64,
    /// Standard InChI identifier (e.g. `"InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"`).
    /// `None` for DirectSmiles dictionary entries.
    pub inchi: Option<String>,
    /// 27-character InChIKey (SHA-256-based hash of the InChI).
    /// `None` for DirectSmiles dictionary entries.
    pub inchi_key: Option<String>,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Resolve any chemical name to SMILES and associated identifiers.
///
/// Pipeline: normalize → dictionary lookup → IUPAC parser.
///
/// # Errors
///
/// Returns [`ResolveError::NotFound`] if the name is not in the dictionary and
/// cannot be parsed as a systematic IUPAC name. Returns [`ResolveError::ParseError`]
/// on structural parse failures.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::resolve;
///
/// let r = resolve("propan-2-one").unwrap();
/// assert_eq!(r.smiles, "CC(=O)C");
/// assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
/// assert!(r.confidence >= 0.85);
/// ```
pub fn resolve(name: &str) -> Result<ResolveResult, ResolveError> {
    let normalized = normalizer::normalize_lowercase(name);
    // Was non-trivial normalization needed?
    // Non-ASCII input (CJK, Greek) or structural normalization (fullwidth, prolonged mark).
    let was_normalized = name.chars().any(|c| !c.is_ascii())
        || name.to_ascii_lowercase() != normalized;

    if let Some(entry) = dict::lookup_synonym(&normalized) {
        match entry {
            dict::DictEntry::CanonicalName(canonical) => {
                let graph = parser::parse_iupac(canonical)?;
                let smiles = parser::smiles::to_smiles(&graph);
                let molecular_formula = Some(parser::formula::molecular_formula(&graph));
                let molecular_weight = Some(parser::formula::molecular_weight(&graph));
                let (inchi, inchi_key) = parser::inchi::mol_to_inchi_pair(&graph);
                let confidence = if was_normalized { 0.90 } else { 0.95 };
                return Ok(ResolveResult {
                    smiles,
                    canonical_name: canonical.to_string(),
                    source: ResolveSource::Dictionary,
                    molecular_formula,
                    molecular_weight,
                    confidence,
                    inchi: Some(inchi),
                    inchi_key: Some(inchi_key),
                });
            }
            dict::DictEntry::DirectSmiles(smiles) => {
                let confidence = if was_normalized { 0.90 } else { 1.00 };
                return Ok(ResolveResult {
                    smiles: smiles.to_string(),
                    canonical_name: normalized.clone(),
                    source: ResolveSource::Dictionary,
                    molecular_formula: None,
                    molecular_weight: None,
                    confidence,
                    inchi: None,
                    inchi_key: None,
                });
            }
        }
    }

    let graph = parser::parse_iupac(&normalized)?;
    let smiles = parser::smiles::to_smiles(&graph);
    let molecular_formula = Some(parser::formula::molecular_formula(&graph));
    let molecular_weight = Some(parser::formula::molecular_weight(&graph));
    let (inchi, inchi_key) = parser::inchi::mol_to_inchi_pair(&graph);
    Ok(ResolveResult {
        smiles,
        canonical_name: normalized,
        source: ResolveSource::Parser,
        molecular_formula,
        molecular_weight,
        confidence: 0.85,
        inchi: Some(inchi),
        inchi_key: Some(inchi_key),
    })
}

/// Resolve a batch of chemical names, returning one result per input in order.
///
/// Each entry is independently resolved; errors in one name do not affect others.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::resolve_batch;
///
/// let results = resolve_batch(&["methane", "ethanol", "water"]);
/// assert_eq!(results.len(), 3);
/// assert_eq!(results[0].as_ref().unwrap().smiles, "C");
/// ```
pub fn resolve_batch(names: &[&str]) -> Vec<Result<ResolveResult, ResolveError>> {
    names.iter().map(|&n| resolve(n)).collect()
}

/// Normalize a chemical name (CJK-safe, zero-copy when already normalized).
///
/// Applies full-width → half-width, katakana prolonged mark → hyphen,
/// Greek letter → ASCII, and whitespace collapse/trim.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::normalize;
///
/// assert_eq!(normalize("α-D-glucose"), "alpha-D-glucose");
/// assert_eq!(normalize("ethanol"), "ethanol"); // zero-copy
/// ```
pub fn normalize(name: &str) -> Cow<'_, str> {
    normalizer::normalize(name)
}

/// Dictionary-only lookup (no IUPAC parsing).
///
/// Returns `None` if the name (after lowercase) is not in any dictionary.
/// Use [`resolve`] for full resolution including systematic names.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::lookup;
///
/// assert!(lookup("benzene").is_some());
/// assert!(lookup("xyzzy").is_none());
/// ```
pub fn lookup(name: &str) -> Option<dict::DictEntry<'static>> {
    dict::lookup_synonym(name)
}

/// Parse a SMILES string into an IUPAC systematic name (straight-chain acyclic only).
///
/// Supports the same functional groups as the IUPAC parser:
/// alkanes, alkenes, alkynes, alcohols, ketones, aldehydes, carboxylic acids,
/// amines, thiols, nitriles, and amides.
///
/// # Errors
///
/// Returns [`ResolveError::ParseError`] for:
/// - Branched or cyclic molecules
/// - Aromatic atoms (lowercase SMILES)
/// - Unsupported elements
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::smiles_to_iupac;
///
/// assert_eq!(smiles_to_iupac("CCCCO").unwrap(), "butan-1-ol");
/// assert_eq!(smiles_to_iupac("CCCC").unwrap(), "butane");
/// assert!(smiles_to_iupac("CC(C)CC").is_err()); // branched → error
/// ```
pub fn smiles_to_iupac(smiles: &str) -> Result<String, ResolveError> {
    let graph = parser::smiles_parser::parse_smiles(smiles)?;
    parser::iupac_namer::mol_to_iupac(&graph)
}

/// Generate a Standard InChI identifier from a SMILES string.
///
/// Scope matches [`smiles_to_iupac`]: non-aromatic acyclic molecules with
/// C, H, N, O, S, F, Cl, Br, I atoms.
///
/// # Errors
///
/// Returns [`ResolveError::ParseError`] if the SMILES cannot be parsed.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::smiles_to_inchi;
///
/// let inchi = smiles_to_inchi("CCO").unwrap();
/// assert!(inchi.starts_with("InChI=1S/C2H6O/"));
/// ```
pub fn smiles_to_inchi(smiles: &str) -> Result<String, ResolveError> {
    let graph = parser::smiles_parser::parse_smiles(smiles)?;
    Ok(parser::inchi::mol_to_inchi(&graph))
}

/// Generate an InChIKey (27-character SHA-256-based hash) from a SMILES string.
///
/// # Errors
///
/// Returns [`ResolveError::ParseError`] if the SMILES cannot be parsed.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::smiles_to_inchikey;
///
/// let key = smiles_to_inchikey("CCO").unwrap();
/// assert_eq!(key.len(), 27);
/// assert!(key.chars().all(|c| c.is_ascii_uppercase() || c == '-'));
/// ```
pub fn smiles_to_inchikey(smiles: &str) -> Result<String, ResolveError> {
    let graph = parser::smiles_parser::parse_smiles(smiles)?;
    let inchi = parser::inchi::mol_to_inchi(&graph);
    Ok(parser::inchi::inchi_to_key(&inchi))
}
