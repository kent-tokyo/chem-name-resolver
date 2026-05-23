mod synonyms;
mod word_map;

use crate::normalizer;

#[derive(Debug, Clone, PartialEq)]
pub enum DictEntry<'a> {
    /// Canonical IUPAC name (may differ from input by synonym expansion).
    CanonicalName(&'a str),
    /// Pre-computed SMILES for names too complex for the IUPAC parser.
    DirectSmiles(&'a str),
}

/// Normalize input, then look up in synonym and katakana maps.
/// Returns None if the name is not in the dictionary (caller should try parser).
pub fn lookup_synonym(name: &str) -> Option<DictEntry<'static>> {
    let normalized = normalizer::normalize_lowercase(name);
    let key: &str = &normalized;

    if let Some(iupac) = synonyms::SYNONYM_TO_IUPAC.get(key) {
        return Some(DictEntry::CanonicalName(iupac));
    }
    if let Some(smiles) = synonyms::SYNONYM_TO_SMILES.get(key) {
        return Some(DictEntry::DirectSmiles(smiles));
    }
    // Katakana word lookup
    if let Some(iupac) = word_map::KATAKANA_TO_IUPAC.get(key) {
        return Some(DictEntry::CanonicalName(iupac));
    }
    // Chinese (Hanzi) lookup — try IUPAC mapping first, then direct SMILES
    if let Some(iupac) = word_map::HANZI_TO_IUPAC.get(key) {
        return Some(DictEntry::CanonicalName(iupac));
    }
    if let Some(smiles) = word_map::HANZI_TO_SMILES.get(key) {
        return Some(DictEntry::DirectSmiles(smiles));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trivial_name_lookup() {
        assert_eq!(lookup_synonym("water"), Some(DictEntry::DirectSmiles("O")));
        assert_eq!(lookup_synonym("WATER"), Some(DictEntry::DirectSmiles("O")));
        assert_eq!(lookup_synonym("acetone"), Some(DictEntry::CanonicalName("propan-2-one")));
    }

    #[test]
    fn direct_smiles_lookup() {
        assert_eq!(lookup_synonym("benzene"), Some(DictEntry::DirectSmiles("c1ccccc1")));
        assert_eq!(lookup_synonym("toluene"), Some(DictEntry::DirectSmiles("Cc1ccccc1")));
    }

    #[test]
    fn katakana_lookup() {
        assert_eq!(lookup_synonym("メタン"), Some(DictEntry::CanonicalName("methane")));
        assert_eq!(lookup_synonym("エタノール"), Some(DictEntry::CanonicalName("ethanol")));
    }

    #[test]
    fn unknown_returns_none() {
        assert_eq!(lookup_synonym("completely_unknown_xyz_123"), None);
    }

    #[test]
    fn hanzi_lookup() {
        assert_eq!(lookup_synonym("甲烷"), Some(DictEntry::CanonicalName("methane")));
        assert_eq!(lookup_synonym("乙醇"), Some(DictEntry::CanonicalName("ethanol")));
        assert_eq!(lookup_synonym("苯"), Some(DictEntry::DirectSmiles("c1ccccc1")));
        assert_eq!(lookup_synonym("丙酮"), Some(DictEntry::CanonicalName("propan-2-one")));
        assert_eq!(lookup_synonym("乙酸"), Some(DictEntry::CanonicalName("ethanoic acid")));
    }
}
