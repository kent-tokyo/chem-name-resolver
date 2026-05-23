use crate::parser::locant::parse_locant_list;

/// A substituent attached at a specific carbon position.
#[derive(Debug, Clone, PartialEq)]
pub struct Substituent {
    pub locants: Vec<u8>,
    pub kind: SubstituentKind,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SubstituentKind {
    Chloro,
    Bromo,
    Fluoro,
    Iodo,
    Methyl,
    Ethyl,
    Propyl,
    Isopropyl,
    Butyl,
    SecButyl,
    IsoButyl,
    TertButyl,
    Pentyl,
    Hexyl,
    Hydroxy,
    Oxo,
    Amino,
    Mercapto,
    Cyano,
    Acetyl,
    Formyl,
}

const SUBSTITUENTS: &[(&str, SubstituentKind)] = &[
    ("mercapto",    SubstituentKind::Mercapto),
    ("chloro",      SubstituentKind::Chloro),
    ("bromo",       SubstituentKind::Bromo),
    ("fluoro",      SubstituentKind::Fluoro),
    ("iodo",        SubstituentKind::Iodo),
    ("methyl",      SubstituentKind::Methyl),
    ("ethyl",       SubstituentKind::Ethyl),
    // Branched variants before n-alkyl to avoid prefix collision (e.g. "isopropyl"
    // must be matched before "propyl", "tert-butyl" before "butyl").
    ("isopropyl",   SubstituentKind::Isopropyl),
    ("tert-butyl",  SubstituentKind::TertButyl),
    ("sec-butyl",   SubstituentKind::SecButyl),
    ("isobutyl",    SubstituentKind::IsoButyl),
    ("hexyl",       SubstituentKind::Hexyl),
    ("pentyl",      SubstituentKind::Pentyl),
    ("butyl",       SubstituentKind::Butyl),
    ("propyl",      SubstituentKind::Propyl),
    ("hydroxy",     SubstituentKind::Hydroxy),
    ("cyano",       SubstituentKind::Cyano),
    ("amino",       SubstituentKind::Amino),
    ("oxo",         SubstituentKind::Oxo),
    ("acetyl",      SubstituentKind::Acetyl),
    ("formyl",      SubstituentKind::Formyl),
];

/// Strip an optional di/tri/tetra multiplier prefix (e.g. "di" from "dichloro").
/// The multiplier count is implicit in the locant list; we only need the root name.
fn strip_substituent_multiplier(s: &str) -> &str {
    if s.starts_with("tetra") { &s[5..] }
    else if s.starts_with("tri") { &s[3..] }
    else if s.starts_with("di") { &s[2..] }
    else { s }
}

/// Parse one substituent prefix like "2-chloro-" from `input`.
/// Returns (Substituent, remaining) on success.
pub fn parse_one_substituent(input: &str) -> Option<(Substituent, &str)> {
    // Must start with a locant.
    let (locants, rest) = parse_locant_list(input)?;

    // Strip optional di/tri/tetra multiplier before the substituent name.
    // E.g. "2,3-dichloro-" → locants=[2,3], rest after strip = "chloro-"
    let name_part = strip_substituent_multiplier(rest);

    for (name, kind) in SUBSTITUENTS {
        if name_part.starts_with(name) {
            let after = &name_part[name.len()..];
            let after = after.strip_prefix('-').unwrap_or(after);
            return Some((Substituent { locants, kind: kind.clone() }, after));
        }
    }
    None
}

/// Parse zero or more substituent prefixes from the start of `input`.
pub fn parse_substituents(input: &str) -> (Vec<Substituent>, &str) {
    let mut subs = Vec::new();
    let mut rest = input;
    while let Some((sub, next)) = parse_one_substituent(rest) {
        subs.push(sub);
        rest = next;
    }
    (subs, rest)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_chloro() {
        let (subs, rest) = parse_substituents("2-chlorobutane");
        assert_eq!(subs.len(), 1);
        assert_eq!(subs[0].kind, SubstituentKind::Chloro);
        assert_eq!(subs[0].locants, vec![2]);
        assert_eq!(rest, "butane");
    }

    #[test]
    fn methyl_substituent() {
        let (subs, rest) = parse_substituents("3-methylpentane");
        assert_eq!(subs.len(), 1);
        assert_eq!(subs[0].kind, SubstituentKind::Methyl);
        assert_eq!(rest, "pentane");
    }

    #[test]
    fn no_substituent() {
        let (subs, rest) = parse_substituents("pentane");
        assert!(subs.is_empty());
        assert_eq!(rest, "pentane");
    }
}
