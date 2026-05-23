/// Functional group suffix parsed from an IUPAC name.
#[derive(Debug, Clone, PartialEq)]
pub enum Suffix {
    Ane,
    Ene,
    Yne,
    Ol,
    One,
    Al,
    OicAcid,
    Amide,
    Amine,
    Thiol,
    Nitrile,
}

/// Optional multiplier prefix before a suffix (di, tri, tetra).
#[derive(Debug, Clone, PartialEq)]
pub enum Multiplier {
    Di,
    Tri,
    Tetra,
}

impl Multiplier {
    pub fn count(&self) -> u8 {
        match self {
            Multiplier::Di => 2,
            Multiplier::Tri => 3,
            Multiplier::Tetra => 4,
        }
    }
}

/// Parsed suffix group: optionally prefixed with a multiplier and locants.
#[derive(Debug, Clone, PartialEq)]
pub struct SuffixGroup {
    pub multiplier: Option<Multiplier>,
    pub locants: Vec<u8>,
    pub suffix: Suffix,
}

/// Parse a suffix group from the beginning of `input`.
/// Returns (SuffixGroup, remaining) on success.
pub fn parse_suffix(input: &str) -> Option<(SuffixGroup, &str)> {
    // Try infix locant first (e.g., "hex-2-ene": locant between stem and suffix)
    let (locants, rest) = if let Some(r) = try_infix_locant(input) {
        r
    } else {
        (vec![], input)
    };

    // Optional multiplier (di/tri/tetra)
    let (multiplier, rest) = parse_multiplier(rest);

    // Core suffix — longer matches first to avoid prefix collisions.
    // "yn"/"en" are elided forms used in compound suffixes (e.g. "hept-3-yn-1-ol").
    let suffixes: &[(&str, Suffix)] = &[
        ("oic acid", Suffix::OicAcid),
        ("nitrile", Suffix::Nitrile),
        ("amide", Suffix::Amide),
        ("amine", Suffix::Amine),
        ("thiol", Suffix::Thiol),
        ("ane", Suffix::Ane),
        ("ene", Suffix::Ene),
        ("yne", Suffix::Yne),
        ("yn", Suffix::Yne),
        ("en", Suffix::Ene),
        ("one", Suffix::One),
        ("ol", Suffix::Ol),
        ("al", Suffix::Al),
    ];

    for (tag, suf) in suffixes {
        if rest.starts_with(tag) {
            let remaining = &rest[tag.len()..];
            return Some((
                SuffixGroup {
                    multiplier,
                    locants,
                    suffix: suf.clone(),
                },
                remaining,
            ));
        }
    }
    None
}

/// Try to consume an infix locant like "-2-" in "hex-2-ene".
fn try_infix_locant(input: &str) -> Option<(Vec<u8>, &str)> {
    let rest = input.strip_prefix('-')?;
    crate::parser::locant::parse_locant_list(rest)
}

fn parse_multiplier(input: &str) -> (Option<Multiplier>, &str) {
    if input.starts_with("tetra") {
        return (Some(Multiplier::Tetra), &input[5..]);
    }
    if input.starts_with("tri") {
        return (Some(Multiplier::Tri), &input[3..]);
    }
    if input.starts_with("di") {
        return (Some(Multiplier::Di), &input[2..]);
    }
    (None, input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_ane() {
        let (sg, rest) = parse_suffix("ane").unwrap();
        assert_eq!(sg.suffix, Suffix::Ane);
        assert_eq!(sg.multiplier, None);
        assert!(sg.locants.is_empty());
        assert_eq!(rest, "");
    }

    #[test]
    fn diol() {
        let (sg, rest) = parse_suffix("diol").unwrap();
        assert_eq!(sg.suffix, Suffix::Ol);
        assert_eq!(sg.multiplier, Some(Multiplier::Di));
        assert_eq!(rest, "");
    }

    #[test]
    fn infix_locant_ene() {
        let (sg, rest) = parse_suffix("-2-ene").unwrap();
        assert_eq!(sg.suffix, Suffix::Ene);
        assert_eq!(sg.locants, vec![2]);
        assert_eq!(rest, "");
    }

    #[test]
    fn oic_acid() {
        let (sg, rest) = parse_suffix("oic acid").unwrap();
        assert_eq!(sg.suffix, Suffix::OicAcid);
        assert_eq!(rest, "");
    }
}
