/// Parse an alkane chain-length stem from the beginning of `input`.
/// Returns (chain_length, remaining_str) on success.
pub fn parse_stem(input: &str) -> Option<(u8, &str)> {
    // IMPORTANT: longer stems must appear before any stem that is a prefix of them.
    // e.g. "pentadec" before "pent", "hexadec" before "hex", etc.
    // Otherwise "pentadecane" would match "pent" (C5) instead of "pentadec" (C15).
    const STEMS: &[(&str, u8)] = &[
        // C11-C20 (all before their conflicting C1-C10 counterparts)
        ("nonadec", 19),
        ("octadec", 18),
        ("heptadec", 17),
        ("hexadec", 16),
        ("pentadec", 15),
        ("tetradec", 14),
        ("tridec", 13),
        ("dodec", 12),
        ("undec", 11),
        ("eicos", 20), // alternative form (eicosane)
        ("icos", 20),  // IUPAC preferred form (icosane)
        // C1-C10
        ("meth", 1),
        ("eth", 2),
        ("prop", 3),
        ("but", 4),
        ("pent", 5),
        ("hex", 6),
        ("hept", 7),
        ("oct", 8),
        ("non", 9),
        ("dec", 10),
    ];

    for (stem, n) in STEMS {
        if input.starts_with(stem) {
            return Some((*n, &input[stem.len()..]));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_stems_c1_c10() {
        assert_eq!(parse_stem("methane"), Some((1, "ane")));
        assert_eq!(parse_stem("ethanol"), Some((2, "anol")));
        assert_eq!(parse_stem("propan-2-one"), Some((3, "an-2-one")));
        assert_eq!(parse_stem("butane"), Some((4, "ane")));
        assert_eq!(parse_stem("pentane"), Some((5, "ane")));
        assert_eq!(parse_stem("hexane"), Some((6, "ane")));
        assert_eq!(parse_stem("heptane"), Some((7, "ane")));
        assert_eq!(parse_stem("octane"), Some((8, "ane")));
        assert_eq!(parse_stem("nonane"), Some((9, "ane")));
        assert_eq!(parse_stem("decane"), Some((10, "ane")));
    }

    #[test]
    fn extended_stems_c11_c20() {
        assert_eq!(parse_stem("undecane"), Some((11, "ane")));
        assert_eq!(parse_stem("dodecane"), Some((12, "ane")));
        assert_eq!(parse_stem("tridecane"), Some((13, "ane")));
        assert_eq!(parse_stem("tetradecane"), Some((14, "ane")));
        assert_eq!(parse_stem("pentadecane"), Some((15, "ane")));
        assert_eq!(parse_stem("hexadecane"), Some((16, "ane")));
        assert_eq!(parse_stem("heptadecane"), Some((17, "ane")));
        assert_eq!(parse_stem("octadecane"), Some((18, "ane")));
        assert_eq!(parse_stem("nonadecane"), Some((19, "ane")));
        assert_eq!(parse_stem("icosane"), Some((20, "ane")));
        assert_eq!(parse_stem("eicosane"), Some((20, "ane")));
    }

    #[test]
    fn no_prefix_collision() {
        // These must NOT match the shorter stem (e.g. "pent" C5 instead of "pentadec" C15)
        assert_eq!(parse_stem("pentadecane"), Some((15, "ane")));
        assert_eq!(parse_stem("hexadecane"), Some((16, "ane")));
        assert_eq!(parse_stem("heptadecane"), Some((17, "ane")));
        assert_eq!(parse_stem("octadecane"), Some((18, "ane")));
        assert_eq!(parse_stem("nonadecane"), Some((19, "ane")));
        // And the shorter stems still work as-is
        assert_eq!(parse_stem("pentane"), Some((5, "ane")));
        assert_eq!(parse_stem("hexane"), Some((6, "ane")));
        assert_eq!(parse_stem("heptane"), Some((7, "ane")));
        assert_eq!(parse_stem("octane"), Some((8, "ane")));
        assert_eq!(parse_stem("nonane"), Some((9, "ane")));
    }

    #[test]
    fn no_match() {
        assert_eq!(parse_stem("benzene"), None);
        assert_eq!(parse_stem("cyclohexane"), None);
    }
}
