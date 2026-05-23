/// Parse a locant list like "2-", "2,4-", "1,3,5-".
/// Returns (locants, remaining) on success; locants is 1-indexed carbon positions.
pub fn parse_locant_list(input: &str) -> Option<(Vec<u8>, &str)> {
    let mut locants = Vec::new();
    let mut rest = input;

    loop {
        // Parse an integer.
        let digit_end = rest
            .char_indices()
            .take_while(|(_, c)| c.is_ascii_digit())
            .last()
            .map(|(i, c)| i + c.len_utf8())?;

        let n: u8 = rest[..digit_end].parse().ok()?;
        locants.push(n);
        rest = &rest[digit_end..];

        if rest.starts_with(',') {
            rest = &rest[1..]; // eat comma, continue
        } else if rest.starts_with('-') {
            rest = &rest[1..]; // eat trailing hyphen, done
            break;
        } else {
            return None; // locant list must end with '-'
        }
    }

    if locants.is_empty() {
        None
    } else {
        Some((locants, rest))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_locant() {
        assert_eq!(parse_locant_list("2-ol"), Some((vec![2], "ol")));
        assert_eq!(parse_locant_list("1-ene"), Some((vec![1], "ene")));
    }

    #[test]
    fn multiple_locants() {
        assert_eq!(parse_locant_list("2,4-dione"), Some((vec![2, 4], "dione")));
        assert_eq!(parse_locant_list("1,3,5-triol"), Some((vec![1, 3, 5], "triol")));
    }

    #[test]
    fn no_locant() {
        assert_eq!(parse_locant_list("ane"), None);
        assert_eq!(parse_locant_list("ol"), None);
    }
}
