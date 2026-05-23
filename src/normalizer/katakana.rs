/// Map Katakana-Hiragana Prolonged Sound Mark (ー U+30FC) to ASCII hyphen.
/// Other katakana normalizations can be added here as needed.
#[inline]
pub(crate) fn map_katakana(c: char) -> Option<char> {
    match c {
        'ー' => Some('-'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prolonged_sound_mark() {
        assert_eq!(map_katakana('ー'), Some('-'));
    }

    #[test]
    fn regular_katakana_unchanged() {
        assert_eq!(map_katakana('ア'), None);
        assert_eq!(map_katakana('メ'), None);
    }
}
