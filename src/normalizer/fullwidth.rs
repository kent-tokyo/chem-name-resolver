/// Map CJK full-width ASCII characters to their half-width equivalents.
/// U+FF01..=U+FF5E → U+0021..=U+007E (linear arithmetic, no table needed).
#[inline]
pub(crate) fn map_fullwidth(c: char) -> Option<char> {
    let u = c as u32;
    if (0xFF01..=0xFF5E).contains(&u) {
        char::from_u32(u - 0xFEE0)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fullwidth_digits() {
        assert_eq!(map_fullwidth('０'), Some('0'));
        assert_eq!(map_fullwidth('９'), Some('9'));
    }

    #[test]
    fn fullwidth_letters() {
        assert_eq!(map_fullwidth('Ａ'), Some('A'));
        assert_eq!(map_fullwidth('ｚ'), Some('z'));
    }

    #[test]
    fn non_fullwidth_unchanged() {
        assert_eq!(map_fullwidth('a'), None);
        assert_eq!(map_fullwidth('あ'), None);
    }
}
