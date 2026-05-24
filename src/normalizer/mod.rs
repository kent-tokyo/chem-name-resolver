//! Chemical name text normalizer (CJK-safe, zero-copy when already normalized).
//!
//! # Transformations (applied in order)
//!
//! 1. **CJK full-width ASCII** → half-width: `ａ` → `a`, `０` → `0`
//! 2. **Katakana prolonged sound mark** → ASCII hyphen: `ー` → `-`
//! 3. **Greek letters** → ASCII equivalents: `α` → `alpha`, `β` → `beta`
//! 4. **Whitespace** → collapsed and trimmed
//!
//! Uses [`std::borrow::Cow`] to avoid allocation when the input is already normalized.

mod fullwidth;
mod greek;
mod katakana;

use std::borrow::Cow;

use fullwidth::map_fullwidth;
use greek::map_greek;
use katakana::map_katakana;

/// Normalize a chemical name. Returns `Borrowed` if no transformation was needed
/// (already normalized input — zero allocation). Allocates at most once when
/// transformations are required.
///
/// Transformations applied in order:
///   1. CJK full-width ASCII → half-width  (ａ→a, ０→0, …)
///   2. Katakana prolonged sound mark ー → ASCII hyphen '-'
///   3. Greek letters → ASCII equivalents  (α→alpha, β→beta, …)
///   4. Collapse runs of ASCII whitespace to single space, trim
pub fn normalize(input: &str) -> Cow<'_, str> {
    if needs_normalization(input) {
        Cow::Owned(normalize_to_string(input))
    } else {
        Cow::Borrowed(input)
    }
}

/// Normalize and convert to ASCII lowercase. Always returns an owned String.
pub fn normalize_lowercase(input: &str) -> String {
    normalize(input).to_ascii_lowercase()
}

fn needs_normalization(input: &str) -> bool {
    let mut prev_was_space = false;
    let mut leading = true;

    for c in input.chars() {
        if map_fullwidth(c).is_some() || map_katakana(c).is_some() || map_greek(c).is_some() {
            return true;
        }
        if c == ' ' || c == '\t' || c == '\n' || c == '\r' {
            if leading || prev_was_space {
                return true;
            }
            prev_was_space = true;
        } else {
            leading = false;
            prev_was_space = false;
        }
    }
    // trailing space check
    if prev_was_space && !input.is_empty() {
        return true;
    }
    false
}

fn normalize_to_string(input: &str) -> String {
    let mut out = String::with_capacity(input.len() + 16);
    let mut prev_was_space = false;
    let mut leading = true;

    for c in input.chars() {
        if let Some(half) = map_fullwidth(c) {
            flush_space(&mut out, &mut prev_was_space, &mut leading);
            out.push(half);
        } else if let Some(repl) = map_katakana(c) {
            flush_space(&mut out, &mut prev_was_space, &mut leading);
            out.push(repl);
        } else if let Some(ascii) = map_greek(c) {
            flush_space(&mut out, &mut prev_was_space, &mut leading);
            out.push_str(ascii);
        } else if c == ' ' || c == '\t' || c == '\n' || c == '\r' {
            if !leading {
                prev_was_space = true;
            }
        } else {
            flush_space(&mut out, &mut prev_was_space, &mut leading);
            out.push(c);
        }
    }
    out
}

#[inline]
fn flush_space(out: &mut String, prev_was_space: &mut bool, leading: &mut bool) {
    if *prev_was_space {
        out.push(' ');
        *prev_was_space = false;
    }
    *leading = false;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn already_normalized_is_borrowed() {
        assert!(matches!(normalize("2-pentanone"), Cow::Borrowed(_)));
        assert!(matches!(normalize("ethanol"), Cow::Borrowed(_)));
    }

    #[test]
    fn fullwidth_digits_and_hyphens() {
        assert_eq!(normalize("２－ペンタノン"), "2-ペンタノン");
    }

    #[test]
    fn prolonged_sound_mark() {
        assert_eq!(normalize("ジエチルエーテル"), "ジエチルエ-テル");
    }

    #[test]
    fn greek_letters() {
        assert_eq!(normalize("α-D-glucose"), "alpha-D-glucose");
        assert_eq!(normalize("β-carotene"), "beta-carotene");
    }

    #[test]
    fn whitespace_collapse() {
        assert_eq!(normalize("  2 - pentanone  "), "2 - pentanone");
        assert_eq!(normalize("a  b"), "a b");
    }

    #[test]
    fn idempotent() {
        let cases = ["α-D-glucose", "ジエチルエーテル", "２－ペンタノン", "propan-2-one"];
        for s in cases {
            let once = normalize(s).into_owned();
            let twice = normalize(&once).into_owned();
            assert_eq!(once, twice, "not idempotent for: {s}");
        }
    }

    #[test]
    fn normalize_lowercase_works() {
        assert_eq!(normalize_lowercase("Ethanol"), "ethanol");
        assert_eq!(normalize_lowercase("α-D-Glucose"), "alpha-d-glucose");
    }
}
