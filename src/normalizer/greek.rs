/// Map Greek letters to their ASCII equivalents used in chemical nomenclature.
/// Returns a &'static str because the replacement is always longer than 1 char.
#[inline]
pub(crate) fn map_greek(c: char) -> Option<&'static str> {
    match c {
        'α' => Some("alpha"),
        'β' => Some("beta"),
        'γ' => Some("gamma"),
        'δ' => Some("delta"),
        'ε' => Some("epsilon"),
        'ζ' => Some("zeta"),
        'η' => Some("eta"),
        'θ' => Some("theta"),
        'κ' => Some("kappa"),
        'λ' => Some("lambda"),
        'μ' => Some("mu"),
        'ν' => Some("nu"),
        'ξ' => Some("xi"),
        'π' => Some("pi"),
        'σ' => Some("sigma"),
        'τ' => Some("tau"),
        'φ' => Some("phi"),
        'χ' => Some("chi"),
        'ψ' => Some("psi"),
        'ω' => Some("omega"),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_greek() {
        assert_eq!(map_greek('α'), Some("alpha"));
        assert_eq!(map_greek('β'), Some("beta"));
        assert_eq!(map_greek('γ'), Some("gamma"));
        assert_eq!(map_greek('δ'), Some("delta"));
    }

    #[test]
    fn non_greek_unchanged() {
        assert_eq!(map_greek('a'), None);
        assert_eq!(map_greek('あ'), None);
    }
}
