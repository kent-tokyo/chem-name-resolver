use std::collections::BTreeMap;

use super::{Element, MolGraph};

/// Standard atomic weights (IUPAC 2021, abridged).
fn atomic_weight(element: &Element) -> f64 {
    match element {
        Element::H => 1.008,
        Element::C => 12.011,
        Element::N => 14.007,
        Element::O => 15.999,
        Element::F => 18.998,
        Element::P => 30.974,
        Element::S => 32.06,
        Element::Cl => 35.45,
        Element::Br => 79.904,
        Element::I => 126.904,
    }
}

/// Compute the molecular formula in Hill notation (C first, H second, rest alphabetical).
pub fn molecular_formula(graph: &MolGraph) -> String {
    let mut counts: BTreeMap<&'static str, u32> = BTreeMap::new();

    for atom in &graph.atoms {
        *counts.entry(atom.element.symbol()).or_insert(0) += 1;
        if atom.implicit_h > 0 {
            *counts.entry("H").or_insert(0) += atom.implicit_h as u32;
        }
    }

    let mut out = String::new();

    // Hill notation: C first
    if let Some(&c) = counts.get("C") {
        out.push('C');
        if c > 1 {
            out.push_str(&c.to_string());
        }
    }
    // H second
    if let Some(&h) = counts.get("H") {
        out.push('H');
        if h > 1 {
            out.push_str(&h.to_string());
        }
    }
    // Remaining elements in alphabetical order (BTreeMap is already sorted)
    for (sym, &count) in &counts {
        if *sym == "C" || *sym == "H" {
            continue;
        }
        out.push_str(sym);
        if count > 1 {
            out.push_str(&count.to_string());
        }
    }

    out
}

/// Compute the molecular weight (g/mol) using standard atomic weights.
pub fn molecular_weight(graph: &MolGraph) -> f64 {
    let mut mw = 0.0;
    for atom in &graph.atoms {
        mw += atomic_weight(&atom.element);
        mw += atom.implicit_h as f64 * 1.008; // implicit H
    }
    mw
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_iupac;

    fn formula(name: &str) -> String {
        molecular_formula(&parse_iupac(name).unwrap())
    }

    fn weight(name: &str) -> f64 {
        molecular_weight(&parse_iupac(name).unwrap())
    }

    #[test]
    fn methane_formula() {
        assert_eq!(formula("methane"), "CH4");
    }

    #[test]
    fn ethanol_formula() {
        assert_eq!(formula("ethanol"), "C2H6O");
    }

    #[test]
    fn propan_2_one_formula() {
        assert_eq!(formula("propan-2-one"), "C3H6O");
    }

    #[test]
    fn two_chlorobutane_formula() {
        assert_eq!(formula("2-chlorobutane"), "C4H9Cl");
    }

    #[test]
    fn ethanoic_acid_formula() {
        assert_eq!(formula("ethanoic acid"), "C2H4O2");
    }

    #[test]
    fn methane_weight() {
        let mw = weight("methane");
        // C(12.011) + 4*H(1.008) = 16.043
        assert!((mw - 16.043).abs() < 0.001, "methane Mw = {mw:.3}");
    }

    #[test]
    fn ethanol_weight() {
        let mw = weight("ethanol");
        // 2*C + 6*H + O = 24.022 + 6.048 + 15.999 = 46.069
        assert!((mw - 46.069).abs() < 0.001, "ethanol Mw = {mw:.3}");
    }

    #[test]
    fn extended_chain_formula() {
        assert_eq!(formula("icosane"), "C20H42");
        assert_eq!(formula("undecane"), "C11H24");
    }
}
