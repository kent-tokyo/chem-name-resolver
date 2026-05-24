//! SMILES string → [`MolGraph`] parser.
//!
//! Supports the organic-subset SMILES used by the existing SMILES generator:
//! - Atoms: `C N O S F Cl Br I` (implicit H), bracketed atoms `[OH]` etc.
//! - Bonds: implicit single, `=` double, `#` triple
//! - Branches: `(...)`
//! - Ring closures: digit `1`-`9` or `%nn`
//!
//! **Not supported** (returns [`ResolveError::ParseError`]):
//! - Aromatic atoms (`c`, `n`, `o`, …)
//! - Stereo (`@`, `/`, `\`)
//! - Isotopes (`[13C]`)

use std::collections::HashMap;

use crate::error::ResolveError;
use super::{Atom, Bond, BondOrder, Element, MolGraph};

/// Parse a SMILES string into a [`MolGraph`].
///
/// # Errors
///
/// Returns [`ResolveError::ParseError`] for unsupported syntax or malformed input.
///
/// # Examples
///
/// ```rust
/// use chem_name_resolver::parser::smiles_parser::parse_smiles;
///
/// let g = parse_smiles("CCO").unwrap();
/// assert_eq!(g.atoms.len(), 3); // 2 C + 1 O
/// ```
pub fn parse_smiles(smiles: &str) -> Result<MolGraph, ResolveError> {
    let mut graph = MolGraph::default();
    // Stack of atom indices for branch tracking.
    let mut stack: Vec<usize> = Vec::new();
    // Ring-closure map: ring_number → (atom_index, bond_order)
    let mut ring_open: HashMap<u32, (usize, BondOrder)> = HashMap::new();
    // Next bond order override (from `=` or `#` before the atom)
    let mut next_bond = BondOrder::Single;
    // Index of the previously added atom (-1 = none yet)
    let mut prev: Option<usize> = None;

    let mut chars = smiles.char_indices().peekable();

    while let Some((pos, c)) = chars.next() {
        match c {
            // ── Bonds ──────────────────────────────────────────────────────
            '=' => { next_bond = BondOrder::Double; continue; }
            '#' => { next_bond = BondOrder::Triple; continue; }
            '-' => { next_bond = BondOrder::Single; continue; }

            // ── Branches ───────────────────────────────────────────────────
            '(' => {
                if let Some(p) = prev {
                    stack.push(p);
                }
                continue;
            }
            ')' => {
                prev = stack.pop();
                next_bond = BondOrder::Single;
                continue;
            }

            // ── Ring closures ──────────────────────────────────────────────
            '%' => {
                // %nn two-digit ring closure
                let d1 = chars.next().ok_or_else(|| err(pos, "expected digit after %"))?.1;
                let d2 = chars.next().ok_or_else(|| err(pos, "expected two digits after %"))?.1;
                if !d1.is_ascii_digit() || !d2.is_ascii_digit() {
                    return Err(err(pos, "non-digit in %nn ring closure"));
                }
                let num = (d1 as u32 - b'0' as u32) * 10 + (d2 as u32 - b'0' as u32);
                close_ring(&mut graph, &mut ring_open, prev, num, next_bond.clone(), pos)?;
                next_bond = BondOrder::Single;
                continue;
            }
            '0'..='9' => {
                let num = c as u32 - b'0' as u32;
                close_ring(&mut graph, &mut ring_open, prev, num, next_bond.clone(), pos)?;
                next_bond = BondOrder::Single;
                continue;
            }

            // ── Bracketed atom [X] ─────────────────────────────────────────
            '[' => {
                let elem = parse_bracket_atom(smiles, pos, &mut chars)?;
                let idx = add_atom_to_graph(&mut graph, &mut prev, &next_bond, elem);
                prev = Some(idx);
                next_bond = BondOrder::Single;
                continue;
            }

            // ── Organic-subset atoms ───────────────────────────────────────
            'C' => {
                // 'C' alone = Carbon, 'Cl' = Chlorine
                let elem = if chars.peek().map(|&(_, x)| x) == Some('l') {
                    chars.next();
                    Element::Cl
                } else {
                    Element::C
                };
                let idx = add_atom_to_graph(&mut graph, &mut prev, &next_bond, elem);
                prev = Some(idx);
                next_bond = BondOrder::Single;
            }
            'N' | 'O' | 'S' | 'P' | 'F' | 'I' => {
                let elem = char_to_element(c, pos)?;
                let idx = add_atom_to_graph(&mut graph, &mut prev, &next_bond, elem);
                prev = Some(idx);
                next_bond = BondOrder::Single;
            }
            'B' => {
                // Could be Br
                if chars.peek().map(|&(_, x)| x) == Some('r') {
                    chars.next();
                    let idx = add_atom_to_graph(&mut graph, &mut prev, &next_bond, Element::Br);
                    prev = Some(idx);
                } else {
                    return Err(err(pos, "bare 'B' atom not supported; use [B] or Br"));
                }
                next_bond = BondOrder::Single;
            }
            'c' | 'n' | 'o' | 's' | 'p' => {
                return Err(err(pos, "aromatic atoms not supported"));
            }
            '.' => {
                // Disconnected component — not supported
                return Err(err(pos, "disconnected SMILES ('.') not supported"));
            }
            '@' | '/' | '\\' => {
                return Err(err(pos, "stereochemistry not supported"));
            }
            ' ' | '\t' | '\n' | '\r' => {
                // Skip whitespace (should not normally appear)
                continue;
            }
            other => {
                return Err(err(pos, &format!("unexpected character {other:?}")));
            }
        }
    }

    if !ring_open.is_empty() {
        return Err(ResolveError::ParseError {
            pos: smiles.len(),
            msg: "unclosed ring closure in SMILES".to_string(),
        });
    }

    graph.fill_implicit_h();
    Ok(graph)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn err(pos: usize, msg: &str) -> ResolveError {
    ResolveError::ParseError { pos, msg: msg.to_string() }
}

fn char_to_element(c: char, pos: usize) -> Result<Element, ResolveError> {
    match c {
        'C' => Ok(Element::C),
        'N' => Ok(Element::N),
        'O' => Ok(Element::O),
        'S' => Ok(Element::S),
        'P' => Ok(Element::P),
        'F' => Ok(Element::F),
        'I' => Ok(Element::I),
        other => Err(err(pos, &format!("unsupported element symbol {other:?}"))),
    }
}

/// Parse a bracketed atom `[...]`, consuming chars up through `]`.
/// Returns the element parsed.
fn parse_bracket_atom(
    _smiles: &str,
    open_pos: usize,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
) -> Result<Element, ResolveError> {
    let mut sym = String::new();

    // Collect up to ']', skipping charge/H specs
    loop {
        match chars.next() {
            None => return Err(err(open_pos, "unclosed '['")),
            Some((_, ']')) => break,
            Some((_, c)) if c.is_ascii_uppercase() || (sym.is_empty() && c == '*') => {
                sym.push(c);
                // Collect lowercase continuation of symbol (e.g. "Cl", "Br")
                while let Some(&(_, nc)) = chars.peek() {
                    if nc.is_ascii_lowercase() && nc != 'h' {
                        sym.push(nc);
                        chars.next();
                    } else {
                        break;
                    }
                }
            }
            Some((_, 'H' | 'h')) => { /* implicit H count, skip */ }
            Some((_, '+' | '-')) => { /* charge, skip */ }
            Some((_, c)) if c.is_ascii_digit() => { /* count, skip */ }
            Some((_, '@')) => { /* stereo, skip */ }
            Some((_, other)) => {
                return Err(err(open_pos, &format!("unexpected char {other:?} in bracket atom")));
            }
        }
    }

    match sym.as_str() {
        "C" => Ok(Element::C),
        "N" => Ok(Element::N),
        "O" => Ok(Element::O),
        "S" => Ok(Element::S),
        "P" => Ok(Element::P),
        "F" => Ok(Element::F),
        "Cl" => Ok(Element::Cl),
        "Br" => Ok(Element::Br),
        "I" => Ok(Element::I),
        "H" => Ok(Element::H),
        other => Err(err(open_pos, &format!("unsupported element in brackets: {other:?}"))),
    }
}

fn add_atom_to_graph(
    graph: &mut MolGraph,
    prev: &mut Option<usize>,
    bond: &BondOrder,
    elem: Element,
) -> usize {
    let idx = graph.atoms.len();
    graph.atoms.push(Atom { element: elem, charge: 0, implicit_h: 0 });
    graph.bonds.push(Vec::new());
    if let Some(p) = *prev {
        graph.bonds[p].push(Bond { to: idx, order: bond.clone() });
        graph.bonds[idx].push(Bond { to: p, order: bond.clone() });
    }
    idx
}

fn close_ring(
    graph: &mut MolGraph,
    ring_open: &mut HashMap<u32, (usize, BondOrder)>,
    current: Option<usize>,
    num: u32,
    bond: BondOrder,
    pos: usize,
) -> Result<(), ResolveError> {
    let cur = current.ok_or_else(|| err(pos, "ring closure with no current atom"))?;
    if let Some((open_idx, open_bond)) = ring_open.remove(&num) {
        // Use the higher-order bond if two are specified
        let order = if bond != BondOrder::Single { bond } else { open_bond };
        graph.bonds[open_idx].push(Bond { to: cur, order: order.clone() });
        graph.bonds[cur].push(Bond { to: open_idx, order });
    } else {
        ring_open.insert(num, (cur, bond));
    }
    Ok(())
}

// ── Cl/Br two-char lookahead is handled in main loop ─────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::smiles::to_smiles;

    fn roundtrip(s: &str) -> String {
        to_smiles(&parse_smiles(s).unwrap_or_else(|e| panic!("parse_smiles({s:?}): {e}")))
    }

    #[test]
    fn simple_alkanes() {
        assert_eq!(roundtrip("C"), "C");
        assert_eq!(roundtrip("CC"), "CC");
        assert_eq!(roundtrip("CCC"), "CCC");
    }

    #[test]
    fn ethanol() {
        // CCO → parse → to_smiles should give "CCO"
        let g = parse_smiles("CCO").unwrap();
        assert_eq!(g.atoms.len(), 3);
        assert_eq!(to_smiles(&g), "CCO");
    }

    #[test]
    fn double_bond() {
        let g = parse_smiles("CC=CC").unwrap();
        assert_eq!(to_smiles(&g), "CC=CC");
    }

    #[test]
    fn ring() {
        let g = parse_smiles("C1CCCCC1").unwrap();
        assert_eq!(to_smiles(&g), "C1CCCCC1");
    }

    #[test]
    fn branch() {
        // isopropanol: CC(C)O
        let g = parse_smiles("CC(C)O").unwrap();
        assert_eq!(to_smiles(&g), "CC(C)O");
    }

    #[test]
    fn chlorine() {
        let g = parse_smiles("CCCl").unwrap();
        assert_eq!(to_smiles(&g), "CCCl");
    }

    #[test]
    fn nitrile() {
        let g = parse_smiles("CC#N").unwrap();
        assert_eq!(to_smiles(&g), "CC#N");
    }

    #[test]
    fn aromatic_error() {
        assert!(parse_smiles("c1ccccc1").is_err());
    }
}
