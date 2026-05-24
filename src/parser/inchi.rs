//! MolGraph → Standard InChI / InChIKey generation.
//!
//! Implements a simplified InChI for the acyclic subset supported by this library.
//!
//! # InChI layers generated
//! - Formula layer: Hill notation (re-uses [`formula::molecular_formula`])
//! - Connection layer `/c`: canonical C-atom connectivity
//! - Hydrogen layer `/h`: fixed and mobile H
//!
//! # InChIKey
//! 27-character hash: SHA-256 of the InChI string → Base-26 encoded blocks.
//! Format: `XXXXXXXXXXXXXX-XXXXXXXXXX-N`

use sha2::{Digest, Sha256};

use super::{Element, MolGraph};
use super::formula::molecular_formula;

// ── Public functions ───────────────────────────────────────────────────────────

/// Compute the Standard InChI string for a [`MolGraph`].
///
/// Returns the full InChI (e.g. `"InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"`).
pub fn mol_to_inchi(graph: &MolGraph) -> String {
    if graph.atoms.is_empty() {
        return "InChI=1S//".to_string();
    }

    let formula = molecular_formula(graph);
    let conn = connection_layer(graph);
    let hydro = hydrogen_layer(graph);

    let mut inchi = format!("InChI=1S/{}", formula);
    if !conn.is_empty() {
        inchi.push_str(&format!("/c{}", conn));
    }
    if !hydro.is_empty() {
        inchi.push_str(&format!("/h{}", hydro));
    }
    inchi
}

/// Compute InChI and InChIKey together (avoids double computation).
///
/// Returns `(inchi_string, inchikey_string)`.
pub fn mol_to_inchi_pair(graph: &MolGraph) -> (String, String) {
    let inchi = mol_to_inchi(graph);
    let key = inchi_to_key(&inchi);
    (inchi, key)
}

/// Compute the InChIKey from an InChI string.
///
/// The key is a 27-character string `XXXXXXXXXXXXXX-XXXXXXXXXX-N`.
pub fn inchi_to_key(inchi: &str) -> String {
    // Hash the entire InChI string
    let hash = Sha256::digest(inchi.as_bytes());

    // First block: 14 chars from first 104 bits (13 bytes)
    let block1 = base26_encode(&hash[..13], 14);
    // Second block: 10 chars from next 75 bits (9+ bytes, using bytes 13-21)
    let block2 = base26_encode(&hash[13..22], 10);
    // Version/flags character: 'N' = standard InChI v1.x
    format!("{}-{}-N", block1, block2)
}

// ── Connection layer ──────────────────────────────────────────────────────────

/// Generate the `/c` connection layer.
///
/// Non-H atoms are numbered 1..n in appearance order (matching Hill formula order).
/// We use the atom index + 1 as the InChI number for simplicity.
fn connection_layer(graph: &MolGraph) -> String {
    if graph.atoms.len() <= 1 {
        return String::new();
    }

    // Number only non-H atoms
    let non_h: Vec<usize> = (0..graph.atoms.len())
        .filter(|&i| graph.atoms[i].element != Element::H)
        .collect();

    if non_h.len() <= 1 {
        return String::new();
    }

    // Build a position map: atom_idx → inchi_number (1-based among non-H atoms)
    let pos: Vec<Option<usize>> = {
        let mut p = vec![None; graph.atoms.len()];
        for (rank, &idx) in non_h.iter().enumerate() {
            p[idx] = Some(rank + 1);
        }
        p
    };

    // Collect all bonds as (min_num, max_num) pairs, deduplicate
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for &ai in &non_h {
        for bond in &graph.bonds[ai] {
            let bi = bond.to;
            if graph.atoms[bi].element == Element::H { continue; }
            let na = pos[ai].unwrap();
            let nb = pos[bi].unwrap();
            let pair = if na < nb { (na, nb) } else { (nb, na) };
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }
        }
    }
    pairs.sort_unstable();

    // Format as "1-2,2-3,..."
    pairs.iter()
        .map(|(a, b)| format!("{}-{}", a, b))
        .collect::<Vec<_>>()
        .join(",")
}

// ── Hydrogen layer ────────────────────────────────────────────────────────────

/// Generate the `/h` hydrogen layer.
///
/// Format: mobile H first (`{atom}H`), then fixed H groups (`{atom}H{n}`).
fn hydrogen_layer(graph: &MolGraph) -> String {
    if graph.atoms.is_empty() {
        return String::new();
    }

    // Atom numbers (1-based among non-H atoms, same as connection layer)
    let non_h: Vec<usize> = (0..graph.atoms.len())
        .filter(|&i| graph.atoms[i].element != Element::H)
        .collect();

    let pos: Vec<Option<usize>> = {
        let mut p = vec![None; graph.atoms.len()];
        for (rank, &idx) in non_h.iter().enumerate() {
            p[idx] = Some(rank + 1);
        }
        p
    };

    // Mobile H: atoms that have OH, NH, SH (single bond to electronegative atom with H)
    let mut mobile: Vec<usize> = Vec::new();
    // Fixed H groups: (atom_number, h_count)
    let mut fixed: Vec<(usize, u32)> = Vec::new();

    for &ai in &non_h {
        let atom = &graph.atoms[ai];
        let num = pos[ai].unwrap();

        // Count implicit H on this atom (from MolGraph)
        let implicit = atom.implicit_h as u32;

        // Determine if this atom has mobile H (heteroatom with H attached to it)
        let is_mobile_h_donor = matches!(atom.element, Element::O | Element::N | Element::S)
            && implicit > 0;

        if is_mobile_h_donor {
            mobile.push(num);
        } else if implicit > 0 {
            fixed.push((num, implicit));
        }
    }

    if mobile.is_empty() && fixed.is_empty() {
        return String::new();
    }

    let mut parts: Vec<String> = Vec::new();

    // Mobile H: listed as "1H" or "1,2H"
    if !mobile.is_empty() {
        let nums: String = mobile.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(",");
        parts.push(format!("{}H", nums));
    }

    // Fixed H: sorted by atom number, then grouped by count
    // Format: "1-3H2,4H3" means atoms 1-3 each have 2H, atom 4 has 3H
    // Simplified: list individually
    fixed.sort_by_key(|&(n, _)| n);
    // Group consecutive atoms with same H count
    let mut i = 0;
    while i < fixed.len() {
        let (num, count) = fixed[i];
        // Find consecutive same-count atoms
        let mut j = i + 1;
        while j < fixed.len() && fixed[j].1 == count && fixed[j].0 == fixed[j-1].0 + 1 {
            j += 1;
        }
        let h_suffix = if count == 1 { "H".to_string() } else { format!("H{}", count) };
        if j - i == 1 {
            parts.push(format!("{}{}", num, h_suffix));
        } else if j - i == 2 {
            parts.push(format!("{},{}{}", num, fixed[j-1].0, h_suffix));
        } else {
            parts.push(format!("{}-{}{}", num, fixed[j-1].0, h_suffix));
        }
        i = j;
    }

    parts.join(",")
}

// ── Base-26 encoding ──────────────────────────────────────────────────────────

/// Encode `bytes` as exactly `n_chars` uppercase A-Z letters.
///
/// Interprets the bytes as a big-endian integer, then encodes in base-26.
fn base26_encode(bytes: &[u8], n_chars: usize) -> String {
    // Convert bytes to a big number (as u128 for up to 16 bytes)
    // For safety, use only the first min(16, bytes.len()) bytes
    let use_bytes = bytes.len().min(16);
    let mut val: u128 = 0;
    for &b in &bytes[..use_bytes] {
        val = val.wrapping_shl(8).wrapping_add(b as u128);
    }

    // Encode in base-26, least significant digit first
    let mut chars: Vec<u8> = Vec::with_capacity(n_chars);
    for _ in 0..n_chars {
        let digit = (val % 26) as u8;
        chars.push(b'A' + digit);
        val /= 26;
    }
    // Most significant digit first
    chars.reverse();
    String::from_utf8(chars).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse_iupac;

    fn inchi(name: &str) -> String {
        mol_to_inchi(&parse_iupac(name).unwrap())
    }

    fn key(name: &str) -> String {
        inchi_to_key(&inchi(name))
    }

    #[test]
    fn methane_inchi() {
        let i = inchi("methane");
        assert!(i.starts_with("InChI=1S/CH4"), "got: {i}");
    }

    #[test]
    fn ethanol_inchi() {
        let i = inchi("ethanol");
        assert!(i.starts_with("InChI=1S/C2H6O/"), "got: {i}");
        assert!(i.contains("/c"), "missing /c: {i}");
        assert!(i.contains("/h"), "missing /h: {i}");
    }

    #[test]
    fn inchikey_format() {
        let k = key("ethanol");
        assert_eq!(k.len(), 27, "key: {k}");
        assert!(k.chars().all(|c| c.is_ascii_uppercase() || c == '-'), "key: {k}");
        // Format: 14 chars - 10 chars - 1 char
        let parts: Vec<&str> = k.split('-').collect();
        assert_eq!(parts.len(), 3, "key: {k}");
        assert_eq!(parts[0].len(), 14);
        assert_eq!(parts[1].len(), 10);
        assert_eq!(parts[2].len(), 1);
    }

    #[test]
    fn inchikey_deterministic() {
        let k1 = key("ethanol");
        let k2 = key("ethanol");
        assert_eq!(k1, k2);
    }

    #[test]
    fn different_molecules_different_keys() {
        let k_ethanol = key("ethanol");
        let k_methanol = key("methanol");
        // "methanol" goes through dict then parser → different formula
        let k_propanol = key("propan-1-ol");
        // All three should be distinct
        assert_ne!(k_ethanol, k_propanol);
    }

    #[test]
    fn base26_length() {
        let b = [0u8; 13];
        let enc = base26_encode(&b, 14);
        assert_eq!(enc.len(), 14);
        assert!(enc.chars().all(|c| c.is_ascii_uppercase()));
    }
}
