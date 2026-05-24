//! MolGraph → IUPAC systematic name generator (straight-chain acyclic scope).
//!
//! Covers the same functional groups as [`parse_iupac`]: alkanes, alkenes, alkynes,
//! alcohols, ketones, aldehydes, carboxylic acids, amines, thiols, nitriles, amides.
//!
//! # Limitations
//! - Branched chains are not supported (returns [`ResolveError::ParseError`]).
//! - Cyclic molecules are not supported.
//! - No stereochemistry.

use crate::error::ResolveError;
use super::{BondOrder, Element, MolGraph};

/// Convert a [`MolGraph`] into a systematic IUPAC name.
///
/// Only straight-chain acyclic molecules are supported.
/// Branched or cyclic structures return [`ResolveError::ParseError`].
pub fn mol_to_iupac(graph: &MolGraph) -> Result<String, ResolveError> {
    if graph.atoms.is_empty() {
        return Err(ResolveError::ParseError { pos: 0, msg: "empty molecule".to_string() });
    }

    // Collect all carbon indices
    let carbons: Vec<usize> = (0..graph.atoms.len())
        .filter(|&i| graph.atoms[i].element == Element::C)
        .collect();

    if carbons.is_empty() {
        return Err(ResolveError::ParseError { pos: 0, msg: "no carbon atoms".to_string() });
    }

    // Verify it's not cyclic: if bond_count >= atom_count it's cyclic
    let bond_count: usize = graph.bonds.iter().map(|b| b.len()).sum::<usize>() / 2;
    if bond_count >= graph.atoms.len() {
        return Err(ResolveError::ParseError {
            pos: 0,
            msg: "cyclic molecules not supported by iupac_namer".to_string(),
        });
    }

    // Find the longest linear carbon chain (handling simple substituents attached to chain C)
    // First, collect the carbon backbone: C atoms connected only to other C atoms or heteroatoms
    // (not to other C branches — i.e., no C with 3+ C neighbours).
    for &ci in &carbons {
        let c_neighbors: usize = graph.bonds[ci]
            .iter()
            .filter(|b| graph.atoms[b.to].element == Element::C)
            .count();
        if c_neighbors > 2 {
            return Err(ResolveError::ParseError {
                pos: 0,
                msg: "branched chains not supported by iupac_namer".to_string(),
            });
        }
    }

    // Find terminal carbons (≤1 carbon neighbour)
    let terminals: Vec<usize> = carbons.iter()
        .copied()
        .filter(|&ci| {
            graph.bonds[ci].iter()
                .filter(|b| graph.atoms[b.to].element == Element::C)
                .count() <= 1
        })
        .collect();

    // Trace chain from one terminal
    let start = terminals[0];
    let chain = trace_chain(graph, start);
    let chain_len = chain.len();

    if chain_len == 0 {
        return Err(ResolveError::ParseError { pos: 0, msg: "could not trace carbon chain".to_string() });
    }

    // Identify the principal characteristic group (PCG) and its IUPAC priority rank.
    // We try numbering from both ends and pick the one giving the lowest locant for the PCG.
    let chain_rev: Vec<usize> = chain.iter().rev().copied().collect();

    let fwd = analyze_chain(graph, &chain)?;
    let rev = analyze_chain(graph, &chain_rev)?;

    // Choose direction: lower locant set for PCG first, then for substituents
    let chosen = choose_direction(fwd, rev);

    assemble_name(graph, &chosen, chain_len)
}

// ── Chain tracing ─────────────────────────────────────────────────────────────

/// DFS from `start` following only C→C bonds, returns ordered chain.
fn trace_chain(graph: &MolGraph, start: usize) -> Vec<usize> {
    let mut chain = Vec::new();
    let mut prev = usize::MAX;
    let mut cur = start;
    loop {
        chain.push(cur);
        let next = graph.bonds[cur].iter()
            .find(|b| graph.atoms[b.to].element == Element::C && b.to != prev)
            .map(|b| b.to);
        match next {
            Some(n) => { prev = cur; cur = n; }
            None => break,
        }
    }
    chain
}

// ── Chain analysis ────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct ChainInfo {
    /// Chain atoms in order (index 0 = C-1)
    #[allow(dead_code)]
    chain: Vec<usize>,
    /// Principal characteristic group suffix (highest-priority)
    #[allow(dead_code)]
    pcg: Option<PcgKind>,
    /// Locant(s) of PCG (1-based)
    pcg_locants: Vec<usize>,
    /// Suffix string for the PCG (e.g. "-ol", "-one")
    suffix: String,
    /// Prefix substituents sorted alphabetically
    prefixes: Vec<PrefixEntry>,
    /// Bond order changes (for ene/yne): (locant, 'e'=double, 'y'=triple)
    unsat: Vec<(usize, char)>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum PcgKind {
    CarboxylicAcid = 8,
    Amide          = 7,
    Nitrile        = 6,
    Aldehyde       = 5,
    Ketone         = 4,
    Alcohol        = 3,
    Amine          = 2,
    Thiol          = 1,
}

#[derive(Debug, Clone)]
struct PrefixEntry {
    /// Alphabetic name used for sorting ("chloro", "methyl", etc.)
    name: String,
    /// 1-based locants
    locants: Vec<usize>,
}

fn analyze_chain(graph: &MolGraph, chain: &[usize]) -> Result<ChainInfo, ResolveError> {
    let n = chain.len();
    let mut pcg: Option<PcgKind> = None;
    let mut pcg_locants: Vec<usize> = Vec::new();
    let mut unsat: Vec<(usize, char)> = Vec::new(); // (locant, 'e' or 'y')
    let mut prefix_map: std::collections::BTreeMap<String, Vec<usize>> = std::collections::BTreeMap::new();

    // Examine each chain carbon
    for (pos0, &ci) in chain.iter().enumerate() {
        let locant = pos0 + 1; // 1-based

        // Check bond to next chain carbon for unsaturation
        if pos0 + 1 < n {
            let next_ci = chain[pos0 + 1];
            let bond_order = graph.bonds[ci]
                .iter()
                .find(|b| b.to == next_ci)
                .map(|b| &b.order);
            match bond_order {
                Some(BondOrder::Double) => unsat.push((locant, 'e')),
                Some(BondOrder::Triple) => unsat.push((locant, 'y')),
                _ => {}
            }
        }

        // Examine heteroatom neighbours
        for bond in &graph.bonds[ci] {
            if graph.atoms[bond.to].element == Element::C {
                continue; // skip carbon neighbours (chain or branches handled elsewhere)
            }
            let neighbour = &graph.atoms[bond.to];
            match (&neighbour.element, &bond.order) {
                (Element::O, BondOrder::Double) => {
                    // =O: could be ketone (if middle) or aldehyde (if terminal with no other O)
                    // or part of carboxylic acid (check for -OH on same carbon)
                    let has_oh = graph.bonds[ci].iter().any(|b| {
                        b.to != bond.to
                            && graph.atoms[b.to].element == Element::O
                            && b.order == BondOrder::Single
                    });
                    if has_oh {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::CarboxylicAcid, locant);
                    } else if locant == 1 || locant == n {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Aldehyde, locant);
                    } else {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Ketone, locant);
                    }
                }
                (Element::O, BondOrder::Single) => {
                    // -OH (single bond O) — but only if not also =O on same C
                    let has_oxo = graph.bonds[ci].iter().any(|b| {
                        b.to != bond.to
                            && graph.atoms[b.to].element == Element::O
                            && b.order == BondOrder::Double
                    });
                    if !has_oxo {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Alcohol, locant);
                    }
                    // else: this -OH belongs to carboxylic acid, already handled above
                }
                (Element::N, BondOrder::Single) => {
                    // Check if there is also =O on same carbon (amide)
                    let has_oxo = graph.bonds[ci].iter().any(|b| {
                        graph.atoms[b.to].element == Element::O
                            && b.order == BondOrder::Double
                    });
                    if has_oxo {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Amide, locant);
                    } else {
                        update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Amine, locant);
                    }
                }
                (Element::N, BondOrder::Triple) => {
                    update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Nitrile, locant);
                }
                (Element::S, BondOrder::Single) => {
                    update_pcg(&mut pcg, &mut pcg_locants, PcgKind::Thiol, locant);
                }
                // Halogens and other substituents → prefix
                (Element::F, _) => { prefix_map.entry("fluoro".to_string()).or_default().push(locant); }
                (Element::Cl, _) => { prefix_map.entry("chloro".to_string()).or_default().push(locant); }
                (Element::Br, _) => { prefix_map.entry("bromo".to_string()).or_default().push(locant); }
                (Element::I, _) => { prefix_map.entry("iodo".to_string()).or_default().push(locant); }
                _ => {
                    // Other heteroatom patterns: ignore for now
                }
            }
        }

        // Check for carbon substituents (methyl, ethyl, etc.) — neighbours that are C
        // but NOT part of the main chain.
        let chain_set: std::collections::HashSet<usize> = chain.iter().copied().collect();
        for bond in &graph.bonds[ci] {
            if graph.atoms[bond.to].element != Element::C { continue; }
            if chain_set.contains(&bond.to) { continue; }
            // It's a carbon substituent. Measure its length.
            let sub_len = measure_branch(graph, bond.to, ci);
            let sub_name = alkyl_prefix(sub_len);
            prefix_map.entry(sub_name).or_default().push(locant);
        }
    }

    // Build suffix string
    let suffix = build_suffix(graph, chain, &pcg, &pcg_locants, n, &unsat)?;

    // Build prefix list
    let mut prefixes: Vec<PrefixEntry> = prefix_map
        .into_iter()
        .map(|(name, mut locs)| {
            locs.sort_unstable();
            PrefixEntry { name, locants: locs }
        })
        .collect();
    prefixes.sort_by(|a, b| a.name.cmp(&b.name));

    Ok(ChainInfo {
        chain: chain.to_vec(),
        pcg,
        pcg_locants: pcg_locants.clone(),
        suffix,
        prefixes,
        unsat,
    })
}

fn update_pcg(pcg: &mut Option<PcgKind>, locants: &mut Vec<usize>, kind: PcgKind, loc: usize) {
    match pcg {
        None => { *pcg = Some(kind); locants.push(loc); }
        Some(existing) => {
            if kind > *existing {
                *pcg = Some(kind);
                locants.clear();
                locants.push(loc);
            } else if kind == *existing {
                locants.push(loc);
            }
        }
    }
}

/// Measure straight-chain length of a substituent branch.
fn measure_branch(graph: &MolGraph, start: usize, parent: usize) -> usize {
    let mut len = 1;
    let mut prev = parent;
    let mut cur = start;
    loop {
        let next = graph.bonds[cur].iter()
            .find(|b| graph.atoms[b.to].element == Element::C && b.to != prev)
            .map(|b| b.to);
        match next {
            Some(n) => { prev = cur; cur = n; len += 1; }
            None => break,
        }
    }
    len
}

fn alkyl_prefix(len: usize) -> String {
    match len {
        1 => "methyl",
        2 => "ethyl",
        3 => "propyl",
        4 => "butyl",
        5 => "pentyl",
        6 => "hexyl",
        _ => "alkyl",
    }.to_string()
}

// ── Suffix builder ────────────────────────────────────────────────────────────
//
// IUPAC suffix pattern after the stem:
//   - Pure alkane (no PCG, no unsat):  "ane"
//   - Pure alkene/yne (no PCG):        "-{loc}-ene" / "-{loc}-yne"
//   - PCG, no unsat:                   "an-{loc}-{term}"  e.g. "an-1-ol"
//   - PCG + unsat:                     "-{uloc}-en-{ploc}-{term}"  e.g. "-2-en-1-ol"
//
// "an" connector only appears when there is a PCG and no unsaturation.
// Unsaturation infix uses "en"/"yn" (without trailing 'e') when followed by PCG.

fn build_suffix(
    _graph: &MolGraph,
    chain: &[usize],
    pcg: &Option<PcgKind>,
    pcg_locants: &[usize],
    n: usize,
    unsat: &[(usize, char)],
) -> Result<String, ResolveError> {
    let _ = chain;
    let suffix = match pcg {
        None => {
            if unsat.is_empty() {
                "ane".to_string()
            } else {
                // Pure ene/yne: "-2-ene" or "-2-yne"
                unsat_pure_suffix(unsat)
            }
        }
        Some(kind) => {
            // PCG present: "an" + optional_unsat_infix + pcg_part
            let unsat_part = if unsat.is_empty() {
                "an".to_string()
            } else {
                unsat_infix_for_pcg(unsat)
            };
            let pcg_part = pcg_suffix_part(kind, pcg_locants, n);
            format!("{}{}", unsat_part, pcg_part)
        }
    };
    let _ = n;
    Ok(suffix)
}

/// For pure unsaturation (no PCG): "-2-ene", "-3-yne", "-1,3-diene"
fn unsat_pure_suffix(unsat: &[(usize, char)]) -> String {
    let enes: Vec<usize> = unsat.iter().filter(|&&(_, t)| t == 'e').map(|&(l, _)| l).collect();
    let ynes: Vec<usize> = unsat.iter().filter(|&&(_, t)| t == 'y').map(|&(l, _)| l).collect();
    let mut parts = Vec::new();
    if !enes.is_empty() {
        let lstr = locs_str(&enes);
        let mul = multiplier_prefix(enes.len());
        parts.push(format!("-{}-{}ene", lstr, mul));
    }
    if !ynes.is_empty() {
        let lstr = locs_str(&ynes);
        let mul = multiplier_prefix(ynes.len());
        parts.push(format!("-{}-{}yne", lstr, mul));
    }
    parts.join("")
}

/// For unsaturation when followed by a PCG: "-2-en" (no trailing 'e')
fn unsat_infix_for_pcg(unsat: &[(usize, char)]) -> String {
    let enes: Vec<usize> = unsat.iter().filter(|&&(_, t)| t == 'e').map(|&(l, _)| l).collect();
    let ynes: Vec<usize> = unsat.iter().filter(|&&(_, t)| t == 'y').map(|&(l, _)| l).collect();
    let mut parts = Vec::new();
    if !enes.is_empty() {
        let lstr = locs_str(&enes);
        let mul = multiplier_prefix(enes.len());
        parts.push(format!("-{}-{}en", lstr, mul));
    }
    if !ynes.is_empty() {
        let lstr = locs_str(&ynes);
        let mul = multiplier_prefix(ynes.len());
        parts.push(format!("-{}-{}yn", lstr, mul));
    }
    parts.join("")
}

/// The PCG part of the suffix after the "an"/unsat connector.
/// e.g. Alcohol locant [1] → "-1-ol", CarboxylicAcid locant [1] n=2 → "oic acid"
fn pcg_suffix_part(kind: &PcgKind, pcg_locants: &[usize], n: usize) -> String {
    match kind {
        PcgKind::CarboxylicAcid => {
            let locs = locant_str_hide_terminal(pcg_locants, n);
            if pcg_locants.len() > 1 {
                format!("{}edioic acid", locs)
            } else {
                format!("{}oic acid", locs)
            }
        }
        PcgKind::Amide => "amide".to_string(),
        PcgKind::Nitrile => "enitrile".to_string(),
        PcgKind::Aldehyde => "al".to_string(),
        PcgKind::Ketone => {
            let locs = locant_str_show(pcg_locants);
            if pcg_locants.len() > 1 {
                format!("{}edione", locs)
            } else {
                format!("{}one", locs)
            }
        }
        PcgKind::Alcohol => {
            let locs = locant_str_show(pcg_locants);
            if pcg_locants.len() > 1 {
                format!("{}ediol", locs)
            } else {
                format!("{}ol", locs)
            }
        }
        PcgKind::Amine => {
            // Amine starts with vowel → elide 'e', no 'e' connector needed.
            let locs = locant_str_show(pcg_locants);
            format!("{}amine", locs) // "amine" or "-1-amine"
        }
        PcgKind::Thiol => {
            // Thiol starts with consonant → keep 'e' connector.
            let locs = locant_str_show(pcg_locants);
            format!("e{}thiol", locs) // "ethiol" or "e-1-thiol"
        }
    }
}

fn locs_str(locs: &[usize]) -> String {
    locs.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",")
}

/// Always show locant (for ketone, alcohol – position matters)
fn locant_str_show(locants: &[usize]) -> String {
    if locants.is_empty() { return String::new(); }
    format!("-{}-", locs_str(locants))
}

/// Hide locant if at terminal C-1 or C-n (for acid, aldehyde)
fn locant_str_hide_terminal(locants: &[usize], n: usize) -> String {
    if locants.is_empty() { return String::new(); }
    if locants.len() == 1 && (locants[0] == 1 || locants[0] == n) {
        return String::new();
    }
    format!("-{}-", locs_str(locants))
}

fn multiplier_prefix(n: usize) -> &'static str {
    match n { 1 => "", 2 => "di", 3 => "tri", 4 => "tetra", _ => "poly" }
}

// ── Direction chooser ─────────────────────────────────────────────────────────

fn choose_direction(fwd: ChainInfo, rev: ChainInfo) -> ChainInfo {
    // Compare locant sets: lower = better (IUPAC rule)
    let fwd_pcg = &fwd.pcg_locants;
    let rev_pcg = &rev.pcg_locants;

    // Compare PCG locants first
    let pcg_cmp = compare_locant_sets(fwd_pcg, rev_pcg);
    if pcg_cmp != std::cmp::Ordering::Equal {
        return if pcg_cmp == std::cmp::Ordering::Less { fwd } else { rev };
    }

    // Compare unsaturation locants
    let fwd_unsat: Vec<usize> = fwd.unsat.iter().map(|&(l, _)| l).collect();
    let rev_unsat: Vec<usize> = rev.unsat.iter().map(|&(l, _)| l).collect();
    let unsat_cmp = compare_locant_sets(&fwd_unsat, &rev_unsat);
    if unsat_cmp != std::cmp::Ordering::Equal {
        return if unsat_cmp == std::cmp::Ordering::Less { fwd } else { rev };
    }

    // Compare substituent locants alphabetically by prefix name then locants
    let fwd_sub_locs: Vec<usize> = fwd.prefixes.iter().flat_map(|p| p.locants.iter().copied()).collect();
    let rev_sub_locs: Vec<usize> = rev.prefixes.iter().flat_map(|p| p.locants.iter().copied()).collect();
    let sub_cmp = compare_locant_sets(&fwd_sub_locs, &rev_sub_locs);
    if sub_cmp == std::cmp::Ordering::Less { fwd } else { rev }
}

fn compare_locant_sets(a: &[usize], b: &[usize]) -> std::cmp::Ordering {
    for (x, y) in a.iter().zip(b.iter()) {
        let c = x.cmp(y);
        if c != std::cmp::Ordering::Equal { return c; }
    }
    a.len().cmp(&b.len())
}

// ── Name assembler ────────────────────────────────────────────────────────────

fn assemble_name(
    _graph: &MolGraph,
    info: &ChainInfo,
    chain_len: usize,
) -> Result<String, ResolveError> {
    let stem = chain_stem(chain_len)?;
    let mut name = String::new();

    // Prefix part (alphabetical substituents)
    for entry in &info.prefixes {
        let mult = multiplier_prefix(entry.locants.len());
        let lstr: String = entry.locants.iter().map(|l| l.to_string()).collect::<Vec<_>>().join(",");
        name.push_str(&format!("{}-{}{}-", lstr, mult, entry.name));
    }

    name.push_str(stem);
    name.push_str(&info.suffix);
    Ok(name)
}

/// Reverse lookup: chain length → IUPAC stem.
fn chain_stem(n: usize) -> Result<&'static str, ResolveError> {
    match n {
        1 => Ok("meth"),
        2 => Ok("eth"),
        3 => Ok("prop"),
        4 => Ok("but"),
        5 => Ok("pent"),
        6 => Ok("hex"),
        7 => Ok("hept"),
        8 => Ok("oct"),
        9 => Ok("non"),
        10 => Ok("dec"),
        11 => Ok("undec"),
        12 => Ok("dodec"),
        13 => Ok("tridec"),
        14 => Ok("tetradec"),
        15 => Ok("pentadec"),
        16 => Ok("hexadec"),
        17 => Ok("heptadec"),
        18 => Ok("octadec"),
        19 => Ok("nonadec"),
        20 => Ok("icos"),
        _ => Err(ResolveError::ParseError {
            pos: 0,
            msg: format!("chain length {n} out of supported range (1-20)"),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::smiles_parser::parse_smiles;

    fn iupac(smiles: &str) -> String {
        let g = parse_smiles(smiles).unwrap_or_else(|e| panic!("parse_smiles({smiles:?}): {e}"));
        mol_to_iupac(&g).unwrap_or_else(|e| panic!("mol_to_iupac({smiles:?}): {e}"))
    }

    #[test]
    fn simple_alkanes() {
        assert_eq!(iupac("C"), "methane");
        assert_eq!(iupac("CC"), "ethane");
        assert_eq!(iupac("CCC"), "propane");
        assert_eq!(iupac("CCCC"), "butane");
        assert_eq!(iupac("CCCCC"), "pentane");
    }

    #[test]
    fn alcohols() {
        assert_eq!(iupac("CCO"), "ethan-1-ol");
        assert_eq!(iupac("CCCO"), "propan-1-ol");
        assert_eq!(iupac("CCCCO"), "butan-1-ol");
    }

    #[test]
    fn ketone() {
        assert_eq!(iupac("CC(=O)C"), "propan-2-one");
    }

    #[test]
    fn aldehyde() {
        assert_eq!(iupac("CC=O"), "ethanal");
    }

    #[test]
    fn carboxylic_acid() {
        assert_eq!(iupac("CC(=O)O"), "ethanoic acid");
    }

    #[test]
    fn alkene() {
        assert_eq!(iupac("CC=CC"), "but-2-ene");
    }

    #[test]
    fn alkyne() {
        assert_eq!(iupac("CC#CC"), "but-2-yne");
    }

    #[test]
    fn nitrile() {
        assert_eq!(iupac("CC#N"), "ethanenitrile");
    }
}
