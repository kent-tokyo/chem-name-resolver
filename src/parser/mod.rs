pub mod formula;
pub mod scanner;
pub mod smiles;

mod alkane;
mod locant;
mod substituent;
mod suffix;

use crate::error::ResolveError;
use suffix::{Suffix, SuffixGroup};
use substituent::parse_substituents;

// ── Molecular graph types ─────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
pub enum Element {
    C,
    H,
    O,
    N,
    S,
    P,
    F,
    Cl,
    Br,
    I,
}

impl Element {
    pub fn symbol(&self) -> &'static str {
        match self {
            Element::C => "C",
            Element::H => "H",
            Element::O => "O",
            Element::N => "N",
            Element::S => "S",
            Element::P => "P",
            Element::F => "F",
            Element::Cl => "Cl",
            Element::Br => "Br",
            Element::I => "I",
        }
    }

    pub fn valence(&self) -> u8 {
        match self {
            Element::C => 4,
            Element::N => 3,
            Element::O | Element::S => 2,
            Element::F | Element::Cl | Element::Br | Element::I | Element::H | Element::P => 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

impl BondOrder {
    pub fn degree(&self) -> u8 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub element: Element,
    pub charge: i8,
    pub implicit_h: u8,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub to: usize,
    pub order: BondOrder,
}

/// Molecular graph: atoms with adjacency lists.
#[derive(Debug, Clone, Default)]
pub struct MolGraph {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Vec<Bond>>,
}

impl MolGraph {
    fn add_atom(&mut self, element: Element) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(Atom { element, charge: 0, implicit_h: 0 });
        self.bonds.push(Vec::new());
        idx
    }

    fn add_bond(&mut self, a: usize, b: usize, order: BondOrder) {
        self.bonds[a].push(Bond { to: b, order: order.clone() });
        self.bonds[b].push(Bond { to: a, order });
    }

    fn fill_implicit_h(&mut self) {
        for i in 0..self.atoms.len() {
            let used: u8 = self.bonds[i].iter().map(|b| b.order.degree()).sum();
            let valence = self.atoms[i].element.valence();
            self.atoms[i].implicit_h = valence.saturating_sub(used);
        }
    }
}

// ── Entry point ───────────────────────────────────────────────────────────────

/// Parse an IUPAC systematic name (already normalized, lowercase) into a MolGraph.
/// MVP scope: straight-chain and monocyclic C1-C20, basic suffixes, simple substituents.
pub fn parse_iupac(name: &str) -> Result<MolGraph, ResolveError> {
    // Strip the "n-" prefix (e.g. "n-butane" = "butane").
    let name = name.trim().strip_prefix("n-").unwrap_or(name.trim());

    // Detect cyclo- prefix (e.g. "cyclohexane" → ring of 6 carbons).
    let (is_cyclic, name) = if let Some(rest) = name.strip_prefix("cyclo") {
        (true, rest)
    } else {
        (false, name)
    };

    // 1. Parse substituent prefix list ("2-chloro-", "3-methyl-", …)
    let (substituents, rest) = parse_substituents(name);

    // 2. Extract prefix locants that belong to the suffix (e.g. "2,4-" in "2,4-pentanedione").
    //    Distinguished from substituent locants because no substituent name follows them.
    let (prefix_suffix_locants, rest) = extract_prefix_suffix_locants(rest);

    // 3. Parse chain stem
    let (chain_len, rest) = alkane::parse_stem(rest).ok_or_else(|| ResolveError::ParseError {
        pos: 0,
        msg: format!("unrecognized chain stem in: {name:?}"),
    })?;

    // 4. Parse one or more suffix groups (handles "an" connector and compound suffixes)
    let (suffix_groups, remaining) =
        parse_suffix_groups(rest).map_err(|_| ResolveError::ParseError {
            pos: name.len() - rest.len(),
            msg: format!("unrecognized suffix in: {name:?}"),
        })?;

    if !remaining.is_empty() {
        return Err(ResolveError::ParseError {
            pos: name.len() - remaining.len(),
            msg: format!("unexpected trailing text: {remaining:?}"),
        });
    }

    build_graph(chain_len, &suffix_groups, &prefix_suffix_locants, &substituents, name, is_cyclic)
}

/// Extract prefix locants that precede the stem (e.g. "2,4-" in "2,4-pentanedione").
fn extract_prefix_suffix_locants<'a>(input: &'a str) -> (Vec<u8>, &'a str) {
    if let Some((locs, rest)) = locant::parse_locant_list(input) {
        if alkane::parse_stem(rest).is_some() {
            return (locs, rest);
        }
    }
    (vec![], input)
}

/// Parse one or more suffix groups from `input`, handling the "an" connector.
fn parse_suffix_groups(input: &str) -> Result<(Vec<SuffixGroup>, &str), ()> {
    let mut groups = Vec::new();
    let mut rest = input;

    loop {
        // Try "an" connector FIRST to avoid "anedione" being parsed as "ane" + "dione".
        // "ane" as a suffix falls through to the direct parse below because stripping "an"
        // from "ane" leaves "e" which strip_elision_e removes to "", and parse_suffix("")
        // returns None, so we fall through.
        if let Some(r) = rest.strip_prefix("an") {
            let r = strip_elision_e(r);
            if let Some((sg, r2)) = suffix::parse_suffix(r) {
                groups.push(sg);
                rest = r2;
                if rest.is_empty() {
                    break;
                }
                continue;
            }
        }

        // Direct parse (for -ane, -ene, -yne without "an" connector, and infix locant forms)
        if let Some((sg, r)) = suffix::parse_suffix(rest) {
            groups.push(sg);
            rest = r;
            if rest.is_empty() {
                break;
            }
            continue;
        }

        break;
    }

    if groups.is_empty() {
        Err(())
    } else {
        Ok((groups, rest))
    }
}

/// Strip the elision 'e' (from "-ane") when the suffix starts with a consonant.
/// E.g. "edione" → "dione", but "ene" is NOT stripped (it is a suffix itself).
fn strip_elision_e(input: &str) -> &str {
    if input.starts_with('e')
        && !input.starts_with("ene")
        && !input.starts_with("en-")
    {
        &input[1..]
    } else {
        input
    }
}

// ── Graph construction ────────────────────────────────────────────────────────

fn build_graph(
    chain_len: u8,
    suffix_groups: &[SuffixGroup],
    prefix_suffix_locants: &[u8],
    substituents: &[substituent::Substituent],
    name: &str,
    is_cyclic: bool,
) -> Result<MolGraph, ResolveError> {
    let mut g = MolGraph::default();

    // Build carbon chain (C-1 = index 0)
    let carbon_indices: Vec<usize> =
        (0..chain_len as usize).map(|_| g.add_atom(Element::C)).collect();

    for i in 0..carbon_indices.len().saturating_sub(1) {
        g.add_bond(carbon_indices[i], carbon_indices[i + 1], BondOrder::Single);
    }

    // For cyclo- compounds, close the ring by bonding C-1 to C-n
    if is_cyclic {
        let n = carbon_indices.len();
        if n >= 3 {
            g.add_bond(carbon_indices[0], carbon_indices[n - 1], BondOrder::Single);
        }
    }

    // Apply suffix groups
    for (sg_idx, sg) in suffix_groups.iter().enumerate() {
        // Use prefix locants only for the first suffix group when it has none of its own
        let effective_locants: Vec<u8> = if sg.locants.is_empty() && sg_idx == 0 && !prefix_suffix_locants.is_empty() {
            prefix_suffix_locants.to_vec()
        } else {
            sg.locants.clone()
        };
        apply_suffix(&mut g, &carbon_indices, sg, &effective_locants, name)?;
    }

    // Apply substituents
    for sub in substituents {
        apply_substituent(&mut g, &carbon_indices, sub, name)?;
    }

    g.fill_implicit_h();
    Ok(g)
}

fn apply_suffix(
    g: &mut MolGraph,
    carbons: &[usize],
    sg: &SuffixGroup,
    effective_locants: &[u8],
    name: &str,
) -> Result<(), ResolveError> {
    let count = sg.multiplier.as_ref().map(|m| m.count()).unwrap_or(1) as usize;

    // Resolve locants to 0-based indices. Default rules when unspecified:
    let locants_0: Vec<usize> = if effective_locants.is_empty() {
        match sg.suffix {
            Suffix::Ane => vec![],
            Suffix::Ene | Suffix::Yne => (0..count).map(|i| i * 2).collect(),
            Suffix::Ol => (0..count).map(|i| i).collect(),
            Suffix::One => {
                if carbons.len() >= 3 {
                    (0..count).map(|i| i + 1).collect()
                } else {
                    vec![0]
                }
            }
            Suffix::Al | Suffix::OicAcid => {
                if count == 1 {
                    vec![0]
                } else {
                    // dioic acid: both terminal carbons
                    vec![0, carbons.len() - 1]
                }
            }
            // These suffixes handle locants internally; default is unused.
            Suffix::Amine | Suffix::Thiol | Suffix::Nitrile | Suffix::Amide => vec![],
        }
    } else {
        effective_locants.iter().map(|&l| l as usize - 1).collect()
    };

    match sg.suffix {
        Suffix::Ane => {}
        Suffix::Ene => {
            for &ci in &locants_0 {
                validate_bond_locant(ci, carbons.len(), "ene", name)?;
                upgrade_bond(g, carbons[ci], carbons[ci + 1], BondOrder::Double);
            }
        }
        Suffix::Yne => {
            for &ci in &locants_0 {
                validate_bond_locant(ci, carbons.len(), "yne", name)?;
                upgrade_bond(g, carbons[ci], carbons[ci + 1], BondOrder::Triple);
            }
        }
        Suffix::Ol => {
            let indices: Vec<usize> = if locants_0.is_empty() {
                vec![0]
            } else {
                locants_0
            };
            for ci in indices {
                validate_atom_locant(ci, carbons.len(), "ol", name)?;
                let oidx = g.add_atom(Element::O);
                g.add_bond(carbons[ci], oidx, BondOrder::Single);
            }
        }
        Suffix::One => {
            for &ci in &locants_0 {
                validate_atom_locant(ci, carbons.len(), "one", name)?;
                let oidx = g.add_atom(Element::O);
                g.add_bond(carbons[ci], oidx, BondOrder::Double);
            }
        }
        Suffix::Al => {
            let oidx = g.add_atom(Element::O);
            g.add_bond(carbons[0], oidx, BondOrder::Double);
        }
        Suffix::OicAcid => {
            let positions: Vec<usize> = if locants_0.is_empty() {
                if count == 1 {
                    vec![0]
                } else {
                    vec![0, carbons.len() - 1]
                }
            } else {
                locants_0
            };
            for ci in positions {
                validate_atom_locant(ci, carbons.len(), "oic acid", name)?;
                let oidx = g.add_atom(Element::O);
                let ohidx = g.add_atom(Element::O);
                g.add_bond(carbons[ci], oidx, BondOrder::Double);
                g.add_bond(carbons[ci], ohidx, BondOrder::Single);
            }
        }
        Suffix::Amine => {
            let indices: Vec<usize> = if locants_0.is_empty() {
                vec![carbons.len() - 1]
            } else {
                locants_0
            };
            for ci in indices {
                validate_atom_locant(ci, carbons.len(), "amine", name)?;
                let nidx = g.add_atom(Element::N);
                g.add_bond(carbons[ci], nidx, BondOrder::Single);
            }
        }
        Suffix::Thiol => {
            let indices: Vec<usize> = if locants_0.is_empty() {
                vec![carbons.len() - 1]
            } else {
                locants_0
            };
            for ci in indices {
                validate_atom_locant(ci, carbons.len(), "thiol", name)?;
                let sidx = g.add_atom(Element::S);
                g.add_bond(carbons[ci], sidx, BondOrder::Single);
            }
        }
        Suffix::Nitrile => {
            // -nitrile: C-1 gains a triple bond to N (the carbon IS part of the chain).
            let cidx = carbons[0];
            let nidx = g.add_atom(Element::N);
            g.add_bond(cidx, nidx, BondOrder::Triple);
        }
        Suffix::Amide => {
            // -amide: C-1 gets -NH2 (single) and =O (double).
            // N added before O so DFS places =O as branch and N as continuation → "C(=O)N".
            let ci = if locants_0.is_empty() { 0 } else { locants_0[0] };
            validate_atom_locant(ci, carbons.len(), "amide", name)?;
            let nidx = g.add_atom(Element::N);
            let oidx = g.add_atom(Element::O);
            g.add_bond(carbons[ci], nidx, BondOrder::Single);
            g.add_bond(carbons[ci], oidx, BondOrder::Double);
        }
    }
    Ok(())
}

fn apply_substituent(
    g: &mut MolGraph,
    carbons: &[usize],
    sub: &substituent::Substituent,
    name: &str,
) -> Result<(), ResolveError> {
    use substituent::SubstituentKind;
    for &loc in &sub.locants {
        let ci = loc as usize - 1;
        validate_atom_locant(ci, carbons.len(), "substituent", name)?;
        let cidx = carbons[ci];
        match &sub.kind {
            SubstituentKind::Oxo => {
                let oidx = g.add_atom(Element::O);
                g.add_bond(cidx, oidx, BondOrder::Double);
            }
            SubstituentKind::Hydroxy => {
                let oidx = g.add_atom(Element::O);
                g.add_bond(cidx, oidx, BondOrder::Single);
            }
            SubstituentKind::Chloro => {
                let x = g.add_atom(Element::Cl);
                g.add_bond(cidx, x, BondOrder::Single);
            }
            SubstituentKind::Bromo => {
                let x = g.add_atom(Element::Br);
                g.add_bond(cidx, x, BondOrder::Single);
            }
            SubstituentKind::Fluoro => {
                let x = g.add_atom(Element::F);
                g.add_bond(cidx, x, BondOrder::Single);
            }
            SubstituentKind::Iodo => {
                let x = g.add_atom(Element::I);
                g.add_bond(cidx, x, BondOrder::Single);
            }
            SubstituentKind::Methyl => {
                let m = g.add_atom(Element::C);
                g.add_bond(cidx, m, BondOrder::Single);
            }
            SubstituentKind::Ethyl => {
                let m1 = g.add_atom(Element::C);
                let m2 = g.add_atom(Element::C);
                g.add_bond(cidx, m1, BondOrder::Single);
                g.add_bond(m1, m2, BondOrder::Single);
            }
            SubstituentKind::Propyl
            | SubstituentKind::Butyl
            | SubstituentKind::Pentyl
            | SubstituentKind::Hexyl => {
                let chain_len = match &sub.kind {
                    SubstituentKind::Propyl => 3,
                    SubstituentKind::Butyl => 4,
                    SubstituentKind::Pentyl => 5,
                    SubstituentKind::Hexyl => 6,
                    _ => unreachable!(),
                };
                let mut prev = cidx;
                for _ in 0..chain_len {
                    let m = g.add_atom(Element::C);
                    g.add_bond(prev, m, BondOrder::Single);
                    prev = m;
                }
            }
            // -CH(CH3)2 : branch carbon + two methyls
            SubstituentKind::Isopropyl => {
                let branch = g.add_atom(Element::C);
                let me1 = g.add_atom(Element::C);
                let me2 = g.add_atom(Element::C);
                g.add_bond(cidx, branch, BondOrder::Single);
                g.add_bond(branch, me1, BondOrder::Single);
                g.add_bond(branch, me2, BondOrder::Single);
            }
            // -C(CH3)3 : quaternary carbon + three methyls
            SubstituentKind::TertButyl => {
                let branch = g.add_atom(Element::C);
                let me1 = g.add_atom(Element::C);
                let me2 = g.add_atom(Element::C);
                let me3 = g.add_atom(Element::C);
                g.add_bond(cidx, branch, BondOrder::Single);
                g.add_bond(branch, me1, BondOrder::Single);
                g.add_bond(branch, me2, BondOrder::Single);
                g.add_bond(branch, me3, BondOrder::Single);
            }
            // -CH(CH3)CH2CH3 : branch carbon + methyl + ethyl
            SubstituentKind::SecButyl => {
                let branch = g.add_atom(Element::C);
                let me = g.add_atom(Element::C);
                let et1 = g.add_atom(Element::C);
                let et2 = g.add_atom(Element::C);
                g.add_bond(cidx, branch, BondOrder::Single);
                g.add_bond(branch, me, BondOrder::Single);
                g.add_bond(branch, et1, BondOrder::Single);
                g.add_bond(et1, et2, BondOrder::Single);
            }
            // -CH2CH(CH3)2 : methylene + isopropyl
            SubstituentKind::IsoButyl => {
                let ch2 = g.add_atom(Element::C);
                let branch = g.add_atom(Element::C);
                let me1 = g.add_atom(Element::C);
                let me2 = g.add_atom(Element::C);
                g.add_bond(cidx, ch2, BondOrder::Single);
                g.add_bond(ch2, branch, BondOrder::Single);
                g.add_bond(branch, me1, BondOrder::Single);
                g.add_bond(branch, me2, BondOrder::Single);
            }
            SubstituentKind::Amino => {
                let nidx = g.add_atom(Element::N);
                g.add_bond(cidx, nidx, BondOrder::Single);
            }
            SubstituentKind::Mercapto => {
                let sidx = g.add_atom(Element::S);
                g.add_bond(cidx, sidx, BondOrder::Single);
            }
            SubstituentKind::Cyano => {
                // cyano- = -C≡N branch attached to chain carbon
                let cbranch = g.add_atom(Element::C);
                let nidx = g.add_atom(Element::N);
                g.add_bond(cidx, cbranch, BondOrder::Single);
                g.add_bond(cbranch, nidx, BondOrder::Triple);
            }
            SubstituentKind::Acetyl => {
                // acetyl- = -C(=O)CH3: add methyl C first so DFS writes "C(=O)C"
                let carbonyl = g.add_atom(Element::C);
                let methyl = g.add_atom(Element::C);
                let o = g.add_atom(Element::O);
                g.add_bond(cidx, carbonyl, BondOrder::Single);
                g.add_bond(carbonyl, methyl, BondOrder::Single);
                g.add_bond(carbonyl, o, BondOrder::Double);
            }
            SubstituentKind::Formyl => {
                // formyl- = -CHO: one carbon double-bonded to O
                let carbonyl = g.add_atom(Element::C);
                let o = g.add_atom(Element::O);
                g.add_bond(cidx, carbonyl, BondOrder::Single);
                g.add_bond(carbonyl, o, BondOrder::Double);
            }
        }
    }
    Ok(())
}

fn validate_bond_locant(ci: usize, len: usize, tag: &str, name: &str) -> Result<(), ResolveError> {
    if ci + 1 >= len {
        Err(ResolveError::ParseError {
            pos: 0,
            msg: format!("{tag} locant {ci} out of range for {len}-carbon chain in {name:?}"),
        })
    } else {
        Ok(())
    }
}

fn validate_atom_locant(ci: usize, len: usize, tag: &str, name: &str) -> Result<(), ResolveError> {
    if ci >= len {
        Err(ResolveError::ParseError {
            pos: 0,
            msg: format!("{tag} locant {ci} out of range for {len}-carbon chain in {name:?}"),
        })
    } else {
        Ok(())
    }
}

fn upgrade_bond(g: &mut MolGraph, a: usize, b: usize, new_order: BondOrder) {
    for bond in &mut g.bonds[a] {
        if bond.to == b {
            bond.order = new_order.clone();
        }
    }
    for bond in &mut g.bonds[b] {
        if bond.to == a {
            bond.order = new_order.clone();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::smiles::to_smiles;

    fn smiles(name: &str) -> String {
        to_smiles(&parse_iupac(name).unwrap_or_else(|e| panic!("{name}: {e}")))
    }

    #[test]
    fn methane() {
        assert_eq!(smiles("methane"), "C");
    }

    #[test]
    fn ethane() {
        assert_eq!(smiles("ethane"), "CC");
    }

    #[test]
    fn propane() {
        assert_eq!(smiles("propane"), "CCC");
    }

    #[test]
    fn ethanol() {
        // DFS from C0(=C1): C0→C1→O → "CCO"
        assert_eq!(smiles("ethanol"), "CCO");
    }

    #[test]
    fn propan_2_one() {
        // DFS: C0→C1 (branch O at higher idx)→C2 → "CC(=O)C"
        assert_eq!(smiles("propan-2-one"), "CC(=O)C");
    }

    #[test]
    fn but_2_yne() {
        assert_eq!(smiles("but-2-yne"), "CC#CC");
    }

    #[test]
    fn two_four_pentanedione() {
        assert_eq!(smiles("2,4-pentanedione"), "CC(=O)CC(=O)C");
    }
}
