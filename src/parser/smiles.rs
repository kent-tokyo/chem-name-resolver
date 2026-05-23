use super::{BondOrder, Element, MolGraph};

/// Convert a MolGraph to a canonical SMILES string.
///
/// For acyclic molecules: canonical DFS — sort children by subtree-signature
/// (ascending), main chain is the child with the lex-largest signature, others
/// are branches.  Start atom is chosen as the terminal that produces the
/// lex-smallest overall string.
///
/// For cyclic molecules (one ring): index-based DFS with ring-closure notation.
pub fn to_smiles(graph: &MolGraph) -> String {
    if graph.atoms.is_empty() {
        return String::new();
    }

    let n = graph.atoms.len();

    // Detect rings: a tree with n atoms has exactly n-1 bonds.
    let bond_count: usize = graph.bonds.iter().map(|b| b.len()).sum::<usize>() / 2;
    let is_cyclic = bond_count >= n;

    if is_cyclic {
        // Ring: use index-based ordering (ring-closure bonds sorted by index).
        let start = find_ring_start(graph);
        let rc = find_ring_closures(graph, start, n);
        let mut out = String::new();
        let mut visited = vec![false; n];
        write_atom_indexed(graph, start, usize::MAX, &mut visited, &mut out, &rc);
        out
    } else {
        // Acyclic: use canonical key-based branch ordering with the same start
        // atom heuristic as before (terminal C with no heteroatom, lowest index).
        // This preserves the natural chain-first traversal while ensuring the
        // branch ordering is deterministic and algorithm-independent.
        let start = find_acyclic_start(graph);
        let rc_none: Vec<Vec<(u8, usize)>> = vec![Vec::new(); n];
        let mut out = String::new();
        let mut visited = vec![false; n];
        write_atom_canonical(graph, start, usize::MAX, &mut visited, &mut out, &rc_none);
        out
    }
}

// ── Ring-closure pre-pass ─────────────────────────────────────────────────────

fn find_ring_closures(graph: &MolGraph, start: usize, n: usize) -> Vec<Vec<(u8, usize)>> {
    let mut rc = vec![Vec::<(u8, usize)>::new(); n];
    let mut visited = vec![false; n];
    let mut counter = 1u8;
    let mut seen: Vec<(usize, usize)> = Vec::new();
    find_rc_dfs(graph, start, usize::MAX, &mut visited, &mut rc, &mut counter, &mut seen);
    rc
}

fn find_rc_dfs(
    graph: &MolGraph,
    idx: usize,
    parent: usize,
    visited: &mut Vec<bool>,
    rc: &mut Vec<Vec<(u8, usize)>>,
    counter: &mut u8,
    seen: &mut Vec<(usize, usize)>,
) {
    visited[idx] = true;
    let mut neighbors: Vec<_> = graph.bonds[idx].iter().filter(|b| b.to != parent).collect();
    neighbors.sort_by(|a, b| b.to.cmp(&a.to));

    for bond in neighbors {
        let to = bond.to;
        let pair = if idx < to { (idx, to) } else { (to, idx) };
        if seen.iter().any(|&p| p == pair) { continue; }
        seen.push(pair);

        if visited[to] {
            let num = *counter;
            *counter += 1;
            rc[to].push((num, idx));
            rc[idx].push((num, to));
        } else {
            find_rc_dfs(graph, to, idx, visited, rc, counter, seen);
        }
    }
}

// ── Start-atom helpers ────────────────────────────────────────────────────────

fn terminal_carbons(graph: &MolGraph) -> Vec<usize> {
    let terms: Vec<usize> = (0..graph.atoms.len())
        .filter(|&i| {
            graph.atoms[i].element == Element::C
                && graph.bonds[i].iter()
                    .filter(|b| graph.atoms[b.to].element == Element::C)
                    .count() <= 1
        })
        .collect();
    if terms.is_empty() { vec![0] } else { terms }
}

/// For acyclic molecules: prefer a terminal carbon with no heteroatom bonds
/// (plain methyl end). Among ties take the lowest atom index.
fn find_acyclic_start(graph: &MolGraph) -> usize {
    let terms = terminal_carbons(graph);

    // Prefer terminal with no direct heteroatom bonds
    for &tc in &terms {
        let has_heteroatom = graph.bonds[tc]
            .iter()
            .any(|b| graph.atoms[b.to].element != Element::C);
        if !has_heteroatom {
            return tc;
        }
    }

    // All terminals have heteroatoms — pick the one with fewest
    terms.into_iter()
        .min_by_key(|&tc| {
            graph.bonds[tc]
                .iter()
                .filter(|b| graph.atoms[b.to].element != Element::C)
                .count()
        })
        .unwrap_or(0)
}

fn find_ring_start(graph: &MolGraph) -> usize {
    // For ring molecules all ring atoms have ≥2 C-neighbors; fall back to 0.
    terminal_carbons(graph).into_iter().next().unwrap_or(0)
}

// ── Acyclic canonical writer ──────────────────────────────────────────────────

fn write_atom_canonical(
    graph: &MolGraph,
    idx: usize,
    parent: usize,
    visited: &mut Vec<bool>,
    out: &mut String,
    rc: &[Vec<(u8, usize)>],
) {
    visited[idx] = true;
    out.push_str(graph.atoms[idx].element.symbol());

    let rc_targets: Vec<usize> = rc[idx].iter().map(|&(_, o)| o).collect();

    let mut children: Vec<_> = graph.bonds[idx]
        .iter()
        .filter(|b| b.to != parent && !visited[b.to] && !rc_targets.contains(&b.to))
        .collect();

    // Sort ascending by subtree signature.  The lex-largest child (last in the
    // sorted list) becomes the main chain (written without parentheses); every
    // other child is a branch in parentheses.
    children.sort_by(|a, b| {
        let ka = subtree_sig(graph, a.to, idx, a.order.clone(), visited);
        let kb = subtree_sig(graph, b.to, idx, b.order.clone(), visited);
        ka.cmp(&kb)
    });

    for (i, bond) in children.iter().enumerate() {
        let is_branch = i < children.len() - 1;
        if is_branch { out.push('('); }
        match bond.order {
            BondOrder::Double => out.push('='),
            BondOrder::Triple => out.push('#'),
            BondOrder::Single => {}
        }
        write_atom_canonical(graph, bond.to, idx, visited, out, rc);
        if is_branch { out.push(')'); }
    }
}

/// Compute the canonical signature for the subtree rooted at `idx`
/// (reached via a bond of `bond_order` from `parent`).
/// The signature is a string that is lexicographically comparable.
fn subtree_sig(
    graph: &MolGraph,
    idx: usize,
    parent: usize,
    bond_order: BondOrder,
    visited: &[bool],
) -> String {
    let bond_prefix = match bond_order {
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Single => "",
    };
    let symbol = graph.atoms[idx].element.symbol();

    let children: Vec<_> = graph.bonds[idx]
        .iter()
        .filter(|b| b.to != parent && !visited[b.to])
        .collect();

    if children.is_empty() {
        return format!("{}{}", bond_prefix, symbol);
    }

    let mut child_sigs: Vec<String> = children.iter()
        .map(|b| subtree_sig(graph, b.to, idx, b.order.clone(), visited))
        .collect();
    child_sigs.sort();

    // Branches in parens, main chain (lex-largest) appended directly — same
    // layout as the SMILES we will actually write.
    let branches: String = child_sigs[..child_sigs.len() - 1]
        .iter()
        .map(|s| format!("({})", s))
        .collect();
    let main = &child_sigs[child_sigs.len() - 1];
    format!("{}{}{}{}", bond_prefix, symbol, branches, main)
}

// ── Ring (index-ordered) writer ───────────────────────────────────────────────

fn write_atom_indexed(
    graph: &MolGraph,
    idx: usize,
    parent: usize,
    visited: &mut Vec<bool>,
    out: &mut String,
    rc: &[Vec<(u8, usize)>],
) {
    visited[idx] = true;
    out.push_str(graph.atoms[idx].element.symbol());

    for &(num, _) in &rc[idx] {
        push_ring_num(out, num);
    }

    let rc_targets: Vec<usize> = rc[idx].iter().map(|&(_, o)| o).collect();

    let mut children: Vec<_> = graph.bonds[idx]
        .iter()
        .filter(|b| b.to != parent && !visited[b.to] && !rc_targets.contains(&b.to))
        .collect();
    children.sort_by(|a, b| b.to.cmp(&a.to)); // descending index — lowest last (main chain)

    for (i, bond) in children.iter().enumerate() {
        let is_branch = i < children.len() - 1;
        if is_branch { out.push('('); }
        match bond.order {
            BondOrder::Double => out.push('='),
            BondOrder::Triple => out.push('#'),
            BondOrder::Single => {}
        }
        write_atom_indexed(graph, bond.to, idx, visited, out, rc);
        if is_branch { out.push(')'); }
    }
}

fn push_ring_num(out: &mut String, num: u8) {
    if num < 10 {
        out.push((b'0' + num) as char);
    } else {
        out.push('%');
        out.push((b'0' + num / 10) as char);
        out.push((b'0' + num % 10) as char);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{Atom, Bond, BondOrder, Element, MolGraph};

    fn make_ethane() -> MolGraph {
        let atoms = vec![
            Atom { element: Element::C, charge: 0, implicit_h: 3 },
            Atom { element: Element::C, charge: 0, implicit_h: 3 },
        ];
        let bonds = vec![
            vec![Bond { to: 1, order: BondOrder::Single }],
            vec![Bond { to: 0, order: BondOrder::Single }],
        ];
        MolGraph { atoms, bonds }
    }

    #[test]
    fn ethane_smiles() {
        assert_eq!(to_smiles(&make_ethane()), "CC");
    }

    #[test]
    fn empty_graph() {
        assert_eq!(to_smiles(&MolGraph::default()), "");
    }
}
