# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

---

## [0.2.0] — 2026-05-24

### Added
- **Confidence score** (`confidence: f64` in `ResolveResult`): 0.0–1.0 signal of resolution quality.
  - `1.00` — exact DirectSmiles dictionary match
  - `0.95` — exact CanonicalName dictionary match
  - `0.90` — match after normalization (CJK input, full-width, Greek, etc.)
  - `0.85` — IUPAC systematic-name parser
- **SMILES → IUPAC reverse conversion** (`smiles_to_iupac(smiles: &str) -> Result<String, ResolveError>`): generates a systematic IUPAC name from a SMILES string. Scope: straight-chain acyclic molecules only; branched chains return `ParseError`.
  - New: `src/parser/smiles_parser.rs` — organic-subset SMILES → `MolGraph`
  - New: `src/parser/iupac_namer.rs` — `MolGraph` → IUPAC name (bidirectional chain scan, lowest-locant-set rule)
- **Standard InChI generation** (`smiles_to_inchi()`, `smiles_to_inchikey()`): formula, connection, and hydrogen layers. InChIKey is a 27-character SHA-256-based hash.
  - New: `src/parser/inchi.rs`
  - `ResolveResult` now includes `inchi: Option<String>` and `inchi_key: Option<String>` (populated for parser and CanonicalName paths; `None` for DirectSmiles entries)
- **sha2 dependency** (`sha2 = "0.10"`, WASM-compatible via `force-soft` feature) for InChIKey generation
- **Branched alkyl substituents**: `isopropyl-`, `tert-butyl-`, `sec-butyl-`, `isobutyl-` substituents in the IUPAC parser (`SubstituentKind` variants + `apply_substituent` graph construction)
- **`cyclo-` prefix**: parser now builds ring-closure bonds for `cycloheptane`–`cyclodecane` (C7–C10); C3–C6 remain in the dictionary
- **Chinese/kanji chemical name dictionary**: `HANZI_TO_IUPAC` (20 systematic names: 甲烷, 乙醇, 丙酮, …) and `HANZI_TO_SMILES` (19 common names: 苯, 水, 氯仿, …) in `src/dict/word_map.rs`
- **Canonical SMILES**: `src/parser/smiles.rs` rewritten with subtree-signature DFS ordering — children sorted ascending by `subtree_sig`; the lex-largest child becomes the main chain (no parentheses); others become branches. Ring molecules use index-descending order with ring-closure notation (`C1CCCCC1` etc.)
- **Python bindings**: `src/py.rs` + `pyproject.toml` via PyO3/Maturin; exports `resolve_to_smiles`, `resolve_full`, `normalize_name`, `resolve_batch`
- **Comprehensive documentation**: all public API items now have `///` doc comments, `# Examples` with doctests, `# Errors`/`# Panics` sections, and `//!` module-level descriptions

### Changed
- `ResolveResult` struct gains three new public fields: `confidence: f64`, `inchi: Option<String>`, `inchi_key: Option<String>`
- SMILES output is now deterministic (canonical) for acyclic molecules. Test expectations updated throughout `tests/iupac_corpus.rs`
  - Alcohols: `propan-2-ol` → `CC(C)O` (was `CC(O)C`)
  - Carboxylic acids: `ethanoic acid` → `CC(=O)O` (was `CC(O)=O`)
  - Heteroatom substituents now appear last: `2-aminobutane` → `CC(CC)N`, `3-mercaptopentane` → `CCC(CC)S`
  - Multi-substituent chains: `2,3-dichlorobutane` → `CC(C(C)Cl)Cl`
- Test count: 74 unit + 39 integration + 10 doctest + 10 proptest = **133 total**

### Fixed
- `HANZI_TO_IUPAC` / `HANZI_TO_SMILES` wired into `dict::lookup_synonym` (two-tier lookup)
- Ring closure pre-pass uses `seen_pairs` to prevent double-registration of the same undirected bond

---

## [0.1.0] — 2026-05-22

### Added
- **Amines / aromatics / cyclic dictionary** (Phase 3.10): methylamine, dimethylamine, trimethylamine, diethylamine, triethylamine, aniline (TEA); phenol, anisole, styrene, o/m/p-xylene, mesitylene; cyclohexane, cyclohexanol, cyclohexanone, cyclopentane/ol, cyclobutane, cyclopropane
- **CLI binary** (`src/bin/chem.rs`): `chem resolve <name>` / `chem resolve --smiles <name>` (Phase 3.9)
- **WASM `resolve_full`** export returning JSON string with smiles + formula + weight + canonical_name (Phase 3.9)
- **Lab abbreviations dictionary** (Phase 3.9): MeOH, EtOH, DCM, DMSO, DMF, THF, MeCN (+ full names)
- **Branched alkane dictionary** (Phase 3.9): isopentane/2-methylbutane, isohexane/2-methylpentane; MEK/butanone, ethyl/methyl acetate, methyl/ethyl formate
- **`-amide` suffix** (Phase 3.8): ethanamide → `CC(=O)N`
- **`acetyl-` / `formyl-` substituents** (Phase 3.8)
- **`-ic acid` trivial names** (Phase 3.8): propionic, butyric, valeric, caproic acid
- **Nitro compounds** (Phase 3.8): nitromethane, nitroethane, nitrobenzene (DirectSmiles)
- **`resolve_batch()` API** (Phase 3.8)
- **iso/sec/tert aliases + halomethane dictionary** (Phase 3.8)
- **Heteroatom substituents** (Phase 3.6): `amino-`, `mercapto-`, `cyano-`
- **Heteroatom suffixes** (Phase 3.6): `-amine`, `-thiol`, `-nitrile`
- **Molecular formula & weight** in `ResolveResult` (Phase 3.5)
- **C11–C20 chain stems** (Phase 3.5)
- **`serde::Serialize`** on `ResolveResult` / `ResolveSource` (Phase 3.7)
- **`n-` prefix stripping** (Phase 3.7)
- **Multiplier prefix stripping** for substituents (`di-`/`tri-`/`tetra-`) (Phase 3.7)
- **Propyl/Butyl/Pentyl/Hexyl substituents** (Phase 3.7)
- **IUPAC parser MVP** (Phase 3): chain stems C1–C10, suffixes (`-ane/-ene/-yne/-ol/-one/-al/-oic acid`), halogen/alkyl substituents, locant lists
- **In-memory dictionary** (Phase 2): PHF synonym tables + katakana → IUPAC
- **CJK normalizer** (Phase 1): fullwidth→halfwidth, `ー`→`-`, Greek letters, whitespace folding
