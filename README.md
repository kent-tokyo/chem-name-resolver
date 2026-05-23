# chem-name-resolver

A pure-Rust library for resolving IUPAC chemical names to SMILES strings and molecular graphs. The Rust equivalent of Java's [OPSIN](https://opsin.ch.cam.ac.uk/), with WebAssembly support.

## Why chem-name-resolver?

Converting an IUPAC name like `"2,4-pentanedione"` to its SMILES representation `"CC(=O)CC(=O)C"` sounds simple, but every existing solution comes with a significant trade-off:

| | [OPSIN](https://opsin.ch.cam.ac.uk/) | [RDKit](https://www.rdkit.org/) | [OpenBabel](https://openbabel.org/) | [CDK](https://cdk.github.io/) | [Indigo](https://lifescience.opensource.epam.com/indigo/) | [PubChem API](https://pubchem.ncbi.nlm.nih.gov/) | [PubChemPy](https://github.com/mcs07/PubChemPy) | [STOUT v2](https://github.com/Kohulan/STOUT) | [ChemCore](https://crates.io/crates/chemcore) | **chem-name-resolver** |
|---|---|---|---|---|---|---|---|---|---|---|
| Language | Java | Python/C++ | C++ | Java | C++ | REST | Python | Python/ML | Rust | **Rust** |
| WASM | ✗ | △ | △ | ✗ | ✓ | ✗ | ✗ | ✗ | △ | **✓** |
| Offline | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | △ | ✓ | **✓** |
| CJK names | ✗ | ✗ | ✗ | ✗ | ✗ | △ | ✗ | ✗ | ✗ | **✓** |
| IUPAC Parser | ✓ (best) | ✗ | ✗ | ✗ | ✗ | Lookup | ✗ | ✓ (neural) | ✗ | **✓** |
| License | MIT | BSD-3 | GPL-2 | LGPL-2.1 | Apache-2 | Public domain | BSD | MIT | MIT | **MIT/Apache-2** |
| Notes | JVM required | ~50 MB; C++ toolchain; rdkit-js WASM is subset | C++ FFI; copyleft; WASM experimental | JVM; IUPAC parsing delegates to OPSIN | Official WASM (npm); structure ops only | Network-dependent; 67M+ compounds | Thin REST wrapper | GPU recommended; non-deterministic; model ~GB | Dormant since 2020; incomplete SMILES | Pure Rust; no native deps |

△ = partial / experimental

**This library fills the gap**: a pure-Rust, WASM-compatible, offline IUPAC→SMILES engine with CJK support. It enables:

- **Browser-side chemistry** — ship a WASM module and resolve names client-side with zero server round-trips
- **Rust-native tooling** — integrate into CLI tools, database indexers (e.g. [Cheminee](https://github.com/rdkit/Cheminee)), or Axum services without pulling in a JVM or C++ build
- **Japanese/Chinese workflows** — normalize katakana and kanji chemical names in the same pipeline, without a separate preprocessing step
- **Lightweight embedding** — the `release` profile produces a small binary (`opt-level = "s"`, LTO enabled) suitable for edge deployments

## Features

- **Pure Rust** — no C/C++ dependencies (no RDKit, no Boost)
- **WASM-compatible** — compiles to `wasm32-unknown-unknown`
- **CJK support** — resolves Japanese katakana names (メタン, エタノール, …)
- **Zero-copy normalization** — returns `Cow::Borrowed` when input needs no changes
- **JSON serialization** — `ResolveResult` implements `serde::Serialize`

## Quick Start

```rust
use chem_name_resolver::resolve;

// Systematic IUPAC name
let r = resolve("propan-2-one").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
assert!((r.molecular_weight.unwrap() - 58.08).abs() < 0.01);

// Trivial name
let r = resolve("acetone").unwrap();
assert_eq!(r.smiles, "CC(=O)C");

// Japanese katakana
let r = resolve("メタン").unwrap();
assert_eq!(r.smiles, "C");

// n- prefix
let r = resolve("n-butane").unwrap();
assert_eq!(r.smiles, "CCCC");

// JSON output
let json = serde_json::to_string(&r).unwrap();
```

## Coverage

### Normalizer

| Input | Output |
|-------|--------|
| Fullwidth ASCII (`２－`) | Halfwidth (`2-`) |
| Katakana prolonged sound mark (`ー`) | Hyphen (`-`) |
| Greek letters (`α`, `β`, `γ`) | ASCII (`alpha`, `beta`, `gamma`) |
| Consecutive whitespace | Single space |
| `n-` prefix | Stripped (`n-butane` → `butane`) |

### Dictionary

| Type | Examples |
|------|---------|
| Trivial → IUPAC | acetone, acetic acid, glycerol, formaldehyde, propionic/butyric/valeric acid, … |
| Trivial → SMILES | water, benzene, toluene, ether, chloroform, aspirin, glucose, caffeine |
| iso/sec/tert aliases | isopropanol, isobutane, tert-butanol, neopentane, sec-butanol, … |
| Branched alkanes | isopentane, isohexane (+ IUPAC systematic aliases) |
| Lab abbreviations | MeOH, EtOH, DCM, DMSO, DMF, THF, MeCN (+ full names) |
| Halomethanes | chloromethane, bromomethane, iodomethane, dibromomethane, … |
| Common reagents | ethyl acetate, methyl acetate, MEK (+ full names) |
| Amines | methylamine, dimethylamine, trimethylamine, aniline, triethylamine, … |
| Phenols / aromatics | phenol, anisole, styrene, o/m/p-xylene, mesitylene, … |
| Cyclic compounds | cyclohexane, cyclohexanol, cyclohexanone, cyclopentane, cyclopropane, … |
| Nitro compounds | nitromethane, nitroethane, nitrobenzene |
| Katakana → IUPAC | メタン–デカン, エタノール, アセトン, ベンゼン, … |

### IUPAC Parser

**Chain stems:** methane–decane (C1–C10), undecane–icosane/eicosane (C11–C20)

**Suffixes:**

| Suffix | Functional group | Example |
|--------|-----------------|---------|
| `-ane` | alkane | ethane → `CC` |
| `-ene` | alkene | hex-1-ene → `C=CCCCC` |
| `-yne` | alkyne | but-2-yne → `CC#CC` |
| `-ol` / `-diol` | alcohol | propan-2-ol → `CC(C)O` |
| `-one` / `-dione` | ketone | propan-2-one → `CC(=O)C` |
| `-al` | aldehyde | pentanal → `CCCCC=O` |
| `-oic acid` / `-dioic acid` | carboxylic acid | ethanoic acid → `CC(=O)O` |
| `-amine` | amine | ethanamine → `CCN` |
| `-amide` | amide | ethanamide → `CC(=O)N` |
| `-thiol` | thiol | ethanethiol → `CCS` |
| `-nitrile` | nitrile | propanenitrile → `CCC#N` |

Multiplier prefixes `di-`, `tri-`, `tetra-` are supported for all suffixes.

**Substituents:**

| Substituent | Atom/group | Example |
|------------|-----------|---------|
| `chloro-`, `bromo-`, `fluoro-`, `iodo-` | halogens | 2-chlorobutane → `CC(CC)Cl` |
| `methyl-`, `ethyl-`, `propyl-`, `butyl-`, `pentyl-`, `hexyl-` | n-alkyl chains | 3-methylpentane → `CCC(C)CC` |
| `hydroxy-` | –OH | — |
| `oxo-` | =O | — |
| `amino-` | –NH₂ | 2-aminobutane → `CC(CC)N` |
| `mercapto-` | –SH | 3-mercaptopentane → `CCC(CC)S` |
| `cyano-` | –C≡N | 2-cyanopentane → `CC(C#N)CCC` |
| `acetyl-` | –C(=O)CH₃ | 3-acetylheptane → `CCC(C(=O)C)CCCC` |
| `formyl-` | –CHO | 3-formylpentane → `CCC(C=O)CC` |

Multiplier prefixes `di-`, `tri-`, `tetra-` are supported (e.g. `2,3-dichlorobutane` → `CC(C(C)Cl)Cl`).

### Output

```rust
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,           // Dictionary | Parser
    pub molecular_formula: Option<String>, // Hill notation, e.g. "C2H6O"
    pub molecular_weight: Option<f64>,     // g/mol
}
```

`molecular_formula` and `molecular_weight` are `None` when resolved via `DirectSmiles` (e.g. benzene).

## Installation

```toml
[dependencies]
chem-name-resolver = "0.1"

# for JSON output
serde_json = "1"
```

## Building & Testing

```bash
# run all 75 tests
cargo test

# verify WASM build
rustup target add wasm32-unknown-unknown
cargo build --features wasm --target wasm32-unknown-unknown

# benchmarks
cargo bench
```

## WASM Usage

```javascript
import init, { resolve_to_smiles, resolve_full, normalize_name } from './chem_name_resolver.js';

await init();
console.log(resolve_to_smiles("propan-2-one")); // "CC(=O)C"
console.log(normalize_name("α-D-glucose"));     // "alpha-d-glucose"

// Full result as JSON string
const json = resolve_full("ethanol");
// '{"smiles":"CCO","canonical_name":"ethanol","source":"Dictionary","molecular_formula":"C2H6O","molecular_weight":46.069}'
```

## CLI Usage

```bash
cargo install chem-name-resolver --features cli

chem resolve ethanol
# {
#   "smiles": "CCO",
#   "canonical_name": "ethanol",
#   "source": "Dictionary",
#   "molecular_formula": "C2H6O",
#   "molecular_weight": 46.069
# }

chem resolve --smiles "propan-2-one"
# CC(=O)C
```

## Known Limitations

- Cyclic and aromatic compounds are not parsed (dictionary lookup only)
- Stereochemistry (R/S, E/Z) is not supported

## Roadmap

- [x] Branched alkyl substituents (isopropyl, tert-butyl, …)
- [x] `cyclo-` prefix (cyclic compounds)
- [x] CLI binary (`chem resolve "ethanol"`)
- [x] Chinese/kanji chemical name dictionary
- [x] Canonical SMILES (subtree-signature DFS ordering)
- [x] Python bindings (PyO3 / Maturin)

## License

MIT OR Apache-2.0
