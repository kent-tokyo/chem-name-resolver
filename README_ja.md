# chem-name-resolver

Pure Rust の IUPAC 化学名 → SMILES 変換ライブラリ。Java の [OPSIN](https://opsin.ch.cam.ac.uk/) に相当する機能を、Pure Rust + WebAssembly 対応で提供します。

## なぜ chem-name-resolver が必要か

`"2,4-pentanedione"` のような IUPAC 名を SMILES `"CC(=O)CC(=O)C"` に変換する機能は単純に見えますが、既存のツールにはそれぞれ大きなトレードオフがあります。

| | [OPSIN](https://opsin.ch.cam.ac.uk/) | [RDKit](https://www.rdkit.org/) | [OpenBabel](https://openbabel.org/) | [CDK](https://cdk.github.io/) | [Indigo](https://lifescience.opensource.epam.com/indigo/) | [PubChem API](https://pubchem.ncbi.nlm.nih.gov/) | [PubChemPy](https://github.com/mcs07/PubChemPy) | [STOUT v2](https://github.com/Kohulan/STOUT) | [ChemCore](https://crates.io/crates/chemcore) | **chem-name-resolver** |
|---|---|---|---|---|---|---|---|---|---|---|
| 言語 | Java | Python/C++ | C++ | Java | C++ | REST | Python | Python/ML | Rust | **Rust** |
| WASM | ✗ | △ | △ | ✗ | ✓ | ✗ | ✗ | ✗ | △ | **✓** |
| オフライン | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | △ | ✓ | **✓** |
| CJK 名 | ✗ | ✗ | ✗ | ✗ | ✗ | △ | ✗ | ✗ | ✗ | **✓** |
| IUPAC パーサー | ✓ (最高水準) | ✗ | ✗ | ✗ | ✗ | 辞書引き | ✗ | ✓ (ニューラル) | ✗ | **✓** |
| ライセンス | MIT | BSD-3 | GPL-2 | LGPL-2.1 | Apache-2 | パブリックドメイン | BSD | MIT | MIT | **MIT/Apache-2** |
| 備考 | JVM 必須 | バイナリ ~50 MB; C++ ビルドチェーン必須; rdkit-js WASM は機能限定 | C++ FFI; コピーレフト; WASM は実験的 | JVM 必須; IUPAC 解析は OPSIN に委譲 | npm 経由の公式 WASM あり; 構造操作のみ | ネットワーク依存; 6700 万化合物以上 | PubChem REST の薄いラッパー | GPU 推奨; 非決定論的; モデル ~GB | 2020 年以降休止; SMILES 不完全 | Pure Rust。ネイティブ依存ゼロ |

△ = 部分対応 / 実験的

**このライブラリが埋めるギャップ**: Pure Rust・WASM 対応・オフライン動作・CJK 対応を同時に満たす IUPAC→SMILES エンジンは現在存在しません。

**具体的なユースケース:**

- **ブラウザサイド化学処理** — WASM モジュールとして配布し、サーバー通信なしにクライアントで化学名を解決
- **Rust ネイティブツール** — JVM や C++ ビルドなしに CLI ツール・データベースインデクサ (例: [Cheminee](https://github.com/rdkit/Cheminee))・Axum サービスへ統合
- **日本語・中国語ワークフロー** — カタカナ・漢字の化学名を前処理なしで同一パイプラインで正規化
- **軽量組み込み** — `opt-level = "s"` + LTO で小さなバイナリを生成。エッジデプロイに適合

## 特徴

- **Pure Rust** — C/C++ 依存なし (RDKit/Boost 不要)
- **WASM 対応** — `wasm32-unknown-unknown` ターゲットでコンパイル可能
- **CJK 対応** — カタカナ化学名 (メタン, エタノール 等) の正規化・解決
- **ゼロコピー正規化** — 既に正規化済みの入力ではアロケーションゼロ
- **JSON 出力** — `ResolveResult` は `serde::Serialize` を実装

## クイックスタート

```rust
use chem_name_resolver::resolve;

// IUPAC 系統名
let r = resolve("propan-2-one").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
assert!((r.molecular_weight.unwrap() - 58.08).abs() < 0.01);

// 慣用名
let r = resolve("acetone").unwrap();
assert_eq!(r.smiles, "CC(=O)C");

// カタカナ名
let r = resolve("メタン").unwrap();
assert_eq!(r.smiles, "C");

// n- 接頭辞付き
let r = resolve("n-butane").unwrap();
assert_eq!(r.smiles, "CCCC");

// JSON シリアライズ
let json = serde_json::to_string(&r).unwrap();
```

## 対応している変換

### Normalizer

| 入力 | 出力 |
|------|------|
| 全角数字・記号 (`２－`) | 半角 (`2-`) |
| 長音符 (`ー`) | ハイフン (`-`) |
| ギリシャ文字 (`α`, `β`, `γ`) | ASCII (`alpha`, `beta`, `gamma`) |
| 連続空白 | 単一スペース |
| `n-` 接頭辞 | 除去 (`n-butane` → `butane`) |

### Dictionary (同義語辞書)

| 種別 | 例 |
|------|----|
| 慣用名 → IUPAC 系統名 | acetone, acetic acid, glycerol, formaldehyde, propionic/butyric/valeric acid, … |
| 慣用名 → 直接 SMILES | water, benzene, toluene, ether, chloroform, aspirin, glucose, caffeine |
| iso/sec/tert 別名 | isopropanol, isobutane, tert-butanol, neopentane, sec-butanol, … |
| 分岐アルカン | isopentane, isohexane (+ IUPAC 系統名別名) |
| 実験室略称 | MeOH, EtOH, DCM, DMSO, DMF, THF, MeCN (正式名も対応) |
| ハロメタン | chloromethane, bromomethane, iodomethane, dibromomethane, … |
| 一般試薬 | ethyl acetate, methyl acetate, MEK (+ 正式名) |
| アミン類 | methylamine, dimethylamine, trimethylamine, aniline, triethylamine, … |
| フェノール・芳香族 | phenol, anisole, styrene, o/m/p-xylene, mesitylene, … |
| 環状化合物 | cyclohexane, cyclohexanol, cyclohexanone, cyclopentane, cyclopropane, … |
| ニトロ化合物 | nitromethane, nitroethane, nitrobenzene |
| カタカナ → IUPAC 系統名 | メタン〜デカン, エタノール, アセトン, ベンゼン, … |

### IUPAC パーサー

**アルカン幹:** methane〜decane (C1–C10)、undecane〜icosane/eicosane (C11–C20)

**接尾辞:**

| 接尾辞 | 官能基 | 例 |
|--------|--------|-----|
| `-ane` | アルカン | ethane → `CC` |
| `-ene` | アルケン | hex-1-ene → `C=CCCCC` |
| `-yne` | アルキン | but-2-yne → `CC#CC` |
| `-ol` / `-diol` | アルコール | propan-2-ol → `CC(C)O` |
| `-one` / `-dione` | ケトン | propan-2-one → `CC(=O)C` |
| `-al` | アルデヒド | pentanal → `CCCCC=O` |
| `-oic acid` / `-dioic acid` | カルボン酸 | ethanoic acid → `CC(=O)O` |
| `-amine` | アミン | ethanamine → `CCN` |
| `-amide` | アミド | ethanamide → `CC(=O)N` |
| `-thiol` | チオール | ethanethiol → `CCS` |
| `-nitrile` | ニトリル | propanenitrile → `CCC#N` |

倍数接頭辞 `di-`, `tri-`, `tetra-` に対応。

**置換基:**

| 置換基 | 元素/基 | 例 |
|--------|---------|-----|
| `chloro-`, `bromo-`, `fluoro-`, `iodo-` | ハロゲン | 2-chlorobutane → `CC(CC)Cl` |
| `methyl-`, `ethyl-`, `propyl-`, `butyl-`, `pentyl-`, `hexyl-` | 直鎖アルキル | 3-methylpentane → `CCC(C)CC` |
| `hydroxy-` | ヒドロキシ基 | — |
| `oxo-` | ケト基 | — |
| `amino-` | アミノ基 | 2-aminobutane → `CC(CC)N` |
| `mercapto-` | メルカプト基 | 3-mercaptopentane → `CCC(CC)S` |
| `cyano-` | シアノ基 | 2-cyanopentane → `CC(C#N)CCC` |
| `acetyl-` | アセチル基 | 3-acetylheptane → `CCC(C(=O)C)CCCC` |
| `formyl-` | ホルミル基 | 3-formylpentane → `CCC(C=O)CC` |

倍数接頭辞 `di-`, `tri-`, `tetra-` に対応 (例: `2,3-dichlorobutane` → `CC(C(C)Cl)Cl`)。

### 出力フィールド

```rust
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,             // Dictionary | Parser
    pub molecular_formula: Option<String>, // Hill 記法 (例: "C2H6O")
    pub molecular_weight: Option<f64>,     // g/mol
}
```

`DirectSmiles` 経由 (benzene 等) の場合、`molecular_formula` / `molecular_weight` は `None`。

## インストール

```toml
[dependencies]
chem-name-resolver = "0.1"

# JSON 出力が必要な場合
serde_json = "1"
```

## ビルド & テスト

```bash
# テスト (75 件)
cargo test

# WASM ビルド確認
rustup target add wasm32-unknown-unknown
cargo build --features wasm --target wasm32-unknown-unknown

# ベンチマーク
cargo bench
```

## WASM 使用例

```javascript
import init, { resolve_to_smiles, resolve_full, normalize_name } from './chem_name_resolver.js';

await init();
console.log(resolve_to_smiles("propan-2-one")); // "CC(=O)C"
console.log(normalize_name("α-D-glucose"));     // "alpha-d-glucose"

// フル結果を JSON 文字列で取得
const json = resolve_full("ethanol");
// '{"smiles":"CCO","canonical_name":"ethanol","source":"Dictionary","molecular_formula":"C2H6O","molecular_weight":46.069}'
```

## CLI 使用例

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

## 既知の制限

- 環状・芳香族化合物はパーサー未対応 (辞書引きのみ)
- 立体化学 (R/S, E/Z) は未対応

## ロードマップ

- [x] 分岐アルキル置換基 (isopropyl, tert-butyl 等)
- [x] `cyclo-` 接頭辞 (環状化合物)
- [x] CLI バイナリ (`chem resolve "ethanol"`)
- [x] 漢字・中国語化学名辞書
- [x] Canonical SMILES (サブツリー署名 DFS 順序)
- [x] Python バインディング (PyO3 / Maturin)

## ライセンス

MIT OR Apache-2.0
