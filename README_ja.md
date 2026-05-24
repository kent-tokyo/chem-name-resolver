# chem-name-resolver

Pure Rust の IUPAC 化学名 → SMILES 変換ライブラリ。Java の [OPSIN](https://opsin.ch.cam.ac.uk/) に相当する機能を、Pure Rust + WebAssembly 対応で提供します。SMILES 文字列、IUPAC 名、InChI 識別子、分子特性の解決をサポートします。

## なぜ chem-name-resolver が必要か

`"2,4-pentanedione"` のような IUPAC 名を SMILES `"CC(=O)CC(=O)C"` に変換する機能は単純に見えますが、既存のツールにはそれぞれ大きなトレードオフがあります。

| | [OPSIN](https://opsin.ch.cam.ac.uk/) | [RDKit](https://www.rdkit.org/) | [OpenBabel](https://openbabel.org/) | [CDK](https://cdk.github.io/) | [Indigo](https://lifescience.opensource.epam.com/indigo/) | [PubChem API](https://pubchem.ncbi.nlm.nih.gov/) | [PubChemPy](https://github.com/mcs07/PubChemPy) | [STOUT v2](https://github.com/Kohulan/STOUT) | [ChemCore](https://crates.io/crates/chemcore) | **chem-name-resolver** |
|---|---|---|---|---|---|---|---|---|---|---|
| 言語 | Java | Python/C++ | C++ | Java | C++ | REST | Python | Python/ML | Rust | **Rust** |
| WASM | ✗ | △ | △ | ✗ | ✓ | ✗ | ✗ | ✗ | △ | **✓** |
| オフライン | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | △ | ✓ | **✓** |
| CJK 名 | ✗ | ✗ | ✗ | ✗ | ✗ | △ | ✗ | ✗ | ✗ | **✓** |
| IUPAC パーサー | ✓ (最高水準) | ✗ | ✗ | ✗ | ✗ | 辞書引き | ✗ | ✓ (ニューラル) | ✗ | **✓** |
| 逆変換 (SMILES→IUPAC) | ✗ | △ | △ | ✗ | ✗ | △ | ✗ | ✓ (ニューラル) | ✗ | **✓** |
| ライセンス | MIT | BSD-3 | GPL-2 | LGPL-2.1 | Apache-2 | パブリックドメイン | BSD | MIT | MIT | **MIT/Apache-2** |
| 備考 | JVM 必須 | バイナリ ~50 MB; C++ ビルドチェーン必須; rdkit-js WASM は機能限定 | C++ FFI; コピーレフト; WASM は実験的 | JVM 必須; IUPAC 解析は OPSIN に委譲 | npm 経由の公式 WASM あり; 構造操作のみ | ネットワーク依存; 6700 万化合物以上 | PubChem REST の薄いラッパー | GPU 推奨; 非決定論的; モデル ~GB | 2020 年以降休止; SMILES 不完全 | Pure Rust。ネイティブ依存ゼロ |

△ = 部分対応 / 実験的

**このライブラリが埋めるギャップ**: Pure Rust・WASM 対応・オフライン動作・CJK 対応を同時に満たす IUPAC↔SMILES エンジンは現在存在しません。

**具体的なユースケース:**

- **ブラウザサイド化学処理** — WASM モジュールとして配布し、サーバー通信なしにクライアントで化学名を解決
- **Rust ネイティブツール** — JVM や C++ ビルドなしに CLI ツール・データベースインデクサ (例: [Cheminee](https://github.com/rdkit/Cheminee))・Axum サービスへ統合
- **日本語・中国語ワークフロー** — カタカナ・漢字の化学名を前処理なしで同一パイプラインで正規化
- **軽量組み込み** — `opt-level = "s"` + LTO で小さなバイナリを生成。エッジデプロイに適合

## 特徴

- **Pure Rust** — C/C++ 依存なし (RDKit/Boost 不要)
- **WASM 対応** — `wasm32-unknown-unknown` ターゲットでコンパイル可能
- **CJK 対応** — カタカナ化学名 (メタン, エタノール 等) および漢字 (甲烷, 乙醇 等) の正規化・解決
- **双方向変換**: IUPAC→SMILES (`resolve`) **および** SMILES→IUPAC (`smiles_to_iupac`)
- **InChI/InChIKey** 生成 (SMILES から)
- **信頼度スコア** — 全結果に `0.0〜1.0` の品質シグナルを付与
- **ゼロコピー正規化** — 既に正規化済みの入力では `Cow::Borrowed` を返却
- **JSON 出力** — `ResolveResult` は `serde::Serialize` を実装

## クイックスタート

```rust
use chem_name_resolver::{resolve, smiles_to_iupac, smiles_to_inchi, smiles_to_inchikey};

// IUPAC 系統名 → SMILES
let r = resolve("propan-2-one").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
assert!((r.molecular_weight.unwrap() - 58.08).abs() < 0.01);
assert_eq!(r.confidence, 0.85);              // IUPAC エンジンで解析
assert!(r.inchi.as_deref().unwrap().starts_with("InChI=1S/"));
assert_eq!(r.inchi_key.as_deref().unwrap().len(), 27);

// 慣用名 (辞書)
let r = resolve("acetone").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.confidence, 0.95);              // 正規名辞書エントリ

// カタカナ名
let r = resolve("メタン").unwrap();
assert_eq!(r.smiles, "C");
assert_eq!(r.confidence, 0.90);             // CJK 入力 → 正規化 → 辞書

// SMILES → IUPAC 名
assert_eq!(smiles_to_iupac("CCCCO").unwrap(), "butan-1-ol");
assert_eq!(smiles_to_iupac("CC=CC").unwrap(), "but-2-ene");

// SMILES から InChI / InChIKey
let inchi = smiles_to_inchi("CCO").unwrap();
assert!(inchi.starts_with("InChI=1S/C2H6O/"));
let key = smiles_to_inchikey("CCO").unwrap();
assert_eq!(key.len(), 27);

// JSON 出力
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
| 漢字 → IUPAC/SMILES | 甲烷, 乙醇, 丙酮, 苯, 水, 氯仿, … |

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

倍数接頭辞 `di-`, `tri-`, `tetra-` に対応 (全接尾辞)。

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

### SMILES → IUPAC (`smiles_to_iupac`)

直鎖非環状分子の逆変換。IUPAC パーサーと同じ官能基スコープに対応。分岐鎖および環状・芳香族 SMILES はエラーを返します。

```rust
smiles_to_iupac("CCCCO")  // → "butan-1-ol"
smiles_to_iupac("CC=O")   // → "ethanal"
smiles_to_iupac("CC#CC")  // → "but-2-yne"
smiles_to_iupac("CC(C)CC") // → Err (分岐鎖)
```

### InChI / InChIKey

```rust
smiles_to_inchi("CCO")    // → "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
smiles_to_inchikey("CCO") // → "XXXXXXXXXXXXXX-XXXXXXXXXX-N" (27 文字)
```

> **注意**: 生成される InChI は簡略化された canonical アルゴリズムを使用しており、すべての分子で IUPAC 標準 InChI と一致するとは限りません。

### 出力フィールド

```rust
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,             // Dictionary | Parser
    pub molecular_formula: Option<String>, // Hill 記法 (例: "C2H6O")
    pub molecular_weight: Option<f64>,     // g/mol
    pub confidence: f64,                   // 0.0..=1.0
    pub inchi: Option<String>,             // 標準 InChI 文字列
    pub inchi_key: Option<String>,         // 27 文字の InChIKey
}
```

`DirectSmiles` 経由 (benzene 等、分子グラフを構築しない場合) では、`molecular_formula`、`molecular_weight`、`inchi`、`inchi_key` はすべて `None` になります。

### 信頼度スコア

| シナリオ | スコア |
|----------|--------|
| 辞書完全一致 (DirectSmiles) | `1.00` |
| 辞書完全一致 (CanonicalName) | `0.95` |
| 正規化後に辞書一致 (CJK・ギリシャ文字等) | `0.90` |
| IUPAC 系統名パーサー | `0.85` |

## インストール

```toml
[dependencies]
chem-name-resolver = "0.1"

# JSON 出力が必要な場合
serde_json = "1"
```

## ビルド & テスト

```bash
# テスト 133 件をすべて実行
cargo test

# doctest のみ
cargo test --doc

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
// '{"smiles":"CCO","canonical_name":"ethanol","source":"Dictionary",
//   "molecular_formula":"C2H6O","molecular_weight":46.069,
//   "confidence":0.95,"inchi":"InChI=1S/C2H6O/...","inchi_key":"..."}'
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
#   "molecular_weight": 46.069,
#   "confidence": 0.95,
#   "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
#   "inchi_key": "XXXXXXXXXXXXXX-XXXXXXXXXX-N"
# }

chem resolve --smiles "propan-2-one"
# CC(=O)C
```

## 既知の制限

- 環状・芳香族化合物はパーサー未対応 (辞書引きのみ)
- 立体化学 (R/S, E/Z) は未対応
- `smiles_to_iupac` は直鎖分子のみ対応 (分岐鎖はエラー)
- 生成される InChI/InChIKey は IUPAC 標準値と異なる場合があります

## ロードマップ

- [x] 分岐アルキル置換基 (isopropyl, tert-butyl 等)
- [x] `cyclo-` 接頭辞 (環状化合物)
- [x] CLI バイナリ (`chem resolve "ethanol"`)
- [x] 漢字・中国語化学名辞書
- [x] Canonical SMILES (サブツリー署名 DFS 順序)
- [x] Python バインディング (PyO3 / Maturin)
- [x] `ResolveResult` に信頼度スコアを追加
- [x] SMILES → IUPAC 逆変換 (`smiles_to_iupac`)
- [x] InChI / InChIKey 生成 (`smiles_to_inchi`, `smiles_to_inchikey`)
- [x] 包括的な API ドキュメント + doctest
- [ ] `smiles_to_iupac` の分岐鎖対応
- [ ] 立体化学 (R/S, E/Z)
- [ ] 大規模同義語辞書 (PubChem/ChEBI 経由 `phf_codegen`)

## ライセンス

MIT OR Apache-2.0
