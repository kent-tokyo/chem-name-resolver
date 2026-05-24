# chem-name-resolver

纯 Rust 编写的 IUPAC 化学名称解析库，可将化学名称解析为 SMILES 字符串、IUPAC 名称、InChI 标识符和分子属性。功能相当于 Java 的 [OPSIN](https://opsin.ch.cam.ac.uk/)，同时支持 Pure Rust + WebAssembly。

## 为何选择 chem-name-resolver？

将 `"2,4-pentanedione"` 这样的 IUPAC 名称转换为 SMILES `"CC(=O)CC(=O)C"` 看似简单，但现有工具各有明显的取舍。

| | [OPSIN](https://opsin.ch.cam.ac.uk/) | [RDKit](https://www.rdkit.org/) | [OpenBabel](https://openbabel.org/) | [CDK](https://cdk.github.io/) | [Indigo](https://lifescience.opensource.epam.com/indigo/) | [PubChem API](https://pubchem.ncbi.nlm.nih.gov/) | [PubChemPy](https://github.com/mcs07/PubChemPy) | [STOUT v2](https://github.com/Kohulan/STOUT) | [ChemCore](https://crates.io/crates/chemcore) | **chem-name-resolver** |
|---|---|---|---|---|---|---|---|---|---|---|
| 语言 | Java | Python/C++ | C++ | Java | C++ | REST | Python | Python/ML | Rust | **Rust** |
| WASM | ✗ | △ | △ | ✗ | ✓ | ✗ | ✗ | ✗ | △ | **✓** |
| 离线 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | △ | ✓ | **✓** |
| CJK 名称 | ✗ | ✗ | ✗ | ✗ | ✗ | △ | ✗ | ✗ | ✗ | **✓** |
| IUPAC 解析器 | ✓ (最佳) | ✗ | ✗ | ✗ | ✗ | 词典查询 | ✗ | ✓ (神经网络) | ✗ | **✓** |
| 反向转换 (SMILES→IUPAC) | ✗ | △ | △ | ✗ | ✗ | △ | ✗ | ✓ (神经网络) | ✗ | **✓** |
| 许可证 | MIT | BSD-3 | GPL-2 | LGPL-2.1 | Apache-2 | 公共领域 | BSD | MIT | MIT | **MIT/Apache-2** |
| 备注 | 需要 JVM | 二进制 ~50 MB；需要 C++ 工具链；rdkit-js WASM 功能受限 | C++ FFI；Copyleft；WASM 实验性 | 需要 JVM；IUPAC 解析委托给 OPSIN | 官方 WASM (npm)；仅结构操作 | 依赖网络；6700 万+ 化合物 | PubChem REST 的薄封装 | 推荐 GPU；非确定性；模型 ~GB | 2020 年起停止维护；SMILES 不完整 | Pure Rust；零原生依赖 |

△ = 部分支持 / 实验性

**本库填补的空缺**：目前尚无同时满足 Pure Rust、WASM 支持、离线运行、CJK 支持的 IUPAC↔SMILES 引擎。

**典型使用场景：**

- **浏览器端化学处理** — 以 WASM 模块分发，无需服务器往返即可在客户端解析化学名称
- **Rust 原生工具** — 无需 JVM 或 C++ 构建，集成到 CLI 工具、数据库索引器（如 [Cheminee](https://github.com/rdkit/Cheminee)）或 Axum 服务
- **中文/日文工作流** — 在同一管道中直接规范化汉字、片假名化学名称，无需单独预处理步骤
- **轻量嵌入** — 通过 `opt-level = "s"` + LTO 生成小体积二进制，适合边缘部署

## 特性

- **Pure Rust** — 无 C/C++ 依赖（无需 RDKit/Boost）
- **WASM 支持** — 可编译为 `wasm32-unknown-unknown` 目标
- **CJK 支持** — 支持片假名化学名称（メタン、エタノール 等）及汉字（甲烷、乙醇 等）的规范化与解析
- **双向转换**：IUPAC→SMILES (`resolve`) **及** SMILES→IUPAC (`smiles_to_iupac`)
- **InChI/InChIKey** 生成（从 SMILES）
- **置信度评分** — 每个结果附带 `0.0～1.0` 的质量信号
- **零拷贝规范化** — 输入已规范化时返回 `Cow::Borrowed`，零内存分配
- **JSON 输出** — `ResolveResult` 实现 `serde::Serialize`

## 快速开始

```rust
use chem_name_resolver::{resolve, smiles_to_iupac, smiles_to_inchi, smiles_to_inchikey};

// IUPAC 系统名 → SMILES
let r = resolve("propan-2-one").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
assert!((r.molecular_weight.unwrap() - 58.08).abs() < 0.01);
assert_eq!(r.confidence, 0.85);              // 由 IUPAC 引擎解析
assert!(r.inchi.as_deref().unwrap().starts_with("InChI=1S/"));
assert_eq!(r.inchi_key.as_deref().unwrap().len(), 27);

// 俗名（词典）
let r = resolve("acetone").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.confidence, 0.95);              // 规范名词典条目

// 片假名名称
let r = resolve("メタン").unwrap();
assert_eq!(r.smiles, "C");
assert_eq!(r.confidence, 0.90);             // CJK 输入 → 规范化 → 词典

// SMILES → IUPAC 名称
assert_eq!(smiles_to_iupac("CCCCO").unwrap(), "butan-1-ol");
assert_eq!(smiles_to_iupac("CC=CC").unwrap(), "but-2-ene");

// 从 SMILES 生成 InChI / InChIKey
let inchi = smiles_to_inchi("CCO").unwrap();
assert!(inchi.starts_with("InChI=1S/C2H6O/"));
let key = smiles_to_inchikey("CCO").unwrap();
assert_eq!(key.len(), 27);

// JSON 输出
let json = serde_json::to_string(&r).unwrap();
```

## 支持的转换

### 规范化器

| 输入 | 输出 |
|------|------|
| 全角数字与符号（`２－`） | 半角（`2-`） |
| 长音符（`ー`） | 连字符（`-`） |
| 希腊字母（`α`、`β`、`γ`） | ASCII（`alpha`、`beta`、`gamma`） |
| 连续空格 | 单个空格 |
| `n-` 前缀 | 去除（`n-butane` → `butane`） |

### 词典（同义词库）

| 类型 | 示例 |
|------|------|
| 俗名 → IUPAC 系统名 | acetone、acetic acid、glycerol、formaldehyde、propionic/butyric/valeric acid、… |
| 俗名 → 直接 SMILES | water、benzene、toluene、ether、chloroform、aspirin、glucose、caffeine |
| iso/sec/tert 别名 | isopropanol、isobutane、tert-butanol、neopentane、sec-butanol、… |
| 支链烷烃 | isopentane、isohexane（+ IUPAC 系统名别名） |
| 实验室缩写 | MeOH、EtOH、DCM、DMSO、DMF、THF、MeCN（含全称） |
| 卤代甲烷 | chloromethane、bromomethane、iodomethane、dibromomethane、… |
| 常用试剂 | ethyl acetate、methyl acetate、MEK（+ 全称） |
| 胺类 | methylamine、dimethylamine、trimethylamine、aniline、triethylamine、… |
| 苯酚/芳香族 | phenol、anisole、styrene、o/m/p-xylene、mesitylene、… |
| 环状化合物 | cyclohexane、cyclohexanol、cyclohexanone、cyclopentane、cyclopropane、… |
| 硝基化合物 | nitromethane、nitroethane、nitrobenzene |
| 片假名 → IUPAC 系统名 | メタン～デカン、エタノール、アセトン、ベンゼン、… |
| 汉字 → IUPAC/SMILES | 甲烷、乙醇、丙酮、苯、水、氯仿、… |

### IUPAC 解析器

**链骨干：** methane～decane（C1–C10）、undecane～icosane/eicosane（C11–C20）

**后缀：**

| 后缀 | 官能团 | 示例 |
|------|--------|------|
| `-ane` | 烷烃 | ethane → `CC` |
| `-ene` | 烯烃 | hex-1-ene → `C=CCCCC` |
| `-yne` | 炔烃 | but-2-yne → `CC#CC` |
| `-ol` / `-diol` | 醇 | propan-2-ol → `CC(C)O` |
| `-one` / `-dione` | 酮 | propan-2-one → `CC(=O)C` |
| `-al` | 醛 | pentanal → `CCCCC=O` |
| `-oic acid` / `-dioic acid` | 羧酸 | ethanoic acid → `CC(=O)O` |
| `-amine` | 胺 | ethanamine → `CCN` |
| `-amide` | 酰胺 | ethanamide → `CC(=O)N` |
| `-thiol` | 硫醇 | ethanethiol → `CCS` |
| `-nitrile` | 腈 | propanenitrile → `CCC#N` |

所有后缀均支持倍数前缀 `di-`、`tri-`、`tetra-`。

**取代基：**

| 取代基 | 原子/基团 | 示例 |
|--------|----------|------|
| `chloro-`、`bromo-`、`fluoro-`、`iodo-` | 卤素 | 2-chlorobutane → `CC(CC)Cl` |
| `methyl-`、`ethyl-`、`propyl-`、`butyl-`、`pentyl-`、`hexyl-` | 直链烷基 | 3-methylpentane → `CCC(C)CC` |
| `hydroxy-` | 羟基 | — |
| `oxo-` | 酮基 | — |
| `amino-` | 氨基 | 2-aminobutane → `CC(CC)N` |
| `mercapto-` | 巯基 | 3-mercaptopentane → `CCC(CC)S` |
| `cyano-` | 氰基 | 2-cyanopentane → `CC(C#N)CCC` |
| `acetyl-` | 乙酰基 | 3-acetylheptane → `CCC(C(=O)C)CCCC` |
| `formyl-` | 甲酰基 | 3-formylpentane → `CCC(C=O)CC` |

支持倍数前缀 `di-`、`tri-`、`tetra-`（例：`2,3-dichlorobutane` → `CC(C(C)Cl)Cl`）。

### SMILES → IUPAC (`smiles_to_iupac`)

直链非环状分子的反向转换，支持与 IUPAC 解析器相同的官能团范围。含支链及环状/芳香族 SMILES 将返回错误。

```rust
smiles_to_iupac("CCCCO")  // → "butan-1-ol"
smiles_to_iupac("CC=O")   // → "ethanal"
smiles_to_iupac("CC#CC")  // → "but-2-yne"
smiles_to_iupac("CC(C)CC") // → Err（支链）
```

### InChI / InChIKey

```rust
smiles_to_inchi("CCO")    // → "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
smiles_to_inchikey("CCO") // → "XXXXXXXXXXXXXX-XXXXXXXXXX-N"（27 个字符）
```

> **注意**：生成的 InChI 使用简化的规范化算法，对于所有分子可能与 IUPAC 标准 InChI 不完全一致。

### 输出字段

```rust
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,             // Dictionary | Parser
    pub molecular_formula: Option<String>, // Hill 表示法（例："C2H6O"）
    pub molecular_weight: Option<f64>,     // g/mol
    pub confidence: f64,                   // 0.0..=1.0
    pub inchi: Option<String>,             // 标准 InChI 字符串
    pub inchi_key: Option<String>,         // 27 字符 InChIKey
}
```

通过 `DirectSmiles` 解析（如 benzene，不构建分子图）时，`molecular_formula`、`molecular_weight`、`inchi`、`inchi_key` 均为 `None`。

### 置信度评分

| 场景 | 评分 |
|------|------|
| 词典精确匹配（DirectSmiles） | `1.00` |
| 词典精确匹配（CanonicalName） | `0.95` |
| 规范化后词典匹配（CJK、希腊字母等） | `0.90` |
| IUPAC 系统名解析器 | `0.85` |

## 安装

```toml
[dependencies]
chem-name-resolver = "0.1"

# 需要 JSON 输出时
serde_json = "1"
```

## 构建与测试

```bash
# 运行全部 133 个测试
cargo test

# 仅运行 doctest
cargo test --doc

# 验证 WASM 构建
rustup target add wasm32-unknown-unknown
cargo build --features wasm --target wasm32-unknown-unknown

# 基准测试
cargo bench
```

## WASM 使用示例

```javascript
import init, { resolve_to_smiles, resolve_full, normalize_name } from './chem_name_resolver.js';

await init();
console.log(resolve_to_smiles("propan-2-one")); // "CC(=O)C"
console.log(normalize_name("α-D-glucose"));     // "alpha-d-glucose"

// 以 JSON 字符串获取完整结果
const json = resolve_full("ethanol");
// '{"smiles":"CCO","canonical_name":"ethanol","source":"Dictionary",
//   "molecular_formula":"C2H6O","molecular_weight":46.069,
//   "confidence":0.95,"inchi":"InChI=1S/C2H6O/...","inchi_key":"..."}'
```

## CLI 使用示例

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

## 已知限制

- 环状与芳香族化合物不支持解析器（仅词典查询）
- 不支持立体化学（R/S、E/Z）
- `smiles_to_iupac` 仅支持直链分子（支链 → 错误）
- 生成的 InChI/InChIKey 可能与 IUPAC 标准值不同

## 路线图

- [x] 支链烷基取代基（isopropyl、tert-butyl 等）
- [x] `cyclo-` 前缀（环状化合物）
- [x] CLI 二进制（`chem resolve "ethanol"`）
- [x] 汉字/中文化学名称词典
- [x] 规范化 SMILES（子树签名 DFS 排序）
- [x] Python 绑定（PyO3 / Maturin）
- [x] `ResolveResult` 置信度评分
- [x] SMILES → IUPAC 反向转换（`smiles_to_iupac`）
- [x] InChI / InChIKey 生成（`smiles_to_inchi`、`smiles_to_inchikey`）
- [x] 完整 API 文档 + doctest
- [ ] `smiles_to_iupac` 支链支持
- [ ] 立体化学（R/S、E/Z）
- [ ] 大规模同义词词典（通过 PubChem/ChEBI `phf_codegen`）

## 许可证

MIT OR Apache-2.0
