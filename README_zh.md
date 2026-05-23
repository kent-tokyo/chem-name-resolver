# chem-name-resolver

纯 Rust 编写的 IUPAC 化学名称 → SMILES 转换库。功能相当于 Java 的 [OPSIN](https://opsin.ch.cam.ac.uk/)，同时支持 Pure Rust + WebAssembly。

## 为何选择 chem-name-resolver？

将 `"2,4-pentanedione"` 这样的 IUPAC 名称转换为 SMILES `"CC(=O)CC(=O)C"` 看似简单，但现有工具各有明显的取舍。

| | [OPSIN](https://opsin.ch.cam.ac.uk/) | [RDKit](https://www.rdkit.org/) | [OpenBabel](https://openbabel.org/) | [CDK](https://cdk.github.io/) | [Indigo](https://lifescience.opensource.epam.com/indigo/) | [PubChem API](https://pubchem.ncbi.nlm.nih.gov/) | [PubChemPy](https://github.com/mcs07/PubChemPy) | [STOUT v2](https://github.com/Kohulan/STOUT) | [ChemCore](https://crates.io/crates/chemcore) | **chem-name-resolver** |
|---|---|---|---|---|---|---|---|---|---|---|
| 语言 | Java | Python/C++ | C++ | Java | C++ | REST | Python | Python/ML | Rust | **Rust** |
| WASM | ✗ | △ | △ | ✗ | ✓ | ✗ | ✗ | ✗ | △ | **✓** |
| 离线 | ✓ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | △ | ✓ | **✓** |
| CJK 名称 | ✗ | ✗ | ✗ | ✗ | ✗ | △ | ✗ | ✗ | ✗ | **✓** |
| IUPAC 解析器 | ✓ (最佳) | ✗ | ✗ | ✗ | ✗ | 词典查询 | ✗ | ✓ (神经网络) | ✗ | **✓** |
| 许可证 | MIT | BSD-3 | GPL-2 | LGPL-2.1 | Apache-2 | 公共领域 | BSD | MIT | MIT | **MIT/Apache-2** |
| 备注 | 需要 JVM | 二进制 ~50 MB；需要 C++ 工具链；rdkit-js WASM 功能受限 | C++ FFI；Copyleft；WASM 实验性 | 需要 JVM；IUPAC 解析委托给 OPSIN | 官方 WASM (npm)；仅结构操作 | 依赖网络；6700 万+ 化合物 | PubChem REST 的薄封装 | 推荐 GPU；非确定性；模型 ~GB | 2020 年起停止维护；SMILES 不完整 | Pure Rust；零原生依赖 |

△ = 部分支持 / 实验性

**本库填补的空缺**：目前尚无同时满足 Pure Rust、WASM 支持、离线运行、CJK 支持的 IUPAC→SMILES 引擎。

**典型使用场景：**

- **浏览器端化学处理** — 以 WASM 模块分发，无需服务器往返即可在客户端解析化学名称
- **Rust 原生工具** — 无需 JVM 或 C++ 构建，集成到 CLI 工具、数据库索引器（如 [Cheminee](https://github.com/rdkit/Cheminee)）或 Axum 服务
- **中文/日文工作流** — 在同一管道中直接规范化汉字、片假名化学名称，无需单独预处理步骤
- **轻量嵌入** — 通过 `opt-level = "s"` + LTO 生成小体积二进制，适合边缘部署

## 特性

- **Pure Rust** — 无 C/C++ 依赖（无需 RDKit/Boost）
- **WASM 支持** — 可编译为 `wasm32-unknown-unknown` 目标
- **CJK 支持** — 支持片假名化学名称（メタン、エタノール 等）的规范化与解析
- **零拷贝规范化** — 输入已规范化时零内存分配
- **JSON 输出** — `ResolveResult` 实现 `serde::Serialize`

## 快速开始

```rust
use chem_name_resolver::resolve;

// IUPAC 系统名
let r = resolve("propan-2-one").unwrap();
assert_eq!(r.smiles, "CC(=O)C");
assert_eq!(r.molecular_formula.as_deref(), Some("C3H6O"));
assert!((r.molecular_weight.unwrap() - 58.08).abs() < 0.01);

// 俗名
let r = resolve("acetone").unwrap();
assert_eq!(r.smiles, "CC(=O)C");

// 片假名名称
let r = resolve("メタン").unwrap();
assert_eq!(r.smiles, "C");

// n- 前缀
let r = resolve("n-butane").unwrap();
assert_eq!(r.smiles, "CCCC");

// JSON 序列化
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

支持倍数前缀 `di-`、`tri-`、`tetra-`。

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

### 输出字段

```rust
pub struct ResolveResult {
    pub smiles: String,
    pub canonical_name: String,
    pub source: ResolveSource,             // Dictionary | Parser
    pub molecular_formula: Option<String>, // Hill 表示法（例："C2H6O"）
    pub molecular_weight: Option<f64>,     // g/mol
}
```

通过 `DirectSmiles` 解析（如 benzene）时，`molecular_formula` / `molecular_weight` 为 `None`。

## 安装

```toml
[dependencies]
chem-name-resolver = "0.1"

# 需要 JSON 输出时
serde_json = "1"
```

## 构建与测试

```bash
# 运行全部测试（75 个）
cargo test

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
// '{"smiles":"CCO","canonical_name":"ethanol","source":"Dictionary","molecular_formula":"C2H6O","molecular_weight":46.069}'
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
#   "molecular_weight": 46.069
# }

chem resolve --smiles "propan-2-one"
# CC(=O)C
```

## 已知限制

- 环状与芳香族化合物不支持解析器（仅词典查询）
- 不支持立体化学（R/S、E/Z）

## 路线图

- [x] 支链烷基取代基（isopropyl、tert-butyl 等）
- [x] `cyclo-` 前缀（环状化合物）
- [x] CLI 二进制（`chem resolve "ethanol"`）
- [x] 汉字/中文化学名称词典
- [x] 规范化 SMILES（子树签名 DFS 排序）
- [x] Python 绑定（PyO3 / Maturin）

## 许可证

MIT OR Apache-2.0
