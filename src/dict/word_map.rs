use phf::phf_map;

/// Maps katakana chemical words to their IUPAC English equivalents.
/// Keys are normalized katakana (ー already replaced with - by normalizer).
pub static KATAKANA_TO_IUPAC: phf::Map<&'static str, &'static str> = phf_map! {
    "メタン"        => "methane",
    "エタン"        => "ethane",
    "プロパン"      => "propane",
    "ブタン"        => "butane",
    "ペンタン"      => "pentane",
    "ヘキサン"      => "hexane",
    "ヘプタン"      => "heptane",
    "オクタン"      => "octane",
    "ノナン"        => "nonane",
    "デカン"        => "decane",
    "エタノ-ル"    => "ethanol",
    "メタノ-ル"    => "methanol",
    "プロパノ-ル"  => "propanol",
    "ブタノ-ル"    => "butanol",
    "アセトン"      => "propan-2-one",
    "ベンゼン"      => "benzene",
    "トルエン"      => "toluene",
    "クロロホルム"  => "trichloromethane",
    "ジエチルエ-テル" => "ethoxyethane",
};

/// Maps Chinese (Simplified) chemical names to their IUPAC English equivalents.
/// Covers common alkanes, alcohols, and laboratory solvents.
pub static HANZI_TO_IUPAC: phf::Map<&'static str, &'static str> = phf_map! {
    // Alkanes C1–C10 (系统命名)
    "甲烷"   => "methane",
    "乙烷"   => "ethane",
    "丙烷"   => "propane",
    "丁烷"   => "butane",
    "戊烷"   => "pentane",
    "己烷"   => "hexane",
    "庚烷"   => "heptane",
    "辛烷"   => "octane",
    "壬烷"   => "nonane",
    "癸烷"   => "decane",
    // Alcohols
    "甲醇"   => "methanol",
    "乙醇"   => "ethanol",
    "丙醇"   => "propanol",
    "丁醇"   => "butanol",
    "异丙醇" => "propan-2-ol",
    // Aldehydes
    "甲醛"   => "methanal",
    "乙醛"   => "ethanal",
    // Ketones
    "丙酮"   => "propan-2-one",
    // Carboxylic acids
    "甲酸"   => "methanoic acid",
    "乙酸"   => "ethanoic acid",
    "醋酸"   => "ethanoic acid",
    "丙酸"   => "propanoic acid",
    "丁酸"   => "butanoic acid",
};

/// Maps Chinese names that resolve directly to SMILES (too complex for the parser).
pub static HANZI_TO_SMILES: phf::Map<&'static str, &'static str> = phf_map! {
    "苯"     => "c1ccccc1",
    "甲苯"   => "Cc1ccccc1",
    "水"     => "O",
    "氯仿"   => "C(Cl)(Cl)Cl",
    "苯酚"   => "Oc1ccccc1",
    "葡萄糖" => "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "咖啡因" => "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "乙醚"   => "CCOCC",
    "二氯甲烷" => "ClCCl",
    "二甲亚砜" => "CS(=O)C",
    "四氢呋喃" => "C1CCOC1",
    "乙腈"   => "CC#N",
    "硝基苯" => "c1ccc(cc1)[N+](=O)[O-]",
    "苯胺"   => "Nc1ccccc1",
    "环己烷" => "C1CCCCC1",
    "环戊烷" => "C1CCCC1",
    "环丙烷" => "C1CC1",
    "异丁烷" => "CC(C)C",
    "新戊烷" => "CC(C)(C)C",
};
