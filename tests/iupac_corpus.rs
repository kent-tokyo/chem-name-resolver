use chem_name_resolver::{resolve, resolve_batch};
use serde_json;

const CORPUS: &[(&str, &str)] = &[
    // Alkanes C1-C10
    ("methane",          "C"),
    ("ethane",           "CC"),
    ("propane",          "CCC"),
    ("butane",           "CCCC"),
    ("pentane",          "CCCCC"),
    ("hexane",           "CCCCCC"),
    // Alcohols
    ("ethanol",          "CCO"),
    ("propan-1-ol",      "CCCO"),
    ("propan-2-ol",      "CC(C)O"),
    // Ketones
    ("propan-2-one",     "CC(=O)C"),
    ("2,4-pentanedione", "CC(=O)CC(=O)C"),
    // Alkenes / alkynes
    ("hex-1-ene",        "C=CCCCC"),
    ("hex-2-ene",        "CC=CCCC"),
    ("but-2-yne",        "CC#CC"),
    // Halides / substituents
    ("2-chlorobutane",   "CC(CC)Cl"),
    ("3-methylpentane",  "CCC(C)CC"),
    // Aldehydes
    ("pentanal",         "CCCCC=O"),
    // Acids
    ("ethanoic acid",    "CC(=O)O"),
    ("butanedioic acid", "C(=O)(CCC(=O)O)O"),
    // Complex compound suffix
    ("hept-3-yn-1-ol",   "CCCC#CCCO"),
];

#[test]
fn iupac_corpus() {
    for (name, expected) in CORPUS {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn extended_chains_c11_c20() {
    let cases = [
        ("undecane",     11usize, "C11H24"),
        ("dodecane",     12,      "C12H26"),
        ("tridecane",    13,      "C13H28"),
        ("tetradecane",  14,      "C14H30"),
        ("pentadecane",  15,      "C15H32"),
        ("hexadecane",   16,      "C16H34"),
        ("heptadecane",  17,      "C17H36"),
        ("octadecane",   18,      "C18H38"),
        ("nonadecane",   19,      "C19H40"),
        ("icosane",      20,      "C20H42"),
        ("eicosane",     20,      "C20H42"),
    ];
    for (name, chain_len, expected_formula) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        // SMILES: straight chain of chain_len carbons
        let expected_smiles = "C".repeat(chain_len);
        assert_eq!(result.smiles, expected_smiles, "SMILES mismatch for: {name}");
        assert_eq!(
            result.molecular_formula.as_deref(),
            Some(expected_formula),
            "formula mismatch for: {name}"
        );
    }
}

#[test]
fn molecular_formula_and_weight() {
    let cases: &[(&str, &str, f64)] = &[
        ("methane",       "CH4",    16.043),
        ("ethane",        "C2H6",   30.070),
        ("ethanol",       "C2H6O",  46.069),
        ("propan-2-one",  "C3H6O",  58.080),
        ("ethanoic acid", "C2H4O2", 60.052),
    ];
    for (name, expected_formula, expected_mw) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(
            result.molecular_formula.as_deref(),
            Some(*expected_formula),
            "formula mismatch for: {name}"
        );
        let mw = result.molecular_weight.unwrap();
        assert!(
            (mw - expected_mw).abs() < 0.01,
            "{name}: Mw = {mw:.3}, expected {expected_mw:.3}"
        );
    }
}

#[test]
fn dictionary_names() {
    let cases = [
        ("water",    "O"),
        ("benzene",  "c1ccccc1"),
        ("toluene",  "Cc1ccccc1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn katakana_names() {
    let cases = [("メタン", "C"), ("エタノール", "CCO")];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn heteroatom_suffixes() {
    let cases = [
        // amines
        ("ethanamine",        "CCN"),
        ("propan-1-amine",    "CCCN"),
        // thiols
        ("ethanethiol",       "CCS"),
        ("propane-1-thiol",   "CCCS"),
        // nitriles (C-1 triple-bonded to N)
        ("ethanenitrile",     "CC#N"),
        ("propanenitrile",    "CCC#N"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn heteroatom_substituents() {
    let cases = [
        ("2-aminobutane",     "CC(CC)N"),
        ("3-mercaptopentane", "CCC(CC)S"),
        ("2-cyanopentane",    "CC(C#N)CCC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn substituent_multiplier_prefix() {
    let cases = [
        ("2,3-dichlorobutane",  "CC(C(C)Cl)Cl"),
        ("2,2-dimethylpropane", "CC(C)(C)C"),
        ("2,4-dimethylpentane", "CC(C)CC(C)C"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn n_prefix_strip() {
    let cases = [
        ("n-butane",  "CCCC"),
        ("n-pentane", "CCCCC"),
        ("n-hexane",  "CCCCCC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn dictionary_bug_fixes() {
    let cases = [
        ("ether",          "CCOCC"),
        ("chloroform",     "C(Cl)(Cl)Cl"),
        ("diethyl ether",  "CCOCC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn longer_alkyl_substituents() {
    let cases = [
        ("3-propylhexane",  "CCC(CCC)CCC"),
        ("3-butylheptane",  "CCC(CCCC)CCCC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn resolve_result_is_serializable() {
    let result = resolve("ethanol").unwrap();
    let json = serde_json::to_string(&result).expect("ResolveResult must be serializable");
    assert!(json.contains("CCO"));
    assert!(json.contains("C2H6O"));
}

#[test]
fn additional_dict_entries() {
    let cases = [
        // Branched alkanes
        ("isopentane",          "CCC(C)C"),
        ("2-methylbutane",      "CCC(C)C"),
        ("isohexane",           "CCCC(C)C"),
        ("2-methylpentane",     "CCCC(C)C"),
        // Common ketones
        ("methyl ethyl ketone", "CC(=O)CC"),
        ("mek",                 "CC(=O)CC"),
        ("butanone",            "CC(=O)CC"),
        // Common esters
        ("ethyl acetate",       "CCOC(=O)C"),
        ("etoac",               "CCOC(=O)C"),
        ("methyl acetate",      "COC(=O)C"),
        ("methyl formate",      "COC=O"),
        ("ethyl formate",       "CCOC=O"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn amine_dict_entries() {
    let cases = [
        ("methylamine",   "CN"),
        ("dimethylamine", "CNC"),
        ("trimethylamine","CN(C)C"),
        ("ethylamine",    "CCN"),
        ("diethylamine",  "CCNCC"),
        ("triethylamine", "CCN(CC)CC"),
        ("TEA",           "CCN(CC)CC"),
        ("aniline",       "Nc1ccccc1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn aromatic_dict_entries() {
    let cases = [
        ("phenol",              "Oc1ccccc1"),
        ("anisole",             "COc1ccccc1"),
        ("methoxybenzene",      "COc1ccccc1"),
        ("styrene",             "C=Cc1ccccc1"),
        ("o-xylene",            "Cc1ccccc1C"),
        ("1,2-dimethylbenzene", "Cc1ccccc1C"),
        ("m-xylene",            "Cc1cccc(C)c1"),
        ("p-xylene",            "Cc1ccc(C)cc1"),
        ("mesitylene",          "Cc1cc(C)cc(C)c1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn cyclic_compound_dict_entries() {
    let cases = [
        ("cyclohexane",   "C1CCCCC1"),
        ("cyclohexanol",  "OC1CCCCC1"),
        ("cyclohexanone", "O=C1CCCCC1"),
        ("cyclopentane",  "C1CCCC1"),
        ("cyclopentanol", "OC1CCCC1"),
        ("cyclobutane",   "C1CCC1"),
        ("cyclopropane",  "C1CC1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn resolve_batch_api() {
    let names = ["methane", "ethanol", "acetone"];
    let results = resolve_batch(&names);
    assert_eq!(results.len(), 3);
    assert_eq!(results[0].as_ref().unwrap().smiles, "C");
    assert_eq!(results[1].as_ref().unwrap().smiles, "CCO");
    assert_eq!(results[2].as_ref().unwrap().smiles, "CC(=O)C");
}

#[test]
fn iso_sec_tert_dict_entries() {
    let cases = [
        ("isopropanol",       "CC(C)O"),
        ("isopropyl alcohol", "CC(C)O"),
        ("sec-butanol",       "CC(CC)O"),
        ("2-butanol",         "CC(CC)O"),
        ("isobutane",         "CC(C)C"),
        ("isobutanol",        "CC(C)CO"),
        ("isobutyl alcohol",  "CC(C)CO"),
        ("tert-butanol",      "CC(C)(C)O"),
        ("t-butanol",         "CC(C)(C)O"),
        ("neopentane",        "CC(C)(C)C"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn lab_abbreviations() {
    let cases = [
        ("MeOH",              "CO"),
        ("EtOH",              "CCO"),
        ("DCM",               "ClCCl"),
        ("methylene chloride","ClCCl"),
        ("dichloromethane",   "ClCCl"),
        ("DMSO",              "CS(=O)C"),
        ("dimethyl sulfoxide","CS(=O)C"),
        ("DMF",               "CN(C)C=O"),
        ("dimethylformamide", "CN(C)C=O"),
        ("THF",               "C1CCOC1"),
        ("tetrahydrofuran",   "C1CCOC1"),
        ("MeCN",              "CC#N"),
        ("acetonitrile",      "CC#N"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn halomethane_dict_entries() {
    let cases = [
        ("chloromethane",  "CCl"),
        ("bromomethane",   "CBr"),
        ("iodomethane",    "CI"),
        ("fluoromethane",  "CF"),
        ("dibromomethane", "BrCBr"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn old_ic_acid_names() {
    let cases = [
        ("propionic acid", "CCC(=O)O"),
        ("butyric acid",   "CCCC(=O)O"),
        ("valeric acid",   "CCCCC(=O)O"),
        ("caproic acid",   "CCCCCC(=O)O"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn nitro_dict_entries() {
    let cases = [
        ("nitromethane", "[N+](=O)[O-]C"),
        ("nitroethane",  "CC[N+](=O)[O-]"),
        ("nitrobenzene", "c1ccc(cc1)[N+](=O)[O-]"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn amide_suffix() {
    let cases = [
        // -amide: C(=O)N at C-1; DFS from far terminal gives C(=O)N suffix
        ("ethanamide",     "CC(=O)N"),
        ("propanamide",    "CCC(=O)N"),
        ("butanamide",     "CCCC(=O)N"),
        ("pentanamide",    "CCCCC(=O)N"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn acyl_substituents() {
    let cases = [
        // acetyl- = -C(=O)CH3; DFS writes branch as "(C(=O)C)"
        ("3-acetylheptane",  "CCC(C(=O)C)CCCC"),
        // formyl- = -CHO; DFS writes branch as "(C=O)"
        ("3-formylpentane",  "CCC(C=O)CC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn cyclo_prefix_parser() {
    let cases = [
        // These ring sizes are not in the dictionary → go through the parser
        ("cycloheptane",  "C1CCCCCC1"),
        ("cyclooctane",   "C1CCCCCCC1"),
        ("cyclononane",   "C1CCCCCCCC1"),
        ("cyclodecane",   "C1CCCCCCCCC1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn branched_alkyl_substituents() {
    let cases = [
        ("2-isopropylpentane",  "CC(C(C)C)CCC"),
        ("2-tert-butylpentane", "CC(C(C)(C)C)CCC"),
        ("3-sec-butylheptane",  "CCC(C(C)CC)CCCC"),
        ("3-isobutylheptane",   "CCC(CC(C)C)CCCC"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}

#[test]
fn chinese_names() {
    let cases = [
        ("甲烷",   "C"),
        ("乙烷",   "CC"),
        ("乙醇",   "CCO"),
        ("丙酮",   "CC(=O)C"),
        ("苯",     "c1ccccc1"),
        ("水",     "O"),
        ("氯仿",   "C(Cl)(Cl)Cl"),
        ("乙酸",   "CC(=O)O"),
        ("环己烷", "C1CCCCC1"),
    ];
    for (name, expected) in cases {
        let result = resolve(name).unwrap_or_else(|e| panic!("resolve({name:?}) failed: {e}"));
        assert_eq!(&result.smiles, expected, "SMILES mismatch for: {name}");
    }
}
