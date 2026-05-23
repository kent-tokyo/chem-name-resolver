use phf::phf_map;

/// Maps trivial/common names to canonical IUPAC names (all keys lowercase, normalized).
pub static SYNONYM_TO_IUPAC: phf::Map<&'static str, &'static str> = phf_map! {
    "ethanol"             => "ethanol",
    "alcohol"             => "ethanol",
    "acetone"             => "propan-2-one",
    "acetic acid"         => "ethanoic acid",
    "formic acid"         => "methanoic acid",
    "glycerol"            => "propane-1,2,3-triol",
    "glycerin"            => "propane-1,2,3-triol",
    "acetaldehyde"        => "ethanal",
    "formaldehyde"        => "methanal",
    // iso/sec/tert aliases that the parser can handle via a canonical IUPAC name
    "isopropanol"         => "propan-2-ol",
    "isopropyl alcohol"   => "propan-2-ol",
    "sec-butanol"         => "butan-2-ol",
    "2-butanol"           => "butan-2-ol",
    // Lab abbreviations → canonical IUPAC (parser resolves these)
    "meoh"                => "methanol",
    "etoh"                => "ethanol",
    // Old "-ic acid" trivial names
    "propionic acid"      => "propanoic acid",
    "butyric acid"        => "butanoic acid",
    "valeric acid"        => "pentanoic acid",
    "caproic acid"        => "hexanoic acid",
    // Common ketone aliases
    "methyl ethyl ketone" => "butan-2-one",
    "mek"                 => "butan-2-one",
    "butanone"            => "butan-2-one",
    // Amines: simple ones resolved via canonical IUPAC
    "ethylamine"          => "ethanamine",
    "propylamine"         => "propan-1-amine",
};

/// Maps trivial names directly to SMILES (for names too complex for the Phase 3 parser).
pub static SYNONYM_TO_SMILES: phf::Map<&'static str, &'static str> = phf_map! {
    "water"               => "O",
    "oxidane"             => "O",
    "benzene"             => "c1ccccc1",
    "toluene"             => "Cc1ccccc1",
    "aspirin"             => "CC(=O)Oc1ccccc1C(=O)O",
    "glucose"             => "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "caffeine"            => "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    // Ether-type: parser cannot handle alkoxy prefixes yet
    "ether"               => "CCOCC",
    "ethoxyethane"        => "CCOCC",
    "diethyl ether"       => "CCOCC",
    // Polyhalomethanes: parser cannot handle tri/tetra substituents without locants
    "chloroform"          => "C(Cl)(Cl)Cl",
    "trichloromethane"    => "C(Cl)(Cl)Cl",
    // Branched alkanes/alcohols: parser cannot handle branched chains yet
    "isobutane"           => "CC(C)C",
    "isobutanol"          => "CC(C)CO",
    "isobutyl alcohol"    => "CC(C)CO",
    "tert-butanol"        => "CC(C)(C)O",
    "t-butanol"           => "CC(C)(C)O",
    "neopentane"          => "CC(C)(C)C",
    // Lab solvents
    "dcm"                 => "ClCCl",
    "methylene chloride"  => "ClCCl",
    "dichloromethane"     => "ClCCl",
    "dmso"                => "CS(=O)C",
    "dimethyl sulfoxide"  => "CS(=O)C",
    "dmf"                 => "CN(C)C=O",
    "dimethylformamide"   => "CN(C)C=O",
    "thf"                 => "C1CCOC1",
    "tetrahydrofuran"     => "C1CCOC1",
    "mecn"                => "CC#N",
    "acetonitrile"        => "CC#N",
    // Halomethanes: parser requires a locant for substituents
    "chloromethane"       => "CCl",
    "bromomethane"        => "CBr",
    "iodomethane"         => "CI",
    "fluoromethane"       => "CF",
    "dibromomethane"      => "BrCBr",
    // Nitro compounds: bracket atoms needed for [N+](=O)[O-]; parser cannot handle yet
    "nitromethane"        => "[N+](=O)[O-]C",
    "nitroethane"         => "CC[N+](=O)[O-]",
    "nitrobenzene"        => "c1ccc(cc1)[N+](=O)[O-]",
    // Amines: too complex for current parser (multiple N bonds or no locant)
    "methylamine"         => "CN",
    "dimethylamine"       => "CNC",
    "trimethylamine"      => "CN(C)C",
    "diethylamine"        => "CCNCC",
    "triethylamine"       => "CCN(CC)CC",
    "tea"                 => "CCN(CC)CC",
    "aniline"             => "Nc1ccccc1",
    // Phenols and aromatic compounds
    "phenol"              => "Oc1ccccc1",
    "anisole"             => "COc1ccccc1",
    "methoxybenzene"      => "COc1ccccc1",
    "styrene"             => "C=Cc1ccccc1",
    "vinylbenzene"        => "C=Cc1ccccc1",
    "o-xylene"            => "Cc1ccccc1C",
    "1,2-dimethylbenzene" => "Cc1ccccc1C",
    "m-xylene"            => "Cc1cccc(C)c1",
    "1,3-dimethylbenzene" => "Cc1cccc(C)c1",
    "p-xylene"            => "Cc1ccc(C)cc1",
    "1,4-dimethylbenzene" => "Cc1ccc(C)cc1",
    "mesitylene"          => "Cc1cc(C)cc(C)c1",
    "1,3,5-trimethylbenzene" => "Cc1cc(C)cc(C)c1",
    // Cyclic compounds: parser cannot handle ring closure yet
    "cyclohexane"         => "C1CCCCC1",
    "cyclohexanol"        => "OC1CCCCC1",
    "cyclohexanone"       => "O=C1CCCCC1",
    "cyclopentane"        => "C1CCCC1",
    "cyclopentanol"       => "OC1CCCC1",
    "cyclopentanone"      => "O=C1CCCC1",
    "cyclobutane"         => "C1CCC1",
    "cyclopropane"        => "C1CC1",
    // Branched alkanes not handled by parser
    "isopentane"          => "CCC(C)C",
    "2-methylbutane"      => "CCC(C)C",
    "isohexane"           => "CCCC(C)C",
    "2-methylpentane"     => "CCCC(C)C",
    // Common esters (parser cannot handle two-chain names yet)
    "ethyl acetate"       => "CCOC(=O)C",
    "etoac"               => "CCOC(=O)C",
    "methyl acetate"      => "COC(=O)C",
    "methyl formate"      => "COC=O",
    "ethyl formate"       => "CCOC=O",
};
