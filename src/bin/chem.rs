use std::process;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    match args.as_slice() {
        [_, cmd, name] if cmd == "resolve" => {
            run_resolve(name, false);
        }
        [_, cmd, flag, name] if cmd == "resolve" && flag == "--smiles" => {
            run_resolve(name, true);
        }
        [_, cmd, flag, name] if cmd == "resolve" && flag == "-s" => {
            run_resolve(name, true);
        }
        _ => {
            eprintln!("Usage:");
            eprintln!("  chem resolve <name>           Resolve a chemical name (JSON output)");
            eprintln!("  chem resolve --smiles <name>  Resolve and print SMILES only");
            process::exit(1);
        }
    }
}

fn run_resolve(name: &str, smiles_only: bool) {
    match chem_name_resolver::resolve(name) {
        Ok(result) => {
            if smiles_only {
                println!("{}", result.smiles);
            } else {
                match serde_json::to_string_pretty(&result) {
                    Ok(json) => println!("{json}"),
                    Err(e) => {
                        eprintln!("Serialization error: {e}");
                        process::exit(1);
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("Error: {e}");
            process::exit(1);
        }
    }
}
