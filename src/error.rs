use thiserror::Error;

#[derive(Debug, Error)]
pub enum ResolveError {
    #[error("name not found in dictionary or parseable by IUPAC rules: {0:?}")]
    NotFound(String),

    #[error("IUPAC parse error at position {pos}: {msg}")]
    ParseError { pos: usize, msg: String },

    #[error("valence violation: atom {atom_idx} ({element}) has {bonds} bonds, max {max}")]
    ValenceError {
        atom_idx: usize,
        element: String,
        bonds: u8,
        max: u8,
    },
}
