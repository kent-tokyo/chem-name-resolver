use thiserror::Error;

/// Errors returned by the resolver and parser.
#[derive(Debug, Error)]
pub enum ResolveError {
    /// The name was not found in any dictionary and could not be parsed as a
    /// systematic IUPAC name.
    #[error("name not found in dictionary or parseable by IUPAC rules: {0:?}")]
    NotFound(String),

    /// A structural or syntactic parse error occurred.
    ///
    /// `pos` is the approximate byte offset in the input where the error was
    /// detected (0 if unavailable). `msg` describes the problem.
    #[error("IUPAC parse error at position {pos}: {msg}")]
    ParseError { pos: usize, msg: String },

    /// A bond-count exceeds the element's maximum valence.
    #[error("valence violation: atom {atom_idx} ({element}) has {bonds} bonds, max {max}")]
    ValenceError {
        atom_idx: usize,
        element: String,
        bonds: u8,
        max: u8,
    },
}
