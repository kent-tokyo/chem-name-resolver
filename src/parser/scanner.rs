/// Zero-copy &str scanner. Always uses char_indices() — never byte indexing.
pub struct Scanner<'a> {
    input: &'a str,
    pos: usize, // byte position (from char_indices), never set manually
}

impl<'a> Scanner<'a> {
    pub fn new(s: &'a str) -> Self {
        Scanner { input: s, pos: 0 }
    }

    pub fn remaining(&self) -> &'a str {
        &self.input[self.pos..]
    }

    pub fn is_empty(&self) -> bool {
        self.pos >= self.input.len()
    }

    pub fn peek_char(&self) -> Option<char> {
        self.remaining().chars().next()
    }

    pub fn eat_char(&mut self) -> Option<char> {
        let mut iter = self.remaining().char_indices();
        let (_, c) = iter.next()?;
        let next_byte = iter.next().map(|(i, _)| self.pos + i).unwrap_or(self.input.len());
        self.pos = next_byte;
        Some(c)
    }

    /// Eat while predicate is true; returns the consumed slice.
    pub fn eat_while(&mut self, pred: impl Fn(char) -> bool) -> &'a str {
        let start = self.pos;
        for (i, c) in self.remaining().char_indices() {
            if !pred(c) {
                self.pos += i;
                return &self.input[start..self.pos];
            }
        }
        let end = self.input.len();
        self.pos = end;
        &self.input[start..end]
    }

    /// Eat a specific string prefix; returns true if matched.
    pub fn eat_str(&mut self, s: &str) -> bool {
        if self.remaining().starts_with(s) {
            self.pos += s.len();
            true
        } else {
            false
        }
    }

    /// Eat ASCII digits; returns the consumed slice (may be empty).
    pub fn eat_digits(&mut self) -> &'a str {
        self.eat_while(|c| c.is_ascii_digit())
    }

    pub fn byte_pos(&self) -> usize {
        self.pos
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eat_while_ascii() {
        let mut s = Scanner::new("2-chlorobutane");
        assert_eq!(s.eat_digits(), "2");
        assert_eq!(s.eat_while(|c| c == '-'), "-");
        assert_eq!(s.remaining(), "chlorobutane");
    }

    #[test]
    fn eat_char_multibyte() {
        let mut s = Scanner::new("αβ");
        assert_eq!(s.eat_char(), Some('α'));
        assert_eq!(s.eat_char(), Some('β'));
        assert_eq!(s.eat_char(), None);
    }
}
