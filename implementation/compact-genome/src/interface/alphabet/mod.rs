//! Alphabets for genome sequences.

use std::convert::{TryFrom, TryInto};

pub mod dna_alphabet;
pub mod dna_alphabet_or_n;

/// A character in an alphabet.
pub trait AlphabetCharacter: Into<u8> + TryFrom<u8> {
    /// The index of this character in the alphabet.
    fn index(&self) -> usize;

    /// Constructs the character from the given index, returning `None` if it is invalid.
    fn from_index(index: usize) -> Option<Self>;
}

/// An alphabet as a subset of the ASCII alphabet.
pub trait Alphabet: Sized {
    /// The internal character type used by the alphabet.
    type CharacterType: AlphabetCharacter;

    /// The amount of characters in the alphabet.
    fn size() -> usize;

    /// Converts the given ASCII character into an alphabet character.
    /// If the ASCII character is not mapped to an alphabet character, then `None` is returned.
    fn ascii_to_character(ascii: u8) -> Option<Self::CharacterType> {
        ascii.try_into().ok()
    }

    /// Converts this alphabet character into an ASCII character.
    fn character_to_ascii(character: Self::CharacterType) -> u8 {
        character.into()
    }
}
