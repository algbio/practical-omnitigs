//! Alphabets for genome sequences.

use std::convert::{TryFrom, TryInto};

pub mod dna_alphabet;
pub mod dna_alphabet_or_n;

/// A character in an alphabet.
pub trait AlphabetCharacter: Into<u8> + TryFrom<u8> {
    /// The amount of characters in the alphabet.
    const ALPHABET_SIZE: usize;

    /// The index of this character in the alphabet.
    fn index(&self) -> usize;

    /// Constructs the character from the given index, returning `None` if it is invalid.
    fn from_index(index: usize) -> Result<Self, AlphabetError>;

    /// Constructs the character from the given index, returning `None` if it is invalid.
    /// This method returns a static reference to the character type, so it can only be implemented via lookup in a static table.
    /// It is required to create an implementation of [std::ops::Index] for genome sequence types that do not store the characters in plain format.
    fn from_index_ref(index: usize) -> Result<&'static Self, AlphabetError>;

    /// Constructs the complement of this character.
    fn complement(&self) -> Self;
}

/// An alphabet as a subset of the ASCII alphabet.
pub trait Alphabet: Sized {
    /// The amount of characters in the alphabet.
    const SIZE: usize = Self::CharacterType::ALPHABET_SIZE;

    /// The internal character type used by the alphabet.
    type CharacterType: AlphabetCharacter + Eq + Ord + Clone + 'static;

    /// Converts the given ASCII character into an alphabet character.
    /// If the ASCII character is not mapped to an alphabet character, then `None` is returned.
    fn ascii_to_character(ascii: u8) -> Result<Self::CharacterType, AlphabetError> {
        ascii
            .try_into()
            .map_err(|_| AlphabetError::AsciiNotPartOfAlphabet { ascii })
    }

    /// Converts this alphabet character into an ASCII character.
    fn character_to_ascii(character: Self::CharacterType) -> u8 {
        character.into()
    }
}

/// An error when dealing with alphabets.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum AlphabetError {
    /// An ascii character was attempted to convert to an alphabet character, but it is not part of the alphabet.
    AsciiNotPartOfAlphabet {
        /// The offending ascii character.
        ascii: u8,
    },

    /// An index was attempted to convert to an alphabet character, but it is not part of the alphabet.
    IndexNotPartOfAlphabet {
        /// The offending index.
        index: usize,
    },
}
