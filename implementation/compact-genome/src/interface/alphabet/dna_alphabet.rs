//! The DNA alphabet, consisting of characters A, C, G and T.

use crate::interface::alphabet::{Alphabet, AlphabetCharacter};
use std::convert::TryFrom;

/// A character of a DNA alphabet: A, C, G or T.
pub struct DnaCharacter {
    character: u8,
}

/// The DNA alphabet, consisting of characters A, C, G and T.
pub struct DnaAlphabet;

static DNA_CHARACTER_TO_ASCII_TABLE: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl From<DnaCharacter> for u8 {
    fn from(character: DnaCharacter) -> u8 {
        // Safety: character is private and cannot be constructed out of range.
        unsafe { *DNA_CHARACTER_TO_ASCII_TABLE.get_unchecked(character.character as usize) }
    }
}

static ASCII_TO_DNA_CHARACTER_TABLE: [u8; 256] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

impl TryFrom<u8> for DnaCharacter {
    type Error = ();

    fn try_from(ascii: u8) -> Result<Self, Self::Error> {
        let character = unsafe { *ASCII_TO_DNA_CHARACTER_TABLE.get_unchecked(ascii as usize) };
        if character >= 4 {
            Err(())
        } else {
            Ok(Self { character })
        }
    }
}

impl AlphabetCharacter for DnaCharacter {
    fn index(&self) -> usize {
        self.character as usize
    }

    fn from_index(index: usize) -> Option<Self> {
        if index < 4 {
            Some(Self {character: index as u8})
        } else {
            None
        }
    }
}

impl Alphabet for DnaAlphabet {
    type CharacterType = DnaCharacter;

    fn size() -> usize {
        4
    }
}

#[cfg(test)]
mod tests {
    use crate::interface::alphabet::dna_alphabet::DnaCharacter;
    use std::convert::TryFrom;

    #[test]
    fn test_dna_alphabet_conversion() {
        for ascii in 0u8..=255u8 {
            if ascii == b'A' || ascii == b'C' || ascii == b'G' || ascii == b'T' {
                assert_eq!(
                    u8::from(DnaCharacter::try_from(ascii).unwrap_or_else(|_| panic!(
                        "character {} was expected to be valid, but is not",
                        ascii
                    ))),
                    ascii
                );
            } else {
                assert!(DnaCharacter::try_from(ascii).is_err());
            }
        }
    }
}
