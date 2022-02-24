//! The DNA alphabet including N, consisting of characters A, C, G, N and T.

use crate::interface::alphabet::{Alphabet, AlphabetCharacter};
use std::convert::TryFrom;

/// A character of a DNA alphabet or N: A, C, G, N or T.
pub struct DnaCharacterOrN {
    character: u8,
}

/// The DNA alphabet, consisting of characters A, C, G and T, or N.
pub struct DnaAlphabetOrN;

static DNA_CHARACTER_OR_N_TO_ASCII_TABLE: [u8; 5] = [b'A', b'C', b'G', b'N', b'T'];

impl From<DnaCharacterOrN> for u8 {
    fn from(character: DnaCharacterOrN) -> u8 {
        // Safety: character is private and cannot be constructed out of range.
        unsafe { *DNA_CHARACTER_OR_N_TO_ASCII_TABLE.get_unchecked(character.character as usize) }
    }
}

static ASCII_TO_DNA_CHARACTER_OR_N_TABLE: [u8; 256] = [
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
];

impl TryFrom<u8> for DnaCharacterOrN {
    type Error = ();

    fn try_from(ascii: u8) -> Result<Self, Self::Error> {
        let character = unsafe { *ASCII_TO_DNA_CHARACTER_OR_N_TABLE.get_unchecked(ascii as usize) };
        if character >= 5 {
            Err(())
        } else {
            Ok(Self { character })
        }
    }
}

impl AlphabetCharacter for DnaCharacterOrN {
    fn index(&self) -> usize {
        self.character as usize
    }

    fn from_index(index: usize) -> Option<Self> {
        if index < 5 {
            Some(Self {character: index as u8})
        } else {
            None
        }
    }
}

impl Alphabet for DnaAlphabetOrN {
    type CharacterType = DnaCharacterOrN;

    fn size() -> usize {
        5
    }
}

#[cfg(test)]
mod tests {
    use crate::interface::alphabet::dna_alphabet_or_n::DnaCharacterOrN;
    use std::convert::TryFrom;

    #[test]
    fn test_dna_alphabet_conversion() {
        for ascii in 0u8..=255u8 {
            if ascii == b'A' || ascii == b'C' || ascii == b'G' || ascii == b'N' || ascii == b'T' {
                assert_eq!(
                    u8::from(DnaCharacterOrN::try_from(ascii).unwrap_or_else(|_| panic!(
                        "character {} was expected to be valid, but is not",
                        ascii
                    ))),
                    ascii
                );
            } else {
                assert!(DnaCharacterOrN::try_from(ascii).is_err());
            }
        }
    }
}
