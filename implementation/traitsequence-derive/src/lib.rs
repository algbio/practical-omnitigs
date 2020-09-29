#![warn(missing_docs)]
//! This crate offers derive macros for easier implementation of sequence types.

#[macro_use]
extern crate syn;

use proc_macro::{TokenStream, TokenTree};
use syn::{DeriveInput, Attribute};
use proc_macro_roids::DeriveInputStructExt;
use proc_macro_roids::FieldExt;
use syn::export::quote::__private::ext::RepToTokensExt;
use syn::parse::Parse;

/// Derive the sequence trait and further related traits for types wrapping a sequence.
///
/// This macro can only derive from a struct and expects the struct to have a field marked with `#[inner_sequence]` that implements `Sequence`.
#[proc_macro_derive(Sequence, attributes(inner_sequence))]
pub fn derive_sequence(item: TokenStream) -> TokenStream {
    let ast = parse_macro_input!(item as DeriveInput);
    let mut relevant_fields = ast.fields().iter().filter(|field| {
        field.attrs
            .iter()
            .map(Attribute::parse_meta)
            .filter_map(Result::ok)
            .filter(|meta| meta.path() == &parse_quote!(inner_sequence)).next().is_some()
    });
    let sequence_field = relevant_fields.next().expect("MISSING a field annotated with #[inner_sequence]. Exactly one such field must be present.");
    assert!(relevant_fields.next().is_none(), "Found MORE THAN ONE field annotated with #[inner_sequence]. Exactly one such field must be present.");

    println!("{:#?}", sequence_field);

    todo!()
}