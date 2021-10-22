use num_traits::{NumCast, PrimInt, ToPrimitive};
use std::hash::Hash;
use std::marker::PhantomData;

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone)]
/// A node index that can be `None`.
/// This is a hack to get a small sized `Option<NodeIndex>` by storing the `None` variant as `IndexType::max_value()`.
/// If Rust ever adds support for integer types with invalid values other than 0, this type becomes obsolete.
pub struct OptionalNodeIndex<IndexType: Sized>(IndexType);
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone)]
/// An edge index that can be `None`.
/// This is a hack to get a small sized `Option<EdgeIndex>` by storing the `None` variant as `IndexType::max_value()`.
/// If Rust ever adds support for integer types with invalid values other than 0, this type becomes obsolete.
pub struct OptionalEdgeIndex<IndexType: Sized>(IndexType);
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone)]
/// A valid node index.
pub struct NodeIndex<IndexType: Sized>(IndexType);
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone)]
/// A valid edge index.
pub struct EdgeIndex<IndexType: Sized>(IndexType);

/// A graph index that can be `None`.
/// This is a hack to get a small sized `Option<GraphIndex>` by storing the `None` variant as `IndexType::max_value()`.
/// If Rust ever adds support for integer types with invalid values other than 0, this trait becomes obsolete.
pub trait OptionalGraphIndex<MirrorGraphIndex: GraphIndex<Self>>:
    Default
    + std::fmt::Debug
    + Eq
    + Ord
    + Hash
    + Copy
    + Sized
    + From<usize>
    + From<Option<usize>>
    + From<MirrorGraphIndex>
    + From<Option<MirrorGraphIndex>>
    + Into<Option<MirrorGraphIndex>>
    + std::ops::Add<usize, Output = Self>
{
    // We don't wanna have OptionalGraphIndex: Into<usize>, to make this type strong, i.e. make it hard to accidentally convert it to a different type.
    /// Get this index as `usize`, but return `None` if this index is marked as invalid.
    fn as_usize(self) -> Option<usize>;

    /// A faster method to get the `usize` value of the index, which does not perform any validity checks.
    fn as_usize_unchecked(self) -> usize;

    /// Returns `true` if the index is valid, i.e. if it is not marked as invalid.
    fn is_valid(self) -> bool {
        self.as_usize().is_some()
    }

    /// Returns `true` if the index is `None`.
    fn is_none(self) -> bool {
        self.as_usize().is_none()
    }

    /// Returns `true` if the index is `Some`.
    fn is_some(self) -> bool {
        self.as_usize().is_some()
    }

    /// Returns a new `OptionalGraphIndex` that is marked as invalid.
    fn new_none() -> Self {
        <Self as From<Option<usize>>>::from(None)
    }

    /// Returns the graph index stored in this optional graph index.
    /// Panics if this optional graph index is `None`.
    fn unwrap(self) -> MirrorGraphIndex {
        self.as_usize().unwrap().into()
    }
}

/// A valid graph index.
pub trait GraphIndex<MirrorOptionalGraphIndex: OptionalGraphIndex<Self>>:
    std::fmt::Debug
    + Eq
    + Ord
    + Hash
    + Copy
    + Sized
    + From<usize>
    + Into<MirrorOptionalGraphIndex>
    + std::ops::Add<usize, Output = Self>
{
    // We don't wanna have GraphIndex: Into<usize>, to make this type strong, i.e. make it hard to accidentally convert it to a different type.
    /// Get this index as `usize`.
    fn as_usize(self) -> usize;
}

macro_rules! impl_graph_index {
    ($GraphIndexType:ident, $OptionalGraphIndexType:ident) => {
        impl<IndexType: PrimInt + Hash> OptionalGraphIndex<$GraphIndexType<IndexType>>
            for $OptionalGraphIndexType<IndexType>
        {
            fn as_usize(self) -> Option<usize> {
                if self.0 != IndexType::max_value() {
                    Some(<usize as NumCast>::from(self.0).unwrap())
                } else {
                    None
                }
            }

            fn as_usize_unchecked(self) -> usize {
                <usize as NumCast>::from(self.0).unwrap()
            }
        }

        impl<IndexType: PrimInt + Hash> GraphIndex<$OptionalGraphIndexType<IndexType>>
            for $GraphIndexType<IndexType>
        {
            fn as_usize(self) -> usize {
                <usize as NumCast>::from(self.0).unwrap()
            }
        }

        impl<IndexType: PrimInt> Default for $OptionalGraphIndexType<IndexType> {
            fn default() -> Self {
                Self(IndexType::max_value())
            }
        }

        impl<IndexType: PrimInt + Hash> std::fmt::Debug for $OptionalGraphIndexType<IndexType> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                if let Some(value) = self.as_usize() {
                    write!(f, "{}", value)
                } else {
                    write!(f, "None")
                }
            }
        }

        impl<IndexType: PrimInt + Hash> std::fmt::Debug for $GraphIndexType<IndexType> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                //debug_assert!(self.0 != IndexType::max_value());
                write!(f, "{}", self.as_usize())
            }
        }

        impl<IndexType: PrimInt> From<Option<usize>> for $OptionalGraphIndexType<IndexType> {
            fn from(source: Option<usize>) -> Self {
                if let Some(source) = source {
                    let source = <IndexType as NumCast>::from(source).unwrap();
                    debug_assert!(source != IndexType::max_value());
                    Self(source)
                } else {
                    Self(IndexType::max_value())
                }
            }
        }

        impl<IndexType: PrimInt> From<usize> for $OptionalGraphIndexType<IndexType> {
            fn from(source: usize) -> Self {
                let source = <IndexType as NumCast>::from(source).unwrap();
                debug_assert!(source != IndexType::max_value());
                Self(source)
            }
        }

        impl<IndexType: PrimInt> From<usize> for $GraphIndexType<IndexType> {
            fn from(source: usize) -> Self {
                let source = <IndexType as NumCast>::from(source).unwrap();
                debug_assert!(source != IndexType::max_value());
                Self(source)
            }
        }

        impl<IndexType: PrimInt + Hash> From<$GraphIndexType<IndexType>>
            for $OptionalGraphIndexType<IndexType>
        {
            fn from(source: $GraphIndexType<IndexType>) -> Self {
                Self::from(source.as_usize())
            }
        }

        impl<IndexType: PrimInt + Hash> From<Option<$GraphIndexType<IndexType>>>
            for $OptionalGraphIndexType<IndexType>
        {
            fn from(source: Option<$GraphIndexType<IndexType>>) -> Self {
                if let Some(source) = source {
                    Self::from(source.as_usize())
                } else {
                    Self::new_none()
                }
            }
        }

        impl<IndexType: PrimInt + Hash> From<$OptionalGraphIndexType<IndexType>>
            for Option<$GraphIndexType<IndexType>>
        {
            fn from(source: $OptionalGraphIndexType<IndexType>) -> Self {
                source.as_usize().map(|source| source.into())
            }
        }

        impl<IndexType: PrimInt + Hash> std::ops::Add<usize>
            for $OptionalGraphIndexType<IndexType>
        {
            type Output = Self;

            fn add(self, rhs: usize) -> Self::Output {
                Self::from(self.as_usize().unwrap() + Self::from(rhs).as_usize().unwrap())
            }
        }

        impl<IndexType: PrimInt + Hash> std::ops::Add<usize> for $GraphIndexType<IndexType> {
            type Output = Self;

            fn add(self, rhs: usize) -> Self::Output {
                Self::from(self.as_usize() + Self::from(rhs).as_usize())
            }
        }
    };
}

impl_graph_index!(NodeIndex, OptionalNodeIndex);
impl_graph_index!(EdgeIndex, OptionalEdgeIndex);

/*impl<IndexType: PrimInt> NodeIndex<IndexType> {
    pub fn invalid() -> Self {
        NodeIndex(IndexType::max_value())
    }

    pub fn is_invalid(&self) -> bool {
        *self == Self::invalid()
    }
}

impl<IndexType: PrimInt> EdgeIndex<IndexType> {
    pub fn invalid() -> Self {
        EdgeIndex(IndexType::max_value())
    }

    pub fn is_invalid(&self) -> bool {
        *self == Self::invalid()
    }
}*/

/*impl<IndexType: PrimInt> From<NodeIndex<IndexType>> for usize {
    fn from(node_index: NodeIndex<IndexType>) -> Self {
        NumCast::from(node_index.0).unwrap()
    }
}

impl<IndexType: PrimInt> From<EdgeIndex<IndexType>> for usize {
    fn from(node_index: EdgeIndex<IndexType>) -> Self {
        NumCast::from(node_index.0).unwrap()
    }
}*/

impl<T, IndexType: PrimInt + Hash> std::ops::Index<NodeIndex<IndexType>> for Vec<T> {
    type Output = T;

    fn index(&self, index: NodeIndex<IndexType>) -> &Self::Output {
        &self[index.as_usize()]
    }
}

impl<T, IndexType: PrimInt + Hash> std::ops::Index<EdgeIndex<IndexType>> for Vec<T> {
    type Output = T;

    fn index(&self, index: EdgeIndex<IndexType>) -> &Self::Output {
        &self[index.as_usize()]
    }
}

impl<T, IndexType: PrimInt + Hash> std::ops::IndexMut<NodeIndex<IndexType>> for Vec<T> {
    fn index_mut(&mut self, index: NodeIndex<IndexType>) -> &mut Self::Output {
        &mut self[index.as_usize()]
    }
}

impl<T, IndexType: PrimInt + Hash> std::ops::IndexMut<EdgeIndex<IndexType>> for Vec<T> {
    fn index_mut(&mut self, index: EdgeIndex<IndexType>) -> &mut Self::Output {
        &mut self[index.as_usize()]
    }
}

// Unstable as of now (slice_index_methods)
/*impl<IndexType: PrimInt> SliceIndex<[NodeIndex<IndexType>]> for NodeIndex<IndexType> {
    type Output = NodeIndex<IndexType>;

    fn get(self, slice: &[NodeIndex<IndexType>]) -> Option<&Self::Output> {
        slice.get(self.0)
    }

    fn get_mut(self, slice: &mut [NodeIndex<IndexType>]) -> Option<&mut Self::Output> {
        slice.get_mut(self.0)
    }

    unsafe fn get_unchecked(self, slice: &[NodeIndex<IndexType>]) -> &Self::Output {
        slice.get_unchecked(self.0)
    }

    unsafe fn get_unchecked_mut(self, slice: &mut [NodeIndex<IndexType>]) -> &mut Self::Output {
        slice.get_unchecked_mut(self.0)
    }

    fn index(self, slice: &[NodeIndex<IndexType>]) -> &Self::Output {
        slice.index(self.0)
    }

    fn index_mut(self, slice: &mut [NodeIndex<IndexType>]) -> &mut Self::Output {
        slice.index_mut(self.0)
    }
}*/

/*impl<IndexType: ToPrimitive> ToPrimitive for NodeIndex<IndexType> {
    fn to_i64(&self) -> Option<i64> {
        self.0.to_i64()
    }

    fn to_u64(&self) -> Option<u64> {
        self.0.to_u64()
    }
}

impl<IndexType: ToPrimitive> ToPrimitive for EdgeIndex<IndexType> {
    fn to_i64(&self) -> Option<i64> {
        self.0.to_i64()
    }

    fn to_u64(&self) -> Option<u64> {
        self.0.to_u64()
    }
}*/

/*impl<IndexType: PrimInt> From<IndexType> for NodeIndex<IndexType> {
    fn from(index: IndexType) -> Self {
        Self(index)
    }
}

impl<IndexType: PrimInt> From<IndexType> for EdgeIndex<IndexType> {
    fn from(index: IndexType) -> Self {
        Self(index)
    }
}*/

/*impl<IndexType: PrimInt> From<NodeIndex<IndexType>> for IndexType {
    fn from(index: NodeIndex<IndexType>) -> Self {
        index.0
    }
}

impl<IndexType: PrimInt> From<EdgeIndex<IndexType>> for IndexType {
    fn from(index: EdgeIndex<IndexType>) -> Self {
        index.0
    }
}*/

/// An iterator over a consecutive sequence of graph indices.
pub struct GraphIndices<IndexType, OptionalIndexType> {
    start: IndexType,
    end: IndexType,
    optional_index_type: PhantomData<OptionalIndexType>,
}

impl<
        RawType: ToPrimitive,
        OptionalIndexType: OptionalGraphIndex<IndexType>,
        IndexType: GraphIndex<OptionalIndexType>,
    > From<(RawType, RawType)> for GraphIndices<IndexType, OptionalIndexType>
{
    fn from(raw: (RawType, RawType)) -> Self {
        Self {
            start: IndexType::from(raw.0.to_usize().unwrap()),
            end: IndexType::from(raw.1.to_usize().unwrap()),
            optional_index_type: Default::default(),
        }
    }
}

impl<
        OptionalIndexType: OptionalGraphIndex<IndexType>,
        IndexType: GraphIndex<OptionalIndexType>,
    > Iterator for GraphIndices<IndexType, OptionalIndexType>
{
    type Item = IndexType;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let result = Some(self.start);
            self.start = self.start + 1;
            result
        } else {
            None
        }
    }
}
