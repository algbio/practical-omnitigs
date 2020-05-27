use num_traits::{NumCast, PrimInt, ToPrimitive};
use std::ops::{Index, IndexMut};

#[derive(Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy, OpaqueTypedef)]
#[opaque_typedef(derive(Display, FromInner))]
pub struct NodeIndex<IndexType: Sized>(IndexType);
#[derive(Default, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Clone, Copy, OpaqueTypedef)]
#[opaque_typedef(derive(Display, FromInner))]
pub struct EdgeIndex<IndexType: Sized>(IndexType);

impl<IndexType: PrimInt> NodeIndex<IndexType> {
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
}

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

impl<T, IndexType: PrimInt> Index<NodeIndex<IndexType>> for Vec<T> {
    type Output = T;

    fn index(&self, index: NodeIndex<IndexType>) -> &Self::Output {
        &self[<usize as NumCast>::from(index.0).unwrap()]
    }
}

impl<T, IndexType: PrimInt> Index<EdgeIndex<IndexType>> for Vec<T> {
    type Output = T;

    fn index(&self, index: EdgeIndex<IndexType>) -> &Self::Output {
        &self[<usize as NumCast>::from(index.0).unwrap()]
    }
}

impl<T, IndexType: PrimInt> IndexMut<NodeIndex<IndexType>> for Vec<T> {
    fn index_mut(&mut self, index: NodeIndex<IndexType>) -> &mut Self::Output {
        &mut self[<usize as NumCast>::from(index.0).unwrap()]
    }
}

impl<T, IndexType: PrimInt> IndexMut<EdgeIndex<IndexType>> for Vec<T> {
    fn index_mut(&mut self, index: EdgeIndex<IndexType>) -> &mut Self::Output {
        &mut self[<usize as NumCast>::from(index.0).unwrap()]
    }
}

impl<IndexType: ToPrimitive> ToPrimitive for NodeIndex<IndexType> {
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
}

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

pub struct NodeIndices<IndexType: PrimInt> {
    start: IndexType,
    end: IndexType,
}

pub struct EdgeIndices<IndexType: PrimInt> {
    start: IndexType,
    end: IndexType,
}

impl<RawType: ToPrimitive, IndexType: PrimInt> From<(RawType, RawType)> for NodeIndices<IndexType> {
    fn from(raw: (RawType, RawType)) -> Self {
        Self {
            start: <IndexType as NumCast>::from(raw.0).unwrap(),
            end: <IndexType as NumCast>::from(raw.1).unwrap(),
        }
    }
}

impl<RawType: PrimInt, IndexType: PrimInt> From<(RawType, RawType)> for EdgeIndices<IndexType> {
    fn from(raw: (RawType, RawType)) -> Self {
        Self {
            start: <IndexType as NumCast>::from(raw.0).unwrap(),
            end: <IndexType as NumCast>::from(raw.1).unwrap(),
        }
    }
}

impl<IndexType: PrimInt> Iterator for NodeIndices<IndexType> {
    type Item = NodeIndex<IndexType>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let result = Some(NodeIndex(self.start));
            self.start = self.start + NumCast::from(1).unwrap();
            result
        } else {
            None
        }
    }
}

impl<IndexType: PrimInt> Iterator for EdgeIndices<IndexType> {
    type Item = EdgeIndex<IndexType>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start < self.end {
            let result = Some(EdgeIndex(self.start));
            self.start = self.start + NumCast::from(1).unwrap();
            result
        } else {
            None
        }
    }
}
