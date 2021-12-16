use crate::algo::dijkstra::DijkstraWeight;

macro_rules! impl_dijkstra_weight {
    ($weight_type:ty) => {
        impl DijkstraWeight for $weight_type {
            #[inline]
            fn infinity() -> Self {
                Self::MAX
            }

            #[inline]
            fn zero() -> Self {
                0
            }
        }
    };
}

impl_dijkstra_weight!(usize);
impl_dijkstra_weight!(isize);
impl_dijkstra_weight!(u8);
impl_dijkstra_weight!(i8);
impl_dijkstra_weight!(u16);
impl_dijkstra_weight!(i16);
impl_dijkstra_weight!(u32);
impl_dijkstra_weight!(i32);
impl_dijkstra_weight!(u64);
impl_dijkstra_weight!(i64);
impl_dijkstra_weight!(u128);
impl_dijkstra_weight!(i128);
