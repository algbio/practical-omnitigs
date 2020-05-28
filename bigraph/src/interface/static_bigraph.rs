use crate::{ImmutableGraphContainer, NavigableGraph, NodeBigraph};
use num_traits::{PrimInt};

pub trait StaticBigraph<NodeData, EdgeData, IndexType: PrimInt>:
ImmutableGraphContainer<NodeData, EdgeData, IndexType>
+ for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>
+ NodeBigraph<NodeData, EdgeData, IndexType>
{
    /**
     * Returns true if the [mirror property] of edges is fulfilled.
     * Assumes that the node pairing is correct (See [verify_node_pairing()](NodeBigraphWrapper::verify_node_pairing))
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn verify_mirror_property(&self) -> bool {
        for from_node in self.node_indices() {
            for to_node in self.out_neighbors(from_node).unwrap() {
                let from_node_partner = self.partner_node(from_node).unwrap();
                let to_node_partner = self.partner_node(to_node.node_id).unwrap();
                if !self.contains_edge(to_node_partner, from_node_partner) {
                    return false
                }
            }
        }

        true
    }

    /**
     * Returns true if all unitigs in the graph have length of at most one node.
     */
    fn verify_unitig_length_is_zero(&self) -> bool {
        unimplemented!()
    }
}
impl<
    NodeData,
    EdgeData,
    IndexType: PrimInt,
    T: ImmutableGraphContainer<NodeData, EdgeData, IndexType>
    + for<'a> NavigableGraph<'a, NodeData, EdgeData, IndexType>
    + NodeBigraph<NodeData, EdgeData, IndexType>,
> StaticBigraph<NodeData, EdgeData, IndexType> for T
{
}
