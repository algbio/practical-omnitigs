use crate::interface::static_bigraph::StaticBigraph;
use crate::{DynamicGraph, StaticBigraphFromDigraph};
use num_traits::PrimInt;

pub trait DynamicBigraph<NodeData, EdgeData, IndexType: PrimInt>:
    DynamicGraph<NodeData, EdgeData, IndexType> + StaticBigraph<NodeData, EdgeData, IndexType>
{
    /**
     * Adds edges such that the graph fulfils the [mirror property].
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn add_mirror_edges(&mut self) {
        unimplemented!()
    }
    /**
     * Adds nodes such that the graph becomes a valid bigraph.
     * The indices of existing nodes are not altered.
     */
    fn add_partner_nodes(&mut self);
}

pub trait DynamicBigraphFromDigraph<NodeData, EdgeData, IndexType: PrimInt>:
    DynamicBigraph<NodeData, EdgeData, IndexType>
    + StaticBigraphFromDigraph<NodeData, EdgeData, IndexType>
where
    Self::Topology: DynamicGraph<NodeData, EdgeData, IndexType>,
{
    /**
     * Converts the given topology into a bigraph.
     * Missing partners are added.
     */
    fn new_with_completed_nodes(
        _topology: Self::Topology,
        _binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self {
        unimplemented!()
    }

    /**
     * Converts the given topology into a bigraph that fulfils the [mirror property].
     * Missing partners and mirror edges are added.
     *
     * [mirror property]: https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md
     */
    fn new_with_completed_nodes_and_mirrored_edges(
        topology: Self::Topology,
        binode_mapping_function: fn(&NodeData) -> NodeData,
    ) -> Self {
        let mut bigraph = Self::new_with_completed_nodes(topology, binode_mapping_function);
        bigraph.add_mirror_edges();
        bigraph
    }
}
