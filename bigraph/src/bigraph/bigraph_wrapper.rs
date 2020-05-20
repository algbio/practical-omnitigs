use crate::bigraph::StaticGraph;

pub struct BigraphWrapper<NodeData, EdgeData, IndexType, T: StaticGraph<NodeData, EdgeData, IndexType>> {
    topology: T,

}