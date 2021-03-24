use crate::omnitigs::{NodeCentricUnivocalExtensionAlgorithm, UnivocalExtensionAlgorithm};
use std::borrow::Borrow;
use traitgraph::interface::StaticGraph;
use traitgraph::walks::{EdgeWalk, NodeWalk};

/// Computes the univocal extension of a walk assuming the graph is strongly connected.
/// Might enter an infinite loop if the graph is not strongly connected.
pub struct SCCUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::EdgeIndex>>>
    UnivocalExtensionAlgorithm<Graph, ResultWalk> for SCCUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::EdgeIndex]) -> ResultWalk {
        EdgeWalk::compute_univocal_extension(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is not strongly connected.
pub struct NonSCCUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::EdgeIndex>>>
    UnivocalExtensionAlgorithm<Graph, ResultWalk> for NonSCCUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::EdgeIndex]) -> ResultWalk {
        EdgeWalk::compute_univocal_extension_non_scc(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is strongly connected.
/// Might enter an infinite loop if the graph is not strongly connected.
pub struct SCCNodeCentricUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::NodeIndex>>>
    NodeCentricUnivocalExtensionAlgorithm<Graph, ResultWalk>
    for SCCNodeCentricUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::NodeIndex]) -> ResultWalk {
        NodeWalk::compute_univocal_extension(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is not strongly connected.
pub struct NonSCCNodeCentricUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::NodeIndex>>>
    NodeCentricUnivocalExtensionAlgorithm<Graph, ResultWalk>
    for NonSCCNodeCentricUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::NodeIndex]) -> ResultWalk {
        NodeWalk::compute_univocal_extension_non_scc(walk.borrow(), graph)
    }
}
