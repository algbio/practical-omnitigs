use crate::omnitigs::{NodeCentricUnivocalExtensionAlgorithm, UnivocalExtensionAlgorithm};
use crate::walks::{EdgeOmnitigLikeExt, NodeOmnitigLikeExt};
use std::borrow::Borrow;
use traitgraph::interface::StaticGraph;

/// Computes the univocal extension of a walk assuming the graph is strongly connected.
/// Might enter an infinite loop if the graph is not strongly connected.
pub struct SccUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::EdgeIndex>>>
    UnivocalExtensionAlgorithm<Graph, ResultWalk> for SccUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::EdgeIndex]) -> ResultWalk {
        EdgeOmnitigLikeExt::compute_univocal_extension(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is not strongly connected.
pub struct NonSccUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::EdgeIndex>>>
    UnivocalExtensionAlgorithm<Graph, ResultWalk> for NonSccUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::EdgeIndex]) -> ResultWalk {
        EdgeOmnitigLikeExt::compute_univocal_extension_non_scc(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is strongly connected.
/// Might enter an infinite loop if the graph is not strongly connected.
pub struct SccNodeCentricUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::NodeIndex>>>
    NodeCentricUnivocalExtensionAlgorithm<Graph, ResultWalk>
    for SccNodeCentricUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::NodeIndex]) -> ResultWalk {
        NodeOmnitigLikeExt::compute_univocal_extension(walk.borrow(), graph)
    }
}

/// Computes the univocal extension of a walk assuming the graph is not strongly connected.
pub struct NonSccNodeCentricUnivocalExtensionStrategy;

impl<Graph: StaticGraph, ResultWalk: From<Vec<Graph::NodeIndex>>>
    NodeCentricUnivocalExtensionAlgorithm<Graph, ResultWalk>
    for NonSccNodeCentricUnivocalExtensionStrategy
{
    fn compute_univocal_extension(graph: &Graph, walk: &[Graph::NodeIndex]) -> ResultWalk {
        NodeOmnitigLikeExt::compute_univocal_extension_non_scc(walk.borrow(), graph)
    }
}
