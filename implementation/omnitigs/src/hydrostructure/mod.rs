/// An incremental implementation of the hydrostructure.
/// Computes and stores the hydrostructure in linear time for a macrotig and all of its subwalks of length at least two.
pub mod incremental_hydrostructure;
/// A static implementation of the hydrostructure for a walk.
pub mod static_hydrostructure;

/// The hydrostructure of a walk `W` as defined in the hydrostructure paper.
pub trait Hydrostructure<NodeIndex: Copy, EdgeIndex: Copy> {
    /// Returns true if the given node is in `R⁺(W)`.
    fn is_node_r_plus(&self, node: NodeIndex) -> bool;
    /// Returns true if the given node is in `R⁻(W)`.
    fn is_node_r_minus(&self, node: NodeIndex) -> bool;
    /// Returns true if the given edge is in `R⁺(W)`.
    fn is_edge_r_plus(&self, edge: EdgeIndex) -> bool;
    /// Returns true if the given edge is in `R⁻(W)`.
    fn is_edge_r_minus(&self, edge: EdgeIndex) -> bool;

    /// Returns true if `W` is bridge-like.
    fn is_bridge_like(&self) -> bool;
    /// Returns true if `W` is avertible.
    fn is_avertible(&self) -> bool {
        !self.is_bridge_like()
    }

    /// Returns true if the given node is in the river.
    fn is_node_river(&self, node: NodeIndex) -> bool {
        !self.is_node_r_plus(node) && !self.is_node_r_minus(node)
    }

    /// Returns true if the given node is in the cloud.
    fn is_node_cloud(&self, node: NodeIndex) -> bool {
        !self.is_node_r_plus(node) && self.is_node_r_minus(node)
    }

    /// Returns true if the given node is in the sea.
    fn is_node_sea(&self, node: NodeIndex) -> bool {
        self.is_node_r_plus(node) && !self.is_node_r_minus(node)
    }

    /// Returns true if the given node is in the vapor.
    fn is_node_vapor(&self, node: NodeIndex) -> bool {
        self.is_node_r_plus(node) && self.is_node_r_minus(node)
    }

    /// Returns true if the given edge is in the river.
    fn is_edge_river(&self, edge: EdgeIndex) -> bool {
        !self.is_edge_r_plus(edge) && !self.is_edge_r_minus(edge)
    }

    /// Returns true if the given edge is in the cloud.
    fn is_edge_cloud(&self, edge: EdgeIndex) -> bool {
        !self.is_edge_r_plus(edge) && self.is_edge_r_minus(edge)
    }

    /// Returns true if the given edge is in the sea.
    fn is_edge_sea(&self, edge: EdgeIndex) -> bool {
        self.is_edge_r_plus(edge) && !self.is_edge_r_minus(edge)
    }

    /// Returns true if the given edge is in the vapor.
    fn is_edge_vapor(&self, edge: EdgeIndex) -> bool {
        self.is_edge_r_plus(edge) && self.is_edge_r_minus(edge)
    }
}
