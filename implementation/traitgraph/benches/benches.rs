use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::collections::LinkedList;
use traitgraph::algo::predefined_graphs::create_binary_tree;
use traitgraph::algo::traversal::{BfsQueueStrategy, ForwardNeighborStrategy, PreOrderTraversal};
use traitgraph::implementation::petgraph_impl;
use traitgraph::interface::MutableGraphContainer;

fn bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_10(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new();
    let root = create_binary_tree(&mut graph, 10).unwrap();
    let n = graph.add_node(());
    graph.add_edge(n, n, ());
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, BfsQueueStrategy, LinkedList<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph linkedlist bintree 10", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

criterion_group!(
    benches,
    bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_10
);
criterion_main!(benches);
