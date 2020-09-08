use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::collections::{LinkedList, VecDeque};
use traitgraph::algo::predefined_graphs::create_binary_tree;
use traitgraph::algo::traversal::{
    BfsQueueStrategy, DfsQueueStrategy, ForwardNeighborStrategy, PreOrderTraversal,
};
use traitgraph::implementation::petgraph_impl;

fn bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_10(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 10).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, BfsQueueStrategy, LinkedList<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_bfs_linkedlist_bintree_10", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_20(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 20).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, BfsQueueStrategy, LinkedList<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_bfs_linkedlist_bintree_20", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_bfs_traversal_vec_deque_bintree_10(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 10).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, BfsQueueStrategy, VecDeque<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_bfs_vecdeque_bintree_10", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_bfs_traversal_vec_deque_bintree_20(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 20).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, BfsQueueStrategy, VecDeque<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_bfs_vecdeque_bintree_20", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_dfs_traversal_linked_list_bintree_10(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 10).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, DfsQueueStrategy, LinkedList<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_dfs_linkedlist_bintree_10", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_dfs_traversal_linked_list_bintree_20(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 20).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, DfsQueueStrategy, LinkedList<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_dfs_linkedlist_bintree_20", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_dfs_traversal_vec_deque_bintree_10(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 10).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, DfsQueueStrategy, VecDeque<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_dfs_vecdeque_bintree_10", |b| {
        b.iter(|| {
            traversal.reset(root);
            for e in &mut traversal {
                black_box(e);
            }
        })
    });
}

fn bench_petgraph_preorder_forward_dfs_traversal_vec_deque_bintree_20(criterion: &mut Criterion) {
    let mut graph = petgraph_impl::new::<(), ()>();
    let root = create_binary_tree(&mut graph, 20).unwrap();
    let mut traversal =
        PreOrderTraversal::<_, ForwardNeighborStrategy, DfsQueueStrategy, VecDeque<_>>::new(
            &graph, root,
        );

    criterion.bench_function("petgraph_dfs_vecdeque_bintree_20", |b| {
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
    bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_10,
    bench_petgraph_preorder_forward_bfs_traversal_linked_list_bintree_20,
    bench_petgraph_preorder_forward_bfs_traversal_vec_deque_bintree_10,
    bench_petgraph_preorder_forward_bfs_traversal_vec_deque_bintree_20,
    bench_petgraph_preorder_forward_dfs_traversal_linked_list_bintree_10,
    bench_petgraph_preorder_forward_dfs_traversal_linked_list_bintree_20,
    bench_petgraph_preorder_forward_dfs_traversal_vec_deque_bintree_10,
    bench_petgraph_preorder_forward_dfs_traversal_vec_deque_bintree_20,
);
criterion_main!(benches);
