"""
We define here the 4 primal operations used in the Blossom algorithm: Grow, Augment, Expand and Expand.
"""


function primal_grow(
    edge::Edge,
    forest::Forest,
    matching::Matching,
    g::Graph,
    w::Dict{Edge,U},
) where {U<:Real}
    # Called if the tight edge (u, v) has labels positive and free. 
    # The tree from the positive node is aquires the free node and its matched node.

    # Get the nodes from the edge
    u = forest.pseudonodes[edge.src]
    v = forest.pseudonodes[edge.dst]
    mate = matching.mate[edge.dst]
    tree = forest.trees[u.tree_id]

    # Add the 2 nodes v, mate and 2 edges (u, v) and (v, mate) to the tree
    backbone = tree.backbone
    psnode_to_backbone = tree.psnode_to_backbone
    backbone_to_psnode = tree.backbone_to_psnode
    add_vertex!(backbone)
    psnode_to_backbone[v.node_id] = nv(backbone)
    backbone_to_psnode[nv(backbone)] = v.node_id

    add_vertex!(backbone)
    psnode_to_backbone[mate] = nv(backbone)
    backbone_to_psnode[nv(backbone)] = mate

    # Add the 2 edges
    add_edge!(backbone, psnode_to_backbone[u.node_id], psnode_to_backbone[v.node_id]) # (u, v)
    add_edge!(backbone, psnode_to_backbone[v.node_id], psnode_to_backbone[mate]) # (v, mate)

    # Update the labels and tree_id of the nodes added
    v.label = negative
    forest.pseudonodes[mate].label = positive
    v.tree_id = u.tree_id
    forest.pseudonodes[mate].tree_id = u.tree_id

end

function primal_augment(
    edge::Edge,
    forest::Forest,
    matching::Matching,
    g::Graph,
    w::Dict{Edge,U},
) where {U<:Real}
    # Called if the tight edge (u, v) has labels positive and positive, and are in different trees.
    # The path root1 -> u -> v -> root2 is augmented, all nodes in the trees become free nodes. 

    # Get the 2 nodes and their trees
    u = forest.pseudonodes[edge.src]
    v = forest.pseudonodes[edge.dst]
    tree1 = forest.trees[u.tree_id]
    tree2 = forest.trees[v.tree_id]

    # Perform changes to the matching : inverse the matching of the path u -> root1 and v -> root2
    augment_path(forest, matching, u, tree1, w)
    augment_path(forest, matching, v, tree2, w)

    # Add (u, v) to the matching
    add_edge_matching!(matching, edge, get_weight(w, edge))

    # Set all nodes in the trees to free
    free_tree!(tree1, forest)
    free_tree!(tree2, forest)


end

function augment_path(
    forest::Forest,
    matching::Matching,
    node::PseudoNode,
    tree::AlternatingTree,
    w::Dict{Edge,U},
) where {U<:Real}
    # Augment the path from the node to the root of the tree
    # The path is defined by the in-neighbors of the node in the backbone of the tree
    if node.node_id == tree.root
        return
    end
    # Get the backbone of the tree
    backbone = tree.backbone
    psnode_to_backbone = tree.psnode_to_backbone
    backbone_to_psnode = tree.backbone_to_psnode

    # Get the node id in the backbone
    bb_node_id = psnode_to_backbone[node.node_id]
    bb_parent = first(inneighbors(backbone, bb_node_id))
    parent_label = negative # The parent of the starting node is always negative
    parent = forest.pseudonodes[backbone_to_psnode[bb_parent]]

    while true # Loop until we reach the root
        if parent_label == negative
            # The parent is negative, remove the edge from the matching
            edge = get_edge(parent.node_id, node.node_id)
            remove_edge_matching!(matching, edge, get_weight(w, edge))
        elseif parent_label == positive
            # The parent is positive, add the edge to the matching
            edge = get_edge(parent.node_id, node.node_id)
            add_edge_matching!(matching, edge, get_weight(w, edge))
        else
            # Should not happen
            error("Encountered free node in an AlternatingTree")
        end
        # Move one step up the tree
        bb_node_id = bb_parent
        node = parent
        bb_parents = inneighbors(backbone, bb_node_id)
        if length(bb_parents) == 0
            break
        end
        bb_parent = first(bb_parents)
        parent = forest.pseudonodes[backbone_to_psnode[bb_parent]]
        parent_label = parent.label
    end

end

function free_tree!(tree::AlternatingTree, forest::Forest)
    # Set all nodes in the tree to free
    psnodes_ids = keys(tree.psnode_to_backbone)
    for psnode_id in psnodes_ids
        node = forest.pseudonodes[psnode_id]
        node.label = free
        node.tree_id = nothing
    end
    # Delete the tree from the forest
    delete!(forest.trees, tree.root)
end
