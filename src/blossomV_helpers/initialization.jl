
function empty_matching(n::T, u::U) where {T<:Integer,U<:Real}
    return Matching(Set{Edge}(), Set{T}(), zeros(T, n), u, zero(T))
end

"""
    get_unit_tree(n::T) where {T<:Integer}

    Create a tree with a single node 'n'.
 """
function get_unit_tree(node::PseudoNode{T,U}) where {T<:Integer,U<:Real}

    backbone = SimpleDiGraph(1)
    psnodes_to_backbone = Dict{T,T}()
    backbone_to_psnode = Dict{T,T}()
    psnodes_to_backbone[node.node_id] = 1
    backbone_to_psnode[1] = node.node_id


    t = AlternatingTree(
        backbone,
        node.node_id,
        psnodes_to_backbone,
        backbone_to_psnode,
        0.0,
        nothing,
        PriorityQueue{PseudoNode,Int}(),
        PriorityQueue{PseudoNode,Int}(),
        PriorityQueue{PseudoNode,Int}(),
    )
    return t
end

"""
Initializes the blossoms_dict pseudonodes_dict, pseudoedges_dict and trees_dict
"""
function make_forest(g::Graph{T}, w::Dict{E,U}) where {U<:Real,E<:Edge,T<:Integer}
    # Create a dict of blossom, where pseudonode_id -> Set(nodes) of pseudonodes in the bossom, used for expanding
    blossoms_dict = Dict{T,Set{T}}()
    # TODO : Make this better using a tree structure

    # Create a dict of trees, where tree_id -> AlternatingTree. Empty.
    trees_dict = Dict{T,AlternatingTree{T,U}}()

    # Create a dict of pseudonodes, where pseudonode_id -> PseudoNode. Initialize with all nodes in the graph
    pseudonodes_dict = Dict{T,PseudoNode{T,U}}()
    for v = 1:nv(g)
        pseudonodes_dict[v] = BaseNode(positive, v, Inf, nothing, nothing)
    end

    # Create a dict of EdgeData, where (src, dst) -> EdgeData. Initialize with all edges in the graph
    edgedata_dict = Dict{Tuple{T,T},EdgeData{T,U}}()
    for e in edges(g)
        _src = min(src(e), dst(e))
        _dst = max(src(e), dst(e))
        edgedata_dict[(_src, _dst)] = EdgeData(_src, _dst, zero(U))
    end

    return Forest(pseudonodes_dict, edgedata_dict, trees_dict, blossoms_dict)

end


"""
    From the graph g and weigh map w, create a pseudonode coresponding to the original nodes in g. 

    To initialize y :
    First, for each node v we go through the incident edges, find the smallest weight and set y_v to the half of this weight. 
    (This guarantees that all edges now have non-negative slacks.) 
    We then go over nodes again, and for each node v greedily increase y_v and choose a matching for v, if possible.

"""
function greedy_initialization(g::Graph{T}, w::Dict{E,U}) where {U<:Real,E<:Edge,T<:Integer}
    # Create empty matching
    match = empty_matching(ne(g), zero(U))
    # Create empty forest
    forest = make_forest(g, w)
    pseudonodes_dict = forest.pseudonodes

    # First pass : set y values conservatively
    for v = 1:nv(g)
        if degree(g, v) == 0
            continue
        end
        temp_y = mapreduce(dest -> get_weight(w, v, dest), min, neighbors(g, v))
        pseudonodes_dict[v].ybar = temp_y / 2
    end

    """ Not sure how this is suposed to work
    # Second pass : increase y values greedily
    for v = 1:nv(g)
        
        slack = Inf
        for n in neighbors(g, v)
            edge_slack = get_weight(w, v, n) - pseudonodes_dict[v].ybar - pseudonodes_dict[n].ybar
            if edge_slack < slack
                slack = edge_slack
            end
        end
        pseudonodes_dict[v].ybar += slack
        println("Node ", v, " slack ", slack, " ybar ", pseudonodes_dict[v].ybar)

    end
    """

    # Third pass : add tight edges to initial matching, and create singular trees for non-matched nodes

    for e in edges(g)
        if (
            get_weight(w, src(e), dst(e)) ==
            pseudonodes_dict[src(e)].ybar + pseudonodes_dict[dst(e)].ybar &&
            pseudonodes_dict[src(e)].label == positive &&
            pseudonodes_dict[dst(e)].label == positive
        )
            # Add edge to matching
            add_edge_matching!(match, e, get_weight(w, src(e), dst(e)))
            #println("----- Adding edge ", e)
            # Update node labels
            pseudonodes_dict[src(e)].label = free
            pseudonodes_dict[dst(e)].label = free

        end
    end

    for v = 1:nv(g)
        # If the node is still positive (not in the initial matching), is added to the forest as a singular tree
        if pseudonodes_dict[v].label == positive
            t = get_unit_tree(pseudonodes_dict[v])
            tree_id = add_tree!(forest, t)
            pseudonodes_dict[v].tree_id = tree_id


        end
    end

    return match, forest
end
