mutable struct Matching{T,U}
    edges::Set{Edge} # edges in the matching
    nodes::Set{T} # nodes in the matching 
    total_weight::U # total weight of the matching
    size::T # number of nodes in the matching
    y::Array{U,1} #dual variable len(y)=nv(g)
end

@enum Labels positive = 1 negative = -1 free = 0

mutable struct PseudoNode{T,U}
    label::Labels
    node_id::T
    ybar::U
    tree_id::Union{T,Nothing} # Tree id if this node is in a tree
    is_blossom::Bool
    parent::Union{T,Nothing} # Parent Blossom if this node is an inner node of a blossom
end

mutable struct EdgeData{T,U}
    # Current src and dst nodes are implicit in the edge data: they are the current src and dst PseudoNodes
    orig_src::T # Original node
    orig_dst::T # Original node
    slack_bar::U # Currently unused
    #TODO :  Lazy Slack value slack_bar = weight - ybar[src] - ybar[dst]. Needs to be updated a node label changes.
end

mutable struct TreeEdge
    pq_pp::PriorityQueue{PseudoNode,Int}
    pq_pm::PriorityQueue{PseudoNode,Int}
    pq_mp::PriorityQueue{PseudoNode,Int}
end

mutable struct AlternatingTree{T,U}
    backbone::SimpleDiGraph{T} # Directed graph -> in-neighbors are parents, out-neighbors are children
    pseudo_nodes_ids::Vector{T}
    epsilon_t::U # Accumulates dual updates on this tree. 
    current_edge::Union{TreeEdge,Nothing}
    pq_pp::PriorityQueue{PseudoNode,Int}
    pq_pz::PriorityQueue{PseudoNode,Int}
    pq_m::PriorityQueue{PseudoNode,Int}
end

mutable struct AuxGraph{E,T,U}
    g::SimpleGraph{T}
    #edge_slack_bar::Dict{Edge,U}
    #edge_tree::Dict{E,TreeEdge}
end

# Forest data structure. Stores all main data structures for the algorithm.
mutable struct Forest{T,U}
    #aux_graph::AuxGraph{E,T,U}
    pseudonodes::Dict{T,PseudoNode{T,U}}
    edgedata::Dict{Tuple{T,T},EdgeData{T,U}}
    trees::Dict{T,AlternatingTree{T,U}}
    blossoms::Dict{T,Set{T}}
end
