using DataStructures
using Graphs
"""
minimum_weight_perfect_matching(g, w::Dict{Edge,Real})
minimum_weight_perfect_matching(g, w::Dict{Edge,Real}, cutoff)

Given a graph `g` and an edgemap `w` containing weights associated to edges,
returns a matching with the mimimum total weight among the ones containing
exactly `nv(g)/2` edges.

Edges in `g` not present in `w` will not be considered for the matching.

This function relies on the BlossomV.jl package, a julia wrapper
around Kolmogorov's BlossomV algorithm.

Eventually a `cutoff` argument can be given, to the reduce computational time
excluding edges with weights higher than the cutoff.

The returned object is of type `MatchingResult`.

In case of error try to change the optional argument `tmaxscale` (default is `tmaxscale=10`).
"""
function minimum_weight_perfect_matching end

function minimum_weight_perfect_matching(
    g::Graph,
    w::Dict{E,U},
    cutoff,
    kws...,
) where {U<:Real,E<:Edge}
    wnew = Dict{E,U}()
    for (e, c) in w
        if c <= cutoff
            wnew[e] = c
        end
    end
    return minimum_weight_perfect_matching(g, wnew; kws...)
end

function minimum_weight_perfect_matching(
    g::Graph,
    w::Dict{E,U};
    tmaxscale = 10.0,
) where {U<:AbstractFloat,E<:Edge}
    wnew = Dict{E,Int32}()
    cmax = maximum(values(w))
    cmin = minimum(values(w))

    tmax = typemax(Int32) / tmaxscale # /10 is kinda arbitrary,
    # hopefully high enough to not occur in overflow problems
    for (e, c) in w
        wnew[e] = round(Int32, (c - cmin) / max(cmax - cmin, 1) * tmax)
    end
    match = minimum_weight_perfect_matching(g, wnew)
    weight = zero(U)
    for i = 1:nv(g)
        j = match.mate[i]
        if j > i
            weight += w[E(i, j)]
        end
    end
    return MatchingResult(weight, match.mate)
end

function minimum_weight_perfect_matching(g::Graph, w::Dict{E,U}) where {U<:Integer,E<:Edge}
    m = BlossomV.Matching(nv(g))
    for (e, c) in w
        BlossomV.add_edge(m, src(e) - 1, dst(e) - 1, c)
    end
    BlossomV.solve(m)

    mate = fill(-1, nv(g))
    totweight = zero(U)
    for i = 1:nv(g)
        j = BlossomV.get_match(m, i - 1) + 1
        mate[i] = j <= 0 ? -1 : j
        if i < j
            totweight += w[Edge(i, j)]
        end
    end
    return MatchingResult(totweight, mate)
end

### blossomv start

# Useful types
include("blossomV_helpers/types.jl")

function add_edge!(m::Matching, e::Edge, w::U) where {U<:Real}
    push!(m.edges, e)
    push!(m.nodes, src(e))
    push!(m.nodes, dst(e))
    m.size += 2
    m.total_weight += w
end

get_weight(w::Dict{E,U}, u::T, v::T) where {U<:Real,E<:Edge,T<:Integer} =
    w[Edge(min(u, v), max(u, v))]
get_weight(w::Dict{E,U}, e::Edge) where {U<:Real,E<:Edge} = w[e]
get_edge(u::T, v::T) where {T<:Integer} = Edge(min(u, v), max(u, v))

function get_node_dual(node_id::T, forest::Forest{T,U}) where {T<:Integer,U<:Real}
    ps_node = forest.pseudonodes[node_id]
    if isnothing(ps_node.tree_id)
        return ps_node.ybar
    end
    epsilon_t = forest.trees[ps_node.tree_id]
    y_bar = ps_node.ybar

    offset_sign = 1
    ps_node.label == Labels.__minus && (offset_sign *= -1)

    return node.ybar + offset_sign * epsilon_t
end

function get_edge_slack(
    edge::Edge,
    forest::Forest{T,U},
    w::Dict{E,U},
) where {U<:Real,E<:Edge,T<:Integer}
    src = pseudonode_dict[edge.src]
    dst = pseudonode_dict[edge.dst]
    slack = get_weight(w, edge) - get_node_dual(src, forest) - get_node_dual(dst, forest)
    return slack
end

function add_tree!(f::Forest, t::AlternatingTree)
    tree_id = t.pseudo_nodes_ids[1]
    f.trees[tree_id] = t
    #add_vertex!(f.aux_graph)
    return tree_id
end

function find_tree(f::Forest)
    #TODO : change to iterate over trees in the forest instead of random selection
    # Select a random tree to avoid getting stuck choosing the same tree
    tree = f.trees[rand(keys(f.trees))]
    root = tree.pseudo_nodes_ids[1]
    return tree, root
end

# Initialization functions
include("blossomV_helpers/initialization.jl")

function primal_operation(
    tree::AlternatingTree,
    root::T,
    w::Dict{E,U},
    g::Graph{T},
    forest::Forest{T,U},
) where {T<:Integer,U<:Real,E<:Edge}
    # Find tight edge
    # TODO : use the priority queues to find a tight edge efficiently
    # For now, iterate over vertices in BFS manner

    pseudonodes_dict = forest.pseudonodes

    queue = Queue{T}()
    push!(queue, root)
    explored = Set{T}()

    while !isempty(queue)
        v = dequeue!(queue)
        if v in explored
            continue
        end
        push!(explored, v)
        for n in neighbors(g, v)
            if n in explored # if the neighbor has already been explored, edge has already been checked
                continue
            end
            edge = get_edge(v, n)
            slack = get_edge_slack(edge, forest, w)
            if get_weight(w, edge) ==
               get_node_dual(tree, pseudonodes_dict[v]) +
               get_node_dual(tree, pseudonodes_dict[n])
                # Tight edge found
                # dispatch to primal operation


                # Add nodes to the queue
                push!(queue, v)
                push!(queue, n)
                break
            end
        end
    end
    return nothing
end


function minimum_weight_perfect_matching_new(
    g::Graph{T},
    w::Dict{E,U},
) where {U<:Real,E<:Edge,T<:Integer}

    ########################################
    # 1 : Greedy initialization
    ########################################

    # Initialize base matching and y values conservatively

    matching, forest = greedy_initialization(g, w)

    # Nodes in the initial matching are "free nodes" with label free. Unmatched nodes are positive singular trees in the forest.

    ########################################
    # 2 : Initialize auxilary graph : TODO : Implement usage of auxilary graph
    ########################################

    # Iterate over edges of the graph, add edges to the auxilary graph if they connect nodes in different trees

    ########################################
    # 3 : Main loop : While there are free nodes in the forest -> primal and dual operations
    ########################################

    while length(forest.trees) > 0
        # Select a free node v in the forest, and the tree it is the root of
        tree, root = find_tree(forest)
        # Iterate over nodes and their neighbors, find a tight edge, perform one of the primal operations
        # If no tight edge is found, perform a dual operation

        # Primal operation
        primal_operation(tree, root, w, g, forest)

        # Dual operation
        #dual_operation(tree, pseudonodes_dict, w)


    end

    return matching
end
