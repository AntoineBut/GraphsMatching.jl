using DataStructures
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
    g::Graph, w::Dict{E,U}, cutoff, kws...
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
    g::Graph, w::Dict{E,U}; tmaxscale=10.0
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
    for i in 1:nv(g)
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
    for i in 1:nv(g)
        j = BlossomV.get_match(m, i - 1) + 1
        mate[i] = j <= 0 ? -1 : j
        if i < j
            totweight += w[Edge(i, j)]
        end
    end
    return MatchingResult(totweight, mate)
end

mutable struct Matching{T}
    adj::Array{Array{Int, 1}, 1};
    edges::Dict{Graphs.SimpleGraphs.SimpleEdge, T};
    y::Array{T,1}; #dual variable len(y)=nv(g)
end
function Matching(g::Graph, weights::Dict{Graphs.SimpleGraphs.SimpleEdge, T}) where {T<:Number}
    return Matching(g.fadjlist, weights, zeros(T, nv(g)));
end
function get_weight(m::Matching, verts::Tuple)
    Edge(verts[1],verts[2]) in keys(m.edges) ? m.edges[Edge(verts[1], verts[2])] : m.edges[Edge(verts[2],verts[1])];
end
function allocate_adj(m::Matching, nv::Int)
    m.adj = [Array{Int, 1}() for i in 1:nv];
end

function get_match(m::Matching{T}, edge_num::Integer) where {T<:Number}
    return #m[i].....
end

function add_edge!(m::Matching{T}, edge_src, edge_dest, weight::T) where {T<:Number}
    m.edges[Edge(edge_src, edge_dest)] = weight;
    return;
end

@enum Labels begin
    __minus = -1;
    __zero = 0;
    __plus = 1;
end
mutable struct PsuedoNode{T}
    label::Labels
    ybar::T
end
label(n::PsuedoNode) = Int(n.label)
get_ybar(n::PsuedoNode) = n.ybar;
function set_ybar(n::PsuedoNode{T}, val::T) where {T<:Number}
    n.ybar = val;
    return;
end

mutable struct PsuedoEdge{T}<:AbstractEdge{T}
    weight::T
    src::Int
    dest::Int
    slack_bar::T
end
src(e::PsuedoEdge) = e.src;
dest(e::PsuedoEdge) = e.dest;
get_weight(e::PsuedoEdge) = e.weight;
get_tmp_slack(e::PsuedoEdge) = e.slack_bar;

mutable struct Tree
     pq_pp::PriorityQueue{PsuedoNode, Int};
    pq_pz::PriorityQueue{PsuedoNode, Int};
    pq_m::PriorityQueue{PsuedoNode, Int};
    current_edge;
end
Tree() = Tree(PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}(), -Inf)

mutable struct TreeEdge
    pq_pp::PriorityQueue{PsuedoNode, Int};
    pq_pm::PriorityQueue{PsuedoNode, Int};
    pq_mp::PriorityQueue{PsuedoNode, Int};
end
TreeEdge() = TreeEdge(PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}())

mutable struct AuxGraph
    nodes::Array{Tree, 1};
    edges::Dict{Tuple{Tree, Tree}, TreeEdge};
    tree_pointers::Dict{Int, Tree}; #vertex=>Tree
end
function AuxGraph()
    return AuxGraph(Array{Tree,1}(), Dict{Tuple{Tree,Tree}, TreeEdge}(), Dict{Int, Tree}());
end
function get_tree(ag::AuxGraph, v::Int)
    return ag.tree_pointers[v];
end
function set_ybar(ag::AuxGraph, v::Int, weight::T) where {T<:Number}
    tree = get_tree(ag, v)
    tree.
end

function greedy_initialization(m::Matching, g::Graph, aux_graph::AuxGraph)
    for v in 1:length(g.fadjlist)
        min_weight = Inf;
        for n in g.fadjlist[v]
            tmp_weight = get_weight(m, (v,n));
            if (tmp_weight < min_weight)
                min_weight = tmp_weight;
            end
        end
        set_ybar(aux_graph, v, min_weight/2)
    end
end

function solve(m::Matching, g::Graph)
    aux_graph = AuxGraph();
    push!(aux_graph.nodes, Tree())
    greedy_initialization(m, g, aux_graph);
end
