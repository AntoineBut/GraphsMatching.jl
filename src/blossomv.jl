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

### blossomv start

@enum Labels begin
    __minus = -1;
    __zero = 0;
    __plus = 1;
end

mutable struct PsuedoNode{T}
    label::Labels
    ybar::T
    vertex_number::Int;
    is_blossom::Bool;
end
function PsuedoNode(T::Type, i::Int)
    return PsuedoNode(__zero, T(0), i, false)
end

mutable struct PsuedoEdge{T}<:AbstractEdge{T}
    weight::T
    src::Int
    dest::Int
    slack_bar::T
end

mutable struct Matching{T}
    adj::Array{Array{Int, 1}, 1};
    edges::Dict{Graphs.SimpleGraphs.SimpleEdge, T};
    nodes::Dict{Int, PsuedoNode{T}}
    y::Array{T,1}; #dual variable len(y)=nv(g)
end
function Matching(g::Graph, weights::Dict{Graphs.SimpleGraphs.SimpleEdge, T}) where {T<:Number}
    d = Dict([i=>PsuedoNode(T,i) for i in 1:length(g.fadjlist)])
    return Matching(g.fadjlist, weights, d, zeros(T, nv(g)));
end

mutable struct Tree
    pq_pp::PriorityQueue{PsuedoNode, Int};
    pq_pz::PriorityQueue{PsuedoNode, Int};
    pq_m::PriorityQueue{PsuedoNode, Int};
    nodes::Array{PsuedoNode, 1}; #Union of pq_pp/pz/m
    parent::Dict{PsuedoNode, PsuedoNode};
    children::Dict{PsuedoNode, Array{PsuedoNode,1}};
    current_edge;
    del;
end
Tree() = Tree(PriorityQueue{PsuedoNode, Int}(),
              PriorityQueue{PsuedoNode, Int}(),
              PriorityQueue{PsuedoNode, Int}(),
              Array{PsuedoNode,1}(),
              Dict{PsuedoNode, PsuedoNode}(),
              Dict{PsuedoNode, Array{PsuedoNode, 1}}(),
              -Inf,
              0.0)

mutable struct TreeEdge
    pq_pp::PriorityQueue{PsuedoNode, Int};
    pq_pm::PriorityQueue{PsuedoNode, Int};
    pq_mp::PriorityQueue{PsuedoNode, Int};
end
TreeEdge() = TreeEdge(PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}(), PriorityQueue{PsuedoNode, Int}())

mutable struct AuxGraph
    nodes::Array{Tree, 1};
    edges::Dict{Tuple{Tree, Tree}, TreeEdge};
    tree_pointers::Dict{Int, Union{Tree, Nothing}}; #vertex=>Tree
    ps_nodes::Dict{Int, PsuedoNode}; #vertex num => PS
end
function AuxGraph()
    return AuxGraph(Array{Tree,1}(), Dict{Tuple{Tree,Tree}, TreeEdge}(), Dict{Int, Tree}(), Dict{Int, PsuedoNode}());
end

get_label(n::PsuedoNode) = n.label
set_label!(n::PsuedoNode, val::Labels) = n.label = val;
get_ybar(n::PsuedoNode) = n.ybar;
get_parent(t::Tree, n::PsuedoNode) = t.parent[n];
set_parent!(t::Tree, n::PsuedoNode, v::PsuedoNode) = (t.parent[n] = v)
get_children(t::Tree, n::PsuedoNode) = t.children[n];
function add_child!(t::Tree, parent::PsuedoNode, child::Int)
    if !(child in get_children(t, parent))
        append!(t.children[parent])
    end
    set_parent!(t, parent, child)
    return;
end
function set_ybar(n::PsuedoNode{T}, val::T) where {T<:Number}
    n.ybar = val;
    return;
end

src(e::PsuedoEdge) = e.src;
dest(e::PsuedoEdge) = e.dest;
get_weight(e::PsuedoEdge) = e.weight;
get_tmp_slack(e::PsuedoEdge) = e.slack_bar;

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

function merge_trees(t1::Tree, t2::Tree)

end
function add_psuedonode_to_tree!(t::Tree, ps::PsuedoNode)
    push!(t.nodes, ps)
    # TODO this should be smarter
    # for e in boundary(ps)
    #     bdy_node = e.src == ps ? e.dst : e.src;
    #     if ((get_tree(ps) === get_tree(bdy_node)) &&
    #         (get_label(ps) == get_label(bdy_node) == __plus))
    #         # add e to tree.pq_pp
    #     elseif (get_label(ps) == __plus && get_label(bdy_node) == __zero)
    #         # add e to tree.pq_pz
    #     elseif (get_label(ps) == __minus)
    #         # add ps to tree.pq_m
    #     end
    # end
    return;
end

function get_slack(m::Matching, verts::Tuple)
    return get_weight(m, verts) - m.y[verts[1]] - m.y[verts[2]]
end
function is_tight(m::Matching, verts::Tuple)
    return get_slack(m, verts) > 0 ? false : true;
end
is_tight(m::Matching, e::Edge) = is_tight(m, (e.src, e.dst))

function update_dual_var!(m::Matching, ag::AuxGraph, n::PsuedoNode)
    if (ag.tree_pointers[n.vertex_number] isa Tree)
        m.y[n.vertex_number] += Int(get_label(n))*ag.tree_pointers[n.vertex_number].del
    end
    return;
end
function get_dual_var(m::Matching, n::PsuedoNode)
    return m.y[n.vertex_number];
end

function get_tree(ag::AuxGraph, v::Int)
    return ag.tree_pointers[v];
end
function set_tree!(ag::AuxGraph, v::Int, t::Union{Tree, Nothing})
    if !(t isa Nothing)
        add_psuedonode_to_tree!(t, ag.ps_nodes[v])
    end
    ag.tree_pointers[v] = t;
    return;
end

function set_ybar(ag::AuxGraph, v::Int, weight::T) where {T<:Number}
    tree = get_tree(ag, v)
    #tree.
end

function greedy_initialization!(m::Matching, g::Graph, aux_graph::AuxGraph)
    #initialize y values conservatively
    y = zeros(length(g.fadjlist));
    for v in 1:length(g.fadjlist)
        min_weight = Inf;
        for n in g.fadjlist[v]
            tmp_weight = get_weight(m, (v,n));
            if (tmp_weight < min_weight)
                min_weight = tmp_weight;
            end
        end
        #set_ybar(aux_graph, v, min_weight/2)
        y[v] = min_weight/2;
    end

    #then increase y values greedily
    for v1 in 1:length(g.fadjlist)
        slack = Inf;
        for v2 in g.fadjlist[v1]
            edge_slack = get_weight(m, (v1,v2)) - y[v1] - y[v2]
            if (edge_slack < slack)
                slack = edge_slack;
            end
        end
        y[v1] += slack;
    end
    m.y = y; #todo: make better

    # set tree(s)
    set_tree!(aux_graph, 1, aux_graph.nodes[1])
end

function primal_grow(edge, m::Matching, aux_graph::AuxGraph)
    set_label!(m.nodes[edge.dst], __minus);
    set_tree!(aux_graph, edge.dst, get_tree(aux_graph, edge.src));
    add_child!(m.nodes[edge.src], edge.dst);
    set_parent!(m.nodes[edge.dst], edge.src);
    #flip signs along branch
end
function primal_augment(edge, m::Matching, aux_graph::AuxGraph)
    src_tree = get_tree(aux_graph, edge.src)
    dst_tree = get_tree(aux_graph, edge.dst)
    #TODO merge_trees(src_tree, dst_tree)

    if (dst_tree isa Tree) #i.e., not just a lonely single node
        for ps in dst_tree.nodes
            set_tree!(aux_graph, ps.vertex_number, src_tree);
        end
    end
    
    for ps in src_tree.nodes
        set_label!(ps, __zero)
    end
    return;
end
function primal_shrink()
end

function primal_update(m::Matching, aux_graph::AuxGraph)
    for edge in keys(m.edges)
        if (is_tight(m, edge))
            # Grow
            if (get_label(m.nodes[edge.src]) == __plus && get_label(m.nodes[edge.dst]) == __zero)
                primal_grow(edge, m, aux_graph);
            elseif (get_label(m.nodes[edge.src]) == get_label(m.nodes[edge.dst]) == __plus)
                if (get_tree(aux_graph, edge.src) !== get_tree(aux_graph, edge.dst))
                    primal_augment(edge, m, aux_graph);
                else
                    #primal_shrink(m, aux_graph);
                end
            end
        end
    end
end

function dual_grow_update(m::Matching, aux_graph::AuxGraph, tree::Tree)
    edges = []
    for u in tree.nodes
        for v in m.adj[u.vertex_number]
            if (get_label(u), get_label(m.nodes[v])) == (__plus, __zero)
                push!(edges, Edge(u.vertex_number, v))
            end
        end
    end
    slacks = [get_slack(m, (edge.src, edge.dst)) for edge in edges];
    return slacks;
end

function dual_augment_update(m::Matching, aux_graph::AuxGraph, tree::Tree)
    edges = []
    for u in tree.nodes
        for v in m.adj[u.vertex_number]
            if ((get_tree(aux_graph, v) !== tree) &&
                (get_label(u), get_label(m.nodes[v])) == (__plus, __plus))
                push!(edges, Edge(u.vertex_number, v))
            end
        end
    end
    slacks = [get_slack(m, (edge.src, edge.dst))-get_tree(aux_graph, v).del for edge in edges];
    return slacks;
end

function dual_shrink_update(m::Matching, aux_graph::AuxGraph, tree::Tree)
    edges = []
    for u in tree.nodes
        for v in m.adj[u.vertex_number]
            if ((get_tree(aux_graph, v) === tree) &&
                (get_label(u), get_label(m.nodes[v])) == (__plus, __plus))
                push!(edges, Edge(u.vertex_number, v))
            end
        end
    end
    slacks = [get_slack(m, (edge.src, edge.dst))/2 for edge in edges];
    return slacks;
end

function dual_expand_update(m::Matching, aux_graph::AuxGraph, tree::Tree)
    slacks = []
    for u in tree.nodes
        if (get_label(u) == __minus && is_blossom(u))
            append!(slacks, get_dual_var(m, u))
        end
    end
    return slacks;
end

function dual_noop_update(m::Matching, aux_graph::AuxGraph, tree::Tree)
    edges = []
    for u in tree.nodes
        for v in m.adj[u.vertex_number]
            if ((get_tree(aux_graph, v) !== tree) &&
                (get_label(u), get_label(m.nodes[v])) == (__plus, __minus))
                push!(edges, Edge(u.vertex_number, v))
            end
        end
    end
    slacks = [get_slack(m, (edge.src, edge.dst))+get_tree(aux_graph, v).del for edge in edges];
    return slacks;
end

function dual_update(m::Matching, aux_graph::AuxGraph)
    for tree in aux_graph.nodes
        tree.del = 0.0; #reset dual updates
    end
    for tree in aux_graph.nodes
        grow_slacks = dual_grow_update(m, aux_graph, tree);
        augment_slacks = dual_augment_update(m, aux_graph, tree);
        shrink_slacks = dual_shrink_update(m, aux_graph, tree);
        expand_slacks = dual_expand_update(m, aux_graph, tree);
        noop_slacks = dual_noop_update(m, aux_graph, tree);
        slacks = cat(grow_slacks,
                     augment_slacks,
                     shrink_slacks,
                     expand_slacks,
                     noop_slacks, dims=1);
        if (length(slacks) > 0)
            tree.del = minimum(slacks);
        end
    end
    for n in keys(m.nodes)
        update_dual_var!(m, aux_graph, m.nodes[n])
    end
end

function solve(m::Matching, g::Graph)
    aux_graph = AuxGraph();
    push!(aux_graph.nodes, Tree())
    for i in 1:length(g.fadjlist)
        aux_graph.ps_nodes[i] = PsuedoNode(typeof(m.edges.vals[1]), i)
    end
    greedy_initialization!(m, g, aux_graph);

    # initialize tree pointers
    set_tree!(aux_graph, 1, aux_graph.nodes[1])
    set_label!(m.nodes[1], __plus);
    for v in 2:length(g.fadjlist)
        pushfirst!(aux_graph.nodes, Tree())
        set_tree!(aux_graph, v, aux_graph.nodes[1])
        set_label!(m.nodes[v], __plus);
    end
    print([m.nodes[v] for v in 1:length(g.fadjlist)])
    println(get_tree(aux_graph, 1))

    primal_update(m, aux_graph)
    dual_update(m, aux_graph)
    primal_update(m, aux_graph)
    dual_update(m, aux_graph)
    primal_update(m, aux_graph)
    dual_update(m, aux_graph)
    primal_update(m, aux_graph)
    return aux_graph
end
