using Graphs, Test
using GraphsMatching
using GraphsMatching:
    make_forest, 
    empty_matching,
    Labels,
    positive,
    greedy_initialization,
    Forest,
    AlternatingTree,
    PseudoNode,
    EdgeData,
    get_unit_tree

@testset "blossomV" begin
    edge_list =
        Edge.([
            (1, 2),
            (2, 3),
            (3, 4),
            (1, 4),
            (3, 5),
            (3, 6),
            (5, 6),
            (5, 7),
            (4, 7),
            (4, 8),
            (7, 8),
        ])
    g = SimpleGraph(edge_list)
    n = nv(g)
    e = ne(g)
    w = Dict{Edge,Float64}(
        Edge(1, 2) => 68,
        Edge(1, 4) => 54,
        Edge(2, 3) => 93,
        Edge(3, 4) => 66,
        Edge(3, 5) => 95,
        Edge(3, 6) => 54,
        Edge(5, 6) => 11,
        Edge(5, 7) => 20,
        Edge(4, 7) => 40,
        Edge(4, 8) => 10,
        Edge(7, 8) => 42,
    )
    @testset "code format and static analysis" begin
        #TODO : add JuliaFormatter, JET and Aqua tests
    end

    @testset "internals" begin

        forest = make_forest(g, w)
        pseudonodes_dict = forest.pseudonodes
        edgedata_dict = forest.edgedata
        blossoms_dict = forest.blossoms
        trees_dict = forest.trees
        @test length(blossoms_dict) == 0
        @test length(pseudonodes_dict) == n
        @test pseudonodes_dict[1].node_id == 1
        @test pseudonodes_dict[1].label == positive
        @test pseudonodes_dict[1].ybar == Inf

        @test length(edgedata_dict) == e
        @test length(trees_dict) == 0

        empty_matching1 = empty_matching(0, 0.0)
        @test length(empty_matching1.edges) == 0
        @test length(empty_matching1.nodes) == 0
        @test empty_matching1.total_weight == 0
        @test empty_matching1.size == 0


        println("Testing greedy initialization")
        match, forest = GraphsMatching.greedy_initialization(g, w)

        # Test all nodes have a ybar value that is no longer Inf
        for i = 1:n
            @test forest.pseudonodes[i].ybar != Inf
        end

        MATCHED_NODES = [4, 5, 6, 8]
        UNMATCHED_NODES = [1, 2, 3, 7]
        #@test (nv(forest.aux_graph)) == length(UNMATCHED_NODES)
        @test length(forest.trees) == length(UNMATCHED_NODES) # all unmatched nodes are in the forest

        #@test length(forest.aux_graph_edges) == 0
        #@test (ne(forest.aux_graph)) == 0 # no edges in the aux graph yet


        @test nv(forest.trees[1].backbone) == 1
        @test nv(forest.trees[2].backbone) == 1
        @test nv(forest.trees[3].backbone) == 1 # The 3 unmatched nodes are in the forest, each in a singular tree
        @test length(forest.trees[1].pseudo_nodes_ids) == 1

        nodes_in_forest = Set{Int}()
        for t in values(forest.trees)
            for n in t.pseudo_nodes_ids
                push!(nodes_in_forest, n)
            end
        end
        @test nodes_in_forest == Set(UNMATCHED_NODES) # all unmatched nodes are in the forest






        println("--------------------")
        println("Testing find_tree")
        # The tree can be any of the 3 singular trees. Test that the root is in UNMATCHED_NODES 10 times
        for i = 1:10
            tree, root = GraphsMatching.find_tree(forest)
            @test root in UNMATCHED_NODES
        end


    end


end
