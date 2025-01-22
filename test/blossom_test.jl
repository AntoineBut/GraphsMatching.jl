using Revise
using Graphs, Test
using GraphsMatching
using GraphsMatching:
    make_forest,
    empty_matching,
    Labels,
    positive,
    negative,
    free,
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
        Edge(1, 2) => 60,
        Edge(1, 4) => 50,
        Edge(2, 3) => 100,
        Edge(3, 4) => 100,
        Edge(3, 5) => 100,
        Edge(3, 6) => 40,
        Edge(5, 6) => 10,
        Edge(5, 7) => 20,
        Edge(4, 7) => 40,
        Edge(4, 8) => 40,
        Edge(7, 8) => 10,
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

        MATCHED_NODES = [5, 6, 7, 8]
        UNMATCHED_NODES = [1, 2, 3, 4]
        #@test (nv(forest.aux_graph)) == length(UNMATCHED_NODES)
        @test length(forest.trees) == length(UNMATCHED_NODES) # all unmatched nodes are in the forest

        #@test length(forest.aux_graph_edges) == 0
        #@test (ne(forest.aux_graph)) == 0 # no edges in the aux graph yet

        @test length(match.edges) == 2
        @test length(match.nodes) == 4
        @test match.total_weight == 20
        @test match.size == 4
        @test match.mate[7] == 8
        @test match.mate[8] == 7
        @test match.mate[5] == 6
        @test match.mate[6] == 5


        @test nv(forest.trees[1].backbone) == 1
        @test nv(forest.trees[2].backbone) == 1
        @test nv(forest.trees[3].backbone) == 1
        @test nv(forest.trees[4].backbone) == 1 # The 4 unmatched nodes are in the forest, each in a singular tree
        @test length(forest.trees[1].psnode_to_backbone) == 1

        nodes_in_forest = Set{Int}()
        for t in values(forest.trees)
            for n in keys(t.psnode_to_backbone)
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
        println("--------------------")
        println("Testing get_node_dual and get_edge_slack")

        # Test get_node_dual
        for i = 1:n
            @test GraphsMatching.get_node_dual(i, forest) == forest.pseudonodes[i].ybar
        end

        # Test get_edge_slack
        for e in edge_list
            @test GraphsMatching.get_edge_slack(e, forest, w) ==
                  w[e] - forest.pseudonodes[e.src].ybar - forest.pseudonodes[e.dst].ybar
        end


    end

    @testset "primal operations" begin
        match, forest = GraphsMatching.greedy_initialization(g, w)

        ### Test primal_grow
        edge = Edge(3, 6)
        grown_tree_id = forest.pseudonodes[3].tree_id

        GraphsMatching.primal_grow(edge, forest, match, g, w)
        @test length(forest.trees) == 4
        @test length(match.edges) == 2
        @test length(match.nodes) == 4
        @test forest.trees[grown_tree_id].root == 3
        @test forest.trees[grown_tree_id].psnode_to_backbone[6] == 2
        @test forest.trees[grown_tree_id].psnode_to_backbone[5] == 3

        @test forest.trees[grown_tree_id].backbone_to_psnode[2] == 6
        @test forest.trees[grown_tree_id].backbone_to_psnode[3] == 5

        new_backbone = forest.trees[grown_tree_id].backbone
        @test nv(new_backbone) == 3
        @test ne(new_backbone) == 2
        @test outneighbors(new_backbone, 1) == [2]
        @test outneighbors(new_backbone, 2) == [3]
        @test forest.pseudonodes[6].label == negative
        @test forest.pseudonodes[5].label == positive

        edge = Edge(4, 8)
        grown_tree_id = forest.pseudonodes[4].tree_id
        GraphsMatching.primal_grow(edge, forest, match, g, w)
        @test forest.trees[grown_tree_id].psnode_to_backbone[8] == 2
        @test forest.trees[grown_tree_id].psnode_to_backbone[7] == 3
        @test forest.trees[grown_tree_id].backbone_to_psnode[2] == 8
        @test forest.trees[grown_tree_id].backbone_to_psnode[3] == 7
        @test nv(forest.trees[grown_tree_id].backbone) == 3
        @test ne(forest.trees[grown_tree_id].backbone) == 2


        ### Test primal_augment
        edge = Edge(5, 7)
        GraphsMatching.primal_augment(edge, forest, match, g, w)
        @test length(match.edges) == 3
        @test length(match.nodes) == 6
        @test match.mate[5] == 7
        @test match.mate[7] == 5
        @test match.mate[6] == 3
        @test match.mate[3] == 6 # Matching has been updated correctly : (5, 7) added, others swapped

        @test length(forest.trees) == 2
        @test forest.trees[1].root == 1
        @test forest.trees[2].root == 2 # 2 singular trees left : 1 and 2

        @test forest.pseudonodes[3].label == free
        @test forest.pseudonodes[5].label == free
        @test forest.pseudonodes[6].label == free
        @test forest.pseudonodes[4].label == free
        @test forest.pseudonodes[7].label == free
        @test forest.pseudonodes[8].label == free # all nodes in augmented trees are free


    end


end
