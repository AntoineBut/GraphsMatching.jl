using Revise, Graphs, GraphsMatching
edge_list = Edge.([(1,2), (2,3), (3,4), (1,4), (3,5), (3,6), (5,6)])
g = SimpleGraph(edge_list)
w = Dict{Edge, Float64}(Edge(1,2) => 68,
                        Edge(1,4) => 54,
                        Edge(2,3) => 93,
                        Edge(3,4) => 66,
                        Edge(3,5) => 95,
                        Edge(3,6) => 54,
                        Edge(5,6) => 11)
m = Matching(g,w)
ag = solve(m,g)
