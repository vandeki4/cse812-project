#CSE 812 Semester Project
#Distributed Minimum Spanning Trees
#Ian Murray
#April 8th, 2020

#singlethread.jl
#This is the single-threaded algorithm outlined in
#http://viswa.engin.umich.edu/wp-content/uploads/sites/169/2016/12/5.pdf
#It is used as the baseline to test agaisnt the performance of the IMP graph

using LightGraphs
using GraphPlot

#---Function reduce---
#Implements one iteration of the MDST algorithm
#\param G: The MST of the origin graph
#\param origin: The origin graph
#\return Bool: Whether or node a change was found and made
function reduce!(G, origin)::Bool
    for u in vertices(G)
        for v in vertices(G)
            for w in vertices(G)
                #Condtion 1: Disjoint nodes
                if u != v && u != w && v != w
                    #Condition 2: Edge between u and v
                    if has_edge(origin, v, w)
                        #Condition 3: deg(U) > Threshold
                        if degree(G, u) >= Δ(G) - log2(nv(G))
                            #Condition 4: deg(v) <= deg(u) - 2, deg(w) <= deg(u) - 2
                            if degree(G, v) <= degree(G, u) - 2 && degree(G, w) <= degree(G, u) - 2
                                #Condtion 5: u lies on the path from v - w
                                path = a_star(G, v, w)
                                onPath = false #Indicator variable
                                uPrime = nothing #A neighbor of vertex u
                                for edge in path
                                    if dst(edge) == u
                                        onPath = true
                                        uPrime = edge
                                        break
                                    end
                                end
                                if onPath #Condition 5 (continued)
                                    #Decrease max degree by inducing a cycle and removing an
                                    #outgoing edge from the higher-degree node
                                    add_edge!(G, v, w)
                                    rem_edge!(G, uPrime)

                                    println("---Successful Reduction---")
                                    println("\tu = " * string(u) * ", v = " * string(v) * ", w = " * string(w))
                                    println("\tAdded Edge: " * string(v) * " => " * string(w))
                                    println("\tRemoved " * string(uPrime))
                                    return true
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("---No Reductions Found---")
    return false
end

#---Phase 0: Start Program Execution---
start = time()
numVertices = 250
numEdges = numVertices * 3

#---Phase 1: Graph Creation---
origin = Graph(numVertices, numEdges)
# add_edge!(origin, 1, 2)
# add_edge!(origin, 2, 3)
# add_edge!(origin, 2, 4)
# add_edge!(origin, 1, 5)
# add_edge!(origin, 5, 6)
# add_edge!(origin, 6, 7)
# add_edge!(origin, 1, 8)
# add_edge!(origin, 8, 9)
# add_edge!(origin, 1, 10)
# add_edge!(origin, 10, 11)
# add_edge!(origin, 6, 8)
#Plot Orginial Graph
#gplot(origin, nodelabel = 1:nv(G))

#---Phase 2: Generate MST---
mstEdges = prim_mst(origin)
G = Graph(nv(origin))
for edge in mstEdges
    add_edge!(G, edge)
end
#Plot MST
#gplot(G, nodelabel = 1:nv(G))

#---Phase 3: Generate MDST---
moreWork = true
while moreWork
    println("\nMax Degree of G: " * string(Δ(G)))

    #If reduce returned false, there weren't any changes made and we are done
    if reduce!(G, origin) == false
        global moreWork = false
    end
end

#---Phase 4: Report Results---
#Report total runtime benchmark
elapsed = time() - start
println("\nElapsed Time: " * string(elapsed))
println("Max Degree Reduced from " * string(Δ(origin)) * " to " * string(Δ(G)))
#Plot MDST
gplot(G, nodelabel = 1:nv(G))
