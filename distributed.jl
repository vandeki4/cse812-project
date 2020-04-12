using LightGraphs
using SparseArrays
using GraphPlot
using Compose

push!(LOAD_PATH, "./")
using SparseEdmondsKarp

function constrained_q_match(bipartite_graph, q::Int, bundle_lookup)
  n::Int = nv(bipartite_graph) / 2
  flow_graph = SimpleDiGraph(3n+2)
  
  I::Vector{Int} = []
  J::Vector{Int} = []
  V::Vector{Int} = []

  for edge in edges(bipartite_graph)
    # If we have already added an edge, don't add it to the sparse flow capacity a second time.
    # add_edge is smart enough to tell us if we already have this edge
    if (add_edge!(flow_graph, edge))
      push!(I, src(edge))
      push!(J, dst(edge))
      push!(V, 1)
    end
    if (add_edge!(flow_graph, dst(edge), 3n+2))
      push!(I, dst(edge))
      push!(J, 3n+2)
      push!(V, q)
    end

    e = (2n + bundle_lookup(src(edge)), src(edge))
    if (add_edge!(flow_graph, e[1], e[2]))
      push!(I, e[1])
      push!(J, e[2])
      push!(V, 1)
    end
    if (add_edge!(flow_graph, 3n+1, e[1]))
      push!(I, 3n+1)
      push!(J, e[1])
      push!(V, q)
    end
  end
  
  capacity_matrix = sparse(I, J, V, 3n+2, 3n+2)
  return SparseEdmondsKarp.edmonds_karp_sparse_impl(flow_graph, 3n+1, 3n+2, capacity_matrix)
end

function improve!(graph, tree, high_degree, low_degree)
  retain_vertex = degree(tree) .< high_degree
  vertices_to_retain = filter(n -> n > 0, retain_vertex .* eachindex(retain_vertex))
  println("Retaining $(size(vertices_to_retain, 1)) vertices")

  # Remove edges adjacent to high degree vertices
  is_high_edge(graph, edge, deg) = degree(graph, src(edge)) >= deg || degree(graph, dst(edge)) >= deg

  branches_graph = SimpleGraph(nv(tree))

  for edge in edges(tree)
    if !is_high_edge(tree, edge, high_degree)
      add_edge!(branches_graph, edge)
    end
  end

  # Fast vertex -> branch lookup
  branch_components = connected_components(branches_graph)
  branch_lookup = zeros(Int, nv(graph))
  for (i, comp) in enumerate(branch_components)
    for vert in comp
      branch_lookup[vert] = i
    end
  end

  # Determine which branches are leaf branches
  n_parents_for_branch::Array{Int, 1} = zeros(size(branch_components, 1))
  parent_for_branch::Array{Int, 1} = zeros(size(branch_components, 1))
  rep_vertex_for_branch::Array{Int, 1} = zeros(size(branch_components, 1))
  for edge in edges(tree)
    # we want to examine the removed edges, not the ones in branches
    if has_edge(branches_graph, edge)
      continue
    end

    # parent is only meaningful for leaf branches,
    # but we don't know what ones are leaf branches yet
    if (branch_lookup[src(edge)] != 0)
      n_parents_for_branch[branch_lookup[src(edge)]] += 1
      parent_for_branch[branch_lookup[src(edge)]] = dst(edge)
      rep_vertex_for_branch[branch_lookup[src(edge)]] = src(edge)
    end

    if (branch_lookup[dst(edge)] != 0)
      n_parents_for_branch[branch_lookup[dst(edge)]] += 1
      parent_for_branch[branch_lookup[dst(edge)]] = src(edge)
      rep_vertex_for_branch[branch_lookup[dst(edge)]] = dst(edge)
    end
  end
  # Helper functions for leaf branches and bundle structures
  is_leafbranch(branch_id) = n_parents_for_branch[branch_id] == 1
  is_vert_of_leafbranch(vertex) = is_leafbranch(branch_lookup[vertex])
  representative_vertex_of_leafbranch(vertex) = rep_vertex_for_branch[branch_lookup[vertex]]
  bundle_id_of_vertex(vertex) = parent_for_branch[branch_lookup[vertex]]

  # Imp graph needs to consider those edges that are not part of the tree, but are part of the graph
  g_sub_tree = Iterators.filter(
    edge -> !has_edge(tree, edge),
    edges(graph)
  )

  imp_graph = SimpleDiGraph(2 * nv(graph))
  
  I::Vector{Int} = []
  J::Vector{Int} = []
  V::Vector{Int} = []

  for edge in g_sub_tree
    @assert !has_edge(tree, edge)
    is_good = !is_high_edge(branches_graph, edge, low_degree) # TODO: is this leaf_branches? Might be tree

    # in different branches
    is_good &= branch_lookup[src(edge)] != branch_lookup[dst(edge)]
    # with at least one being a leaf branch
    is_good &= is_vert_of_leafbranch(src(edge)) || is_vert_of_leafbranch(dst(edge))
    if (is_good)
      # for each of the two vertices of this good edge, check to see if it is in a leaf branch.
      # if it is, it contributes to the imp graph. Note that this also decides incomming/outgoing edges.
      if (is_vert_of_leafbranch(src(edge)))
        if add_edge!(imp_graph, branch_lookup[src(edge)], dst(edge) + nv(graph))
          push!(I, branch_lookup[src(edge)])
          push!(J, dst(edge))
          push!(V, src(edge))
        end
      end
      if (is_vert_of_leafbranch(dst(edge)))
        if add_edge!(imp_graph, branch_lookup[dst(edge)], src(edge) + nv(graph))
          push!(I, branch_lookup[dst(edge)])
          push!(J, src(edge))
          push!(V, dst(edge))
        end
      end
    end
  end
  # Reverse lookup, so we can go from an imp graph edge back to the edge that originally created it
  original_edge_from_imp = sparse(I,J,V, nv(graph), nv(graph))

  # Todo: change to the distributed flow alg. described in the paper
  total, flow = constrained_q_match(imp_graph, 1, bundle_id_of_vertex)
  for (i,j,f) in zip(findnz(flow)...)
    if f <= 0
      continue
    end

    if i <= nv(graph)
      dest = j - nv(graph)
      # part of the constrained q matching
      source = original_edge_from_imp[i, dest]

      # TODO: entering & exiting
      rem_edge!(tree, bundle_id_of_vertex(source), representative_vertex_of_leafbranch(source))
      if !(add_edge!(tree, source, dest))
        println("Broken! tried to add already existing edge!")
      end
    end
  end
end

g = erdos_renyi(12, 0.3)

#g = Graph(2048, 4096) 

n_components = size(connected_components(g), 1)

println("Original Graph Components: $n_components")

tree = SimpleGraphFromIterator(kruskal_mst(g))

locs_x, locs_y = spring_layout(g)
labels = nv(g) <= 32 ? (1:nv(g)) : (nothing)
draw(SVG("graph.svg", 16cm, 16cm), gplot(g, locs_x, locs_y, nodelabel=labels))
draw(SVG("tree-0.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))

high_degree = maximum(degree(tree))

#q = 2^(floor(Int, log2(high_degree) - 2)
#h = 1

for i in 1:10 # q >= 0:
  high_degree = maximum(degree(tree))

  low_degree = high_degree - 2 * 1;
  
  improve!(g, tree, high_degree, low_degree)
  print("$high_degree -> ")

  draw(SVG("tree-$i.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))

  n_components = size(connected_components(tree), 1)
  if n_components != 1
    println("FAULT: MDST is no longer a tree!")
    break
  end
end