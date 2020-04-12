using LightGraphs
using SparseArrays
using GraphPlot
using Compose

push!(LOAD_PATH, "./")
using SparseEdmondsKarp

function approximate_q_matching(bipartite_graph, q::Int, bundle_lookup)
  n::Int = nv(bipartite_graph) / 2
  
  I::Vector{Int} = []
  J::Vector{Int} = []
  V::Vector{Int} = []

  for edge in edges(bipartite_graph)
    push!(I, src(edge))
    push!(J, dst(edge))
    push!(V, 1)
  end
  capacity_matrix = sparse(I, J, V, 2n, 2n)
  flow = convert(SparseMatrixCSC{Float64,Int64}, capacity_matrix)
  fill!(flow.nzval, 1 / ne(bipartite_graph))
  
  capacity = zeros(Int, nv(bipartite_graph))
  for vert in vertices(bipartite_graph)
    if vert > n
      # in Q, right side
      capacity[vert] = max(indegree(bipartite_graph, vert), q)
    else
      # in U, left side
      capacity[vert] = max(1, outdegree(bipartite_graph, vert))
    end
  end
  bundle_capacity = zeros(Int, n)
  for v in 1:n
    if bundle_lookup(v) != 0
      bundle_capacity[bundle_lookup(v)] += 1
    end
  end
  for v in 1:n
    if bundle_lookup(v) != 0
      bundle_capacity[bundle_lookup(v)] = max(bundle_capacity[bundle_lookup(v)], q)
    end
  end

  doubled = true
  while doubled
    flow_through_vert = zeros(Float64, nv(bipartite_graph))
    for e in edges(bipartite_graph)
      flow_through_vert[src(e)] += flow[src(e), dst(e)]
      flow_through_vert[dst(e)] += flow[src(e), dst(e)]
    end

    # Whether U nodes have been requested to double
    request_double::Vector{Tuple{Int, Int}} = []
    
    for e in edges(bipartite_graph)
      # if not full
      if flow_through_vert[dst(e)] < (1/8) * capacity[dst(e)]
        push!(request_double, (src(e), dst(e)))
      end
    end

    # flow through a bundle. uses representative parent id, so need n
    bundle_flow = zeros(Float64, n)
    # U nodes
    for v in 1:n
      if bundle_lookup(v) != 0
        bundle_flow[bundle_lookup(v)] += flow_through_vert[v]
      end
    end
    doubled = false
    for edge in request_double
      if bundle_lookup(edge[1]) == 0
        continue
      end
      if bundle_flow[bundle_lookup(edge[1])] < (1/8) * bundle_capacity[bundle_lookup(edge[1])]
        doubled = true
        flow[edge[1], edge[2]] *= 2
      end
    end
  end

  S = Set(Iterators.filter(
    e -> flow[src(e), dst(e)] > rand(),
    edges(bipartite_graph)
  ))

  # remove extras...
  # do not need n, can get away with fewer since there are fewer branches
  bundle_incidence = zeros(Int, n)
  uq_incidence = zeros(Int, 2n)

  result::Vector{Edge} = []
  
  for edge in S
    if uq_incidence[src(edge)] < 1 && 
      uq_incidence[dst(edge)] < q && 
      bundle_incidence[bundle_lookup(src(edge))] < q
      push!(result, edge)
    end

    uq_incidence[src(edge)] += 1
    uq_incidence[dst(edge)] += 1
    bundle_incidence[bundle_lookup(src(edge))] += 1
  end

  println(result)

  return result
end

function improve!(graph, tree, q, high_degree, low_degree)
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
  matches = approximate_q_matching(imp_graph, q, bundle_id_of_vertex)
  for edge in matches
    dest = dst(edge) - nv(graph)
    # part of the constrained q matching
    source = original_edge_from_imp[src(edge), dest]

    # TODO: entering & exiting
    rem_edge!(tree, bundle_id_of_vertex(source), representative_vertex_of_leafbranch(source))
    if !(add_edge!(tree, source, dest))
      println("Broken! tried to add already existing edge!")
    end
  end
end

function prog!(tree, q)
  k = Δ(tree)
  b(j) = k + 1 - j * q

  low_degree = k - 2 * 1;
  
  improve!(g, tree, q, k, low_degree)
end

function mdst(graph, save_graphs=false)
  tree = SimpleGraphFromIterator(prim_mst(g))

  q = 2^max(0, floor(Int, log2(Δ(tree)) - 2))
  h = 1
  i = 0

  while q >= 1
    i += 1

    prog!(tree, q)

    if save_graphs
      draw(SVG("tree-$i.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))
    end

    n_components = size(connected_components(tree), 1)
    if n_components != 1
      println("FAULT: MDST is no longer a tree!")
      # break
    end

    q = Int(floor(q/2))
  end

  return tree
end

g = erdos_renyi(512, 0.1)

#g = Graph(2048, 4096) 

n_components = size(connected_components(g), 1)
println("Original Graph Components: $n_components")

tree = mdst(g)

locs_x, locs_y = spring_layout(g)
labels = nv(g) <= 32 ? (1:nv(g)) : (nothing)
draw(SVG("graph.svg", 16cm, 16cm), gplot(g, locs_x, locs_y, nodelabel=labels))
draw(SVG("tree-final.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))



