module distributed

using LightGraphs
using SparseArrays
using GraphPlot
using Compose
using BenchmarkTools

push!(LOAD_PATH, "./")
using SparseEdmondsKarp

function approximate_q_matching(bipartite_graph, q::Int, bundle_lookup)
  n::Int = nv(bipartite_graph) / 2
  
  I::Array{Int, 1} = zeros(ne(bipartite_graph))
  J::Array{Int, 1} = zeros(ne(bipartite_graph))
  V::Array{Int, 1} = zeros(ne(bipartite_graph))

  for (i ,edge) in enumerate(edges(bipartite_graph))
    I[i] = src(edge)
    J[i] = dst(edge)
    V[i] = 1
  end
  
  capacity_matrix = sparse(I, J, V, 2n, 2n)
  flow = convert(SparseMatrixCSC{Float64,Int64}, capacity_matrix)
  fill!(flow.nzval, 1 / ne(bipartite_graph))
  
  capacity = zeros(Int, nv(bipartite_graph))
  Threads.@threads for vert in vertices(bipartite_graph)
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

    # cannot be done parallel without syncing

    for j in 1:size(flow, 2)
      for i_ind in flow.colptr[j]:(flow.colptr[j+1]-1)
        i = flow.rowval[i_ind]
        v = flow.nzval[i_ind]
        flow_through_vert[i] += v 
        flow_through_vert[j] += v
      end
    end

    #for (i,j,v) in zip(findnz(flow)...)
    #  flow_through_vert[i] += v 
    #  flow_through_vert[j] += v 
    #end
    #for e in edges(bipartite_graph)
    #  flow_through_vert[src(e)] += flow[src(e), dst(e)]
    #  flow_through_vert[dst(e)] += flow[src(e), dst(e)]
    #end

    # Whether U nodes have been requested to double
    request_double = Array{Tuple{Int, Int}, 1}(undef, length(flow.nzval))
    request_count = 0

    for j in 1:size(flow, 2)
      for i_ind in flow.colptr[j]:(flow.colptr[j+1]-1)
        i = flow.rowval[i_ind]
        v = flow.nzval[i_ind]
        if flow_through_vert[j] < (1/8) * capacity[j]
          request_count += 1
          request_double[request_count] = (i, j)
        end
      end
    end

    # for (i,j,v) in zip(findnz(flow)...)
    #   if flow_through_vert[j] < (1/8) * capacity[j]
    #     request_count += 1
    #     request_double[request_count] = (i, j)
    #   end
    # end

    # flow through a bundle. uses representative parent id, so need n
    bundle_flow = zeros(Float64, n)
    # U nodes
    for v in 1:n
      if bundle_lookup(v) != 0
        bundle_flow[bundle_lookup(v)] += flow_through_vert[v]
      end
    end

    doubled = false
    Threads.@threads for i in 1:request_count
      edge = request_double[i]
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
    if bundle_lookup(src(edge)) == 0
      continue
    end
    if uq_incidence[src(edge)] < 1 && 
      uq_incidence[dst(edge)] < q && 
      bundle_incidence[bundle_lookup(src(edge))] < q
      push!(result, edge)
    end

    uq_incidence[src(edge)] += 1
    uq_incidence[dst(edge)] += 1
    bundle_incidence[bundle_lookup(src(edge))] += 1
  end

  return result
end

function improve!(graph, tree, q, γ)
  high_degree = γ
  low_degree = γ - 2q
  
  retain_vertex = degree(tree) .< high_degree
  vertices_to_retain = filter(n -> n > 0, retain_vertex .* eachindex(retain_vertex))

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
  Threads.@threads for i in eachindex(branch_components)
    comp = branch_components[i]
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
    is_good = !is_high_edge(branches_graph, edge, low_degree) # TODO: is this leaf_branches? Might be tree

    # in different branches
    is_good &= branch_lookup[src(edge)] != branch_lookup[dst(edge)]
    # with at least one being a leaf branch
    is_good &= is_vert_of_leafbranch(src(edge)) || is_vert_of_leafbranch(dst(edge))


    # TODO: add_edge! is actually fairly slow.
    # Ideally this would be parallel.
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

  matches = approximate_q_matching(imp_graph, q, bundle_id_of_vertex)

  if size(matches, 1) == 0
    return 0
  end

  filtered_match::Vector{Edge} = []
  while size(filtered_match, 1) < (size(matches, 1) / 8)
    filtered_match = []
    leaf_coins = rand(size(branch_components, 1)) .> 0.5

    for match in matches
      if leaf_coins[branch_lookup[src(match)]] == false ||
        leaf_coins[branch_lookup[dst(match) - nv(graph)]] == true
        push!(filtered_match, match)
      end
    end
  end

  for edge in filtered_match
    dest = dst(edge) - nv(graph)
    # part of the constrained q matching
    source = original_edge_from_imp[src(edge), dest]

    @assert is_vert_of_leafbranch(source)

    if !has_edge(tree, source, dest) && has_edge(tree, bundle_id_of_vertex(source), representative_vertex_of_leafbranch(source))
      add_edge!(tree, source, dest)
      rem_edge!(tree, bundle_id_of_vertex(source), representative_vertex_of_leafbranch(source))
    end
  end

  return size(matches, 1)
end

function prog!(graph, tree, q, h_hat, i, save_graphs, func)
  k = Δ(tree)

  δ = 1 - 1/(log(nv(tree))^0.25)
  c = 1.0001
  τ = 2 / (1 - δ)
  b(j, k, q) = k + 1 - j * q

  t = ceil((k + 1 - h_hat) / q) + 1
  if t <= 0
    t = 1
  end
  
  j = -1
  first = true

  staleness = 0
  last_val = 0

  while first || b(j, k, q) > h_hat
    first = false

    j = argmax([sum(degree(tree) .>= b(s, k, q)) / (τ^s) for s in 1:t])
    if j > floor(k/q - 2)
      j = floor(k/q - 2)
    end
    
    k = Δ(tree)
    # print("$k -> ")
    imp_size = improve!(graph, tree, q, b(j, k, q))

    if (k >= last_val)
      last_val = k
      staleness += 1
    else
      i += 1
      if save_graphs
        func(i, tree)
      end
      staleness = 0
      last_val = Δ(tree)
    end

    Π = δ * q * sum(degree(tree) .>= b(j, k, q)) / (8c)
    # Π += staleness / 30;
    
    if (imp_size <= Π)
      h_hat = max(h_hat, b(j, k, q))
    end

  end
  return (h_hat, i)
end

function mdst(graph, save_graphs=false, func=nothing)
  tree = SimpleGraphFromIterator(prim_mst(graph))

  if save_graphs
    func(0, tree)
  end

  q = 2^max(0, floor(Int, log2(Δ(tree)) - 2))
  h_hat = 1
  i = 0

  while q >= 1
    h_hat, i = prog!(graph, tree, q, h_hat, i, save_graphs, func)

    q = Int(floor(q/2))

    n_components = size(connected_components(tree), 1)
    if n_components != 1
      # println("Fault: n_components = $n_components")
      # break
    end
    if ne(tree) != nv(tree) - 1
      println("Fault: $(ne(tree)) != $(nv(tree) - 1)")
    end
  end

  return tree
end


function main()
  g = erdos_renyi(16, 0.5)

  #g = Graph(2048, 4096)

  n_components = size(connected_components(g), 1)
  println("Original Graph Components: $n_components")

  locs_x = zeros(nv(g))
  locs_y = zeros(nv(g))
  labels = nothing
  if nv(g) <= 24
    locs_x, locs_y = circular_layout(g)
    labels = 1:nv(g)
  else
    locs_x, locs_y = spring_layout(g)
  end
  draw(SVG("graph.svg", 16cm, 16cm), gplot(g, locs_x, locs_y, nodelabel=labels))

  function save(i, tree)
    draw(SVG("tree-$i.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))
  end

  tree = mdst(g, false, save)
  println(Δ(tree))

  draw(SVG("tree-final.svg", 16cm, 16cm), gplot(tree, locs_x, locs_y, nodelabel=labels))
end

function benchmark_test()
  t64 = @benchmark mdst(erdos_renyi(64, 1000))
  t128 = @benchmark mdst(erdos_renyi(128, 4000))
  t256 = @benchmark mdst(erdos_renyi(256, 16000))
  t512 = @benchmark mdst(erdos_renyi(512, 64000))
  t2048 = @benchmark mdst(erdos_renyi(512, 128000))

  return (t64, t128, t256, t512, t2048)
end

# results = benchmark()
# println("Threads: $(Threads.nthreads())")
# t64 = @benchmark mdst(erdos_renyi(512, 128000))
# println(time(t64))

# using ProfileView

# mdst(erdos_renyi(64, 1000))
# @profview mdst(erdos_renyi(512, 128000))

end