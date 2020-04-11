using LightGraphs
using Plots

function improve(graph, tree, high_degree, low_degree)
  retain_vertex = degree(tree) .< high_degree
  vertices_to_retain = filter(n -> n > 0, retain_vertex .* eachindex(retain_vertex))

  # Remove edges adjacent to high degree vertices
  is_high_edge(graph, edge, deg) = degree(graph, src(edge)) >= deg || degree(graph, dst(edge)) >= deg

  branches_graph = SimpleGraphFromIterator(
    Iterators.filter(
      edge -> !is_high_edge(tree, edge, high_degree),
      edges(tree)
    )
  )

  # Fast vertex -> branch lookup
  branch_components = connected_components(branches_graph)
  branch_lookup = collect(1:size(graph, 1))
  for (index, comp) in enumerate(branch_components)
    for vert in comp
      branch_lookup[vert] = index
    end
  end

  # Determine which branches are leaf branches
  n_parents_for_branch::Array{Int, 1} = zeros(size(branch_components, 1))
  parent_for_branch::Array{Int, 1} = zeros(size(branch_components, 1))
  for edge in edges(tree)
    if has_edge(branches_graph, edge)
      continue
    end
    n_parents_for_branch[branch_lookup[src(edge)]] += 1
    n_parents_for_branch[branch_lookup[dst(edge)]] += 1

    # parent is only meaningful for leaf branches,
    # but we don't know what ones are leaf branches yet
    parent_for_branch[branch_lookup[src(edge)]] = dst(edge)
    parent_for_branch[branch_lookup[dst(edge)]] = src(edge)
  end
  is_leafbranch(branch_id) = n_parents_for_branch[branch_id] == 1
  is_vert_of_leafbranch(vertex) = is_leafbranch(branch_lookup[vertex])

  g_sub_tree = Iterators.filter(
    edge -> !has_edge(tree, edge),
    edges(graph)
  )

  good = 0
  bad = 0
  for edge in g_sub_tree
    is_good = !is_high_edge(branches_graph, edge, low_degree) # TODO: is this leaf_branches? Might be tree

    # in different branches
    is_good &= branch_lookup[src(edge)] != branch_lookup[dst(edge)]
    # with at least one being a leaf branch
    is_good &= is_vert_of_leafbranch(src(edge)) || is_vert_of_leafbranch(dst(edge))
    if (is_good)
      good += 1
    else
      bad += 1
    end
  end

  for (id, branch) in enumerate(branch_components)
    if !is_leafbranch(id)
      continue
    end
    println(parent_for_branch[id])
  end

  println("Leaf branch count: $(size(branch_components, 1))")
end

g = erdos_renyi(512, 0.1)
n_components = size(connected_components(g), 1)

println("Original Graph Components: $n_components")

tree = SimpleGraphFromIterator(kruskal_mst(g))

high_degree = maximum(degree(tree))
low_degree = high_degree - 2 * 1;

improve(g, tree, high_degree, low_degree)