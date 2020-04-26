using LightGraphs
using BenchmarkTools

push!(LOAD_PATH, "./")
using distributed

verts = parse(Int, ARGS[1])
edges = parse(Int, ARGS[2])


function make_connected_graph(verts, edges, max_tries = 100)
    g = erdos_renyi(verts, edges)
    tries = 0
    while !is_connected(g) && tries < max_tries
        tries += 1
        g = erdos_renyi(verts, edges)
    end
    if tries == max_tries
        throw("Failed to make a random connected graph ($verts, $edges)!")
    end
    return g
end

function run_test(g)
    ds_time = @benchmark distributed.mdst(g)

    print("Vertices,")
    print("Edges,")
    print("Threads,")
    print("Multi-Threaded Time\n")


    print("$(nv(g)),")
    print("$(ne(g)),")
    print("$(Threads.nthreads()),")
    print("$(time(ds_time))\n")
end

g = make_connected_graph(verts, edges)
run_test(g)