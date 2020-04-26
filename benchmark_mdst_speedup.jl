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

function run_test(verts, edges)
    n = 25
    times = []
    for i in 1:n
        # g = make_connected_graph(verts, edges)
        ds_time = @benchmark distributed.mdst(g) setup=(g = make_connected_graph($verts, $edges))
 
        push!(times, time(ds_time))
    end

    print("Vertices,")
    print("Edges,")
    print("Threads,")
    print("Multi-Threaded Time\n")


    print("$verts,")
    print("$edges,")
    print("$(Threads.nthreads()),")
    print("$(sum(times) / length(times))\n")
end

run_test(verts, edges)