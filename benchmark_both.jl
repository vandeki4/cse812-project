using LightGraphs
using BenchmarkTools

push!(LOAD_PATH, "./")
using singlethread
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
    st_opt = singlethread.mdst(g)

    end_degrees = []
    n = 25
    times_correct = 0
    for i in 1:n
        ds_opt = distributed.mdst(g)
        if is_connected(ds_opt)
            times_correct += 1
            push!(end_degrees, Δ(ds_opt))
        end
    end
    correct_rate = times_correct / n
    distributed_mean = sum(end_degrees) / length(end_degrees)
    distributed_var = sum((end_degrees .- distributed_mean).^2) / length(end_degrees)

    st_time = @benchmark singlethread.mdst(g)
    ds_time = @benchmark distributed.mdst(g)

    print("Single-Threaded Time,")
    print("Multi-Threaded Time,")
    
    print("Start Degree,")
    print("Single-Threaded End Degree,")
    print("Multi-Threaded End Degree (avg),")
    print("Multi-Threaded End Degree (var),")
    print("Multi-Threaded Correct Rate\n")

    print("$(time(st_time)),")
    print("$(time(ds_time)),")
    
    print("$(Δ(g)),")
    print("$(Δ(st_opt)),")
    print("$distributed_mean,")
    print("$distributed_var,")
    print("$correct_rate\n")
end

g = make_connected_graph(verts, edges)
run_test(g)