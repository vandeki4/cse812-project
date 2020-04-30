# CSE 812 Semester Project
### Distributed Minimum Spanning Trees
Ian Murray, Brenden Vandekieft, Zach McCullough

## Requirements
- [Julia](https://julialang.org/). We have tested on 1.0.1, 1.3.1, and 1.4.3.
- [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)
- [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl)
- [GraphPlot](https://github.com/JuliaGraphs/GraphPlot.jl)
- [Compose](https://github.com/GiovineItalia/Compose.jl)

After installing Julia, the packages can be installed with

```julia
] add BenchmarkTools,LightGraphs,GraphPlot,Compose
```

The tests can then be run as:

```bash
julia benchmark_both.jl N_VERTICES N_EDGES
```

```bash
JULIA_NUM_THREADS=n julia benchmark_mdst_speedup.jl N_VERTICES N_EDGES
```

They each output a single row of CSV output, including the header. 

## File List

- [`distributed.jl`](distributed.jl): Distributed MDST algorithm immplemetation (Brendan, Zac)
- [`singlethread.jl`](singlethread.jl): Single-threaded MDST algorithm implementation, [outlined here](http://viswa.engin.umich.edu/wp-content/uploads/sites/169/2016/12/5.pdf) (Ian)
- [`results.xlsx`](results.xlsx): Various preliminary test results. We do not use these in the final report, they are here for the sake of completeness.
- [`benchmark_mdst_speedup.csv`](benchmark_mdst_speedup.csv): Comparison of singlethreaded and distributed algorithm, run on HPCC
- [`benchmark_both.csv`](benchmark_both.csv): Results of distributed algorithm scaling on MSU HPCC

## Note
Lastly a link to this Github repo is accessible at [https://github.com/vandeki4/cse812-project](https://github.com/vandeki4/cse812-project)
