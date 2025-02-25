using Pkg
Pkg.activate("CNDP")
# importing the package. 
using Printf
using LinearAlgebra
using Graphs, SparseArrays
using DataFrames, XLSX, CSV
using MutableNamedTuples
using Random
using Optim
using DelimitedFiles

include("../src/variables.jl")
include("../src/misc.jl")
include("../src/utils.jl")
include("../src/optimizer.jl")
include("../src/dataloader.jl")
include("experiments_SiouxFalls.jl")
include("NRMFD_EBA_SiouxFalls.jl")

for jam in 0.8:0.1:1.2
    SiouxFalls_initial(jam=jam)
    HJ_SiouxFalls(jam = jam)
    EDO_SiouxFalls(jam = jam)
    AL_SiouxFalls(jam = jam)
    DC_SiouxFalls(jam = jam)
    NRMFD_EBA_SF(; jam = jam, verbose=0)
    comparison_SiouxFall(jam = jam)
end

for method in ["HJ", "EDO", "AL", "DC", "EBA"]
    post_SiouxFall(; method=method)
end
