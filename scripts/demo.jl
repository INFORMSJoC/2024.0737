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
using JuMP, Ipopt
using Distributions

include("../src/variables.jl")
include("../src/misc.jl")
include("../src/utils.jl")
include("../src/optimizer.jl")
include("../src/dataloader.jl")
include("experiments_SiouxFalls.jl")

DC_SiouxFalls(jam = 1)
