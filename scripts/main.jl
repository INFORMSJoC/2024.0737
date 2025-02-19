using Pkg
Pkg.activate("CNDP")
# importing the package. 
using Printf
using LinearAlgebra
using Graphs, SparseArrays
using DataFrames, XLSX
using MutableNamedTuples
using Random
using Optim
using ArgParse
using DelimitedFiles
using Distributions

include("../src/variables.jl")
include("../src/misc.jl")
include("../src/utils.jl")
include("../src/optimizer.jl")
include("../src/dataloader.jl")
include("experiments.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dataset", "-d"
            help = "choose dataset"
            arg_type = String
            default = "Anaheim"
        "--jam", "-j"
            help = "set the jam rate"
            arg_type = Int
            default = 1
        "--per", "-p"
            help = "set the percentage of expanding roads"
            arg_type = Float64
            default = 0.1
        "--post", "-t"
            help = "post processing"
            arg_type = Bool
            default = false
    end

    return parse_args(s)
end


function main()

    args = parse_commandline()
    seed = 42

    network_name = args["dataset"]
    per = args["per"]
    jam = args["jam"]
    post = args["post"]

    root = "../results/CNDP/"

    methods = ["DC", "AL", "EDO", "HJ"]
    methods2 = ["DC", "AL", "EDO", "HJ"]

    println("network:", network_name)

    @printf("jam ratio is : %d \n", jam)

    @printf("expansion ratio is : %3.1f \n", per)

    if network_name == "Anaheim"
        dmean, dstd = 2., 4.
        dc_para = MutableNamedTuple(EPS = 1e-4, rho = 1., beta = 1e-8, tau = 2)
        EPS_EDO = 1e-2
        hj_para = MutableNamedTuple(eps = 125, alpha = 0.9, delta = 1000)
        al_para = MutableNamedTuple(mu = 0., rho = 1, beta = 100, alpha1 = 10, alpha2 = 100, EPS_AL = 1e-4, EPS_IN = 1e-2) 
        due_acc = 1e-3
        due_acc2 = 1e-2
        if jam == 1
            due_acc = 1e-5
            due_acc2 = 1e-4
        elseif jam == 3
            due_acc = 1e-4
            due_acc2 = 1e-3
        elseif jam == 5
            due_acc = 1e-3
            due_acc2 = 1e-2
        end  
    elseif network_name == "Barcelona"
        dmean, dstd = 2e6, 4e6
        dc_para = MutableNamedTuple(EPS = 5e-4, rho = 1., beta = 1e-8, tau = 2)
        EPS_EDO = 1e-2
        hj_para = MutableNamedTuple(eps = 0.125, alpha = 0.9, delta = 1.)
        al_para = MutableNamedTuple(mu = 0., rho = 1, beta = 100, alpha1 = 10, alpha2 = 200, EPS_AL = 1e-4, EPS_IN = 1e-2) 
        due_acc = 3e-3
        due_acc2 = 3e-2
        if jam == 1
            due_acc = 1e-4
            due_acc2 = 1e-3
        elseif jam == 3
            due_acc = 1e-2
            due_acc2 = 1e-1
        elseif jam == 5
            due_acc = 1e-2
            due_acc2 = 1e-1
        end
    elseif network_name == "Chicago-Sketch"
        dmean, dstd = 2e2, 4e2
        dc_para = MutableNamedTuple(EPS = 1e-4, rho = 1., beta = 1e-8, tau = 2) 
        EPS_EDO = 1e-2
        hj_para = MutableNamedTuple(eps = 500, alpha = 0.9, delta = 4000)
        al_para = MutableNamedTuple(mu = 0., rho = 1, beta = 100, alpha1 = 10, alpha2 = 100, EPS_AL = 1e-4, EPS_IN = 1e-2) 
        due_acc = 1e-2
        due_acc2 = 1e-1
        if jam == 1
            due_acc = 1e-5
            due_acc2 = 1e-3
        elseif jam == 3
            due_acc = 1e-3
            due_acc2 = 1e-1
        elseif jam == 5
            due_acc = 1e-3
            due_acc2 = 1e-1
        end
    else
        println("experiments on this network not implemented")
        return
    end

    if "DC" in methods
    DC_tntp(root, network_name, 1., per, seed, dmean, dstd, dc_para, jam, due_acc)
    @printf("DC\t")
    end
    
    if "EDO" in methods
    EDO_tntp(root, network_name, 1., per, seed, dmean, dstd, EPS_EDO, jam, due_acc, true)
    @printf("EDO\t")
    end 

    if "HJ" in methods
    HJ_tntp(root, network_name, 1., per, seed, dmean, dstd, hj_para, jam, due_acc) 
    @printf("HJ\t")
    end 
    
    if "AL" in methods
    AL_tntp(root, network_name, 1., per, 42, dmean, dstd, al_para, jam, due_acc) 
    @printf("AL\t")
    end

    if post == true
        @printf("post processing...")
        post_tntp(root, network_name, 42, dmean, dstd, due_acc / 10, 1, jam, methods2)
        @printf("done\n")
    end

end

main()
