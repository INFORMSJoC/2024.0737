# ==============================================================================
# This file defines several mutable structures used for CNDP
#
# The following structures are defined:
#
# 1. OD_sol:
#    - Records the Origin-Destination (OD) solution.
#    - Contains the origin and destination node IDs, number of paths, 
#      and associated paths and flows.
#
# 2. Enhancement_plan:
#    - Represents an enhancement plan for optimizing network performance.
#    - Stores optimization variables y (decision variables), d (coefficients), 
#      and inds (indices of enhanced links).
#
# 3. TA_Data:
#    - Stores key transportation network data and parameters for traffic 
#      assignment modeling.
#    - Includes details about the network (nodes, links, zones), 
#      link attributes (capacity, travel time, speed limit, tolls), 
#      travel demand data, and optimization parameters.
#
# These data structures are essential for running traffic assignment 
# algorithms and optimization procedures in transportation network analysis.
#
# Author: Haian Yin
# Date: 2022-09-06
# ==============================================================================

mutable struct OD_sol
    O::Int
    D::Int 
    npath::Int 
    paths::Array{Array{Int32}}
    flows::Array{Float64}
end

mutable struct Enhancement_plan 
    y::Array{Float64, 1}
    d::Array{Float64, 1}
    inds::Array{Int32, 1}
end

mutable struct TA_Data
    network_name::String

    number_of_zones::Int
    number_of_nodes::Int
    first_thru_node::Int
    number_of_links::Int

    init_node::Array{Int,1}
    term_node::Array{Int,1}
    capacity::Array{Float64,1}
    link_length::Array{Float64,1}
    free_flow_time::Array{Float64,1}
    b::Array{Float64,1}
    power::Array{Float64,1}
    speed_limit::Array{Float64,1}
    toll::Array{Float64,1}
    link_type::Array{Int64,1}

    total_od_flow::Float64

    travel_demand::Array{Float64,2}
    od_pairs::Array{Tuple{Int64,Int64},1}

    toll_factor::Float64
    distance_factor::Float64

    best_objective::Float64
end