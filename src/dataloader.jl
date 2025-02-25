# ==============================================================================
# Traffic Assignment Network Data Loading and Enhancement
#
# This Julia script is used for loading and processing transportation network data, specifically for solving traffic assignment and capacity expansion problems.
#
# Variables:
#   - SiouxFalls_dir: Directory path for the SiouxFalls dataset.
#   - tntp_dir: Directory path for the generic transportation network datasets.
#
# Functions:
#   - read_ta: Reads network and trip data, and returns them in a structured format.
#   - load_SiouxFalls: Loads a specific dataset (SiouxFalls) for the capacity expansion problem.
#   - load_SiouxFalls2: Loads a specific dataset (SiouxFalls) from the TNTP dataset.
#   - load_tntp: Loads a generic network dataset for traffic assignment with the possibility of expansion planning.
#   - load_ta_network: Loads and processes network and trip data from a specified source.
#
# Author: Haian Yin
# Date: 2024-08-26  # Date when the code was written or last modified
# ==============================================================================

SiouxFalls_dir = "../data/SiouxFalls/SiouxFalls"   
tntp_dir = "../data/TransportationNetworks"

# define functions for loading data
function read_ta(netname)
    m_net = XLSX.readxlsx(netname * "_net" * ".xlsx")[1][:][2:end, :]
    m_trips = XLSX.readxlsx(netname * "_trips" * ".xlsx")[1][:]

    ta_data = TA_Data(
        netname, # network_name::String
        0, # number_of_zones::Int
        0, # number_of_nodes::Int
        0, # first_thru_node::Int
        size(m_net)[1], #  number_of_links::Int

        m_net[:, 2], # init_node::Array{Int,1}
        m_net[:, 3], # term_node::Array{Int,1}
        m_net[:, 6], # capacity::Array{Float64,1}
        [0], # link_length::Array{Float64,1}
        m_net[:, 4], # free_flow_time::Array{Float64,1}
        m_net[:, 5], # b::Array{Float64,1}
        m_net[:, 7], # power::Array{Float64,1}
        [0], # speed_limit::Array{Float64,1}
        [0], # toll::Array{Float64,1}
        [0], # link_type::Array{Int64,1}
        
        sum(m_trips), # total_od_flow::Float64
        m_trips, # travel_demand::Array{Float64,2}
        [(0, 0)], # od_pairs::Array{Tuple{Int64,Int64},1}

        0, # toll_factor::Float64
        0, # distance_factor::Float64

        -1 # best_objective::Float64
        )

    return ta_data
end


function load_SiouxFalls(jam = 1)
    # this funciton is used for load the dataset of Sioux Falls capacity expansion problem 
    # td = read_ta("SiouxFalls")
    td = read_ta(SiouxFalls_dir)
    td.travel_demand *= jam

    inds = [16, 17, 19, 20, 25, 26, 29, 39, 48, 74]
    d =  0.001 * [26., 40., 26., 40., 25., 25., 48., 34., 48., 34.]
    y = zeros(lastindex(inds))
    ep = Enhancement_plan(y, d, inds)
    flow = zeros(td.number_of_links)
    return td, flow, ep 
end

include("utils.jl")

function load_tntp(name, per, seed, dL, dU, jam = 1)
    # name:: name of dataset 
    # per:: percent of expansion roads, 0<per<=1
    # seed:: seed of random generator 
    td = load_ta_network(name) 
    td.travel_demand *= jam
    
    rng = MersenneTwister(seed)
    inds = randsubseq(rng, 1:td.number_of_links, per)
    tmp1 = sum(td.capacity[inds] .* td.capacity[inds])

    y = zeros(lastindex(inds))
    d = zeros(lastindex(inds))
    ep = Enhancement_plan(y, d, inds)

    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:length(td.init_node)))

    qls_opt = (lambda = 1., delta = 0.5, gamma = 1e-8)
    bi_para = MutableNamedTuple(rho = 1., beta = 1e-8, tau = 100)
    # m_od_f = Array{Any}(nothing, size(td.travel_demand)[1], size(td.travel_demand)[2])
    m_od_f = Array{Any}(nothing, size(td.travel_demand)[1])
    for i in 1:size(td.travel_demand)[1]
        m_od_f[i] = Array{Any}(nothing, 0)
    end
    flow = zeros(td.number_of_links)

    all_or_nothing!(m_od_f, flow, td, graph, link_dic);
    # for i in 1:10
    nii = 0
    # while rgap(flow, td, graph, link_dic, ep) > due_acc
    for _ in 1:10
        nii += 1
        flow = ISMO_f!(m_od_f, flow, td, ep, graph, link_dic, qls_opt)
        if nii > 10; break; end
        # break
    end

    tmp = 0.0
    for i in 1:td.number_of_links
        if flow[i] > 0
            tmp += td.free_flow_time[i] * flow[i] * (1 + td.b[i] .* (flow[i]/(td.capacity[i])) ^ (td.power[i]))
        end
    end

    # dmean = 10*tmp1 / tmp 
    # dstd  = dmean / 2 
    # d = 1e-3 .* (dmean .+ dstd .* randn(rng, length(inds)))
    # d = max.(d, 0)
    d = 1e-3 .* (dL .+ (dU - dL) .* rand(rng, length(inds)))

    y = zeros(lastindex(inds))
    ep = Enhancement_plan(y, d, inds)
    flow = zeros(td.number_of_links)

    return td, flow, ep 
end

function load_SiouxFalls2(jam=1)
    # this funciton is used for load the dataset of Sioux Falls capacity expansion problem 
    td = load_ta_network("SiouxFalls")
    td.capacity *= 0.001
    td.travel_demand *= 0.0011
    td.travel_demand *= jam
    td.free_flow_time *= 0.01
    inds = [16, 17, 19, 20, 25, 26, 29, 39, 48, 74]
    d =  0.001*[26., 40., 26., 40., 25., 25., 48., 34., 48., 34.]
    y = zeros(lastindex(inds))
    ep = Enhancement_plan(y, d, inds)
    flow = zeros(td.number_of_links)
    return td, flow, ep 
end

function read_ta_network(network_name)
  network_dir = joinpath(tntp_dir, network_name)
  println(pwd())
  println(tntp_dir)
  @assert ispath(network_dir)

  network_data_file = ""
  trip_table_file = ""

  for f in readdir(network_dir)
    if occursin(".zip", lowercase(f))
      zipfile = joinpath(network_dir, f)
      run(unpack_cmd(zipfile, network_dir, ".zip", ""))
      rm(zipfile)
    end
  end

  for f in readdir(network_dir)
    if occursin("_net", lowercase(f)) && occursin(".tntp", lowercase(f))
      network_data_file = joinpath(network_dir, f)
    elseif occursin("_trips", lowercase(f)) && occursin(".tntp", lowercase(f))
      trip_table_file = joinpath(network_dir, f)
    end
  end

  @assert network_data_file != ""
  @assert trip_table_file != ""

  return network_data_file, trip_table_file
end


# Traffic Assignment Data structure
using BinDeps, DataFrames, OrderedCollections
using Distributed, Printf, LinearAlgebra, SparseArrays

search_sc(s, c) = something(findfirst(isequal(c), s), 0)

function load_ta_network(network_name; best_objective=-1.0, toll_factor=0.0, distance_factor=0.0)
  network_data_file, trip_table_file = read_ta_network(network_name)

  load_ta_network(network_name, network_data_file, trip_table_file, best_objective=best_objective, toll_factor=toll_factor, distance_factor=distance_factor)
end


function load_ta_network(network_name, network_data_file, trip_table_file; best_objective=-1.0, toll_factor=0.0, distance_factor=0.0)

  @assert ispath(network_data_file)
  @assert ispath(trip_table_file)

  ##################################################
  # Network Data
  ##################################################


  number_of_zones = 0
  number_of_links = 0
  number_of_nodes = 0
  first_thru_node = 0

  n = open(network_data_file, "r")

  while (line = readline(n)) != ""
    if occursin("<NUMBER OF ZONES>", line)
      number_of_zones = parse(Int, line[search_sc(line, '>')+1:end])
    elseif occursin("<NUMBER OF NODES>", line)
      number_of_nodes = parse(Int, line[search_sc(line, '>')+1:end])
    elseif occursin("<FIRST THRU NODE>", line)
      first_thru_node = parse(Int, line[search_sc(line, '>')+1:end])
    elseif occursin("<NUMBER OF LINKS>", line)
      number_of_links = parse(Int, line[search_sc(line, '>')+1:end])
    elseif occursin("<END OF METADATA>", line)
      break
    end
  end

  @assert number_of_links > 0

  init_node = Array{Int64}(undef, number_of_links)
  term_node = Array{Int64}(undef, number_of_links)
  capacity = zeros(number_of_links)
  link_length = zeros(number_of_links)
  free_flow_time = zeros(number_of_links)
  b = zeros(number_of_links)
  power = zeros(number_of_links)
  speed_limit = zeros(number_of_links)
  toll = zeros(number_of_links)
  link_type = Array{Int64}(undef, number_of_links)

  idx = 1
  while !eof(n)
    line = readline(n)
    if occursin("~", line) || line == ""
      continue
    end

    if occursin(";", line)
      line = strip(line, [' ', '\n', ';'])
      line = replace(line, ";" => "")

      numbers = split(line)
      init_node[idx] = parse(Int64, numbers[1])
      term_node[idx] = parse(Int64, numbers[2])
      capacity[idx] = parse(Float64, numbers[3])
      link_length[idx] = parse(Float64, numbers[4])
      free_flow_time[idx] = parse(Float64, numbers[5])
      b[idx] = parse(Float64, numbers[6])
      power[idx] = parse(Float64, numbers[7])
      speed_limit[idx] = parse(Float64, numbers[8])
      toll[idx] = parse(Float64, numbers[9])
      link_type[idx] = Int(round(parse(Float64, numbers[10])))

      idx = idx + 1
    end
  end

  ##################################################
  # Trip Table
  ##################################################

  number_of_zones_trip = 0
  total_od_flow = 0

  f = open(trip_table_file, "r")

  while (line = readline(f)) != ""
    if occursin("<NUMBER OF ZONES>", line)
      number_of_zones_trip = parse(Int, line[search_sc(line, '>')+1:end])
    elseif occursin("<TOTAL OD FLOW>", line)
      total_od_flow = parse(Float64, line[search_sc(line, '>')+1:end])
    elseif occursin("<END OF METADATA>", line)
      break
    end
  end

  @assert number_of_zones_trip == number_of_zones # Check if number_of_zone is same in both txt files
  @assert total_od_flow > 0

  travel_demand = zeros(number_of_zones, number_of_zones)
  od_pairs = Array{Tuple{Int64,Int64}}(undef, 0)

  origin = -1

  while !eof(f)
    line = readline(f)

    if line == ""
      origin = -1
      continue
    elseif occursin("Origin", line)
      origin = parse(Int, split(line)[2])
    elseif occursin(";", line)
      pairs = split(line, ";")
      for i = 1:size(pairs)[1]
        if occursin(":", pairs[i])
          pair = split(pairs[i], ":")
          destination = parse(Int64, strip(pair[1]))
          od_flow = parse(Float64, strip(pair[2]))

          # println("origin=$origin, destination=$destination, flow=$od_flow")

          travel_demand[origin, destination] = od_flow
          push!(od_pairs, (origin, destination))
        end
      end
    end
  end

  # Preparing data to return
  ta_data = TA_Data(
    network_name,
    number_of_zones,
    number_of_nodes,
    first_thru_node,
    number_of_links,
    init_node,
    term_node,
    capacity,
    link_length,
    free_flow_time,
    b,
    power,
    speed_limit,
    toll,
    link_type,
    total_od_flow,
    travel_demand,
    od_pairs,
    toll_factor,
    distance_factor,
    best_objective)

  return ta_data

end # end of load_network function