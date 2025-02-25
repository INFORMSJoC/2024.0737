# ==============================================================================
# This file contains functions for computing shortest paths in a graph using
# different algorithms: Dijkstra's Algorithm and Bellman-Ford Algorithm,  
# which are useful for finding the feasible solution of CNDP.
# 
# This implementation is primarily derived from the repository:
# https://github.com/chkwon/TrafficAssignment.jl/blob/master/src/misc.jl
#
# The following functions are implemented:
# 1. TA_dijkstra_shortest_paths - Computes shortest paths using Dijkstra's Algorithm.
# 2. TA_bf_shortest_paths - Computes shortest paths using the Bellman-Ford Algorithm.
# 3. create_graph - Creates a directed graph from a set of initial and terminal nodes.
# 4. get_vector - Traces the shortest path from origin to destination and returns a vector.
# 5. add_demand_vector! - Adds demand to the path vector based on a given demand value.
# 6. get_path - Traces the shortest path from origin to destination and returns the sorted path.
#
#
# Author: Haian Yin
# Date: 2022-08-29  # Date when the code was written or last modified
# ==============================================================================

function TA_dijkstra_shortest_paths(graph, travel_time, origin, init_node, term_node, first_thru_node)
    no_node = nv(graph)
    no_arc = ne(graph)

    distmx = Inf*ones(no_node, no_node)
    for i in 1:no_arc
      if term_node[i] >= first_thru_node
          distmx[init_node[i], term_node[i]] = travel_time[i]
      end
    end

    state = dijkstra_shortest_paths(graph, origin, distmx)
    return state
end

function TA_dijkstra_shortest_paths(graph, travel_time, origin, init_node, term_node)
    no_node = nv(graph)
    no_arc = ne(graph)

    distmx = Inf*ones(no_node, no_node)
    for i in 1:no_arc
      distmx[init_node[i], term_node[i]] = travel_time[i]
    end

    state = dijkstra_shortest_paths(graph, origin, distmx)
    return state
end

function TA_bf_shortest_paths(graph, travel_time, origin, init_node, term_node, first_thru_node)
    no_node = nv(graph)
    no_arc = ne(graph)

    distmx = Inf*ones(no_node, no_node)
    for i in 1:no_arc
      if term_node[i] >= first_thru_node
          distmx[init_node[i], term_node[i]] = travel_time[i]
      end
    end

    state = bellman_ford_shortest_paths(graph, origin, distmx)
    return state
end

function TA_bf_shortest_paths(graph, travel_time, origin, init_node, term_node)
    no_node = nv(graph)
    no_arc = ne(graph)

    distmx = Inf*ones(no_node, no_node)
    for i in 1:no_arc
      distmx[init_node[i], term_node[i]] = travel_time[i]
    end

    state = bellman_ford_shortest_paths(graph, origin, distmx)
    return state
end

function create_graph(init_node, term_node)
    @assert Base.length(init_node)==Base.length(term_node)

    no_node = max(maximum(init_node), maximum(term_node))
    no_arc = Base.length(init_node)

    graph = DiGraph(no_node)
    for i=1:no_arc
        add_edge!(graph, init_node[i], term_node[i])
    end
    return graph
end

function get_vector(state, origin, destination, link_dic)
    current = destination
    parent = -1
    x = zeros(Int, maximum(link_dic))

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]

        # println("origin=$origin, destination=$destination, parent=$parent, current=$current")

        if parent != 0
            link_idx = link_dic[parent,current]
            if link_idx != 0
                x[link_idx] = 1
            end
        end

        current = parent
    end

    return x
end

function add_demand_vector!(x, demand, state, origin, destination, link_dic)
  current = destination
  parent = -1

  while parent != origin && origin != destination && current != 0
      parent = state.parents[current]

      if parent != 0
          link_idx = link_dic[parent, current]
          if link_idx != 0
              x[link_idx] += demand
          end
      end

      current = parent
  end
end

function get_path(state, origin, destination, link_dic)
    current = destination
    parent = -1
    path = Any[]

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]

        # println("origin=$origin, destination=$destination, parent=$parent, current=$current")

        if parent != 0
            link_idx = link_dic[parent,current]
            if link_idx != 0
                push!(path, link_idx)
            end
        end

        current = parent
    end

    return sort(path)
end
