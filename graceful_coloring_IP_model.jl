"""
Title: Graceful Graph Coloring - Integer Programming Model
Description: Julia implementation of an Integer Programming (MIP) model 
             to solve the graceful graph coloring problem.
Author: Atilio Gomes Luiz
Date: August 23, 2025
Dependencies:
    - JuMP (optimization modeling)
    - CPLEX (integer programming solver)
    - Graphs.jl (graph manipulation)
    - DelimitedFiles (file reading/writing)
    - Dates (date handling)
    - MathOptInterface (optimization interface)
Usage: Define the graph as a SimpleGraph object and use the MIP model
       to find the corresponding graceful coloring.
"""
# Importing the necessary packages
using DelimitedFiles
using JuMP, CPLEX
using Graphs
using Dates
using MathOptInterface
const MOI = MathOptInterface

# Defining a constant that controls debug messages
const DEBUG = false  # false to deactivate

# The @debug macro for wrapping any code block
macro debug(block)
    return :(if DEBUG
                $(esc(block))
            end)
end

# Function to read a graph from a file
function read_simple_graph(filepath::String)::SimpleGraph
    """
    Reads an edge list from a file, representing an undirected simple graph.
    Normalizes the vertices to be consecutive 0-based indices
    (internally converted to 1-based for `SimpleGraph`).
    Discards self-loops and multiple edges.
    """
    # read the edges from the file
    edges = open(filepath, "r") do file
        [parse.(Int, split(line)) for line in eachline(file)]
    end

    # Remove loops and multiple edges 
    unique_edges = Set{Tuple{Int,Int}}()
    for (u, v) in edges
        if u != v
            push!(unique_edges, u < v ? (u, v) : (v, u))
        end
    end

    # Create a sorted list of all the present vertices
    all_vertices = sort(collect(unique(vcat(collect.(unique_edges)...))))

    # Creates a normalization map: original -> consecutive index (1-based)
    vertex_map = Dict(v => i for (i, v) in enumerate(all_vertices))

    # Number of normalized vertices
    n = length(all_vertices)

    # Creates a simple graph (1-based)
    g = SimpleGraph(n)

    # Add the normalized edges 
    for (u, v) in unique_edges
        u_norm = vertex_map[u] 
        v_norm = vertex_map[v]
        add_edge!(g, u_norm, v_norm)
    end

    return g
end

# Finds a feasible solution for the input graph,
# that is, finds a valid graceful coloring for the graph
function viable_solution(graph::SimpleGraph)::Dict{Int,Int}
    colors = Dict{Int,Int}()

    adj = Dict()
    for v in vertices(graph)
        adj[v] = neighbors(graph, v)
    end

    temp_graph = deepcopy(adj)
    current_color = 1

    while length(temp_graph) > 0
        two_packing = Vector{Int}()
        auxg = deepcopy(temp_graph)
        for w in keys(temp_graph)
            if w in keys(auxg)
                push!(two_packing, w)
                for z in adj[w]
                    for x in adj[z]
                        if x != w && x in keys(auxg)
                            # delete vertex x and its incident edges
                            for u in auxg[x]
                                idx = findfirst(==(x), auxg[u])
                                if idx !== nothing
                                    auxg[u][idx] = auxg[u][end] 
                                    pop!(auxg[u])               
                                end
                            end
                            delete!(auxg, x)
                        end
                    end
                    if z in keys(auxg)
                        # delete vertex z and its incident edges
                        for u in auxg[z]
                            idx = findfirst(==(z), auxg[u])
                            if idx !== nothing
                                auxg[u][idx] = auxg[u][end]   
                                pop!(auxg[u])                 
                            end
                        end
                        delete!(auxg, z)
                    end
                end
                if w in keys(auxg)
                    # delete vertex w and its incident edges
                    for u in auxg[w]
                        idx = findfirst(==(w), auxg[u])
                        if idx !== nothing
                            auxg[u][idx] = auxg[u][end]  
                            pop!(auxg[u])                
                        end
                    end
                    delete!(auxg, w)
                end
            end
        end
        vertices_to_be_colored = Set{Int}()
        for w in two_packing
            safe_vertex = true
            for z in adj[w] 
                for y in adj[z]
                    if haskey(colors, z) && haskey(colors, y) && w != y && current_color == 2*colors[z]-colors[y]
                        safe_vertex = false
                    end
                end
            end
            if safe_vertex == true 
                push!(vertices_to_be_colored, w)
            end
        end
        for v in vertices_to_be_colored
            colors[v] = current_color
            # delete vertex v and its incident edges
            for u in keys(temp_graph)
                idx = findfirst(==(v), temp_graph[u])
                if idx !== nothing
                    temp_graph[u][idx] = temp_graph[u][end]  
                    pop!(temp_graph[u])                      
                end
            end
            delete!(temp_graph, v)
        end
        current_color += 1
    end
    return colors
end


# Function that computes a graceful coloring of the input graph
function graceful_coloring(graph::SimpleGraph, timelimitMinutes::Int, viable_solution::Dict{Int,Int}=Dict())
    """
    This function takes as input the simple graph, the maximum time (in minutes)
    allowed for the solver to run on this graph, and a viable solution. If the solver does not 
    find the optimal solution within this time, the best solution found will be returned.
    """
    # Preparing the optimization model
    model = Model(CPLEX.Optimizer)

    set_optimizer_attribute(model, MOI.Silent(), true) # activate log if false

    # Setting the maximum time (in seconds) for CPLEX running time
    set_optimizer_attribute(model, MOI.TimeLimitSec(), timelimitMinutes*60)

    max_degree = maximum(degree(graph, v) for v in vertices(graph))
    Big_M = 2*(max_degree^2) - max_degree + 1

    # Defining decision variables
    @variable(model, 1 <= x[v in vertices(graph)], Int)
    @variable(model, 1 <= z <= Big_M, Int)

    pairs1 = [(i, j) for i in vertices(graph) for j in neighbors(graph, i) if i < j]
    @variable(model, b[pairs1], Bin)

    # constraint: all vertices must be assigned a label at most z (x_i <= z)
    for i in vertices(graph)
        @constraint(model, x[i] <= z)
    end

    # constraint: any two adjacent vertices must have distinct colors (x_i != x_j)
    for i in vertices(graph)
        for j in neighbors(graph, i)
            if i < j 
                @constraint(model, x[i]-x[j]>=1-Big_M*(1-b[(i,j)]))
                @constraint(model, x[j]-x[i]>=1-Big_M*b[(i,j)])
            end
        end
    end

    # constraint: any two vertices at distance 2 must have distinct colors (x_i != x_j)
    pairs2 = [(i, j) for v in vertices(graph) for i in neighbors(graph, v) for j in neighbors(graph, v) if i < j && !(i in neighbors(graph, j))]
    pairs2 = unique(pairs2)
    @variable(model, c[pairs2], Bin)
    for v in vertices(graph)
        for i in neighbors(graph, v)
            for j in neighbors(graph, v)
                if i < j && !(i in neighbors(graph, j))
                    @constraint(model, x[i]-x[j]>=1-Big_M*(1-c[(i,j)]))
                    @constraint(model, x[j]-x[i]>=1-Big_M*c[(i,j)])
                end
            end
        end
    end

    # constraint: adjacent edges must have distinct colors
    triples = [(i, j, k) for j in vertices(graph) for i in neighbors(graph, j) for k in neighbors(graph, j) if i < k]
    triples = unique(triples)
    @variable(model, d[triples], Bin)
    MAX = 4*(max_degree^2) - 2*max_degree
    for j in vertices(graph)
        for i in neighbors(graph, j)
            for k in neighbors(graph, j)
                if i < k && !((i,j,k) in d)
                    @constraint(model, (x[i]+x[k]-2*x[j])>=1-MAX*(1-d[(i,j,k)]))
                    @constraint(model, (2*x[j]-(x[i]+x[k]))>=1-MAX*d[(i,j,k)])
                end
            end
        end
    end

    # Setting the objective function
    @objective(model, Min, z)

    # Printing the model
    @debug println(model) 

    # Setting an initial viable solution
    for (i, val) in viable_solution
        set_start_value(x[i], val)
    end
    set_start_value(z, maximum(values(viable_solution)))

    start_time = now()
    JuMP.optimize!(model)
    end_time = now()

    @debug begin
        println(stderr, "Termination status: ", termination_status(model))
        println(stderr, "Primal status: ", primal_status(model))
        println(stderr, "Dual status: ", dual_status(model))
    end

    status = termination_status(model)

    if has_values(model)
        obj = JuMP.objective_value(model)
        x_star = JuMP.value.(x)
    else
        println(stderr, "Error: solution not found")
        exit(1)
    end

    return x_star, obj, Dates.value(end_time - start_time), status
end

# -------------------------------------------------------------------
# Main Function
# -------------------------------------------------------------------
function main()
    # Reads the file path from the terminal
    # Checks if the user provided any arguments
    if length(ARGS) < 1
        println(stderr, "Usage: julia roman.jl <arquivo.txt>")
        exit(1)
    end

    # Takes the first argument as the file path
    filepath = ARGS[1]

    # Verify whether the file exists
    if !isfile(filepath)
        println(stderr, "Error: file not found -> ", filepath)
        exit(1)
    end

    # Splits the path and the file extension
    filename, ext = splitext(filepath)

    # Read the graph and runs the solver 
    g = read_simple_graph(filepath)
    colors = viable_solution(g)
    greedy_span = maximum(values(colors))
    max_degree = maximum(degree(g, v) for v in vertices(g))
    min_degree = minimum(degree(g, v) for v in vertices(g))
    n = nv(g)  
    m = ne(g)  
    density = 2.0*m / (n*(n-1))
    
    @debug println("colors of greedy coloring: ", colors)

    TIME_LIMIT_MINUTES = 15 

    x, opt, elapsed_time, status = graceful_coloring(g, TIME_LIMIT_MINUTES, colors)

    @debug begin
        println("opt = ", opt)
        for v in vertices(g)
            println("label($v) = $(x[v])")
        end
    end

    # Creates the 'output' directory if it does not exist
    output_dir = "output"
    isdir(output_dir) || mkdir(output_dir)

    # Full path of the output CSV file
    result_file_path = joinpath(output_dir, basename(filename) * ".csv")

    # Writing the result in a CSV file
    open(result_file_path, "w") do io
        # Header
        println(io, "G,|V|,|E|,density,maxDegree,minDegree,greedy,span,time(milliseconds),status")
        # Row with results
        println(io, "$(basename(filename)),$(nv(g)),$(ne(g)),$density,$max_degree,$min_degree,$greedy_span,$opt,$elapsed_time,$status")
    end

    println(stdout, "Results saved in: $result_file_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
