# Solve small define-combine bilevel MIPs using 
# MibS solver from COIN-OR via BilevelJuMP interface. 

using BilevelJuMP
using MibS_jll

using Printf
using Test

ENV["TMPDIR"] = joinpath(pwd(), "temporary")
mkdir(ENV["TMPDIR"])
mktempdir()

"""
Convert index of a grid cell to (row, col) pair. 
Index starts at 0 in top-left corner and increases 
left-to-right, top-to-bottom. 
"""
function index_to_rowcol(index, num_cols)
    return (index ÷ num_cols, index % num_cols)
end

"""
Convert index of a grid cell to (row, col) pair. 
Index starts at 1 in top-left corner and increases 
left-to-right, top-to-bottom. 
Row and column also start at 1, 
so the top-left corner is (1, 1). 
"""
function index_to_rowcol_onebased(index, num_cols)
    rowcol_zerobased = index_to_rowcol(index - 1, num_cols)
    return (rowcol_zerobased[1] + 1, rowcol_zerobased[2] + 1)
end

"""
Convert (row, col) of a grid cell to index. 
Index starts at 0 in top-left corner and increases 
left-to-right, top-to-bottom. 
"""
function rowcol_to_index(row, col, num_cols)
    return row * num_cols + col
end

"""
Return all ordered pairs of edges, 
using indices starting at 0 in the top-left corner and 
moving left-to-right, top-to-bottom, 
in a grid graph with the given number of rows and columns. 
"""
function get_directed_grid_edges(num_rows, num_cols)
    undir_edges = Vector{Tuple{Int, Int}}()

    for i = 0:num_rows-1
        for j = 0:num_cols-1
            current = rowcol_to_index(i, j, num_cols)

            # Right edge
            right_neighbor = current + 1
            if j + 1 < num_cols
                push!(undir_edges, (current + 1, right_neighbor + 1))
            end

            # Down edge
            down_neighbor = current + num_cols
            if i + 1 < num_rows
                push!(undir_edges, (current + 1, down_neighbor + 1))
            end
        end
    end

    dir_edges = Vector{Tuple{Int, Int}}()
    for i = 1:length(undir_edges)
        u, v = undir_edges[i]
        push!(dir_edges, (u, v))
        push!(dir_edges, (v, u))
    end

    return dir_edges
end

"""
Assuming a `num_rows` by `num_cols` grid of units, 
populate a matrix (2D array) representing all pairwise 
squared Euclidean distances between unit centers. 
"""
function populate_squared_distance_matrix(num_rows, num_cols)
    num_units = num_rows * num_cols
    D = zeros(num_units, num_units)

    for i = 1:num_units
        row_i, col_i = index_to_rowcol(i - 1, num_cols)
        for j = 1:num_units
            row_j, col_j = index_to_rowcol(j - 1, num_cols)
            D[i, j] = (row_i - row_j)^2 + (col_i - col_j)^2
        end
    end
    
    return D
end

"""
For L-fixing. Carefully constructs an ordering 
to put a large set B of indices at the end 
such that none of them can serve as centers. 

TODO: Should also take parameter num_districts, 
or L (lower bound on district size), 
since the maximum allowed size of a connected component in G[B]
depends on the number of units per district. 

Returns new positions (pos[i] is the new position of unit i) 
and the size of the set B. 
"""
function shuffle_indices(num_rows, num_cols)
    supported_dimensions = [(20, 20), (10, 10), (8, 8), (4, 4)]
    
    num_districts = 4 # TODO: Allow other numbers of districts
    
    if !((num_rows, num_cols) in supported_dimensions)
        println("shuffle_indices is not yet implemented for general grid dimensions.\nCurrently supported: ", supported_dimensions)
        return 1:(num_rows * num_cols), 0 # default ordering
    end

    num_units = num_rows * num_cols
    
    rows_in_B = nothing
    rows_not_in_B = nothing
    
    if num_rows == 20
        rows_in_B = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]
        rows_not_in_B = [6, 13]
    elseif num_rows == 10
        rows_in_B = [0, 1, 2, 3, 5, 6, 7, 8]
        rows_not_in_B = [4, 9]
    elseif num_rows == 8
        rows_in_B = [0, 1, 2, 4, 5, 6]
        rows_not_in_B = [3, 7]
    elseif num_rows == 4
        rows_in_B = [0, 2]
        rows_not_in_B = [1, 3]
    end
    println([r for r in rows_in_B])
    cols_in_B = rows_in_B # By symmetry; would need to change to support asymmetric grid dimensions
    cols_not_in_B = rows_not_in_B

    pairs_in_B = [(i, j) for i in rows_in_B for j in cols_in_B]
    pairs_not_in_B = []
    for i = 0:num_rows-1
        for j = 0:num_cols-1
            if i in rows_not_in_B || j in cols_not_in_B
                pairs_not_in_B = vcat(pairs_not_in_B, (i, j))
            end
        end
    end

    all_pairs_sorted = vcat(pairs_not_in_B, pairs_in_B)
    
    pos = zeros(Int, num_units)
    for i = 1:num_units
        row, col = index_to_rowcol(i - 1, num_cols)
        pos[i] = findfirst(x -> x == (row, col), all_pairs_sorted)
    end

    return pos, length(pairs_in_B)
end

"""
Build a bilevel MILP model for the define-combine procedure 
on a grid graph with the given number of rows, columns, and districts. 

Then, solve the model using MibS and return the solution. 

Note: The `voter_grid` parameter is a 2D float array of units'
definer fractional vote-shares (between 0.0 and 1.0). 
"""
function build_and_solve_grid_model(num_rows, num_cols, num_districts, voter_grid)
    num_units = num_rows * num_cols 
    num_units_per_district = num_units ÷ num_districts
    big_M = num_units_per_district + 1
    epsilon = 1.0e-8 # Needs to be smaller than smallest possible fractional part of a district's net votes
    
    # Create model, x (definer assignment), and y (combiner assignment) variables
    m = BilevelModel()
    @variable(Upper(m), x[1:num_units, 1:num_units], Bin)
    @variable(Lower(m), y[1:num_units, 1:num_units], Bin)

    # Add objectives: definer's net utility, combiner's net utility. 
    # Need new variables, a_j for definer win and b_j for combiner win. 
    # @variable(Lower(m), a[1:num_units], Bin)
    # @variable(Lower(m), b[1:num_units], Bin)
    # for i = 1:num_units
    #     @constraint(Lower(m), a[i] + b[i] <= 1)
    # end
    # @objective(Upper(m), Max, sum(a) - sum(b))
    # @objective(Lower(m), Max, sum(b) - sum(a))
    @objective(Upper(m), Min, sum(y[j, j] for j in 1:num_units))
    @objective(Lower(m), Max, sum(y[j, j] for j in 1:num_units))

    # Add big-M constraints for defining a_j and b_j using voter_grid and y_ij
    # for j = 1:num_units
    #     @constraint(Lower(m), sum((voter_grid[index_to_rowcol_onebased(i, num_cols)...] - 0.5) * y[i, j] for i in 1:num_units) 
    #         <= big_M * a[j])
    #     @constraint(Lower(m), sum((voter_grid[index_to_rowcol_onebased(i, num_cols)...] - 0.5) * y[i, j] for i in 1:num_units) 
    #         >= epsilon - big_M * (1 - a[j]))
    # end

    # for j = 1:num_units
    #     @constraint(Lower(m), sum((0.5 - voter_grid[index_to_rowcol_onebased(i, num_cols)...]) * y[i, j] for i in 1:num_units) 
    #         <= big_M * b[j])
    #     @constraint(Lower(m), sum((0.5 - voter_grid[index_to_rowcol_onebased(i, num_cols)...]) * y[i, j] for i in 1:num_units) 
    #         >= epsilon - big_M * (1 - b[j]))
    # end

    # Add sum(x_jj) <= 2 * num_districts (only 2N subdistrict centers)
    @constraint(Upper(m), sum(x[j, j] for j in 1:num_units) <= 2 * num_districts)
    
    # Add y_jj <= x_jj constraints (y centers must also be x centers)
    for j = 1:num_units
        @constraint(Lower(m), y[j, j] <= x[j, j])
    end

    # Add assignment constraints: every unit assigned exactly once
    for i = 1:num_units
        @constraint(Upper(m), sum(x[i, j] for j in 1:num_units) == 1)
        @constraint(Lower(m), sum(y[i, j] for j in 1:num_units) == 1)
    end

    # Add population balance constraints
    popbal_tol = 1.0e-4
    x_ideal_district_pop = num_units ÷ (num_districts * 2)
    x_popbal_lb = x_ideal_district_pop * (1.0 - popbal_tol)
    x_popbal_ub = x_ideal_district_pop * (1.0 + popbal_tol)
    y_ideal_district_pop = num_units ÷ num_districts
    y_popbal_lb = y_ideal_district_pop * (1.0 - popbal_tol)
    y_popbal_ub = y_ideal_district_pop * (1.0 + popbal_tol)
    
    for j = 1:num_units
        continue
        # @constraint(Upper(m), sum(x[i, j] for i in 1:num_units) <= x_popbal_ub)
        # @constraint(Upper(m), sum(x[i, j] for i in 1:num_units) >= x_popbal_lb)
        # @constraint(Lower(m), sum(y[i, j] for i in 1:num_units) <= y_popbal_ub)
        # @constraint(Lower(m), sum(y[i, j] for i in 1:num_units) >= y_popbal_lb)
    end

    # Add some symmetry-breaking constraints. 
    # TODO: Add L-fixing/pos-based symmetry breaking. 
    for i = 1:num_units
        for j = i+1:num_units
            @constraint(Upper(m), x[i, j] <= 0)
        end
    end
    
    for i = 1:num_units
        for j = i+1:num_units
            @constraint(Lower(m), y[i, j] <= 0)
        end
    end

    # Add constraints that y assignments follow x assignments: 
    # If x_ij = 1, then y_ik <= y_jk for all units i, j, k. 
    # for k = 1:num_units
    #     for i = 1:num_units
    #         for j = 1:num_units
    #             if (i < j) && (i < k) && (j < k) # TODO: Use positions (result from L-fixing) rather than raw indices
    #                 @constraint(Lower(m), y[i, k] <= y[j, k] + (1 - x[i, j]))
    #             end
    #         end
    #     end
    # end

    # Get directed edge set
    dir_edge_list = get_directed_grid_edges(num_rows, num_cols)
    num_dir_edges = length(dir_edge_list)
    # Get undirected edges (u, v) in canonical order (u < v)
    undir_edge_list = [(e[1], e[2]) for e in dir_edge_list if e[1] < e[2]]
    num_undir_edges = length(undir_edge_list)
    @printf "Directed edges:  \t %i\n" num_dir_edges
    @printf "Undirected edges:\t %i\n" num_undir_edges
    # display(undir_edge_list)
    println()

    # Add Shirabe flow variables for contiguity
    # TODO: It seems these have to be integral for MibS to be able to solve the bilevel MILP... 
    @variable(Lower(m), f[1:num_units, 1:num_dir_edges] >= 0, Int)

    # Add cut edge variables with quadratic constraints to define them
    # TODO: May have to manually linearize these quadratic constraints if MibS is unhappy
    # @variable(Upper(m), c[1:num_undir_edges, 1:num_units, 1:num_units], Bin)
    # for (index, e) in enumerate(undir_edge_list)
    #     for k = 1:num_units
    #         for ell = 1:num_units
    #             @constraint(Upper(m), c[index, k, ell] == x[e[1], k] * x[e[2], ell])
    #         end
    #     end
    # end

    # Add constraints linking y with x via c
    # for k = 1:num_units
    #     for ell = 1:num_units
    #         @constraint(Lower(m), y[k, ell] <= sum(c[e, k, ell] for e in 1:num_undir_edges))
    #     end
    # end
    
    return m
end



index = 11
num_rows = 4
num_cols = 4
num_units = num_rows * num_cols
num_districts = 4
voter_grid = zeros(Float64, 4, 4)
voter_grid[1:2, 1:4] .= 1.0

# d2 = populate_squared_distance_matrix(num_rows, num_cols)
# display(d2)

row, col = index_to_rowcol(index, num_cols)
@printf "\nrow = %i, col = %i\n" row col

pos, size_of_B  = shuffle_indices(num_rows, num_cols)
@printf "\nL-fixing: B contains %i units.\n" size_of_B

model = build_and_solve_grid_model(num_rows, num_cols, num_districts, voter_grid)
display(model)

solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs, verbose_results=true)
# println(solution)
@printf "Status: %s\n" solution.status
@printf "Definer opt obj.: %.1f\n" solution.objective

asst = Dict() # y, or Combiner, assignment

for i in 1:num_rows
    for j in 1:num_cols
        curr_index = rowcol_to_index(i - 1, j - 1, num_cols) + 1
        for center_index in 1:num_units
            varname = string("y[", curr_index, ",", center_index, "]")
            if solution.all_lower[varname] ≈ 1
                asst[curr_index] = center_index
            end
        end
    end
end

for i in 1:num_rows
    for j in 1:num_cols
        curr_index = rowcol_to_index(i - 1, j - 1, num_cols) + 1
        @printf "%3i" asst[curr_index]
    end
    println()
end
