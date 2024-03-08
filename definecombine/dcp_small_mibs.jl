# Solve small define-combine bilevel MIPs using 
# MibS solver from COIN-OR via BilevelJuMP interface. 

using BilevelJuMP
using MibS_jll

using Printf
using Test

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
    @objective(Upper(m), Min, 4)#Max, sum(x[j, j] for j in 1:num_units))
    @objective(Lower(m), Min, 17)#Max, sum(y[j, j] for j in 1:num_units))

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
    @constraint(Upper(m), c_num_subdist_ub, sum(x[j, j] for j in 1:num_units) <= 2 * num_districts)
    @constraint(Upper(m), c_num_subdist_lb, sum(x[j, j] for j in 1:num_units) >= 2 * num_districts)
    println(c_num_subdist_ub)
    println(c_num_subdist_lb)
    
    # Add y_jj <= x_jj constraints (y centers must also be x centers)
    # for j = 1:num_units
    #     @constraint(Lower(m), y[j, j] <= x[j, j])
    # end

    # Add sum(y_jj) == num_districts (exactly N district centers)
    @constraint(Lower(m), c_num_dist_ub, sum(y[j, j] for j in 1:num_units) <= 2) # num_districts
    @constraint(Lower(m), c_num_dist_lb, sum(y[j, j] for j in 1:num_units) >= 2) # num_districts
    println(c_num_dist_ub)
    println(c_num_dist_lb)
    
    println("Are these constraints being ignored?!?")
    # @constraint(Lower(m), c_center3, y[3, 3] <= 0)
    # println(c_center3)

    i1 = 2
    i2 = 5
    @constraint(Lower(m), c_only_centers, y[i1, i1] + y[i2, i2] >= 2)
    println(c_only_centers) 
    
    # @constraint(Lower(m), c_center3_2, y[3, 3] >= y[5, 2])
    # println(c_center3_2)
    
    # @constraint(Lower(m), y[5, 5] <= 0)
    # @constraint(Lower(m), y[7, 7] <= 0)
    # @constraint(Lower(m), y[9, 9] <= 0)
    # @constraint(Lower(m), y[11, 11] >= 1)
    println("No, but the indices seem to be changing...")

    # @constraint(Lower(m), y[5, 2] <= 0) # No noticeable effect for these four
    # @constraint(Lower(m), y[10, 13] <= 0)
    # @constraint(Lower(m), y[1, 11] <= 0)
    # @constraint(Lower(m), y[16, 9] <= 0)
    # @constraint(Lower(m), y[4, 13] <= 0)

    for i = 1:num_units
        if i != i1 && i != i2
            @constraint(Lower(m), y[i, i] <= 0)
        end
    end

    # TEMP (TODO: delete): add constraints forcing all variables to zero except (sub)district centers
    # for i = 1:num_units
    #     for j = 1:num_units
    #         if i != j
    #             @constraint(Upper(m), x[i, j] <= 0)
    #             @constraint(Lower(m), y[i, j] <= 0)
    #         end
    #     end
    # end

    # Add assignment constraints: every unit assigned exactly once
    # c_x_asst_ub = [
    #     @constraint(Upper(m), sum(x[i, j] for j in 1:num_units) <= 1)
    #     for i = 1:num_units
    # ]
    # c_x_asst_lb = [
    #     @constraint(Upper(m), sum(x[i, j] for j in 1:num_units) >= 1)
    #     for i = 1:num_units
    # ]
    # c_y_asst_ub = [
    #     @constraint(Lower(m), sum(y[i, j] for j in 1:num_units) <= 1)
    #     for i = 1:num_units
    # ]
    # c_y_asst_lb = [
    #     @constraint(Lower(m), sum(y[i, j] for j in 1:num_units) >= 1)
    #     for i = 1:num_units
    # ]
    # for i = 1:num_units
    #     set_name(c_x_asst_ub[i], "c_x_asst_ub[$(i)]")
    #     set_name(c_x_asst_lb[i], "c_x_asst_lb[$(i)]")
    #     set_name(c_y_asst_ub[i], "c_y_asst_ub[$(i)]")
    #     set_name(c_y_asst_lb[i], "c_y_asst_lb[$(i)]")
    # end

    # @constraint(Lower(m), c_y7_atmost_one, y[7, 7] + y[7, 11] + y[7, 14] <= 1)

    # @constraint(Lower(m), sum(y[i, j] for i in 1:num_units, j in 1:num_units) >= 16) # TEST, TODO: remove

    # Add population balance constraints
    # popbal_tol = 1.0e-4
    # x_ideal_district_pop = num_units ÷ (num_districts * 2)
    # x_popbal_lb = x_ideal_district_pop * (1.0 - popbal_tol)
    # x_popbal_ub = x_ideal_district_pop * (1.0 + popbal_tol)
    # y_ideal_district_pop = num_units ÷ num_districts
    # y_popbal_lb = y_ideal_district_pop * (1.0 - popbal_tol)
    # y_popbal_ub = y_ideal_district_pop * (1.0 + popbal_tol)
    # @printf "subdistrict pop bounds: %.4f to %.4f\n" x_popbal_lb x_popbal_ub
    # @printf "District    pop bounds: %.4f to %.4f\n" y_popbal_lb y_popbal_ub
    
    
    # c_x_popbal_ub = [
    #     @constraint(Upper(m), sum(x[i, j] for i in 1:num_units) <= x_popbal_ub * x[j, j])
    #     for j = 1:num_units
    # ]
    # c_x_popbal_lb = [
    #     @constraint(Upper(m), sum(x[i, j] for i in 1:num_units) >= x_popbal_lb * x[j, j])
    #     for j = 1:num_units
    # ]
    # c_y_popbal_ub = [
    #     @constraint(Lower(m), sum(y[i, j] for i in 1:num_units) <= y_popbal_ub * y[j, j])
    #     for j = 1:num_units
    # ]
    # c_y_popbal_lb = [
    #     @constraint(Lower(m), sum(y[i, j] for i in 1:num_units) >= y_popbal_lb * y[j, j])
    #     for j = 1:num_units
    # ]
    # for j = 1:num_units
    #     set_name(c_x_popbal_ub[j], "c_x_popbal_ub[$(j)]")
    #     set_name(c_x_popbal_lb[j], "c_x_popbal_lb[$(j)]")
    #     set_name(c_y_popbal_ub[j], "c_y_popbal_ub[$(j)]")
    #     set_name(c_y_popbal_lb[j], "c_y_popbal_lb[$(j)]")
    # end
    

    # Add some symmetry-breaking constraints. 
    # TODO: Add L-fixing/pos-based symmetry breaking. 
    # for i = 1:num_units
    #     for j = i+1:num_units
    #         @constraint(Upper(m), x[i, j] <= 0)
    #     end
    # end
    
    # for i = 1:num_units
    #     for j = i+1:num_units
    #         @constraint(Lower(m), y[i, j] <= 0)
    #     end
    # end

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
    # dir_edge_list = get_directed_grid_edges(num_rows, num_cols)
    # num_dir_edges = length(dir_edge_list)
    # # Get undirected edges (u, v) in canonical order (u < v)
    # undir_edge_list = [(e[1], e[2]) for e in dir_edge_list if e[1] < e[2]]
    # num_undir_edges = length(undir_edge_list)
    # @printf "Directed edges:  \t %i\n" num_dir_edges
    # @printf "Undirected edges:\t %i\n" num_undir_edges
    # # display(undir_edge_list)
    # println()

    # Add Shirabe flow variables for contiguity
    # TODO: It seems these have to be integral for MibS to be able to solve the bilevel MILP... 
    # @variable(Lower(m), f[1:num_units, 1:num_dir_edges] >= 0, Int)

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

x_asst = Dict() # x, or Definer, assignment
y_asst = Dict() # y, or Combiner, assignment

dist_centers = Set()
subdist_centers = Set()

x_sum = 0.0
y_sum = 0.0
y7_sum = 0.0

for i in 1:num_units
    for j in 1:num_units
        y_varname = string("y[", i, ",", j, "]")
        if solution.all_lower[y_varname] ≈ 1
            push!(dist_centers, j)
        end
    end
end

# TODO: I think the bug is actually in here, where I try to convert/recover the indices for assignments
# for i in 1:num_rows
#     for j in 1:num_cols
#         curr_index = rowcol_to_index(i - 1, j - 1, num_cols) + 1
#         for center_index in 1:num_units
#             x_varname = string("x[", curr_index, ",", center_index, "]")
#             if solution.all_upper[x_varname] ≈ 1
#                 x_asst[curr_index] = center_index
#                 if curr_index == center_index
#                     push!(subdist_centers, curr_index)
#                 end
#             end

#             y_varname = string("y[", curr_index, ",", center_index, "]")
#             global y_sum += solution.all_lower[y_varname]
#             if curr_index == 7
#                 global y7_sum += solution.all_lower[y_varname]
#             end
            
#             if solution.all_lower[y_varname] > 0 # ≈ 1
#                 @printf "%3i is assigned to %3i by combiner: (%3.2f)\n" curr_index center_index solution.all_lower[y_varname]
#                 y_asst[curr_index] = center_index
#                 if curr_index == center_index
#                     push!(dist_centers, curr_index)
#                 end
#             end
#         end
#     end
# end

@printf "Sum of y variables: %.1f\n" y_sum
@printf "Sum of y[7, j] variables: %.1f\n" y7_sum

println("\n*** Definer centers ***")
for unit in subdist_centers
    @printf "%2i, a.k.a. (%2i, %2i), is a subdistrict center.\n" unit index_to_rowcol_onebased(unit, num_cols)...
end

println("\n*** Definer assignment ***")
# for i in 1:num_rows
#     for j in 1:num_cols
#         curr_index = rowcol_to_index(i - 1, j - 1, num_cols) + 1
#         @printf "%3i" x_asst[curr_index]
#     end
#     println()
# end

println("\n*** Combiner centers ***")
for unit in dist_centers
    @printf "%2i, a.k.a. (%2i, %2i), is a district center.\n" unit index_to_rowcol_onebased(unit, num_cols)...
end

println("\n*** Combiner assignment ***")
# for i in 1:num_rows
#     for j in 1:num_cols
#         curr_index = rowcol_to_index(i - 1, j - 1, num_cols) + 1
#         assigned_center = haskey(y_asst, curr_index) ? string(y_asst[curr_index]) : "?"
#         @printf "%3s" assigned_center
#     end
#     println()
# end

# println(solution)
nz_lower = solution.nonzero_lower
println(nz_lower)
println(length(nz_lower))

for low_name in keys(solution.all_lower)
    # @printf "%s = %.2f\n" low_name solution.all_lower[low_name]
    if solution.all_lower[low_name] > 0
        @printf "%s is a nonzero lower var: %.2f\n" low_name solution.all_lower[low_name]
    end
end

# for nz_low_key in keys(nz_lower)
#     println(nz_low_key)
#     println(solution.all_lower[nz_low_key])
# end

