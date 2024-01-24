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
Convert (row, col) of a grid cell to index. 
Index starts at 0 in top-left corner and increases 
left-to-right, top-to-bottom. 
"""
function rowcol_to_index(row, col, num_cols)
    return row * num_cols + col
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
    supported_dimensions = [(20, 20), (10, 10), (8, 8)] # , (4, 4)
    
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




index = 11
num_rows = 8
num_cols = 8

d2 = populate_squared_distance_matrix(num_rows, num_cols)
display(d2)

row, col = index_to_rowcol(index, num_cols)
@printf "\nrow = %i, col = %i\n" row col

pos, size_of_B  = shuffle_indices(num_rows, num_cols)
@printf "\nL-fixing: B contains %i units.\n" size_of_B