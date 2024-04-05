# Demonstration of variable arrays bug in BilevelJuMP with MibS solver. 

# Model of the problem
# First level
# ```math
# \max_{{\bf x}} \sum_{i,j} x_{i,j},\\
# \notag s.t.\\
# x \in \{0, 1\}^{n \times n},\\
# ```
# Second level
# ```math
# \max_{{\bf y}} \sum_{i,j} y_{i,j},\\
# \notag s.t.\\
# \sum_{j \in [n]} y_{j,j} == 1,\\
# y \in \{0, 1\}^{n \times n}.\\
# ```

using BilevelJuMP
using MibS_jll
using Printf

function build_model(num_units, center_index)
    # Create model and assignment variables
    m = BilevelModel()
    @variable(Upper(m), x[1:num_units, 1:num_units], Bin)
    @variable(Lower(m), y[1:num_units, 1:num_units], Bin)

    # Add max-sum objectives
    @objective(Upper(m), Max, sum(x[j, j] for j in 1:num_units))
    @objective(Lower(m), Max, sum(y[j, j] for j in 1:num_units))

    # Add sum(y_jj) == 1, i.e., one unit is chosen as the center
    @constraint(Lower(m), c_num_dist_ub, sum(y[j, j] for j in 1:num_units) <= 1)
    @constraint(Lower(m), c_num_dist_lb, sum(y[j, j] for j in 1:num_units) >= 1)

    # Force the unit labeled center_index to be the center
    @constraint(Lower(m), c_only_center, y[center_index, center_index] >= 1)

    # Force all other units to NOT be the center
    for i = 1:num_units
        if i != center_index
            @constraint(Lower(m), y[i, i] <= 0)
        end
    end

    return m
end


### For each possible center_index, build and solve model using MibS
num_units = 16

for c_i = 1:num_units
    model = build_model(num_units, c_i)
    # display(model)
    solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
    for low_name in keys(solution.all_lower)
        if solution.all_lower[low_name] > 0
            @printf "When center_index is %i, the only nonzero lower variable is %s = %.2f. Expected: y[%i, %i]\n" c_i low_name solution.all_lower[low_name] c_i c_i
        end
    end
end

### Output: 
#=
When center_index is 1, the only nonzero lower variable is y[13,4] = 1.00. Expected: y[1, 1]
When center_index is 2, the only nonzero lower variable is y[5,16] = 1.00. Expected: y[2, 2]
When center_index is 3, the only nonzero lower variable is y[10,13] = 1.00. Expected: y[3, 3]
When center_index is 4, the only nonzero lower variable is y[6,1] = 1.00. Expected: y[4, 4]
When center_index is 5, the only nonzero lower variable is y[15,9] = 1.00. Expected: y[5, 5]
When center_index is 6, the only nonzero lower variable is y[12,10] = 1.00. Expected: y[6, 6]
When center_index is 7, the only nonzero lower variable is y[5,2] = 1.00. Expected: y[7, 7]
When center_index is 8, the only nonzero lower variable is y[12,2] = 1.00. Expected: y[8, 8]
When center_index is 9, the only nonzero lower variable is y[1,11] = 1.00. Expected: y[9, 9]
When center_index is 10, the only nonzero lower variable is y[4,3] = 1.00. Expected: y[10, 10]
When center_index is 11, the only nonzero lower variable is y[13,7] = 1.00. Expected: y[11, 11]
When center_index is 12, the only nonzero lower variable is y[7,10] = 1.00. Expected: y[12, 12]
When center_index is 13, the only nonzero lower variable is y[1,9] = 1.00. Expected: y[13, 13]
When center_index is 14, the only nonzero lower variable is y[10,3] = 1.00. Expected: y[14, 14]
When center_index is 15, the only nonzero lower variable is y[13,3] = 1.00. Expected: y[15, 15]
When center_index is 16, the only nonzero lower variable is y[15,5] = 1.00. Expected: y[16, 16]
=#
