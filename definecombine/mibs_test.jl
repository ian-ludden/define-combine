# # MibS: Example 1 (Experimental feature)

# Model of the problem
# First level
# ```math
# \min_{x} -3x -7y,\\
# \notag s.t.\\
# -3x + 2y \leq 12,\\ 
# x + 2y \leq 20,\\
# x \leq 10,\\
# x \in \mathbb{Z},\\
# ```
# Second level
# ```math
# \min_{y} y,\\
# \notag s.t.\\
# 2x - y <= 7,\\
# -2x + 4y <= 16,\\
# y <= 5\\
# y \in \mathbb{Z}\\
# ```

using BilevelJuMP
using Printf
using Test
using MibS_jll

function test_example_one()
    # Taken from https://github.com/joaquimg/BilevelJuMP.jl/blob/master/docs/src/examples/MibS_example1.jl
    model = BilevelModel()

    # First we need to create all of the variables in the upper and lower problems:
    
    # Upper level variables
    @variable(Upper(model), x, Int)
    
    #Lower level variables
    @variable(Lower(model), y, Int)
    
    # Then we can add the objective and constraints of the upper problem:
    
    # Upper level objecive function
    @objective(Upper(model), Min, -3x - 7y)
    
    # Upper constraints
    @constraints(Upper(model), begin
        u1, -3x + 2y <= 12
        u2, x + 2y <= 20
        u3, x <= 10
    end)
    
    # Followed by the objective and constraints of the lower problem:
    
    # Lower objective function
    @objective(Lower(model), Min, y)
    
    # Lower constraints
    @constraint(Lower(model), l1, 2x - y <= 7)
    @constraint(Lower(model), l2, -2x + 4y <= 16)
    @constraint(Lower(model), l3, y <= 5)
    
    # Using MibS Solver
    solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
    @printf "Status: %s\n" solution.status
    @printf "Opt obj.: %.1f\n" solution.objective
    @printf "x value: %.1f\n" solution.all_upper["x"]
    @printf "y value: %.1f\n" solution.all_lower["y"]
    
    @printf "\nUpper variables:\n"
    for var ∈ solution.all_upper
        @printf "var: %s\n" var
    end
    
    @printf "\nLower variables:\n"
    for var ∈ solution.all_lower
        @printf "var: %s\n" var
    end
    
    # Auto testing
    @test solution.status == true
    @test solution.objective ≈ -53
    @test solution.nonzero_upper == Dict(0 => 6.0)
    @test solution.nonzero_lower == Dict(0 => 5.0)
    @test solution.all_upper["x"] == 6.0
    @test solution.all_lower["y"] == 5.0
    return 0
end

function test_example_two()
    # Same problem as in build_sample_model_two of pyomo_pao_example.py, 
    # Example 1 of Moore + Bard (1990)
    
    model = BilevelModel()
    # Upper level variables
    @variable(Upper(model), x, Int)
    
    #Lower level variables
    @variable(Lower(model), y, Int)
    
    # Upper level objecive function
    @objective(Upper(model), Max, x + 10y)
    
    # Upper constraints
    @constraints(Upper(model), begin
        u1, x >= 0
    end)
    
    # Lower objective function
    @objective(Lower(model), Max, -y)
    
    # Lower constraints
    @constraint(Lower(model), l1, -25x +20y <= 30)
    @constraint(Lower(model), l2, x + 2y <= 10)
    @constraint(Lower(model), l3, 2x - y <= 15)
    @constraint(Lower(model), l4, 2x + 10y >= 15)
    @constraint(Lower(model), l5, y >= 0)
    
    # Using MibS Solver
    solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
    # println(solution)
    @printf "Status: %s\n" solution.status
    @printf "Opt obj.: %.1f\n" solution.objective
    @printf "x value: %.1f\n" solution.all_upper["x"]
    @printf "y value: %.1f\n" solution.all_lower["y"]
    
    @printf "\nUpper variables:\n"
    for var ∈ solution.all_upper
        @printf "var: %s\n" var
    end
    
    @printf "\nLower variables:\n"
    for var ∈ solution.all_lower
        @printf "var: %s\n" var
    end
    
    # Auto testing
    @test solution.status == true
    @test solution.objective ≈ -22.000000000001
    # @test solution.nonzero_upper == Dict(0 => 6.0)
    # @test solution.nonzero_lower == Dict(0 => 5.0)
    @test solution.all_upper["x"] == 2.0
    @test solution.all_lower["y"] == 2.0
    return 0
end


function test_example_three()
    # Same problem as in build_sample_model_three of pyomo_pao_example.py, 
    # from PAO examples (Bard 1998, 5.1.1): 
    # https://pao.readthedocs.io/en/latest/examples.html
    
    model = BilevelModel()
    # Upper level variables
    @variable(Upper(model), x, Int) # Default type is Real
    
    #Lower level variables
    @variable(Lower(model), y, Int) # Default type is Real
    
    # Upper level objecive function
    @objective(Upper(model), Min, x - 4y)
    
    # Upper constraints
    @constraints(Upper(model), begin
        u1, x >= 0
    end)
    
    # Lower objective function
    @objective(Lower(model), Min, y)
    
    # Lower constraints
    @constraint(Lower(model), l1, -x -y <= -3)
    @constraint(Lower(model), l2, -2x + y <= 0)
    @constraint(Lower(model), l3, 2x + y <= 12)
    @constraint(Lower(model), l4, 3x - 2y <= 4)
    @constraint(Lower(model), l5, y >= 0)
    
    # Using MibS Solver
    solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
    # println(solution)
    @printf "Status: %s\n" solution.status
    @printf "Opt obj.: %.1f\n" solution.objective
    @printf "x value: %.1f\n" solution.all_upper["x"]
    @printf "y value: %.1f\n" solution.all_lower["y"]
    
    @printf "\nUpper variables:\n"
    for var ∈ solution.all_upper
        @printf "var: %s\n" var
    end
    
    @printf "\nLower variables:\n"
    for var ∈ solution.all_lower
        @printf "var: %s\n" var
    end
    
    # Auto testing
    @test solution.status == true
    @test solution.objective ≈ -12.000000000001
    # @test solution.nonzero_upper == Dict(0 => 6.0)
    # @test solution.nonzero_lower == Dict(0 => 5.0)
    @test solution.all_upper["x"] == 4.0
    @test solution.all_lower["y"] == 4.0
    return 0
end

test_example_one()
println("\n")
test_example_two()
println("\n")
test_example_three()
