# Solve single-level redistricting MIP using Gurobi. 
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import random
import re

NUM_ROWS = 10
NUM_COLS = 20
NUM_DISTRICTS = 10

NUM_UNITS = NUM_ROWS * NUM_COLS
UNITS_PER_DISTRICT = NUM_UNITS // NUM_DISTRICTS

NUM_DEFINER_SUPPORTERS = int(NUM_UNITS * 0.55)
NUM_OTHER_SUPPORTERS = NUM_UNITS - NUM_DEFINER_SUPPORTERS

def index_to_row_col_tuple(index):
    return(index // NUM_COLS, index % NUM_COLS)

if __name__ == '__main__':
    try:
        random.seed(2021) # For reproducibility
        assert(NUM_UNITS // NUM_DISTRICTS == NUM_UNITS / NUM_DISTRICTS) # Require exact population balance

        unit_weights = [-1] * NUM_OTHER_SUPPORTERS + [1] * NUM_DEFINER_SUPPORTERS
        random.shuffle(unit_weights)

        # Create MIP model
        m = gp.Model("mip_redist")

        # Create variables x_{ij}, indicators for whether 
        # unit i is assigned to (the district centered at) unit j
        # unit_ids = [(i, j) for i in range(NUM_ROWS) for j in range(NUM_COLS)]
        unit_ids = [i for i in range(NUM_UNITS)]

        # Create (squared) distance matrix
        d = np.zeros((len(unit_ids), len(unit_ids)))
        for i in unit_ids:
            x1, y1 = index_to_row_col_tuple(i)
            for j in unit_ids:
                x2, y2 = index_to_row_col_tuple(j)
                d[i, j] = (x1 - x2)**2 + (y1 - y2)**2

        x = m.addVars(unit_ids, unit_ids, name="x", vtype=GRB.BINARY, obj=d)

        # Create edges for grid graph
        # x_edges = set()
        # for i, j in unit_ids:
        #     for u, v in unit_ids:
        #         if abs(i - u) + abs(j - v) == 1:
        #             x_edges.add((i,j), (u,v))
     
    
        # Set objective: maximize compactness in terms of Hess objective, squared distances
        # m.setObjective(d @ x, GRB.MAXIMIZE)
        # m.setObjective(47, GRB.MINIMIZE)

        # Add constraint: sum of x_{jj} = NUM_DISTRICTS
        m.addConstr(sum(x[i, i] for i in unit_ids) == NUM_DISTRICTS, "numdistricts")

        # # Add constraints: every unit fully assigned
        m.addConstrs(
            (x.sum(i, '*') == 1 
                for i in unit_ids), "fullyassigned")

        # # # Add constraints: every district has exactly NUM_UNITS / NUM_DISTRICTS units
        for j in unit_ids:
            m.addConstr(sum(x[i, j] for i in unit_ids) == UNITS_PER_DISTRICT * x[j, j], "districtsize %s" % (j))

        # Add constraints: can only assign to centers
        for i in unit_ids:
            for j in unit_ids:
                if i == j: continue; # Only need this constraint for distinct i, j
                m.addConstr(x[i, j] - x[j, j] <= 0, "onlycenters[%s, %s]" % (i, j))
        
        # Optimize model
        m.optimize()
        print('Obj: %g' % m.objVal)

        # Check which units are centers, and extract partition
        partition_dict = {}
        assignment_dict = {}
        tol = m.Params.IntFeasTol

        for v in m.getVars():
            var_name_parts = re.split('\[|,|\]', v.varName)
            i = int(var_name_parts[-3])
            j = int(var_name_parts[-2])

            if v.x >= 1 - tol:
                if j not in partition_dict:
                    j_row, j_col = index_to_row_col_tuple(j)
                    print('Unit %g at (%g, %g) is a center.' % (j, j_row, j_col))
                    partition_dict[j] = set()
                
                partition_dict[j].add(i)
                assignment_dict[i] = j
           
        print('\nPartition dict:')
        print(partition_dict)
        print()

        print('\nAssignment dict:')
        print(assignment_dict)
        print()

        # Print visualization of solution (district plan) using letters for districts
        centers = sorted([key for key in partition_dict])
        center_to_char = {centers[i] : chr(65 + i) for i in range(len(centers))}

        index = 0
        while index < NUM_UNITS:
            if index % NUM_COLS == 0:
                print()
            
            i, j = index_to_row_col_tuple(index)
            print(center_to_char[assignment_dict[index]], end='')

            index += 1
        
        print()

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError as e:
        print('Encountered an attribute error.')
