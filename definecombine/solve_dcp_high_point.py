# Solve define-combine high-point problem (leader problem with follower feasibility constraints). 
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import random
import re

NUM_ROWS = 4
NUM_COLS = 6
NUM_DISTRICTS = 4
NUM_SUBDISTRICTS = 2 * NUM_DISTRICTS

NUM_UNITS = NUM_ROWS * NUM_COLS
UNITS_PER_DISTRICT = NUM_UNITS // NUM_DISTRICTS
UNITS_PER_SUBDISTRICT = NUM_UNITS // NUM_SUBDISTRICTS

NUM_DEFINER_SUPPORTERS = int(NUM_UNITS * 0.50)
NUM_OTHER_SUPPORTERS = NUM_UNITS - NUM_DEFINER_SUPPORTERS

def index_to_row_col_tuple(index):
    return(index // NUM_COLS, index % NUM_COLS)

if __name__ == '__main__':
    try:
        random.seed(2021) # For reproducibility
        assert(NUM_UNITS // NUM_SUBDISTRICTS == NUM_UNITS / NUM_SUBDISTRICTS) # Require exact population balance

        print("Definer has", NUM_DEFINER_SUPPORTERS, "supporters.")

        unit_weights = [-1] * NUM_OTHER_SUPPORTERS + [1] * NUM_DEFINER_SUPPORTERS
        random.shuffle(unit_weights)
        unit_ids = [i for i in range(NUM_UNITS)]

        # Print voter distribution
        index = 0
        while index < NUM_UNITS:
            if index % NUM_COLS == 0:
                print()
            
            weight = unit_weights[index]
            output_char = 1 if weight == 1 else 0 # Use 0 to represent -1 since it's a single character
            print(output_char, end='')

            index += 1
        
        print()

        # Create MIP model
        m = gp.Model("mip_redist")

        # Create variables x_{ij}, indicators for whether 
        # unit i is assigned to (the subdistrict centered at) unit j
        x = m.addVars(unit_ids, unit_ids, name="x", vtype=GRB.BINARY)

        # Create variables y_{ij}, indicators for whether 
        # unit i is assigned to (the district centered at) unit j 
        y = m.addVars(unit_ids, unit_ids, name="y", vtype=GRB.BINARY)

        # Create variables a_{j} and b_{j}, indicators for whether 
        # the district centered at j is a win for party a (resp., party b)
        a = m.addVars(unit_ids, name="a", vtype=GRB.BINARY, obj=-1)
        b = m.addVars(unit_ids, name="b", vtype=GRB.BINARY, obj=1)

        # Create edges for grid graph
        x_edges = set()
        for u in unit_ids:
            row_u, col_u = index_to_row_col_tuple(u)
            for v in unit_ids:
                row_v, col_v = index_to_row_col_tuple(v)
                if abs(row_u - row_v) + abs(col_u - col_v) == 1:
                    x_edges.add((u, v))
                    x_edges.add((v, u))

        edges_in = []
        edges_out = []
        for u in unit_ids:
            edges_in_u = []
            edges_out_u = []
            for x_edge in x_edges:
                if x_edge[1] == u:
                    edges_in_u.append(x_edge)
                if x_edge[0] == u:
                    edges_out_u.append(x_edge)
            
            edges_in.append(edges_in_u)
            edges_out.append(edges_out_u)
        
        all_directed_edges_possible_duplicates = []
        for u in unit_ids:
            all_directed_edges_possible_duplicates.extend(edges_in[u])
            all_directed_edges_possible_duplicates.extend(edges_out[u])

        all_directed_edges = list(set(all_directed_edges_possible_duplicates))

        # Create flow variables for SHIR contiguity
        f = m.addVars(unit_ids, all_directed_edges, name="f", vtype=GRB.CONTINUOUS) # use default lower bound of 0

        # Create cut-edges variables
        c = m.addVars(all_directed_edges, unit_ids, unit_ids, name="c", vtype=GRB.BINARY)

        # Create both-units-are-subdistrict-centers variables
        z = m.addVars(unit_ids, unit_ids, name="z", vtype=GRB.BINARY)
    

        # Add constraint: sum of x_{jj} = NUM_SUBDISTRICTS
        m.addConstr(sum(x[i, i] for i in unit_ids) == NUM_SUBDISTRICTS, "num_subdistricts")

        # # Add constraints: every unit fully assigned
        m.addConstrs(
            (x.sum(i, '*') == 1 
                for i in unit_ids), "fullyassigned")

        # # # Add constraints: every subdistrict has exactly NUM_UNITS / NUM_SUBDISTRICTS units
        for j in unit_ids:
            m.addConstr(sum(x[i, j] for i in unit_ids) == UNITS_PER_SUBDISTRICT * x[j, j], "districtsize %s" % (j))

        # Add constraints: can only assign to centers
        for i in unit_ids:
            for j in unit_ids:
                if i == j: continue; # Only need this constraint for distinct i, j
                m.addConstr(x[i, j] - x[j, j] <= 0, "onlycenters[%s, %s]" % (i, j))
        
        # Add constraints: 
        # a_{j} is 1 if and only if definer, party a, wins district centered at j;
        # b_{j} is 1 if and only if combiner, barty b, wins district centered at j. 
        for j in unit_ids:
            m.addConstr(sum(y[i, j] * unit_weights[i] 
                for i in unit_ids) <= UNITS_PER_DISTRICT * a[j], 
                "a def1 %s" % (j))
        
            m.addConstr(sum(y[i, j] * unit_weights[i]
                for i in unit_ids) >= 1 - (UNITS_PER_DISTRICT + 1) * (1 - a[j]), 
                "a def2 %s" % (j))

            m.addConstr(sum(y[i, j] * unit_weights[i] * (-1)
                for i in unit_ids) <= UNITS_PER_DISTRICT * b[j], 
                "b def1 %s" % (j))

            m.addConstr(sum(y[i, j] * unit_weights[i] * (-1)
                for i in unit_ids) >= 1 - (UNITS_PER_DISTRICT + 1) * (1 - b[j]), 
                "b def2 %s" % (j))

        
        # Add flow-based (SHIR) contiguity constraints (alternative is separator-based, CUT, from Lykhovyd, Validi & Buchanan 2021)
        for j in unit_ids:
            for u in unit_ids:
                if u == j: continue
                # 1. Flow out minus flow in is equal to assignment
                m.addConstr(sum(f[j, uv[0], uv[1]] for uv in edges_out[u]) 
                            - sum(f[j, vu[0], vu[1]] for vu in edges_in[u]) 
                            == x[u, j], "flow is assignment %s %s" % (j, u))
                # 2. Flow is at most UNITS_PER_DISTRICT times assignment
                m.addConstr(sum(f[j, uv[0], uv[1]] for uv in edges_out[u]) <= (UNITS_PER_SUBDISTRICT) * x[u, j], "flow capacity %s %s" % (j, u))
            
            # 3. Flow leaving center to itself is zero
            m.addConstr(sum(f[j, jv[0], jv[1]] for jv in edges_out[j]) == 0, "flow zero %s" % (j))


        # Add constraints for cut-edge variables
        for k in unit_ids:
            for ell in unit_ids:
                for i, j in all_directed_edges:
                    m.addConstr(c[i, j, k, ell] <= x[i, k], "c def 1 %s %s %s %s" % (i, j, k, ell))
                    m.addConstr(c[i, j, k, ell] <= x[j, ell], "c def 2 %s %s %s %s" % (i, j, k, ell))
                    m.addConstr(c[i, j, k, ell] >= x[i, k] + x[j, ell] - 1, "c def 3 %s %s %s %s" % (i, j, k, ell))

        # Add constraints for both-units-are-subdistrict-centers variables
        for k in unit_ids:
            for ell in unit_ids:
                if ell == k: continue
                m.addConstr(z[k, ell] <= x[k, k], "z def 1 %s %s" % (k, ell))
                m.addConstr(z[k, ell] <= x[ell, ell], "z def 2 %s %s" % (k, ell))
                m.addConstr(z[k, ell] >= x[k, k] + x[ell, ell] - 1, "z def 3 %s %s" % (k, ell))

        # Add constraints: y centers follow x centers
        for j in unit_ids:
            m.addConstr(y[j, j] <= x[j, j], "y center only if x center %s" % (j))
        
        # Add constraint: NUM_DISTRICTS total y centers
        m.addConstr(sum(y[j, j] for j in unit_ids) == NUM_DISTRICTS, "num_districts")
        
        # Add constraints: y assignments only to centers
        for j in unit_ids:
            for i in unit_ids:
                m.addConstr(y[i, j] <= y[j, j], "y assign only to centers %s %s" % (i, j))

        # Add constraints: each unit assigned in y
        for i in unit_ids:
            m.addConstr(sum(y[i, j] for j in unit_ids) == 1, "y assign unit %s to exactly one center" % (i))
        
        # Add constraints: if x_{ij} = 1, then y_{ik} <= y_{jk} for all k (only need for k not equal to i nor j)
        for i in unit_ids:
            for j in unit_ids:
                for k in unit_ids:
                    if k == i or k == j: continue
                    m.addConstr(y[i, k] <= y[j, k] + (1 - x[i, j]), "subdistricts follow centers %s %s %s" % (i, j, k))
        
        # Add constraints: subdistrict center assigned to subdistrict center requires at least one cut-edge between subdistricts
        for k in unit_ids:
            for ell in unit_ids:
                if ell == k: continue
                m.addConstr(y[k, ell] <= (sum(c[i,j,k,ell] for i,j in all_directed_edges)) + (1 - z[k, ell]), "adjacent subdistricts %s %s" % (k, ell))

        # Optimize model
        m.optimize()
        print('Obj: %g' % m.objVal)

        # Check which units are subdistrict centers, and extract subdistrict partition
        partition_dict = {}
        assignment_dict = {}
        tol = m.Params.IntFeasTol

        for v in m.getVars():
            var_name_parts = re.split('\[|,|\]', v.varName)
            if not (var_name_parts[0] == 'x'): continue
            i = int(var_name_parts[-3])
            j = int(var_name_parts[-2])

            if v.x >= 1 - tol:
                if j not in partition_dict:
                    j_row, j_col = index_to_row_col_tuple(j)
                    print('Unit %g at (%g, %g) is a subdistrict center.' % (j, j_row, j_col))
                    partition_dict[j] = set()
                
                partition_dict[j].add(i)
                assignment_dict[i] = j
           
        print('\nPartition dict:')
        print(partition_dict)
        print()


        y_partition_dict = {}
        y_assignment_dict = {}

        for v in m.getVars():
            var_name_parts = re.split('\[|,|\]', v.varName)
            if not (var_name_parts[0] == 'y'): continue
            i = int(var_name_parts[-3])
            j = int(var_name_parts[-2])

            if v.x >= 1 - tol:
                if j not in y_partition_dict:
                    j_row, j_col = index_to_row_col_tuple(j)
                    print('Unit %g at (%g, %g) is a district center.' % (j, j_row, j_col))
                    y_partition_dict[j] = set()
                
                y_partition_dict[j].add(i)
                y_assignment_dict[i] = j
           
        print('\ny partition dict:')
        print(y_partition_dict)
        print()


        count_a_wins = 0
        count_b_wins = 0
        for key in y_partition_dict:
            district_weight = 0
            for unit_id in y_partition_dict[key]:
                district_weight += unit_weights[unit_id]
            print('District centered at unit', key, 'has total weight', district_weight)
            if district_weight > 0:
                count_a_wins += 1
            if district_weight < 0:
                count_b_wins += 1
            
        print('Wins for a:\t', count_a_wins)
        print('Wins for b:\t', count_b_wins)
        print('Ties:      \t', NUM_DISTRICTS - count_a_wins - count_b_wins)

        print('\ny assignment dict:')
        print(y_assignment_dict)
        print()

        # Print visualization of solution (district plan) using letters for districts
        centers = sorted([key for key in y_partition_dict])
        center_to_char = {centers[i] : chr(65 + i) for i in range(len(centers))}

        index = 0
        while index < NUM_UNITS:
            if index % NUM_COLS == 0:
                print()
            
            print(center_to_char[y_assignment_dict[index]], end='')

            index += 1
        
        print()

        # Also print visualization of subdistrict plan, using letters for subdistricts
        centers = sorted([key for key in partition_dict])
        center_to_char = {centers[i] : chr(65 + i) for i in range(len(centers))}

        index = 0
        while index < NUM_UNITS:
            if index % NUM_COLS == 0:
                print()
            
            print(center_to_char[assignment_dict[index]], end='')

            index += 1
        
        print()

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError as e:
        print('Encountered an attribute error.')
