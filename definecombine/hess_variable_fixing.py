# Fix variables in Hess redistricting model, using techniques from 
# Validi and Buchanan (2022), "Political districting to minimize cut edges"
import gurobipy as gp
from gurobipy import GRB
import numpy as np

TOL = 1e-6

def populate_squared_distance_matrix(num_rows, num_cols):
    """
    Assuming a num_rows by num_cols grid of units, 
    populate a matrix representing all pairwise 
    squared Euclidean distances between unit centers. 
    """
    num_units = num_rows * num_cols
    D = np.zeros((num_units, num_units))

    for i in range(num_units):
        row_i, col_i = index_to_rowcol(i, num_cols)
        for j in range(num_units):
            row_j, col_j = index_to_rowcol(j, num_cols)
            D[i, j] = (row_i - row_j)**2 + (col_i - col_j)**2
    
    return D

def index_to_rowcol(index, num_cols):
    return (index // num_cols, index % num_cols)

def rowcol_to_index(row, col, num_cols):
    return row * num_cols + col

def shuffle_indices(num_rows, num_cols):
    """
    For L-fixing. Carefully constructs an ordering 
    to put a large set B of indices at the end 
    such that none of them can serve as centers. 

    Returns new positions (pos[i] is the new position of unit i) 
    and the size of the set B. 
    """
    if num_rows != 20 or num_cols != 20:
        print("shuffle_indices is not yet implemented for general grid dimensions.")
        return np.arange(num_rows * num_cols), 0 # default ordering
    
    num_units = num_rows * num_cols

    pairs_in_B = []
    rows_in_B = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]
    cols_in_B = rows_in_B # By symmetry
    for i in rows_in_B:
        for j in cols_in_B:
            pairs_in_B.append((i, j))

    rows_not_in_B = [6, 13]
    cols_not_in_B = rows_not_in_B
    pairs_not_in_B = []
    for i in range(num_rows):
        for j in range(num_cols):
            if i in rows_not_in_B or j in cols_not_in_B:
                pairs_not_in_B.append((i, j))

    all_pairs_sorted = pairs_not_in_B + pairs_in_B

    pos = []
    
    for i in range(num_units):
        row, col = index_to_rowcol(i, num_cols)
        pos.append(all_pairs_sorted.index((row, col)))

    return pos, len(pairs_in_B)

def get_directed_grid_edges(num_rows, num_cols):
    """
    Return all ordered pairs of edges, 
    using indices starting at 0 in the top-left corner and 
    moving left-to-right, top-to-bottom, 
    in a grid graph with the given number of rows and columns. 
    """
    edges = [] 

    for i in range(num_rows):
        for j in range(num_cols):
            current = i * num_cols + j

            # Right edge
            right_neighbor = i * num_cols + (j + 1)
            if j + 1 < num_cols:
                edges.append((current, right_neighbor))

            # Down edge
            down_neighbor = (i + 1) * num_cols + j
            if i + 1 < num_rows:
                edges.append((current, down_neighbor))

    return edges


if __name__ == '__main__':
    try:
        # Create a new model
        m = gp.Model("hess_mip")

        num_rows = 6
        num_cols = 6
        num_districts = 2
        num_units = num_rows * num_cols
        max_units_per_district = num_units / num_districts

        # For now, assume population of each unit is 1
        p = np.ones(num_units)

        # Set up weights
        d2 = populate_squared_distance_matrix(num_rows, num_cols)
        print(d2)

        # Get directed edge set
        dir_edge_set = get_directed_grid_edges(num_rows, num_cols)
        num_dir_edges = len(dir_edge_set)

        # Create variables
        x = m.addVars(num_units, num_units, vtype=GRB.BINARY, name="x")
        f = m.addVars(num_units, num_dir_edges, vtype=GRB.CONTINUOUS, lb=0.0, name="f")

        # Set objective
        hess_obj = gp.quicksum((p[i] * d2[i, j]) * x[i, j]
                               for i in range(num_units) 
                               for j in range(num_units))
        m.setObjective(hess_obj, GRB.MINIMIZE)

        # Add constraint: each unit is assigned to a district
        for i in range(num_units):
            m.addConstr(gp.quicksum(x[i, j] for j in range(num_units)) == 1, f"assign_{i}")

        # Add constraint: exactly num_districts centers
        m.addConstr(gp.quicksum(x[j, j] for j in range(num_units)) == num_districts, f"count_centers")

        # Add constraint: population balance
        L = num_units / num_districts
        U = num_units / num_districts
        for j in range(num_units):
            m.addConstr(L * x[j, j] <= gp.quicksum(p[i] * x[i, j] for i in range(num_units)), f"pop_lower_bound_{j}")
            m.addConstr(U * x[j, j] >= gp.quicksum(p[i] * x[i, j] for i in range(num_units)), f"pop_upper_bound_{j}")

        # Add constraint: can only assign units to centers
        for i in range(num_units):
            for j in range(num_units):
                m.addConstr(x[i, j] <= x[j, j], f"assignonlyifcenter_{i}_{j}")

        # Diagonal fixing (TODO: improve ordering; see Section 6.1 and later)
        # for j in range(num_units):
        #     for i in range(j + 1, num_units):
        #         m.addConstr(x[i, j] == 0, f"fixzero_x[{i},{j}]")

        # L-fixing (specific to 20x20 grid into 10 districts, currently)
        pos, size_of_B = shuffle_indices(num_rows, num_cols)
        print(pos[:10], size_of_B)
        
        count_L_fixings = 0

        for i in range(num_units):
            for j in range(num_units):
                if pos[i] < pos[j]:
                    count_L_fixings += 1
                    m.addConstr(x[i, j] == 0, f"fixzero_x[{i},{j}]")

        for j in range(num_units):
            if pos[j] >= num_units - size_of_B:
                count_L_fixings += 1
                m.addConstr(x[j, j] == 0, f"fixnotcenter_{j}")

        print(f"L-fixing: fixed {count_L_fixings} variable(s) to zero.")


        # Add SHIR flow constraints for contiguity
        
        for j in range(num_units):
            edges_into_j_indices = []
            for dir_edge_index, dir_edge in enumerate(dir_edge_set):
                if dir_edge[1] == j:
                    edges_into_j_indices.append(dir_edge_index)

            # (2d) No incoming flow of type j for unit j
            m.addConstr(gp.quicksum(f[j, e] for e in edges_into_j_indices) == 0, f"noflowcirculation_{j}")

            for i in range(num_units):
                if i == j:
                    continue
                    
                edges_into_i_indices = []
                edges_out_i_indices = []
                for dir_edge_index, dir_edge in enumerate(dir_edge_set):
                    if dir_edge[1] == i:
                        edges_into_i_indices.append(dir_edge_index)
                    if dir_edge[0] == i:
                        edges_out_i_indices.append(dir_edge_index)
                
                # (2c) Unit can only receive flow from its assigned center
                m.addConstr(gp.quicksum(f[j, e] for e in edges_into_i_indices) 
                            <= max_units_per_district * x[i, j], f"flowonlyifassigned_({i},{j})")
                
                # (2b) If i assigned to j, then i receives exactly one unit of j flow
                m.addConstr(gp.quicksum(f[j, e_in] for e_in in edges_into_i_indices)
                            - gp.quicksum(f[j, e_out] for e_out in edges_out_i_indices)
                            == x[i, j], f"absorboneunitflow_({i},{j})")
                

        # Optimize model
        m.optimize()

        # for v in m.getVars():
        #     print('%s %g' % (v.varName, v.x))

        # Visualize results
        centers = []
        for j in range(num_units):
            if x[j, j].x >= 1 - TOL:
                centers.append(j)

        print(centers)
        print(sorted(centers))
        
        for a in range(num_rows):
            for b in range(num_cols):
                index = rowcol_to_index(a, b, num_cols)

                for j in range(num_units):
                    if x[index, j].x >= 1 - TOL:
                        district_number = 1 + centers.index(j)
                        if index in centers:
                            print(f"C{district_number:<2d}", end='')
                        else:
                            print(f" {district_number:<2d}", end='')    
            print()

        print('Obj: %g' % m.objVal)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')



