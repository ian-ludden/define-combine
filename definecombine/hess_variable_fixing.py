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

    TODO: Should also take parameter num_districts, or L (lower bound on district size), 
    since the maximum allowed size of a connected component in G[B]
    depends on the number of units per district. 

    Returns new positions (pos[i] is the new position of unit i) 
    and the size of the set B. 
    """
    supported_dimensions = [(20, 20), (10, 10), (8, 8)]

    if (num_rows, num_cols) not in supported_dimensions:
        print("shuffle_indices is not yet implemented for general grid dimensions.\nCurrently supported:", supported_dimensions)
        return np.arange(num_rows * num_cols), 0 # default ordering
    
    num_units = num_rows * num_cols

    pairs_in_B = []
    
    rows_in_B = [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19] if num_rows == 20 \
        else ([0, 1, 2, 3, 5, 6, 7, 8] if num_rows == 10 \
              else [0, 1, 2, 4, 5, 6]) # num_rows == 8
    cols_in_B = rows_in_B # By symmetry; would need to change this to support asymmetric dimensions
    for i in rows_in_B:
        for j in cols_in_B:
            pairs_in_B.append((i, j))

    rows_not_in_B = [6, 13] if num_rows == 20 else ([4, 9] if num_rows == 10 \
                                                    else [3, 7]) # num_rows == 8
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
    undir_edges = [] 

    for i in range(num_rows):
        for j in range(num_cols):
            current = i * num_cols + j

            # Right edge
            right_neighbor = i * num_cols + (j + 1)
            if j + 1 < num_cols:
                undir_edges.append((current, right_neighbor))

            # Down edge
            down_neighbor = (i + 1) * num_cols + j
            if i + 1 < num_rows:
                undir_edges.append((current, down_neighbor))

    dir_edges = []
    for u, v in undir_edges:
        dir_edges.append((u, v))
        dir_edges.append((v, u))

    return dir_edges


if __name__ == '__main__':
    try:
        # Create a new model
        m = gp.Model("hess_mip")

        num_rows = 8
        num_cols = 8
        num_districts = 4
        num_units = num_rows * num_cols
        max_units_per_district = num_units / num_districts

        # For now, assume population of each unit is 1
        p = np.ones(num_units)

        # Set up weights
        d2 = populate_squared_distance_matrix(num_rows, num_cols)
        print(d2)

        # Get directed edge set
        dir_edge_list = get_directed_grid_edges(num_rows, num_cols)
        num_dir_edges = len(dir_edge_list)
        # Get undirected edges (u, v) in canonical order (u < v)
        undir_edge_list = [(u, v) for u, v in dir_edge_list if u < v]
        num_undir_edges = len(undir_edge_list)
        print("Directed edges:", num_dir_edges, "\nUndirected edges:", num_undir_edges)

        # Create variables
        ## Unit-to-unit assignment variables
        x = m.addVars(num_units, num_units, vtype=GRB.BINARY, name="x")
        ## Center-to-unit flow variables
        f = m.addVars(num_units, num_dir_edges, vtype=GRB.CONTINUOUS, lb=0.0, name="f")
        ## Cut-edge indicator variables
        ## TODO: Change to undirected edges to reduce number of variables
        y = m.addVars(num_dir_edges, vtype=GRB.BINARY, name="y") 
        ## Cut-edge indicator variables for extended formulation
        z = m.addVars(num_undir_edges, num_units, vtype=GRB.CONTINUOUS, lb=0.0, name="z")

        # Set objective
        ## Hess objective: sum of districts' moments of inertia
        hess_obj = gp.quicksum((p[i] * d2[i, j]) * x[i, j]
                               for i in range(num_units) 
                               for j in range(num_units))
        
        ## Min-cut objective: how many edges are cut (see Validi and Buchanan 2022, (2a)-(2b))
        min_cut_obj = gp.quicksum(y[e] for e in range(num_dir_edges))

        # m.setObjective(hess_obj, GRB.MINIMIZE)
        m.setObjective(min_cut_obj, GRB.MINIMIZE)

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

        print(f"L-fixing: fixed {count_L_fixings} variable(s) to zero with distinct i \\neq j.")

        for j in range(num_units):
            if pos[j] >= num_units - size_of_B:
                count_L_fixings += 1
                m.addConstr(x[j, j] == 0, f"fixnotcenter_{j}")

        print(f"L-fixing: fixed {count_L_fixings} variable(s) to zero total.")

        # Add SHIR flow constraints for contiguity
        for j in range(num_units):
            edges_into_j_indices = []
            for dir_edge_index, dir_edge in enumerate(dir_edge_list):
                if dir_edge[1] == j:
                    edges_into_j_indices.append(dir_edge_index)

            # (2d) No incoming flow of type j for unit j
            m.addConstr(gp.quicksum(f[j, e] for e in edges_into_j_indices) == 0, f"noflowcirculation_{j}")

            for i in range(num_units):
                if i == j:
                    continue
                    
                edges_into_i_indices = []
                edges_out_i_indices = []
                for dir_edge_index, dir_edge in enumerate(dir_edge_list): # TODO: This part is very slow, need to speed it up (use dictionaries?)
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
                

        # TODO: Delete, old approach that is replaced by z-variable-based approach 
        # Add constraints for cut-edge indicator variables
        # for j in range(num_units):
        #     for e_index in range(num_dir_edges):
        #         e = dir_edge_list[e_index]
        #         m.addConstr(x[e[0], j] - x[e[1], j] <= y[e_index], f"cutedge_y_({e[0]}, {e[1]})_j{j}_e{e_index}")


        # Add constraints for extended formulation cut-edge indicator variables (z)
        for j in range(num_units):
            for e_index in range(num_undir_edges):
                u, v = undir_edge_list[e_index]
                m.addConstr(x[u, j] - x[v, j] <= z[e_index, j], f"cutedge_z_({u}, {v})_j{j}_e{e_index}")


        # Add constraints relating y and z
        for e_index in range(num_dir_edges):
            e = dir_edge_list[e_index]
            
            if e in undir_edge_list: # i.e., u < v
                undir_index = undir_edge_list.index(e)
                m.addConstr(y[e_index] == gp.quicksum([z[undir_index, j] for j in range(num_units)]), f"z{undir_index}_to_y{e_index}")
            else:
                m.addConstr(y[e_index] == 0, f"fixyzero_{e_index}")


        # Z-fixing: if x_uj is fixed to zero, or x_vj is fixed to 1, then we can fix z_{uv}^{j} to zero
        m.update() # To make sure constraint names have been indexed
        for e_index in range(num_undir_edges):
            u, v = undir_edge_list[e_index]
            for j in range(num_units):
                try:
                    c = m.getConstrByName(f"fixzero_x[{u},{j}]")
                    if c:
                        m.addConstr(z[(u,v), j] == 0.0, f"fixzero_z[({u},{v})_{j}]")
                except Exception as e:
                    pass # constraint doesn't exist, can't z-fix


        # Warm-start: pick district centers
        if (num_rows, num_cols) == (8, 8):
            x[3, 3].Start = 1.0
            x[7, 7].Start = 1.0
            x[35, 35].Start = 1.0
            x[39, 39].Start = 1.0
        if (num_rows, num_cols) == (10, 10):
            x[4, 4].Start = 1.0
            x[9, 9].Start = 1.0
            x[54, 54].Start = 1.0
            x[59, 59].Start = 1.0
        if (num_rows, num_cols) == (20, 20):
            # x[6, 6].Start = 1.0
            # x[13, 13].Start = 1.0
            # x[120, 120].Start = 1.0
            # x[135, 135].Start = 1.0
            # x[166, 166].Start = 1.0
            # x[173, 173].Start = 1.0
            # x[260, 260].Start = 1.0
            # x[275, 275].Start = 1.0
            # x[326, 326].Start = 1.0
            # x[366, 366].Start = 1.0
            for r in range(16):
                for c in range(20):
                    i = r * 20 + c
                    if c <= 4 and r <= 7:
                        x[i, 120].Start = 1.0
                    elif c <= 9 and r <= 7:
                        x[i, 6].Start = 1.0
                    elif c <= 14 and r <= 7:
                        x[i, 13].Start = 1.0
                    elif r <= 7:
                        x[i, 135].Start = 1.0
                    elif c <= 4 and r <= 15:
                        x[i, 260].Start = 1.0
                    elif c <= 9 and r <= 15:
                        x[i, 166].Start = 1.0
                    elif c <= 14 and r <= 15:
                        x[i, 173].Start = 1.0
                    elif r <= 15:
                        x[i, 275].Start = 1.0
                    elif r <= 17:
                        x[i, 326].Start = 1.0
                    else:
                        x[i, 366].Start = 1.0




        # Optimize model
        m.optimize()

        # Check infeasibility
        if m.Status == GRB.INFEASIBLE:
            m.computeIIS()
            m.write('iismodel.ilp')

            # Print out the IIS constraints and variables
            print('\nThe following constraints and variables are in the IIS:')
            for c in m.getConstrs():
                if c.IISConstr: print(f'\t{c.constrname}: {m.getRow(c)} {c.Sense} {c.RHS}')

            for v in m.getVars():
                if v.IISLB: print(f'\t{v.varname} ≥ {v.LB}')
                if v.IISUB: print(f'\t{v.varname} ≤ {v.UB}')

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



