# Solve define-combine greedy heuristic: maximize subdistrict wins
from cgi import print_arguments
import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import numpy as np
import random
import re

NUM_ROWS = 4
NUM_COLS = 6
NUM_DISTRICTS = 4
NUM_SUBDISTRICTS = 2 * NUM_DISTRICTS

NUM_UNITS = NUM_ROWS * NUM_COLS
# UNITS_PER_DISTRICT = NUM_UNITS // NUM_DISTRICTS
UNITS_PER_SUBDISTRICT = NUM_UNITS // NUM_SUBDISTRICTS

NUM_DEFINER_SUPPORTERS = int(NUM_UNITS * 0.50)
NUM_OTHER_SUPPORTERS = NUM_UNITS - NUM_DEFINER_SUPPORTERS

BIGM = UNITS_PER_SUBDISTRICT + 1 # at least (max vote-delta per unit) X (max units per subdistrict)
EPSILON = 0.5 # smaller than resolution of vote counts


"""
Find the optimal combiner response (max-weight perfect matching) 
given the definer's subdistrict partition. 

Parameters:
===============
G - networkx.Graph representing unit adjacency graph
subdist_assignment_dict - dictionary where keys are unit indices, 
                          values are assigned subdistrict indices
subdist_partition_dict - dictionary where keys are subdistrict indices, 
                          values are sets of units in each subdistrict

Returns:
===============
M - the combiner's best response, 
    a matching on the subdistrict indices 
    to form districts
"""
def optimal_combiner_response(G, subdist_assignment_dict, subdist_partition_dict): 
    node_list = list(G.nodes)
    grid_names_to_unit_ids = {node_list[i]: i for i in unit_ids}
    subdist_edges = set()

    for edge in G.edges:
        p = grid_names_to_unit_ids[edge[0]]
        q = grid_names_to_unit_ids[edge[1]]
        u = min(subdist_assignment_dict[p], subdist_assignment_dict[q])
        v = max(subdist_assignment_dict[p], subdist_assignment_dict[q])
        if u < v:
            subdist_edges.add((u, v))
    
    subdistG = nx.Graph(subdist_edges)

    edge_weights = []
    for u, v in subdist_edges:
        net_votes = sum([G.nodes[node_list[i]]['net_votes'] for i in subdist_partition_dict[u]])\
                  + sum([G.nodes[node_list[i]]['net_votes'] for i in subdist_partition_dict[v]])
        edge_weight = None
        if net_votes < 0:
            edge_weight = 1 # Flip to view from combiner's perspective
        elif net_votes > 0:
            edge_weight = -1
        else:
            edge_weight = 0
        edge_weights.append(edge_weight)

    edge_list = list(subdistG.edges)
    nx.set_edge_attributes(subdistG, {edge_list[i]: edge_weights[i] for i in range(len(edge_list))}, name='weight')
    return nx.max_weight_matching(subdistG, maxcardinality=True)


def build_model(G):
    try:      
        # Create MIP model
        m = gp.Model("msip")

        unit_ids = np.arange(len(G.nodes))
        node_list = list(G.nodes)
        grid_names_to_unit_ids = {node_list[i]: i for i in unit_ids}

        ### VARIABLES ###
        # Create variables x_{ij}, indicators for whether 
        # unit i is assigned to (the subdistrict centered at) unit j
        x = m.addVars(unit_ids, unit_ids, name="x", vtype=GRB.BINARY)

        # Create variables w_{v} and l_{v}, indicators for whether 
        # the subdistrict centered at v 
        # is a win for party A (resp., party B)
        w = m.addVars(unit_ids, name="w", vtype=GRB.BINARY)
        ell = m.addVars(unit_ids, name="ell", vtype=GRB.BINARY)

        edges_in = []
        edges_out = []
        for u in unit_ids:
            edges_in_u = []
            edges_out_u = []
            for x_edge in G.edges:
                p = grid_names_to_unit_ids[x_edge[0]]
                q = grid_names_to_unit_ids[x_edge[1]]
                if q == u:
                    edges_in_u.append((p, q))
                if p == u:
                    edges_out_u.append((p, q))
            
            edges_in.append(edges_in_u)
            edges_out.append(edges_out_u)
        
        all_directed_edges = []
        for edge in G.edges:
            p = grid_names_to_unit_ids[edge[0]]
            q = grid_names_to_unit_ids[edge[1]]
            all_directed_edges.extend([(p, q), (q, p)])

        # Create flow variables for SHIR contiguity
        f = m.addVars(unit_ids, all_directed_edges, name="f", vtype=GRB.CONTINUOUS) # use default lower bound of 0


        ### CONSTRAINTS ###
        # Add constraint: force all assignments in x to go from lower units to higher units
        # (breaks some symmetry)
        for i in unit_ids:
            for j in unit_ids:
                if j < i:
                    m.addConstr(x[i, j] == 0, "forced_zero_x_{}_{}".format(i, j))
        

        # Add constraint: sum of x_{jj} = NUM_SUBDISTRICTS
        m.addConstr(sum(x[i, i] for i in unit_ids) == NUM_SUBDISTRICTS, "num_subdistricts")


        # # Add constraints: every unit fully assigned
        m.addConstrs(
            (x.sum(i, '*') == 1 
                for i in unit_ids), "fullyassigned")


        # Add constraints: every subdistrict has exactly NUM_UNITS / NUM_SUBDISTRICTS units (perfect population balance)
        for j in unit_ids:
            m.addConstr(sum(x[i, j] for i in unit_ids) == UNITS_PER_SUBDISTRICT * x[j, j], "subdistrict size %s" % (j))


        # # Add constraints: every district has exactly NUM_UNITS / NUM_DISTRICTS units (perfect population balance)
        # for j in unit_ids:
        #     m.addConstr(sum(y[i, j] for i in unit_ids) == UNITS_PER_DISTRICT * y[j, j], "district size %s" % (j))


        # Add constraints: can only assign to centers
        for i in unit_ids:
            for j in unit_ids:
                if i == j: continue; # Only need this constraint for distinct i, j
                m.addConstr(x[i, j] - x[j, j] <= 0, "onlycenters[%s, %s]" % (i, j))

        # Add constraints: 
        # w_{v} is 1 if and only if definer wins or ties 
        # subdistrict centered at v;
        # ell_{v} is 1 if and only if definer loses or ties
        # subdistrict centered at v.
        for v in unit_ids:
            m.addConstr(sum(x[i, v] * unit_weights[i] 
                for i in unit_ids) <= BIGM * w[v] - EPSILON, 
                "w def1 %s" % (v))
        
            m.addConstr(sum(x[i, v] * unit_weights[i]
                for i in unit_ids) >= - BIGM * (1 - w[v]), 
                "w def2 %s" % (v))

            m.addConstr(sum(x[i, v] * unit_weights[i]
                for i in unit_ids) >= EPSILON - BIGM * ell[v], 
                "ell def1 %s" % (v))

            m.addConstr(sum(x[i, v] * unit_weights[i]
                for i in unit_ids) <= BIGM * (1 - ell[v]), 
                "ell def2 %s" % (v))

        
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


        # Set objective: maximize sum a_uv minus b_uv
        m.setObjective(sum(w[v] - ell[v] for v in unit_ids), sense=GRB.MAXIMIZE)
        
        return m

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError as e:
        print('Encountered an attribute error.')



def summarize_model_results(m, G):
    print('Obj: %g' % m.objVal)

    nodelist = list(G.nodes)

    # Check which units are subdistrict centers, and extract subdistrict partition
    partition_dict = {}
    assignment_dict = {}
    tol = m.Params.IntFeasTol

    x_vals = np.zeros((NUM_UNITS, NUM_UNITS))

    for v in m.getVars():
        var_name_parts = re.split('\[|,|\]', v.varName)
        if not (var_name_parts[0] == 'x'): continue
        i = int(var_name_parts[-3])
        j = int(var_name_parts[-2])

        x_vals[i, j] = v.x

    # Assign each unit to the unit with max weight in x
    assignments = np.argmax(x_vals, axis=1)
    print(assignments)
    print(np.unique(assignments))

    assignment_dict = {unit_id: assignments[unit_id] for unit_id in unit_ids}
    for subdistrict_center in np.unique(assignments):
        units_in_subdistrict = set()
        for unit_id in unit_ids:
            if assignment_dict[unit_id] == subdistrict_center:
                units_in_subdistrict.add(unit_id)
        partition_dict[subdistrict_center] = units_in_subdistrict

        net_votes = sum([G.nodes[nodelist[i]]['net_votes'] for i in partition_dict[subdistrict_center]])
        print('Subdistrict centered at', subdistrict_center, 'has', net_votes, 'net votes for definer.')
    print()

    # Determine best response from combiner 
    # by building and solving perfect matching instance
    M = optimal_combiner_response(G, assignment_dict, partition_dict)
    subdist_to_dist = {}
    for edge in M:
        subdist_to_dist[edge[0]] = min(edge)
        subdist_to_dist[edge[1]] = min(edge)

    y_assignment_dict = {}
    for unit_id in unit_ids:
        subdist_index = assignment_dict[unit_id]
        dist_index = subdist_to_dist[subdist_index]
        y_assignment_dict[unit_id] = dist_index

    y_partition_dict = {}
    for district_center in np.unique(list(y_assignment_dict.values())):
        units_in_district = set()
        for unit_id in unit_ids:
            if y_assignment_dict[unit_id] == district_center:
                units_in_district.add(unit_id)
        y_partition_dict[district_center] = units_in_district

    count_a_wins = 0
    count_b_wins = 0
    for key in y_partition_dict:
        district_weight = 0
        for unit_id in y_partition_dict[key]:
            district_weight += G.nodes[nodelist[unit_id]]['net_votes']
        print('District centered at unit', key, 'has total weight', district_weight)
        if district_weight > 0:
            count_a_wins += 1
        if district_weight < 0:
            count_b_wins += 1
        
    print('Wins for definer:\t', count_a_wins)
    print('Wins for combiner:\t', count_b_wins)
    print('Net for definer:\t', count_a_wins - count_b_wins)

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



if __name__ == '__main__':
    random.seed(2022) # For reproducibility
    assert(NUM_UNITS // NUM_SUBDISTRICTS == NUM_UNITS / NUM_SUBDISTRICTS) # Require exact population balance

    print("Definer has", NUM_DEFINER_SUPPORTERS, "supporters.")

    unit_weights = [-1] * NUM_OTHER_SUPPORTERS + [1] * NUM_DEFINER_SUPPORTERS
    random.shuffle(unit_weights)
    unit_ids = [i for i in range(NUM_UNITS)]

    # Build unit adjacency graph
    G = nx.grid_graph(dim=(NUM_COLS, NUM_ROWS))
    print(G)

    nodelist = list(G.nodes)
    nx.set_node_attributes(G, {nodelist[i]: unit_weights[i] for i in range(len(G.nodes))}, 'net_votes')

    # for i in G.nodes:
    #     print(i, G.nodes[i]['net_votes'])


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
    
    try:
        # Build model
        m = build_model(G)

        # Optimize model
        m.optimize()
        summarize_model_results(m, G)

        # print("\n*** RESET, RELAX, and RESOLVE ***\n")
        # m.reset()
        # m = m.relax()
        # # m.addConstr(m.getVarByName("x[0,6]") == 0, "forced_zero_x_{}_{}".format(0, 6))
        
        # m.optimize()
        # summarize_model_results(m)

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError as e:
        print('Encountered an attribute error.')
