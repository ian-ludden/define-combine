import csv
import networkx as nx
import numpy as np
import sys
from tqdm import tqdm

DEBUG = False

"""
Solve DCP on a grid graph by exhaustive search 
over all possible definer strategies (subdistrict partitions) 
and finding the best combiner response (perfect matching) for each. 

Parameters:
===============
    num_rows - number of rows in grid graph
    num_cols - number of columns in grid graph
    num_districts - number of districts
    voter_grid - 2-D NumPy array with shape (num_rows, num_cols)
                 with entry [i,j] == 1 (0) if 
                 unit (i,j) is for definer (combiner)

Returns:
===============
    best_assignment_grid - best subdistrict partition (optimal definer strategy)
    best_definer_utility - maximum definer utility under optimal play
    best_combiner_matching - combiner's best response to the subdistrict partition 
                             represented by best_assignment_grid
    best_subdistrict_votes - combiner vote counts for each subdistrict in best_assignment_grid
"""
def solve_dcp_grid(num_rows, num_cols, num_districts, voter_grid):
    num_subdistricts = 2 * num_districts
    units_per_subdistrict = (num_rows * num_cols) // num_subdistricts
    enumeration_filename = f"enumerations\\enum_[{num_cols},{num_rows}]_[{units_per_subdistrict}]_{num_subdistricts}_rc.txt"
    
    # Keep track of best definer map so far
    best_assignment_grid = None
    best_definer_utility = -1 * num_districts
    best_combiner_matching = []
    best_subdistrict_votes = None
    num_partitions = 0

    with open(enumeration_filename, 'r') as enum_file:
        csv_enum_file = csv.reader(enum_file)
        for line in tqdm(csv_enum_file, desc="Trying all subdistrict partitions"):
            num_partitions += 1
            assignment_grid = np.array(line, dtype=int).reshape(num_rows, num_cols)

            # Tally combiner votes in each subdistrict
            subdistrict_votes = np.zeros((num_subdistricts,), dtype=int)
            for i in range(num_rows):
                for j in range(num_cols):
                    subdistrict_index = assignment_grid[i, j] - 1
                    if voter_grid[i, j] <= 0:
                        subdistrict_votes[subdistrict_index] += 1

            # Compute edges
            edge_set = set()
            for i in range(num_rows):
                for j in range(num_cols):
                    subdistrict = assignment_grid[i, j]
                    if i < num_rows - 1:
                        subdistrict_below = assignment_grid[i + 1, j]
                        edge_set.add(tuple(sorted([subdistrict, subdistrict_below])))
                    if j < num_cols - 1:
                        subdistrict_right = assignment_grid[i, j + 1]
                        edge_set.add(tuple(sorted([subdistrict, subdistrict_right])))
            
            # Remove self-loops
            for subdistrict in range(1, 2 * num_districts):
                edge_set.remove((subdistrict, subdistrict))

            edge_list = list(edge_set)

            # Use networkx to solve the combiner's max-weight perfect matching problem
            G = nx.Graph(edge_list)
            for edge in G.edges:
                u, v = edge
                edge_weight = 0
                dist_sum = subdistrict_votes[u - 1] + subdistrict_votes[v - 1]
                if dist_sum > units_per_subdistrict:
                    edge_weight = 1
                elif dist_sum < units_per_subdistrict:
                    edge_weight = -1
                else:
                    pass
                G[u][v]['weight'] = edge_weight

            combiner_matching = nx.max_weight_matching(G, maxcardinality=True)
            definer_utility = - 1 * sum(G[u][v]['weight'] for u, v in combiner_matching)

            if definer_utility > best_definer_utility:
                best_definer_utility = definer_utility
                best_assignment_grid = assignment_grid
                best_combiner_matching = combiner_matching
                best_subdistrict_votes = subdistrict_votes

    if DEBUG: print("Examined all", num_partitions, "available definer strategies.\n")

    return best_assignment_grid, best_definer_utility, best_combiner_matching, best_subdistrict_votes
    



if __name__ == '__main__':
    if len(sys.argv) < 5:
        raise Exception("Missing argument(s). Usage: python enum_dcp.py [num_rows] [num_cols] [num_districts] [voter_distribution_filename] [first_player (optional, \"R\" or \"D\", default is \"D\")].")
        
    num_rows = int(sys.argv[1])
    num_cols = int(sys.argv[2])
    num_districts = int(sys.argv[3])
    voter_distr_fname = str(sys.argv[4])

    voter_grid = np.genfromtxt(voter_distr_fname, dtype=int, delimiter=1)

    first_player = "D"
    if len(sys.argv) >= 6:
        first_player = str(sys.argv[5])
    
    if first_player == "R": # Invert 1s and 0s in voter distribution
        voter_grid = 1 - voter_grid

    best_definer_strategy, \
        best_definer_utility, \
        best_combiner_response, \
        subdistrict_votes = \
            solve_dcp_grid(num_rows=num_rows, num_cols=num_cols, \
            num_districts=num_districts, voter_grid=voter_grid)

    definer_wins_for_max_utility = 0
    for edge in best_combiner_response:
        combiner_votes = subdistrict_votes[edge[0] - 1] + subdistrict_votes[edge[1] - 1]
        definer_votes = (num_rows * num_cols) // num_districts - combiner_votes
        definer_wins_for_max_utility += (definer_votes > combiner_votes)
    
    
    # Print CSV summary of instance with results
    print(first_player, num_rows, num_cols, num_districts, voter_distr_fname, best_definer_utility, definer_wins_for_max_utility, sep=',')

    if DEBUG:
        print("Maximum definer utility:", best_definer_utility, "from the following subdistrict plan:\n")
        print(best_definer_strategy)
        print("with combiner response (perfect matching):\n")
        for edge in best_combiner_response:
            combiner_votes = subdistrict_votes[edge[0] - 1] + subdistrict_votes[edge[1] - 1]
            definer_votes = (num_rows * num_cols) // num_districts - combiner_votes
            print("{}\t{} ({}) votes for definer (combiner)".format(edge, int(definer_votes), int(combiner_votes)))
