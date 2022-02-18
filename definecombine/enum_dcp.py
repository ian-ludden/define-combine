import networkx as nx
import numpy as np
import sys

# TODO: Enumerative approach to solving DCP on grid instances, 
# similar to optimal recursive bisection on grid graphs

def solve_dcp_grid(num_rows, num_cols, num_districts, voter_grid):
    num_subdistricts = 2 * num_districts
    units_per_subdistrict = (num_rows * num_cols) // num_subdistricts
    enumeration_filename = "../enumerations/enum_[{num_cols},{num_rows}]_[{units_per_subdistrict}]_{num_subdistricts}_rc.txt"
    
    with open(enumeration_filename, 'r') as enum_file:
        pass
    
    pass



if __name__ == '__main__':
    if len(sys.argv) < 5:
        raise Exception("Missing argument(s). Usage: python enum_dcp.py [num_rows] [num_cols] [num_districts] [voter_distribution_filename].")
        
    num_rows = int(sys.argv[1])
    num_cols = int(sys.argv[2])
    num_districts = int(sys.argv[3])
    voter_distr_fname = str(sys.argv[4])

    voter_grid = np.genfromtxt(voter_distr_fname, dtype=int, delimiter=1)

    solve_dcp_grid(num_rows=num_rows, num_cols=num_cols, \
                   num_districts=num_districts, voter_grid=voter_grid)
