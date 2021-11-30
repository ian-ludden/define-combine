from functools import partial
import time
import networkx as nx
import numpy as np
from random import randint

from gerrychain import accept, constraints
from gerrychain.chain import MarkovChain
from gerrychain.grid import Grid
from gerrychain.proposals import recom
from gerrychain.updaters import cut_edges, Tally

from utils.custom_updaters import aux_graph, definer_utility

num_parts = 10 # also referred to as N
num_cols = num_parts * 8
num_rows = num_parts * 4

votes_per_unit = 1
population_per_unit = 1

unit_ids = [(i, j) for i in range(num_cols) for j in range(num_rows)]
# Create 2N subparts as two rows of N
init_assignment_dict = {(i, j): (j // (num_rows // 2))*num_parts + (i // (num_cols // num_parts)) + 1 for i,j in unit_ids}

grid = Grid(
    dimensions=(num_cols, num_rows), 
    assignment=init_assignment_dict, 
    updaters={
        "cut_edges": cut_edges, 
        "votes_definer": Tally("votes_definer"), 
        "votes_total": Tally("votes_total"), 
        "population": Tally("population"), 
        "aux_graph": aux_graph, 
        "definer_utility": definer_utility}
    )

nx.set_node_attributes(grid.graph, init_assignment_dict, 'init_assignment')

nx.set_node_attributes(grid.graph, votes_per_unit, 'votes_total')
nx.set_node_attributes(grid.graph, population_per_unit, 'population')

if votes_per_unit > 3:
    votes_definer_dict = {node: randint(2, votes_per_unit - 2) for node in grid.graph.nodes}
else:
    votes_definer_dict = {node: randint(0, 1) for node in grid.graph.nodes}

nx.set_node_attributes(grid.graph, votes_definer_dict, 'votes_definer')

if num_parts < 5: # Otherwise, subparts with indices 10+ mess up the string representation
    print(grid)

print(grid.graph)
print("Number of cut edges:", len(grid.cut_edges))
print("Auxiliary graph edges:", len(grid.aux_graph.edges), grid.aux_graph.edges,'\n')

print("\nDefiner votes per part:")
print(grid["votes_definer"])

print("\nDefiner votes per unit:\n")
for j in range(num_rows):
    for i in range(num_cols):    
        print(grid.graph.nodes[(i,j)]["votes_definer"], end="")
    print()
print("\n")

ideal_population = sum(grid["population"].values()) / len(grid)

proposal = partial(recom, 
                   pop_col="population", 
                   pop_target=ideal_population,
                   epsilon=0.02,
                   node_repeats=2
                  )

compactness_bound = constraints.UpperBound(
    lambda p: len(p["cut_edges"]), 
    2 * len(grid["cut_edges"])
)

pop_constraint=constraints.within_percent_of_ideal_population(grid, 0.02)

chain = MarkovChain(
    proposal=proposal, 
    constraints=[
        pop_constraint, 
        compactness_bound
    ], 
    accept=accept.always_accept, 
    initial_state=grid, 
    total_steps=1000
)

H = grid.aux_graph
oldH = grid.aux_graph
index = 0
def_util_arr = np.zeros((chain.total_steps,))

# Set up timer
tic = time.perf_counter()

for partition in chain:
    index += 1
    if (index == 1) or (index % 100 == 0): print("*** Partition", index, "***")
    def_util_arr[index - 1] = partition["definer_utility"]

toc = time.perf_counter()
elapsed_time = toc - tic
average_time = elapsed_time / chain.total_steps
print(f"The ReCom chain ran in {elapsed_time:0.4f} seconds total, for an average of {average_time:0.4f} seconds per step.")

if num_parts < 5:
    print(partition)
print(partition.graph)
print("Number of cut edges:", len(partition.cut_edges))
print("Auxiliary graph edges:", len(H.edges), H.edges,'\n')

argmax_def_util = np.argmax(def_util_arr)
max_def_util = def_util_arr[argmax_def_util]

print(f"The best possible definer utility is {max_def_util}, achieved by map #{argmax_def_util + 1}.")
print("\nDefiner utilities:\n", def_util_arr)
