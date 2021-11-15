from functools import partial
import networkx as nx
import numpy as np
from random import randint

from gerrychain import accept, constraints
from gerrychain.chain import MarkovChain
from gerrychain.grid import Grid
from gerrychain.proposals import recom
from gerrychain.updaters import cut_edges, Tally

from utils.custom_updaters import aux_graph, definer_utility

num_cols = 40
num_rows = 20
votes_per_unit = 10
population_per_unit = 10

unit_ids = [(i, j) for i in range(num_cols) for j in range(num_rows)]
init_assignment_dict = {(i, j): (j // (num_rows // 2))*4 + (i // (num_cols // 4)) + 1 for i,j in unit_ids}

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

votes_definer_dict = {node: randint(2, votes_per_unit - 2) for node in grid.graph.nodes}
nx.set_node_attributes(grid.graph, votes_definer_dict, 'votes_definer')

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
    total_steps=10
)

H = grid.aux_graph
oldH = grid.aux_graph
index = 0
def_util_arr = np.zeros((chain.total_steps,))
for partition in chain:
    index += 1
    if (index == 1) or (index % 100 == 0): print("*** Partition", index, "***")
    def_util_arr[index - 1] = partition["definer_utility"]

print(partition)
print(partition.graph)
print("Number of cut edges:", len(partition.cut_edges))
print("Auxiliary graph edges:", len(H.edges), H.edges,'\n')
print("Definer utilities:\n", def_util_arr)

argmax_def_util = np.argmax(def_util_arr)
max_def_util = def_util_arr[argmax_def_util]

print(f"The best possible definer utility is {max_def_util}, achieved by map #{argmax_def_util + 1}")
