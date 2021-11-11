from functools import partial
import networkx as nx
from random import randint

from gerrychain import accept, constraints
from gerrychain.chain import MarkovChain
from gerrychain.grid import Grid
from gerrychain.proposals import recom
from gerrychain.updaters import cut_edges, Tally

from utils.aux_graph_updater import aux_graph

num_cols = 20
num_rows = 10
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
        "aux_graph": aux_graph}
    )

nx.set_node_attributes(grid.graph, init_assignment_dict, 'init_assignment')

nx.set_node_attributes(grid.graph, votes_per_unit, 'votes_total')
nx.set_node_attributes(grid.graph, population_per_unit, 'population')

votes_definer_dict = {node: randint(1, votes_per_unit) for node in grid.graph.nodes}
nx.set_node_attributes(grid.graph, votes_definer_dict, 'votes_definer')

print(grid)

print(grid.graph)
print("Number of cut edges:", len(grid.cut_edges))
print("Auxiliary graph edges:", len(grid.aux_graph.edges), grid.aux_graph.edges,'\n')

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
for partition in chain:
    index += 1
    print("*** Partition", index, "***")
    H = partition.aux_graph
    delta_edges = (H.edges - oldH.edges) | (oldH.edges - H.edges)
    oldH = H.copy()
    print(delta_edges, "edge(s) changed.")
    print()

print(partition)
print(partition.graph)
print("Number of cut edges:", len(partition.cut_edges))
print("Auxiliary graph edges:", len(H.edges), H.edges,'\n')