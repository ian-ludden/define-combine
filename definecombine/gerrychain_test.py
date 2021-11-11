# Test Gerrychain packages
from gerrychain.grid import Grid
from gerrychain.updaters import cut_edges, Tally
import networkx as nx
from random import randint

from utils.aux_graph_updater import aux_graph

num_cols = 10
num_rows = 10
votes_per_unit = 10

init_assignment_values = [i for i in range(1, num_cols + 1)] * num_rows
unit_ids = [(i, j) for i in range(num_rows) for j in range(num_cols)]
init_assignment_dict = {unit_ids[i]: init_assignment_values[i] for i in range(num_rows * num_cols)}

grid = Grid(
    dimensions=(num_rows, num_cols), 
    assignment=init_assignment_dict, 
    updaters={
        "cut_edges": cut_edges, 
        "votes_definer": Tally("votes_definer"), 
        "votes_total": Tally("votes_total"), 
        "aux_graph": aux_graph}
    )

nx.set_node_attributes(grid.graph, init_assignment_dict, 'init_assignment')

nx.set_node_attributes(grid.graph, votes_per_unit, 'votes_total')

votes_definer_dict = {node: randint(1, votes_per_unit) for node in grid.graph.nodes}
nx.set_node_attributes(grid.graph, votes_definer_dict, 'votes_definer')

# for unit_id in grid.graph.nodes:
#         print(unit_id, "is in part", grid.assignment[unit_id])
#         print(grid.graph.nodes[unit_id])

print(grid.graph)
print("Number of cut edges:", len(grid.cut_edges))
print("Auxiliary graph edges:", grid.aux_graph.edges)