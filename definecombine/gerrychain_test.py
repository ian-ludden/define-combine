# Test Gerrychain packages
from gerrychain import Partition, Graph
from gerrychain.grid import Grid
from gerrychain.updaters import cut_edges, Tally
from random import randint

from .utils import aux_graph

num_cols = 10
num_rows = 10
votes_per_unit = 10

graph = Graph()
for i in range(num_rows):
    for j in range(num_cols):
        unit_id = i * num_cols + j
        graph.add_node(
            unit_id, 
            init_assignment=(i // 5) * 2 + (j // 5) + 1, 
            votes_total=votes_per_unit, 
            votes_definer=randint(1, votes_per_unit)
            )

for node in graph.nodes:
    print(node)

partition = Partition(
    graph=graph, 
    assignment="init_assignment", 
    updaters={"cut_edges": cut_edges, 
              "votes_definer": Tally("votes_definer"), 
              "votes_total": Tally("votes_total"), 
              "aux_graph": aux_graph}
    )

print(partition)
        
print(partition.graph)

for node in partition.graph:
    print("Node", node, "is in part", partition.assignment[node])

print("Number of cut edges:", len(partition.cut_edges))
print("Auxiliary graph edges:", partition.aux_graph.edges)
