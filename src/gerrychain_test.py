# Test Gerrychain packages
import geopandas
from gerrychain import Partition, Graph
from gerrychain.updaters import cut_edges

graph = Graph.from_json("./data/PA_VTDs.json")
partition = Partition(graph, "CD_2011", {"cut_edges": cut_edges})

print(partition.graph)
