from gerrychain import Partition
import networkx as nx

def aux_graph(partition):
    parent = partition.parent

    if not parent:
        cut_edges_set = partition.cut_edges
        
        G = nx.Graph()