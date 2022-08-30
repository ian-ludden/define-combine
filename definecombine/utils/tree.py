import random
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def uniform_random_spanning_tree(graph, choice=random.choice):
    """ Builds a spanning tree chosen uniformly from the space of all
        spanning trees of the graph.

        This implementation is from 
        GerryChain's tree.py:
        https://github.com/mggg/GerryChain/blob/main/gerrychain/tree.py

        :param graph: Networkx Graph
        :param choice: :func:`random.choice`
    """
    node_indices = np.arange(graph.number_of_nodes())
    node_names = list(graph.nodes)
    index_to_name = {i : node_names[i] for i in node_indices}
    name_to_index = {node_names[i] : i for i in node_indices}
    root = choice(list(node_indices))
    tree_nodes = set([root])
    next_node = {root: None}

    for node in node_indices:
        u = node
        while u not in tree_nodes:
            next_node[u] = name_to_index[choice(list(graph.neighbors(index_to_name[u])))]
            u = next_node[u]

        u = node
        while u not in tree_nodes:
            tree_nodes.add(u)
            u = next_node[u]

    G = nx.Graph()
    for node in tree_nodes:
        if next_node[node] is not None:
            G.add_edge(node, next_node[node])

    return G


# TODO: Implement Mulvey-Davis algorithm for 
# determining whether a tree is 2-splittable 
# within an error epsilon


if __name__ == '__main__':
    print("NetworkX version:", nx.__version__)
    
    graph = nx.grid_2d_graph(m=16, n=16)
    node_names = list(graph.nodes)

    nx.draw_networkx(graph, 
        pos={node: node for node in graph.nodes}, 
        with_labels=False, 
        node_size=10)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    ITERATIONS = 10

    for i in range(ITERATIONS):
        span_tree = uniform_random_spanning_tree(graph)
        assert nx.is_tree(span_tree)

        nx.draw_networkx(span_tree, 
            pos={node: node_names[node] 
                for node in span_tree.nodes}, 
            with_labels=False, 
            node_size=10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()   
