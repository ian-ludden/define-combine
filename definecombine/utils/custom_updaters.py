import networkx as nx

def aux_graph(partition):
    parent = partition.parent

    if parent: # Parent partition exists
        pass # Do not yet know a way to efficiently update aux_graph
        
    cut_edges_set = partition.cut_edges

    num_parts = len(partition.parts)
    
    H = nx.Graph()
    H.add_nodes_from([i for i in range(1, num_parts + 1)])

    for u, v in cut_edges_set:
        u_part = partition.assignment[u]
        v_part = partition.assignment[v]

        aux_edge = tuple(sorted((u_part, v_part)))

        if aux_edge in H.edges:
            continue

        u_part_votes_definer = partition.votes_definer[u_part]
        v_part_votes_definer = partition.votes_definer[v_part]
        
        u_part_total_votes = partition.votes_total[u_part]
        v_part_total_votes = partition.votes_total[v_part]

        frac_votes_definer = (u_part_votes_definer + v_part_votes_definer) \
            * 1. / (u_part_total_votes + v_part_total_votes)

        edge_weight = 1
        tol = 1e-6
        if (frac_votes_definer > 0.5 + tol): edge_weight = 0
        if (abs(frac_votes_definer - 0.5) <= tol): edge_weight = 0.5 

        H.add_edge(aux_edge[0], aux_edge[1], weight=edge_weight)
    
    return H


def definer_utility(partition):
    H = partition.aux_graph
    mwpm = nx.max_weight_matching(H, maxcardinality=True)
    mwpm_val = 0
    for u,v in mwpm:
        mwpm_val += H[u][v]['weight']
    
    return (len(H.nodes) // 2) - mwpm_val
