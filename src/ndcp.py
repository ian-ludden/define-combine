from math import floor, ceil
import networkx as nx
import numpy as np

def max_pack(pA, N, D):
    # Creates as many packed subdistricts as possible, then puts remainder all together
    pArem = pA
    subdistA = np.zeros((2*N,))
    for i in range(2 * N):
        if (pArem >= D):
            subdistA[2 * N - 1 - i] = D
            pArem -= D
        else:
            subdistA[2 * N - 1 - i] = pArem
            break
    
    return subdistA

def max_crack(pA, N, D):
    # Creates as many floor(D/2) + 1 subdistricts as possible, then puts remainder all together
    max_subdist = floor(D / 2.) + 1

    pArem = pA
    subdistA = np.zeros((2*N,))
    for i in range(2 * N):
        if (pArem >= max_subdist):
            subdistA[2 * N - 1 - i] = max_subdist
            pArem -= max_subdist
        else:
            subdistA[2 * N - 1 - i] = pArem
            pArem = 0
            break

    # Evenly distribute any remainder
    i = 0
    while (pArem > 0):
        subdistA[i] += 1
        pArem -= 1
        i += 1
        i = i % (2 * N)
    
    return subdistA

def max_pack_minus_one(pA, N, D):
    # Creates as many (D-1) subdistricts as possible while maintaining that 
    # the rest contain at least 2
    max_subdist = D - 1

    pArem = pA
    subdistA = np.zeros((2*N,))
    for i in range(2 * N):
        if (pArem >= max_subdist + 2 * (2 * N - (i + 1))):
            subdistA[2 * N - 1 - i] = max_subdist
            pArem -= max_subdist
        else:
            subdistA[2 * N - 1 - i] = 2
            pArem -= 2

    # Evenly distribute any remainder
    i = 0
    while (pArem > 0):
        subdistA[i] += 1
        pArem -= 1
        i += 1
        i = i % (2 * N)
    
    return subdistA


""" 
    Solves the nongeometric define-combine procedure (NDCP). 

    Parameters:
        pA  - number of voters for player A
        N   - number of districts
        D   - size of each subdistrict

    Implied are:
        P = N * 2D  - total population
        pB = P - pA - number of voters for player B
    
    Returns:
        The utility of player A (adding 1 for wins, 1/2 for ties) 
        under optimal play
"""
def solve_ndcp(pA, N, D):
    min_mwpm_val = D
    # TEST: Find min-weight perfect matching for each of the three definer strategies considered
    for definer_strategy in [max_pack, max_crack, max_pack_minus_one]:
        subdistA = definer_strategy(pA, N, D)

        mwpm, G = combine_optimally(subdistA, N, D)   
              
        mwpm_val = 0
        for u, v in mwpm:
            mwpm_val += G[u][v]['weight']
        
        if mwpm_val < min_mwpm_val:
            min_mwpm_val = mwpm_val

        print(str(definer_strategy))
        print(subdistA)
        print("Definer: ", N - mwpm_val)
        print()

    return N - min_mwpm_val


def combine_optimally(subdistricts, N, D):
    # Given definer partition of pA across 2N subdistricts, 
    # return best perfect matching for the combiner. 

    # Compute edge weights (from B's perspective)
    edge_weight = np.zeros((2*N, 2*N))
    edge_list = []

    for i in range(2 * N):
        for j in range(i + 1, 2 * N):
            edge_list.append((i, j))

            dist_sum = subdistricts[i] + subdistricts[j]
            if dist_sum < D:
                edge_weight[i, j] = 1
            elif dist_sum == D:
                edge_weight[i, j] = 0.5
            else:
                pass

    # Use networkx to solve the min-weight perfect matching problem
    mwpm = set()
    G = nx.complete_graph(2 * N)
    for edge in G.edges:
        u, v = edge
        G[u][v]['weight'] = edge_weight[u, v]

    return (nx.max_weight_matching(G, maxcardinality=True), G)


if __name__ == '__main__':
    pA = 3*99 + 9*2
    N = 6
    D = 100

    print(pA, N, D)

    uA = solve_ndcp(pA, N, D)
    print("Definer: %.1f" % uA)
    print("Combiner: %.1f" % (N-uA))

    print(sum(max_pack(pA, N, D)), max_pack(pA, N, D))
    print(sum(max_crack(pA, N, D)), max_crack(pA, N, D))
    print(sum(max_pack_minus_one(pA, N, D)), max_pack_minus_one(pA, N, D))