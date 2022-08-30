from math import floor, ceil
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import networkx as nx
import numpy as np
import random
import sys

# Conjectured thresholds for U_A >= q
def threshold_i_integral(N, D, q):
    return 2 * q * D

def threshold_ii_integral(N, D, q):
    return (N - 1 + q) * (floor(D / 2.) + 1) + ceil(D / 2.)

def threshold_iii_integral(N, D, q):
    return (2 * q - 1) * (D - 1) + (2 * N - 2 * q + 1) * 2


# Conjectured thresholds for U_A >= q - 1/2
def threshold_I_half_integral(N, D, q):
    return (2 * q - 1) * D

def threshold_II_half_integral(N, D, q):
    return (N - 3 + q) * (floor(D / 2.) + 1) + ceil(D / 2.) + D

def threshold_III_half_integral(N, D, q):
    return (2 * q - 1) * (D - 1) + (2 * N - 2 * q) * 2 + 1


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



def t21(pA, N, D):
    q = floor((pA // D - 1) / 2.)
    subdistA = np.zeros((2*N,))
    for i in range(2 * q - 1):
        subdistA[i] = D
    
    pArem = pA - (2 * q - 1) * D

    i = 2 * q - 1
    while pArem > 0:
        if subdistA[i] < D:
            subdistA[i] += 1
            pArem -= 1
        i += 1
        i = i % (2 * N)
    return subdistA


def t22(pA, N, D):
    # Make two ceil(D/2) subdistricts, then as many floor(D/2) + 1 as possible
    subdistA = np.zeros((2*N,))
    subdistA[0] = ceil(D / 2)
    subdistA[1] = ceil(D / 2)

    pArem = pA - sum(subdistA)

    i = 2
    while pArem >= floor(D / 2 + 1):
        subdistA[i] += floor(D / 2) + 1
        pArem -= subdistA[i]
        i += 1
        i = i % (2 * N)
    
    subdistA[i] += pArem # Allocate any leftovers

    return subdistA


def t23(pA, N, D):
    q = 2
    while (2 * q - 1) * (D - 1) + (2 * N - 2 * q) * 2 + 1 <= pA:
        q += 1
    q -= 1

    # Make 2q - 1 subdistricts with D - 1, then one with 1, then 2 * N - 2 * q with 2.
    subdistA = np.zeros((2*N,))
    for i in range(2 * q - 1):
        subdistA[i] = D - 1

    subdistA[2 * q - 1] = 1

    for i in range(2 * q, 2 * N):
        subdistA[i] = 2

    pArem = pA - sum(subdistA)

    i = 0
    while pArem > 0: # Evenly distribute remainder
        subdistA[i] += 1
        pArem -= 1
        i += 1
        i = i % (2 * N)
    return subdistA


def random_subdistricts(pA, N, D):
    # Randomly allocates pA across the subdistricts
    pArem = pA
    subdistA = np.zeros((2*N,))
    for i in range(2 * N):
        if pArem <= 0:
            break
        max_subdist = min(D, pArem)
        subdistA[i] = random.randint(0, max_subdist)
        pArem -= subdistA[i]

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
    min_mwpm_val = N
    # TEST: Find min-weight perfect matching for each of the three definer strategies considered
    for definer_strategy in [max_pack, max_crack, max_pack_minus_one, random_subdistricts]:
        subdistA = definer_strategy(pA, N, D)

        mwpm, G = combine_optimally(subdistA, N, D)   

        mwpm_val = 0.0
        for u, v in mwpm:
            mwpm_val += G[u][v]['weight']
        
        if mwpm_val < min_mwpm_val:
            min_mwpm_val = mwpm_val

        # print(str(definer_strategy))
        # print(subdistA)
        # print("Definer: ", N - mwpm_val)
        # print()

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



def random_instance(Nmax, Dmax):
    N = random.randint(4, Nmax)
    D = random.randint(6, Dmax)

    # Adversarially choose pA to be one less than some threshold; 
    # if we're going to find a counterexample, this is a good place to look
    q = random.randint(2, N)
    is_integral_uA = random.random() > 0.5

    if is_integral_uA:
        pA = -1 + min([threshold_i_integral(N, D, q), 
                       threshold_ii_integral(N, D, q), 
                       threshold_i_integral(N, D, q)]
                     )
    else:
        pA = -1 + min([threshold_I_half_integral(N, D, q),
                       threshold_II_half_integral(N, D, q),
                       threshold_III_half_integral(N, D, q)]
                     )

    return (pA, N, D)


"""
Compute district vote-shares for definer in 
nongeometric define-combine procedure (NDCP). 

Parameters:
===============
    N - number of districts (int)
    D - units per subdistrict (int)
    definer_votes - total votes (units) for the definer (int)

Returns:
===============
    list of definer vote-shares in each district, sorted ascending
"""
def district_vote_shares(N, D, definer_votes):
    best_definer_utility = 0
    best_district_vote_shares = None
    # Find min-weight perfect matching for each of the three definer strategies considered
    subdistricts = []
    for index, definer_strategy in enumerate([max_pack, max_crack, max_pack_minus_one, t21, t22, t23]):
        subdistricts.append(definer_strategy(definer_votes, N, D))

        mwpm, G = combine_optimally(subdistricts[index], N, D)

        mwpm_val = 0.0
        for u, v in mwpm:
            mwpm_val += G[u][v]['weight']

        definer_utility = N - mwpm_val
        
        if definer_utility >= best_definer_utility:
            best_definer_utility = definer_utility
            best_district_vote_shares = [subdistricts[index][u] + subdistricts[index][v] for u, v in mwpm]
    
    if not best_district_vote_shares:
        raise Exception("No best definer strategy found.")
    return [dvs * 1. / (2 * D) for dvs in best_district_vote_shares]


"""
Compute seats-votes curve for 
the nongeometric define-combine procedure (NDCP).

Parameters:
===============
    N - number of districts (int)
    D - units per subdistrict (int)
    definer_votes - total votes (units) for the definer (int)

Returns:
===============
    x - list of x-coordinates of points on seats-votes curve
    y - list of y-coordinates of points on seats-votes curve
"""
def seats_votes_curve(N, D, definer_votes):
    vs = district_vote_shares(N, D, definer_votes)
    print('District vote-shares from DCP({},{},{}): {}'.format(N, D, definer_votes, vs))
    unique_vals, counts = np.unique(vs, return_counts=True)
    x = [0]
    y = [0]

    for index, val in enumerate(unique_vals):
        if val < 0.5:
            # Determine how much must be added to total vote-share to make these districts wins
            diff = 0.5 - val
            vs_hypothetical = [min(v + diff, 1.0) for v in vs]
            x_hypothetical = sum(vs_hypothetical) * D * 2
            y_hypothetical_2 = sum(map(lambda v : v >= 0.5, vs_hypothetical))
            y_hypothetical_1 = y_hypothetical_2 - counts[index]
        else: # val >= 0.5
            # Determine how much must be removed from total vote-share to make these losses
            diff = val - 0.5
            vs_hypothetical = [max(v - diff, 0.0) for v in vs]
            x_hypothetical = sum(vs_hypothetical) * D * 2
            y_hypothetical_2 = sum(map(lambda v : v >= 0.5, vs_hypothetical))
            y_hypothetical_1 = y_hypothetical_2 - counts[index]

        x.extend([x_hypothetical, x_hypothetical])
        y.extend([y_hypothetical_1, y_hypothetical_2])

    x.append(2 * N * D)
    y.append(N)

    x = sorted(x)
    y = sorted(y)

    return sorted(x), sorted(y)


def plot_utility_curve(N, D, output_filename=None):
    plt.figure()

    y = np.linspace(0, N, 2*N+1)
    x = np.zeros(y.shape)

    for i in range(1, len(x)):
        def_util = y[i]
        q = ceil(def_util)

        if def_util == int(def_util):
            min_threshold = min(threshold_i_integral(N, D, q), 
                         threshold_ii_integral(N, D, q),
                         threshold_iii_integral(N, D, q))
        else:
            min_threshold = min(threshold_I_half_integral(N, D, q),
                                threshold_II_half_integral(N, D, q), 
                                threshold_III_half_integral(N, D, q))
        
        x[i] = min_threshold

    # Double all intermediate values to achieve step function, and add top right point
    # y = np.repeat(y, 2)
    # x = np.repeat(x, 2)[1:]
    # x = np.append(x, [2 * N * D])

    # plt.plot(x, y)
    # plt.xlabel("Definer Vote-Share")
    # plt.ylabel("Definer Utility")

    # plt.title(f"Definer Utility Curve for $N$ = {N}, $D$ = {D}")
    # plt.xticks(np.linspace(0, 2 * N * D, 2 * N + 1))
    # if N > 10:
    #     plt.xticks(np.linspace(0, 2 * N * D, 13))
    
    # Plot conjectured asymptotic utility curve
    x2 = [0, 2 / 3. * N * D, N * D, 2 * N * D]
    x2 = [xi / (2. * N * D) for xi in x2]
    y2 = [0, N / 3., N, N]
    y2 = [yi / N for yi in y2]
    plt.plot(x2, y2, linestyle='solid')
    plt.legend(['Exact', 'Conjectured Asymptotic'])

    # ICYF
    x3 = np.linspace(0, 1, 10001)
    y3 = np.zeros(x3.shape)
    for i, xi in enumerate(x3):
        if xi <= 0.5:
            y3[i] = 2 * (xi ** 2)
        else:
            y3[i] = 1 - 2 * (1 - xi) ** 2
    plt.plot(x3, y3, linestyle='dashdot')

    # Bisection
    x4 = [0, 1]
    y4 = [0, 1]
    plt.plot(x4, y4, linestyle='dotted')

    # One party draws all (One-sided)
    x5 = [0, 0.5, 1]
    y5 = [0, 1, 1]
    plt.plot(x5, y5, linestyle='dashed')

    plt.xticks(np.linspace(0, 1, 11))

    plt.xlabel("Player 1 Fractional Vote-Share")
    plt.ylabel("Player 1 Utility")
    plt.legend(['DCP', 'I-cut-you-freeze', 'Bisection', 'Player 1 only'])
    plt.title('Limiting Nongeometric Utility Curves for Redistricting Protocols')

    if output_filename:
        plt.savefig(output_filename, dpi=200)
    
    plt.show()
    return


if __name__ == '__main__':
    SEARCH_FOR_COUNTEREXAMPLES = False
    PLOT_UTILITY = True
    PRINT_THRESHOLDS = False

    if PRINT_THRESHOLDS:
        if len(sys.argv) < 3:
            raise Exception("Missing argument(s). Usage: python ndcp.py [N] [D].")
        
        N = int(sys.argv[1])
        D = int(sys.argv[2])

        for q in range(1, N + 1):
            for uD in [q - 0.5, q]:
                print("Definer utility:", uD)
                if uD == q:
                    min_threshold = min(threshold_i_integral(N, D, q), threshold_ii_integral(N, D, q), threshold_iii_integral(N, D, q))
                else:
                    min_threshold = min(threshold_I_half_integral(N, D, q), threshold_II_half_integral(N, D, q), threshold_III_half_integral(N, D, q))
                    if uD == 3.5:
                        print("\t", threshold_I_half_integral(N, D, q))
                        print("\t", threshold_II_half_integral(N, D, q))
                        print("\t", threshold_III_half_integral(N, D, q))
                        
                print("\tMinimum vote-share:", min_threshold)
        

    if PLOT_UTILITY:
        if len(sys.argv) < 3:
            raise Exception("Missing argument(s). Usage: python ndcp.py [N] [D] [(optional) output filename].")
        
        N = int(sys.argv[1])
        D = int(sys.argv[2])

        output_filename = None

        if len(sys.argv) >= 4:
            output_filename = sys.argv[3]
        
        plot_utility_curve(N, D, output_filename)


