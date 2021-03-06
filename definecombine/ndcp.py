from math import floor, ceil
import matplotlib.pyplot as plt
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
    y = np.repeat(y, 2)[1:]
    y = np.append(y, [N])
    x = np.repeat(x, 2)[:-1]
    x = np.append(x, [2 * N * D])

    plt.plot(x, y)
    plt.xlabel("Definer Support ($P^D$)")
    plt.ylabel("Definer Utility ($U_D$)")
    plt.title(f"Definer Utility Curve for $N$ = {N}, $D$ = {D}")
    plt.xticks(np.linspace(0, 2 * N * D, N + 1))
    if N > 10:
        plt.xticks(np.linspace(0, 2 * N * D, 13))
    
    # Plot conjectured asymptotic utility curve
    x2 = [0, 2 / 3. * N * D, N * D, 2 * N * D]
    y2 = [0, N / 3., N, N]
    plt.plot(x2, y2, '--')

    plt.legend(['Exact', 'Conjectured Asymptotic'])

    if output_filename:
        plt.savefig(output_filename, dpi=200)
    
    plt.show()
    return


if __name__ == '__main__':
    SEARCH_FOR_COUNTEREXAMPLES = False
    PLOT_UTILITY = False
    PRINT_THRESHOLDS = True

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


    if SEARCH_FOR_COUNTEREXAMPLES:
        Nmax = 25
        Dmax = 100
        num_iterations = 1000000

        num_counterexamples = 0
        for iter_index in range(1, num_iterations + 1):
            if iter_index % 1000 == 0:
                print("Iteration %d." % iter_index)

            pA, N, D = random_instance(Nmax, Dmax)
            
            uA = solve_ndcp(pA, N, D)

            # print(sum(max_pack(pA, N, D)), max_pack(pA, N, D))
            # print(sum(max_crack(pA, N, D)), max_crack(pA, N, D))
            # print(sum(max_pack_minus_one(pA, N, D)), max_pack_minus_one(pA, N, D))

            q = ceil(uA)
            is_integral_uA = abs(uA - int(uA) < 0.1)

            threshold_fns = []
            if is_integral_uA:
                threshold_fns = [threshold_i_integral, threshold_ii_integral, threshold_iii_integral]
            else:
                threshold_fns = [threshold_I_half_integral, threshold_II_half_integral, threshold_III_half_integral]
            
            thresholds = [fn(N, D, q) for fn in threshold_fns]
            is_counterexample = True
            for index, fn in enumerate(threshold_fns):
                if pA >= thresholds[index]:
                    is_counterexample = False
                    break
            if is_counterexample:
                print('***COUNTEREXAMPLE***')
                print("Iteration %d" % iter_index)
                print(pA, N, D)
                print("Definer: %.1f" % uA)
                print("Combiner: %.1f" % (N-uA))
                print()
                num_counterexamples += 1
        print('Finished %d iterations. Found %d counterexamples.' % (num_iterations, num_counterexamples))
