from math import floor, ceil
import matplotlib.pyplot as plt
import numpy as np
import sys

from ndcp import threshold_i_integral, threshold_ii_integral, threshold_iii_integral
from ndcp import threshold_I_half_integral, threshold_II_half_integral, threshold_III_half_integral

marker_styles = ["o", "+", "x"]
colors = ["xkcd:blue", "xkcd:orange", "xkcd:purple"]

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Missing arguments. Usage: ")
        print("\tpython ndcp_threshold_plots.py [number of districts] [voters per subdistrict]")
        exit(0)
    
    N = int(sys.argv[1])
    D = int(sys.argv[2])

    while N:

        if N < 4:
            print("The number of districts must be at least 4.")
            exit(0)
        if D < 6:
            print("The number of voters per subdistrict must be at least 6.")
            exit(0)

        P = 2 * N * D
        floorD = floor(D/2.0)
        ceilD = ceil(D/2.0)

        whole_thresholds = [threshold_i_integral, threshold_ii_integral, threshold_iii_integral]
        half_thresholds = [threshold_I_half_integral, threshold_II_half_integral, threshold_III_half_integral]

        fig = plt.figure()

        x_whole = [i for i in range(N + 1)]
        y_whole = np.zeros((len(whole_thresholds), len(x_whole)))

        for fn_index, fn in enumerate(whole_thresholds):
            y_one = min(2*D, D + 2*N - 1, N * (floorD + 1) + ceilD)
            y_whole_fn = [0, y_one]

            for q in range(2, N + 1):
                threshold = fn(N, D, q)
                threshold = min(threshold, P)
                y_whole_fn.append(threshold)
            
            y_whole[fn_index, :] = y_whole_fn
            
            print(fn)
            print(y_whole[fn_index, :])

            plt.scatter(x_whole, y_whole[fn_index, :], marker=marker_styles[fn_index], c=colors[fn_index])

        x_half = [i - 0.5 for i in range(1, N + 1)]
        y_half = np.zeros((len(half_thresholds), len(x_half)))

        for fn_index, fn in enumerate(half_thresholds):
            y_half_fn = [D] # Takes D to guarantee a tie

            for q in range(2, N + 1):
                threshold = fn(N, D, q)
                threshold = min(threshold, P)
                y_half_fn.append(threshold)
            
            y_half[fn_index, :] = y_half_fn
            
            print(fn)
            print(y_half[fn_index, :])
            
            plt.scatter(x_half, y_half[fn_index, :], marker=marker_styles[fn_index], c=colors[fn_index])

        y_whole_min_index = np.argmin(y_whole, axis=0)
        for i, fn_index in enumerate(y_whole_min_index):
            plt.bar(i, y_whole[fn_index, i], width=0.5, color=colors[fn_index], alpha=0.5)
        
        y_half_min_index = np.argmin(y_half, axis=0)
        for i, fn_index in enumerate(y_half_min_index):
            plt.bar(i + 0.5, y_half[fn_index, i], width=0.5, color=colors[fn_index], alpha=0.5)

        plt.legend(["(i) or (I)", "(ii) or (II)", "(iii) or (III)"])
        plt.title(f"N = {N}, D={D}")
        plt.xlabel("Seats for Definer")
        plt.ylabel("Votes for Definer")
        plt.savefig(f"../out/thresholds_N{N}_D{D}.png")
        plt.show()
    
        repeat = input("Run again? (y/n)")
        if not (repeat.lower() == "y"):
            exit(0)
        
        N_str = input("N:")
        D_str = input("D:")

        try:
            N = int(N_str)
            D = int(D_str)
        except Exception as e:
            print(e)
            print("Invalid N or D provided.")
            exit(0)



