import os
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns

def main():
    
    filename = "results/all_bounds_iGap_series_002_run03.txt"
    # filename = "all_bounds_GDELT_GCC_2010_weekly.txt"
    # filename = "all_bounds_GDELT_GCC_2011_weekly.txt"

    # sol = np.loadtxt(filename.replace('all_bounds', 'solution'))

    sol = np.loadtxt(filename.replace('all_bounds', 'root_solution'))

    sns.set_style("white")
    
    # plt.xkcd()
    fig = plt.figure()

    import itertools
    marker = itertools.cycle(('*', 's', 'o', 'v')) 
    # for n in y:
    #         plt.plot(x,n, marker = marker.next(), linestyle='')

    count = 0
    for row in sol:
        count += 1
        print(str(row[0])+" "+str(row[1]))
        plot_root_pruned(np.loadtxt(filename), row[0], row[1], next(marker))

        if count == 4:
            break;

    plt.xlabel("Time Period", fontsize=16) 
    plt.ylabel("Percentage Pruned for Time Period", fontsize=16)
    plt.legend(loc='lower left')

    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='16')

    # plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().tight_layout()
    plt.show()
    
    filename = filename.replace(".txt",".pdf")
    filename = filename.replace("results/","plots/")
    fig.savefig(filename, format='pdf', dpi=600)

def all_files():
    
    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if f.endswith(".txt"):
            l = np.loadtxt(f)
            # simple_plot([row[2] for row in l], f)
            if f.startswith("all_"):
                analyze_bounds(l, f)


def plot_root_pruned(all_bounds, alpha, solution, m):

    ROOT = 1.0/alpha

    # name = 'Bounds norm for '+str(alpha)
    actual_size = all_bounds[len(all_bounds)-1][1]
    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    total_count  = 0.0
    pruned_count = 0.0

    for row in all_bounds:
        total_count += 1.0
        time_length = row[1] - row[0]
        norm = row[2] / mt.pow(time_length + 1, ROOT)
        if (norm > solution):
            pruned_count +=1.0
            time_values[time_length].append(1)
        else:
            time_values[time_length].append(0)

    total_percentage_pruned = (pruned_count / total_count) * 100
    pruned = []
    for key in time_values.keys():
        pruned.append((time_values[key].count(1) / len(time_values[key])) * 100)

    plt.plot(pruned, marker=m, label=r'$\alpha$ ' +str(alpha)+" Pruned: "+str(int(total_percentage_pruned))+"%")
    

if __name__ == "__main__":
    main()
