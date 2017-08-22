import os
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import seaborn as sns


def main():
    
    filename = "com_d_moreVertices_time_0100_run000.txt"

    sol = np.loadtxt("solution.txt")
    # Normalizing by time span
    time_length = 10
    sol = sol * mt.pow(time_length, 0.2)
    print(sol)

    plot_pruned(filename, sol)


def plot_pruned(com_bounds, solution):
    mpl.rc('xtick', labelsize=20) 
    mpl.rc('ytick', labelsize=20) 
    # mpl.rcParams.update({'font.size': 22})
    fig, ax =  plt.subplots()

    # fig =  sns.plt.figure()
    # sns.set_style("whitegrid")
    
    # sns.set_context({"figure.figsize": (10, 5)})


    all_bounds = com_bounds.replace("com","all")
    
    m_all = np.matrix(pruned(np.loadtxt(all_bounds), solution))
    m_com = np.matrix(pruned(np.loadtxt(com_bounds), solution))
    
    m_grp = np.matrix(np.loadtxt('pruned_by_grouping.txt'))

    # sns.heatmap(np.loadtxt(com_bounds))
    # plt.show()
    
    m = np.add(m_all, m_com)
    m = np.add(m, m_grp)


    # Colors
    c1 = plt.cm.Reds
    c2 = plt.cm.Blues
    c3 = plt.cm.GnBu
    c4 = plt.cm.Greys

    plt.pcolor(m.tolist(), cmap=c2)

    ax.set_ylabel("Length", fontsize=20)
    ax.set_xlabel('Starting Time', fontsize=20)

    
    darkb = '#091d4c'
    medb = '#2569ad'
    litb = '#92c0db'

    grouped_patch = mpatches.Patch(color=darkb, linewidth=0)
    compos_patch = mpatches.Patch(color=medb, linewidth=0)
    all_patch = mpatches.Patch(color=litb, linewidth=0)

    # plt.legend((p,), ('Entry',))
    ax.legend((grouped_patch, compos_patch, all_patch), ('Grouping', 'Composite','Full'), fontsize=26)

    

    # sns.heatmap(m, cmap="YlGnBu")
    # ax = sns.heatmap(m, xticklabels=False, yticklabels=20)

    plt.gcf().subplots_adjust(bottom=0.15)
    plt.show()

    # plt.xlabel("Time Periods")
    # plt.ylabel("Bound Size")
    fig.savefig("heat.pdf", format='pdf', dpi=1200)

def pruned(bounds, solution):
    for x in np.nditer(bounds, op_flags=['readwrite']):
        # print(x, solution)
        if x in (None, 0.0) or x is False:
            x[...] = 0
        elif x > solution:
            x[...] = 1
        else:
            x[...] = 0.1
    
    return bounds

if __name__ == "__main__":
    main()
