import os
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import seaborn as sns

def main():
    
    f = "results/all_bounds_iGap_series_002_run03.txt"
    # filename = "results/all_bounds_GDELT_GCC_2010_weekly.txt"
    # filename = "all_bounds_GDELT_GCC_2011_weekly.txt"

    # sol = np.loadtxt(filename.replace('all_bounds', 'solution'))

    sol = np.loadtxt(f.replace('all_bounds', 'root_solution'))

    # plt.xkcd()

    # fig, axes  = plt.subplots(4, 4, sharex='col', sharey='row')
    fig = plt.figure(figsize=(10,10), dpi=100)
    
    count = 1
    
    for row in sol:
        print(str(row[0])+" "+str(row[1]))
        plot_pruned(f, row[0], row[1], fig, count)
        count += 1
        if count == 17:
            break

    fig.patch.set_visible(False)
    # axes.axis('off')

    # colorbar legend
    # test = plt.pcolor([[0,1],[0.2,2]], cmap=plt.cm.YlGnBu)
    # fig.subplots_adjust(right=0.8)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    # cbar = plt.colorbar(test, cax=cbar_ax)
    
    # ax = fig.add_subplot()
    # test = plt.pcolor([[0,1],[0.2,2]], cmap=plt.cm.YlGnBu)
    # lblue_patch = mpl.patches.Patch(color = 'cornflowerblue', label='Composite All')
    # dblue_patch = mpl.patches.Patch(color = 'midnightblue', label='All Bounds')
    # plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol = 2, handles=[lblue_patch, dblue_patch])



    plt.show()

    filename = f.replace(".txt", ".pdf").replace("results/", "plots/")
    # dir_name = str(f).replace(".txt","") + "/"
    # if not os.path.exists(dir_name):
    #     os.makedirs(dir_name)
    # filename = dir_name+"combined16.pdf"
    filename = "heat_"+filename;
    fig.savefig(filename, format='pdf', dpi=1000)

def plot_pruned(all_bounds, alpha, solution, fig, i):

    com_bounds = all_bounds.replace("all","com")
    
    m_all = np.matrix(pruned(np.loadtxt(all_bounds), alpha, solution))
    m_com = np.matrix(pruned(np.loadtxt(com_bounds), alpha, solution))

    c1 = plt.cm.Reds
    c2 = plt.cm.Blues

    c3 = plt.cm.GnBu

    ax = fig.add_subplot(4,4,i)
   
    s = len(m_all)

    m = np.add(m_all, m_com)

    ax.pcolor(m.tolist(), cmap = c2, label=str(alpha))
    
    # ax.pcolor(m_all, cmap =c1)
    leg = ax.legend(loc='upper right', handlelength=1, handletextpad=1, fancybox=True)
    for item in leg.legendHandles:
        item.set_visible(False)

    if i != 13:
        ax.axis('off')
    else:
        plt.xlabel("Time Periods")
        plt.ylabel("Bound Size")

def pruned(bounds, alpha, solution):
    ROOT = 1.0/alpha

    # name = 'Bounds norm for '+str(alpha)
    actual_size = bounds[len(bounds)-1][1]
    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    for row in bounds:
        time_length = row[1] - row[0]
        norm = row[2] / mt.pow(time_length + 1, ROOT)
        if (norm > solution):
            time_values[time_length].append(0.1)
        else:
            time_values[time_length].append(1)

    s = int(actual_size) + 1
    m = []

    for i in range(0, s):
        m.append(time_values[i])

    for i in range(0, s):
        l = m[i]
        for x in range(0, i):
            l.append(0)
    
    return m

if __name__ == "__main__":
    main()
