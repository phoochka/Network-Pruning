import os
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def main():
    
    filename = "results/all_bounds_iGap_series_002_run03.txt"
    # filename = "all_bounds_GDELT_GCC_2010_weekly.txt"
    # filename = "all_bounds_GDELT_GCC_2011_weekly.txt"

    # sol = np.loadtxt(filename.replace('all_bounds', 'solution'))

    sol = np.loadtxt(filename.replace('all_bounds', 'root_solution'))

    for row in sol:
        print(str(row[0])+" "+str(row[1]))
        plot_root_pruned(np.loadtxt(filename), row[0], row[1], filename)

    # single_file(filename, -0.01, 0.003442)
    # single_file(filename, -0.03, 0.003307)
    # single_file(filename, -0.05, 0.003177)
    # single_file(filename, -0.06, 0.002973)
    # single_file(filename, -0.09, 0.002202)
    # single_file(filename, -0.13, 0.001476)


    # filename = "all_bounds_iGap_series_002_run03.txt"
    # single_file(filename, -0.005, 0.008808) 
    # single_file(filename, -0.04, 0.008213) 
    # plot_cube_pruned(np.loadtxt(filename), 0.002502, filename)
    # all_files()
 
def single_file(f, alpha, solution):
    all_bounds = np.loadtxt(f)
    analyze_bounds(all_bounds, f, alpha, solution)


def all_files():
    
    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if f.endswith(".txt"):
            l = np.loadtxt(f)
            # simple_plot([row[2] for row in l], f)
            if f.startswith("all_"):
                analyze_bounds(l, f)

def analyze_bounds(raw_bounds, f, alpha = None, solution = None):

    # alphas = [-0.001, -0.005, -0.01, -0.1]
    fig = sns.plt.figure()

    # for alpha in alphas:
        # plot_heat(raw_bounds, alpha)

    if (alpha):
        plot_alpha_pruned(raw_bounds, alpha, solution)
    else:
        plot_heat_bounds(raw_bounds)

    # sns.plt.show()
    
    #
    # plt.xlabel('Time Periods')
    # plt.ylabel('Conductance')
    # lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc =2, borderaxespad=0.)
    #
    # plt.tight_layout()
    # art = []
    # art.append(lgd)
    #
    # fig.savefig(filename, format='pdf', dpi=1200, additional_artists = art, bbox_inches="tight")

    filename = f.replace(".txt", ".png")
    filename = "alpha_test/heat_"+str(alpha)+"_"+filename
    fig.savefig(filename, format='png', dpi=600)


def plot_root_pruned(all_bounds, alpha, solution, f):

    ROOT = 1.0/alpha

    # name = 'Bounds norm for '+str(alpha)
    actual_size = all_bounds[len(all_bounds)-1][1]
    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    for row in all_bounds:
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


    fig, ax = plt.subplots()
    # sns.heatmap(m)
    heatmap = ax.pcolor(m, cmap=plt.cm.Reds)
    # sns.plt.show()
    # filename = f.replace(".txt", ".pdf")
    dir_name = f.replace(".txt","") + "/"

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    filename = dir_name+"all_"+str(alpha)+".pdf"
    print("saving at "+filename)
    fig.savefig(filename, format='pdf', dpi=1200)
    
def plot_alpha_pruned(all_bounds, alpha, solution):

    name = 'Bounds norm for '+str(alpha)
    actual_size = all_bounds[len(all_bounds)-1][1]
    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    for row in all_bounds:
        time_length = row[1] - row[0]
        norm = row[2] * mt.exp(alpha * time_length + 1)
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


    sns.heatmap(m)

def plot_heat_bounds(all_bounds):

    name = 'All Bounds'
    actual_size = all_bounds[len(all_bounds)-1][1]
    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    for row in all_bounds:
        time_length = row[1] - row[0]
        # norm = row[2] * m.exp(alpha * time_length + 1)
        norm = row[2]
        time_values[time_length].append(norm)


    s = int(actual_size)
    m = []

    for i in range(0, s):
        m.append(time_values[i])

    for i in range(0, s):
        l = m[i]
        for x in range(0, i):
            l.append(0)


    sns.heatmap(m)

if __name__ == "__main__":
    main()
