import os
import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

CUBE = 1.0/3.0


# allb = [(row[2] / m.pow( (row[1] - row[0] + 1), CUBE)) for row in all_bounds]

def main():
    
    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if f.endswith(".txt"):
            l = np.loadtxt(f)
            # simple_plot([row[2] for row in l], f)
            if f.startswith("all_"):
                analyze_bounds(l, f)




    plt.show()

def analyze_bounds(raw_bounds, f):

    print(f)
    
    alphas = [1, 2, 3]
    fig = plt.figure()

    for alpha in alphas:
        plot_bounds(raw_bounds, alpha, f)

    filename = f.replace(".txt", ".pdf")
    filename = "plot_"+filename

    plt.xlabel('Time Periods')
    plt.ylabel('Conductance')
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc =2, borderaxespad=0.)

    plt.tight_layout()
    art = []
    art.append(lgd)

    fig.savefig(filename, format='pdf', dpi=1200, additional_artists = art, bbox_inches="tight")
    

def plot_bounds(all_bounds, alpha, f):

    actual_size = all_bounds[len(all_bounds)-1][1]

    time_values = {}

    for i in range(0,int(actual_size)+1):
        time_values[i] = []

    for row in all_bounds:
        time_length = row[1] - row[0]
        norm = row[2] / m.pow(time_length + 1, 1)
        time_values[time_length].append(norm)

    means = []
    stds = []

    maxes = []
    mins = []

    for length in time_values:
        total = 0
        t_min = min(time_values[length]) 
        t_max = max(time_values[length]) 
        
        for i in time_values[length]:
            total = total + i

        stds.append(np.std(time_values[length]))
        mean = total/len(time_values[length])
        means.append(mean)

        maxes.append(t_max - mean)
        mins.append(mean - t_min)

    maxes = np.array(maxes)
    mins = np.array(mins)
    assym_error = [maxes, mins]

    print(assym_error)
    print(len(assym_error))

    # plt.errorbar(means, stds, xerr=0.2, yerr=0.4)
    # plt.show()

    times = list(range(0,len(means)))
    name = 'Bounds norm for '+str(alpha)

    # plt.errorbar(times, means, yerr=stds, label=name)
    plt.errorbar(times, means, yerr=assym_error, label=name)


if __name__ == "__main__":
    main()
