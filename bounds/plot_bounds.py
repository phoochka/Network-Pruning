import os
import numpy as np
from matplotlib import pyplot as plt


def main():

    # all_bounds = np.loadtxt("all_bounds_GDELT_GCC_2010_weekly.txt", unpack=False)
    # simple_plot([row[2] for row in all_bounds], 'all_bounds.pdf')

    # com_bounds = np.loadtxt("com_bounds_GDELT_GCC_2010_weekly.txt", unpack=False)
    # simple_plot([row[2] for row in com_bounds], 'com_bounds.pdf')

    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:
        if f.endswith(".txt"):
            l = np.loadtxt(f)
            simple_plot([row[2] for row in l], f)


def simple_plot(l, name):
    name = name.replace(".txt", ".pdf")
    print(name)
    fig = plt.figure()
    plt.plot(l)
    fig.savefig(name, format='pdf', dpi=1200)

if __name__ == "__main__":
    main()
