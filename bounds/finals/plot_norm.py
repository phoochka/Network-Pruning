import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    data = np.loadtxt("data_norm.csv", delimiter=",", skiprows=1)
    for row in data:
        print(row)

    fig =  sns.plt.figure()
    sns.set_style("whitegrid")
    
    sns.set_context({"figure.figsize": (20, 5)})
    
    # Define colors
    # blue = sns.xkcd_rgb["denim blue"]
    # red = sns.xkcd_rgb["pale red"]
    # green = sns.xkcd_rgb["medium green"]
    pink = '#feb1c0'
    blue = '#00006c'

    # Overlay
    sns.barplot(x = data[:,0], y =(data[:,3]/data[:,3] * 100) , color = 'white')

    # MixedGroup
    grouped  = sns.barplot(x = data[:,0], y = (data[:,2]/data[:,3] * 100), color = blue)

    # Composite
    composite = sns.barplot(x = data[:,0], y = (data[:,1]/data[:,3] * 100), color = pink)



    group_bar = plt.Rectangle((0,0),1,1,fc=blue, edgecolor = 'none')
    com_bar = plt.Rectangle((0,0),1,1,fc=pink,  edgecolor = 'none')

    # l = plt.legend([group_bar, com_bar], ['Composite', 'Grouped'],bbox_to_anchor=(0., 1., 1., .1), loc=3, ncol=2, mode="expand", prop={'size':17})
    l = plt.legend([group_bar, com_bar], ['Composite', 'Grouped'], loc='lower right', ncol=2, prop={'size':25})
    l.draw_frame(True)

    #Optional code - Make plot look nicer
    sns.despine(left=True)
    composite.set_ylabel("Pruning %")
    composite.set_xlabel(r'$\alpha$ ')

    #Set fonts to consistent 16pt size
    for item in ([composite.xaxis.label, composite.yaxis.label] +
                         composite.get_xticklabels() + composite.get_yticklabels()):
            item.set_fontsize(23)

    plt.gcf().subplots_adjust(bottom=0.15)

    plt.show()

    fig.savefig("norm_stack.pdf", format='pdf', dpi=1800)
    

if __name__ == "__main__":
    main()
