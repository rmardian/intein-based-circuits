import matplotlib.pyplot as plt
import seaborn as sns


def plot_heatmap(heatmaps, xlabel, ylabel, cmap='cividis'):
    '''generating 2-dimensional induction-level plot'''
    
    f, axs = plt.subplots(2, 3, sharex=False, sharey=False, figsize=(18, 6))
    axr = axs.ravel()
    for i, ax in enumerate(axr):
        if i < len(heatmaps):
            sns.heatmap(heatmaps[i][1], annot=True, fmt='.2f', cmap=cmap, ax=ax)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_title(heatmaps[i][0])
        else:
            ax.set_visible(False)
    plt.tight_layout()
    plt.show()