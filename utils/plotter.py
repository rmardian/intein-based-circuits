import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_responses(data, promoters, inducers, row=2, col=2, xlabel='RPU', ylabel='RPU'):
    
    f, axs = plt.subplots(row, col, sharex=False, sharey=False, figsize=(col*6, row*2))
    axr = axs.ravel()
    for i, ax in enumerate(axr):
        if i < len(promoters):
            ax.scatter(inducers[i], data[filter(lambda x: x.startswith(promoters[i]), data.index)])
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.loglog()
            ax.set_title(promoters[i])
        else:
            ax.set_visible(False)
    plt.tight_layout()
    sns.despine()

def plot_heatmap(heatmaps, row=2, col=3, xlabel='RPU', ylabel='RPU', cmap='cividis'):
    '''generating 2-dimensional induction-level plot'''
    
    f, axs = plt.subplots(row, col, sharex=False, sharey=False, figsize=(col*6, row*3))
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

def plot_trace(params, param_name='parameter', num_cols=2):
    
    num_rows = int(len(params)/num_cols)
    f, ax = plt.subplots(num_rows, num_cols*2, figsize=(7*num_cols, 2*num_rows))
    
    for i in tqdm(range(num_rows)):
        for j in range(num_cols):
            
            param = params[i*num_cols+j]
        
            mean = np.mean(param)
            median = np.median(param)
            cred_min, cred_max = np.percentile(param, 2.5), np.percentile(param, 97.5)

            ax[i, j*2].plot(param)
            ax[i, j*2].axhline(mean, color='r', lw=2, linestyle='--')
            ax[i, j*2].axhline(median, color='c', lw=2, linestyle='--')
            ax[i, j*2].axhline(cred_min, linestyle=':', color='k')
            ax[i, j*2].axhline(cred_max, linestyle=':', color='k')

            ax[i, j*2+1].hist(param, 30, density=True);
            sns.kdeplot(param, shade=True, ax=ax[i, j*2+1])
            ax[i, j*2+1].axvline(mean, color='r', lw=2, linestyle='--',label='mean')
            ax[i, j*2+1].axvline(median, color='c', lw=2, linestyle='--',label='median')
            ax[i, j*2+1].axvline(cred_min, linestyle=':', color='k', label='95% CI')
            ax[i, j*2+1].axvline(cred_max, linestyle=':', color='k')
            ax[i, j*2+1].set_ylabel(None)
            ax[i, j*2+1].legend()
    
    plt.suptitle(param_name)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.show()