import numpy as np
from multiprocessing import Pool

import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class WorkerPool(object):
    def __init__(self, cores):
        self.pool = Pool(processes=cores)
        self.size = cores
    def map(self, function, tasks):
        return self.pool.map(function, tasks)

def priorTransformFromJs(p,js):
    pt = np.zeros(len(p))
    for idx in range(len(p)):
        lb = js['Distributions'][idx]['Minimum']
        ub = js['Distributions'][idx]['Maximum']

        pt[idx] = lb+p[idx]*(ub-lb)
    return pt

def resample_equal_with_idx(samples, weights, rstate=None):

    if rstate is None:
        rstate = np.random

    if abs(np.sum(weights) - 1.) > 1e-9:  # same tol as in np.random.choice.
        # Guarantee that the weights will sum to 1.
        warnings.warn("Weights do not sum to 1 and have been renormalized.")
        weights = np.array(weights) / np.sum(weights)

    # Make N subdivisions and choose positions with a consistent random offset.
    nsamples = len(weights)
    positions = (rstate.random() + np.arange(nsamples)) / nsamples

    # Resample the data.
    idx = np.zeros(nsamples, dtype=np.int)
    cumulative_sum = np.cumsum(weights)
    i, j = 0, 0
    while i < nsamples:
        if positions[i] < cumulative_sum[j]:
            idx[i] = j
            i += 1
        else:
            j += 1

    return samples[idx], idx


def getPosteriorFromResult(result):
    
    weights = np.exp(result.logwt - result.logz[-1]) # normalized weights
    samples, idx = resample_equal_with_idx(result.samples, weights)
    return samples, idx


# Plot histogram of sampes in diagonal
def plot_histogram(ax, theta):
  dim = theta.shape[1]
  num_bins = 30

  for i in range(dim):

    if (dim == 1):
      ax_loc = ax
    else:
      ax_loc = ax[i, i]

    hist, bins, _ = ax_loc.hist(
        theta[:, i], num_bins, density=True, color='lightgreen', ec='black')

    if i == 0:

      # Rescale hist to scale of theta -> get correct axis titles
      widths = np.diff(bins)
      if (dim > 1):
        hist = hist / np.max(hist) * (
            ax_loc.get_xlim()[1] - ax_loc.get_xlim()[0])
        bottom = ax_loc.get_xlim()[0]
        ax_loc.cla()
        ax_loc.bar(
            bins[:-1],
            hist,
            widths,
            color='lightgreen',
            ec='black',
            bottom=bottom)
        ax_loc.set_ylim(ax_loc.get_xlim())
        ax_loc.set_xticklabels([])
      else:
        ax_loc.cla()
        ax_loc.bar(bins[:-1], hist, widths, color='lightgreen', ec='black')

    elif i == theta.shape[1] - 1:
      ax_loc.set_yticklabels([])

    else:
      ax_loc.set_xticklabels([])
      ax_loc.set_yticklabels([])
    ax_loc.tick_params(axis='both', which='both', length=0)


#Plot scatter plot in upper triangle of figure
def plot_upper_triangle(ax, theta, lik):
  dim = theta.shape[1]
  if (dim == 1):
    return

  for i in range(dim):
    for j in range(i + 1, dim):
      if lik:
        ax[i, j].scatter(
            theta[:, j], theta[:, i], marker='o', s=3, alpha=0.5, c=lik)
      else:
        ax[i, j].plot(theta[:, j], theta[:, i], marker='.', s=1, alpha=0.5)
      ax[i, j].set_xticklabels([])
      ax[i, j].set_yticklabels([])
      ax[i, j].grid(b=True, which='both')


#Plot 2d histogram in lower triangle of figure
def plot_lower_triangle(ax, theta):
  dim = theta.shape[1]
  if (dim == 1):
    return

  for i in range(dim):
    for j in range(i):
      # returns bin values, bin edges and bin edges
      H, xe, ye = np.histogram2d(theta[:, j], theta[:, i], 10, density=True)
      # plot and interpolate data
      ax[i, j].imshow(
          H.T,
          aspect="auto",
          interpolation='spline16',
          origin='lower',
          extent=np.hstack((ax[j, j].get_xlim(), ax[i, i].get_xlim())),
          cmap=plt.get_cmap('jet'))

      if i < theta.shape[1] - 1:
        ax[i, j].set_xticklabels([])
      if j > 0:
        ax[i, j].set_yticklabels([])

def plotNetsedResult(result, savepath=None):

    samples, idx  = getPosteriorFromResult(result)
 
    numdim     = len(samples[0])
    numentries = len(samples)
    
    llk     = (result.logl[idx]).tolist()

    fig, ax = plt.subplots(numdim, numdim, figsize=(8, 8))
    samplesTmp = np.reshape(samples, (numentries, numdim))
    plt.suptitle(
      '{0} Plotter - \nNumber of Samples {1}'.format(
          str("Nested Sampling"), str(numentries)).strip(),
          fontweight='bold',
          fontsize=12)

    plot_histogram(ax, samplesTmp)
    plot_upper_triangle(ax, samplesTmp, llk)
    plot_lower_triangle(ax, samplesTmp)

#    for i in range(numdim):
#      ax[i, 0].set_ylabel(genList[idx]['Variables'][i]['Name'])
#      ax[-1, i].set_xlabel(genList[idx]['Variables'][i]['Name'])
    
    if (savepath==None):
        print("TEST")
        plt.show()
    else:
        plt.savefig(savepath)



