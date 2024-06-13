# Name:         plotter.py
# Purpose:      Plot histograms of variant datasets
#
# Author:       Tristan J. Hayeck, Christopher J. Sottolano, and Melissa A. Gilbert
#
# Created:      2024
# Copyright:    (c) Tristan J. Hayeck 2024
# Licence:      <your licence>
# -------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from itertools import cycle
matplotlib.use('pdf')


def getKernelDensityEstimation(values, x, bandwidth, kernel):
    """
    Calculate kernel density estimate
    :param values: List containing variant score means across replicates sets
    :param x: x values corresponding to the generated distribution
    :param bandwidth: Bandwidth of the kernel
    :param kernel: numpy array of kernel density estimate
    :return:
    """
    model = KernelDensity(kernel=kernel, bandwidth=bandwidth)
    model.fit(values[:, np.newaxis])
    log_density = model.score_samples(x)
    return np.exp(log_density)


def CeilXlim(xMin):
    """
    Return new xMin list where each element (x) is the smallest integer i, such that i >= x
    :param xMin: Vector containing minimum or maximum of variant score means across replicates
    :return: xMin
    """
    if xMin < 0:
        xMin = np.ceil(xMin) - 1
    else:
        xMin = np.ceil(xMin)
    return xMin


def plotSeperateHist(dataRef, datalib, LABEL, COLORs, N_bin, xRefLine, outputFile, fixedYAxis=False, yMaxOffset=1.2):
    """
    Plot histograms
    :param dataRef: Variant score means across replicates containing variant control set
    :param datalib: Variant score means across replicates containing library set
    :param LABEL: List containing legend labels for each plot
    :param COLORs: List colors for each plot
    :param N_bin: Number of bins used for plotting histograms
    :param xRefLine: Value used to represent the pathogenic threshold verrical line
    :param outputFile: Filename of output pdf
    :param fixedYAxis: (bool) If True, sets Y-axis as fixed the maximum histogram value
    :param yMaxOffset: Scalar multiplier of the histogram y-axis
    :return:
    """
    if len(COLORs) == 0:
        COLORs = ['grey']
    colorCycler = cycle(COLORs)
    # Determine the bins using the min and max of the data
    xMin = []
    xMax = []
    DATA = [dataRef, datalib]
    for iplt in range(len(DATA)):
        if len(DATA[iplt]) == 0:
            print(f'List {iplt} empty. Skipping.')
        else:
            for itern in range(len(DATA[iplt])):
                xMin.append(min(DATA[iplt][itern]))
                xMax.append(max(DATA[iplt][itern]))
    xMin = CeilXlim(min(xMin))
    xMax = CeilXlim(max(xMax))

    ## Plot Histogram
    HIST_BINS = np.linspace(xMin, xMax, N_bin)
    fig, axs = plt.subplots(iplt+1, 1, figsize=(8, 6))
    plt.subplots_adjust(hspace=.0)
    yAxLabelPos = -1
    yMax = []
    for iplt in range(len(DATA)):
        if len(DATA[iplt]) == 0:
            yAxLabelPos = iplt
            axs[iplt].xaxis.set_visible(False)
            axs[iplt].yaxis.set_visible(False)
            axs[iplt].spines['top'].set_visible(False)
            axs[iplt].spines['right'].set_visible(False)
            axs[iplt].spines['bottom'].set_visible(False)
            axs[iplt].spines['left'].set_visible(False)
        else:
            ySubMax = []
            for itern in range(len(DATA[iplt])):
                data = DATA[iplt][itern]
                color = next(colorCycler)
                axs[iplt].hist(data, bins=HIST_BINS, alpha=0.6, color='white', label=LABEL[itern] + '\nN=' + str(len(data)))
                axs[iplt].hist(data, bins=HIST_BINS, alpha=0.6, color=color, edgecolor='black', linewidth=0.5)
                axs[iplt].legend(loc='upper right', fontsize=12, frameon=False)
                y1, y2 = axs[iplt].get_ylim()
                ySubMax.append(y2)

                if iplt == 0:
                    # The x values corresponding to the generated distribution
                    x = np.linspace(data.min(), data.max(), data.shape[0])
                    # Get kernel density estimation for current DAF table
                    kernel = pd.DataFrame(getKernelDensityEstimation(data.to_numpy(), x.reshape(-1, 1), 0.05, 'gaussian'))
                    # kernel rescaled based on y-axis of histogram
                    kernel = kernel/max(kernel[0])*y2
                    axs[iplt].plot(x, kernel, c=color)
                    # Hide x-axis labels from all but bottom graph
                    if iplt != len(DATA):
                        axs[iplt].xaxis.set_visible(False)
            yMax.append(max(ySubMax))
            if fixedYAxis:
                axs[iplt].set(ylim=[0, max(yMax)])
            else:
                for iplta in range(len(yMax)):
                    axs[iplt].set(ylim=[0, yMax[iplta]*yMaxOffset])
                    # set x-axis upper limit
                    axs[iplt].set_xlim(left=0)
                    axs[iplt].set_xlim(right=1.25)

    # Add x axis label
    axs[iplt].set_xlabel('Abundance score', fontsize=16)
    # add vertical line to library
    axs[iplt].axvline(x=xRefLine, linewidth=2, color='black', linestyle=':')
    # Add y axis label
    axs[iplt].set_ylabel('Number of variants', fontsize=16)
    if yAxLabelPos == -1:
        axs[iplt].yaxis.set_label_coords(-0.08, 1.02)
    elif yAxLabelPos == 0:
        axs[iplt].set_position([0.125, 0.3, 0.8, 0.4])
    plt.savefig(outputFile, dpi=200)
