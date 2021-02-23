import cooler
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool
import pandas as pd
import numpy

threshold = 0.8


# Replace all values in rows and columns with NaN if there are more than >= threshold zeros:

def zeros_to_nan(arr):
    arr = arr.astype(float)
    n = arr.shape[0]
    for i in range(len(arr)):
        if ((arr[i] == 0).sum(0) / n) >= threshold:
            arr[i] = numpy.nan
            arr[:, i] = numpy.nan
    return arr


# Normalization to Obs/Exp matrix:

def distribution_at_dist(arr, d):
    n = arr.shape[0]
    return np.array([arr[i, j] for i, j in zip(range(0, n - d), range(d, n))])


def normalize_intra(arr):
    n = arr.shape[0]
    averages_at_dist = [np.nanmean(distribution_at_dist(arr, d)) for d in range(0, n)]
    ans = np.zeros_like(arr, dtype='float64')
    for i in range(n):
        for j in range(n):
            ans[i, j] = arr[i, j] / averages_at_dist[abs(i - j)]

    return ans