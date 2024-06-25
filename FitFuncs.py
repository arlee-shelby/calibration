import numpy as np

def get_pixels(nab_file):
    return np.unique(nab_file.singleWaves().headers()['pixel'].to_numpy())

def get_peak(results,x1,x2):
    hist, bins = np.histogram(results.data()['energy'], bins = np.arange(x1, x2))

    return max(hist), list(hist).index(max(hist))+x1

# def get_xray(results):
#     histogram_lower, bin_edges_lower = np.histogram(results.data()['energy'], bins = np.arange(27, 50))
    
#     return max(histogram_lower), list(histogram_lower).index(max(histogram_lower))+27




