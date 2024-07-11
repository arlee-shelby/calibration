import numpy as np
from scipy.optimize import curve_fit

def get_peak(results,x1,x2):
	hist, bins = np.histogram(results.data()['energy'], bins = np.arange(x1, x2))

	return max(hist), list(hist).index(max(hist))+x1

def get_counts(results,x1,x2):
	hist, bins = np.histogram(results.data()['energy'], bins = np.arange(x1, x2))
	counts = len(np.nonzero(hist)[0])

	return counts

def results(run, n, rise, length, decay):
	run.singleWaves().resetCuts()
	run.singleWaves().defineCut('pixel', '=', n)
	results = run.singleWaves().determineEnergyTiming(method='trap', params=[rise, length, decay])
	return results
