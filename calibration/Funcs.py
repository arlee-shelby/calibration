import numpy as np

#get list of all pixels
def get_pixels(nab_file):
	return np.unique(nab_file.singleWaves().headers()['pixel'].to_numpy())

#get peak value and index within range x1 to x2
def get_peak(results,x1,x2):
	hist, bins = np.histogram(results.data()['energy'], bins = np.arange(x1, x2))

	return max(hist), list(hist).index(max(hist))+x1

#define threshold gaussian
def threshold(x,pars):
	return pars[0]*np.exp(-1.0 * (x/pars[2])**2.0)

#define background polynomial for low lying peaks
def poly(x,pars):
	return pars[0] + pars[1]*x +pars[2]*x**2

#define gaussian
def gaussian(x,pars):
	return pars[0]*np.exp(-1.0*((x-pars[1])/pars[2])**2.0)

#define higher lying peak background one
def line1(x,pars):
	return pars[1]*x + pars[0]

#define higher lying peak background 2
def line2(x,pars):
	return -pars[0]*x + pars[1]

#define double gaussian with peak and apmlitude ratio options
def double_gaus(x,pars,amp=1.0,peak=1.0):
	return pars[0]*np.exp(-1.0 * ((x-pars[1])/np.abs(pars[2]))**2.0) + (pars[0]*amp)*np.exp(-1.0 * ((x-(pars[1]*peak))/np.abs(pars[2]))**2.0)

def get_counts(results,x1,x2):
	hist, bins = np.histogram(results.data()['energy'], bins = np.arange(x1, x2))
	counts = len(np.nonzero(hist)[0])

	return counts

def results(run, n, rise, length, decay):
	run.singleWaves().resetCuts()
	run.singleWaves().defineCut('pixel', '=', n)
	results = run.singleWaves().determineEnergyTiming(method='trap', params=[rise, length, decay])
	return results
