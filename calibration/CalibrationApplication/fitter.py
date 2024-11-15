import numpy as np

class Fitter():
    def threshold(x,pars):
	    return pars[0]*np.exp(-1.0 * (x/pars[1])**2.0)

    def poly(x,pars):
        return pars[0] + pars[1]*x +pars[2]*x**2

    def gaussian(x,pars):
        return pars[0]*np.exp(-1.0*((x-pars[1])/pars[2])**2.0)

    def line1(x,pars):
        return pars[1]*x + pars[0]

    def line2(x,pars):
        return pars[0]*x + pars[1]

    def double_gaus(x,pars,amp=1.0,peak=1.0):
        return pars[0]*np.exp(-1.0 * ((x-pars[1])/np.abs(pars[2]))**2.0) + (pars[0]*amp)*np.exp(-1.0 * ((x-(pars[1]*peak))/np.abs(pars[2]))**2.0)