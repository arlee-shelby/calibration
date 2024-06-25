import numpy as np


def poly(x,*pars):
    return (pars[0] + pars[1]*x +pars[2]*x**2)
def gaussian(x,*pars):
    return pars[0]*np.exp(-1.0 * ((x-pars[1])/pars[2])**2.0)

def threshold(x,*pars):
    return pars[0]*np.exp(-1.0 * ((x)/pars[2])**2.0)

def double_gaus(x,*pars):
    return pars[0]*np.exp(-1.0 * ((x-pars[1])/np.abs(pars[2]))**2.0) + (pars[0]*16.05/79.8)*np.exp(-1.0 * ((x-(pars[1]*27.3523/22.59))/np.abs(pars[2]))**2.0) 

def multi_gaus(x, *pars):
    num_gaus = int((len(pars))/3)
    # print(pars)
    if num_gaus==2:
        g = threshold(x,*pars[0:3]) + double_gaus(x,*pars[3:6])
    if num_gaus==3:
        g = threshold(x,*pars[0:3]) + double_gaus(x,*pars[3:6]) + gaussian(x,*pars[6:9])
    if num_gaus==4:
        g = threshold(x,*pars[0:3]) + double_gaus(x,*pars[3:6]) + gaussian(x,*pars[6:9]) + gaussian(x,*pars[9:12])
    
    offset = pars[-1]
    return g + offset

def gaus_poly(x,*pars):
    list_pars = pars[0]
    return multi_gaus(x,*pars[0:-3])+poly(x,*pars[-3:])