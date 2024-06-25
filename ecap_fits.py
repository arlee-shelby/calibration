import numpy as np

def gaussian(x,*pars):
    return pars[0]*np.exp(-1.0 * ((x-pars[1])/pars[2])**2.0)
def double_gaus(x,*pars):
    return pars[0]*np.exp(-1.0 * ((x-pars[1])/np.abs(pars[2]))**2.0) + (pars[0]*1.137/5.60)*np.exp(-1.0 * ((x-(pars[1]*390.872/387.461))/np.abs(pars[2]))**2.0) 

def line1(x,*pars):
    return pars[0]*x + pars[1]
def line2(x,*pars):
    return -pars[0]*x + pars[1]
def gaus(x,*pars):
    offset = pars[-1]
    
    y = np.zeros(x.shape)
    
    gaus_loc = np.logical_and(x>pars[1]-pars[2],x<(pars[4]*390.872/387.461)+pars[5])

    y[x<pars[1]-pars[2]] = line1(x[x<pars[1]-pars[2]],*pars[6:8])

    y[x>(pars[4]*390.872/387.461)+pars[5]] = line2(x[x>(pars[4]*390.872/387.461)+pars[5]],*pars[7:9])

    y[gaus_loc] = (line2(x[gaus_loc],*pars[7:9])-line1(x[gaus_loc],*pars[6:8]))/len(gaus_loc)

    gaus = double_gaus(x,*pars[3:6])+ gaussian(x,*pars[0:3]) + offset
    
    
    return gaus +y