import numpy as np
from config import conf
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
import FitFuncs

class Fitter(Source):
    def __init__:
		self.capture = Source.capture
		self.xray = Source.xray
        
		if self.xray=='OFF':
			self.amp = Source.CEamp
			self.peak = Source.CEpeak

		if self.capture=='OFF':
			self.amp = Source.Xamp
			self.peak = Source.Xpeak
    
	def get_bckgrd(self,x,pars,reg=None):
		y = np.zeros(x.shape)
        
		if self.capture == 'zero':
			return FitFuncs.line1(x,pars[:2])
        
		if self.capture !='OFF':
			extrap_reg = np.logical_and(x>reg[0],x<reg[1])
			y[extrap_reg] = (FitFuncs.line2(x[extrap_reg],pars[1:])-FitFuncs.line1(x[extrap_reg],pars[:2]))/len(extrap_reg)
			y[x<reg[0]] = FitFuncs.line1(x[x<reg[0]],pars[:2])
			y[x>reg[1]] = FitFuncs.line2(x[x>reg[1]],pars[1:])
			return UnivariateSpline(x, y, k=1)(x)
        
		elif self.xray!='OFF':
			return FitFuncs.poly(x,pars)
        
        
	def get_CEfit(self,x,*pars):
		reg = {}
		if self.capture =='three':
			reg[0] = pars[1]-pars[2]
			reg[1] = (pars[4]*self.CEpeak)+pars[5]
			return FitFuncs.gaussian(x,pars[0:3]) + FitFuncs.double_gaus(x,pars[3:6],amp=self.amp,peak=self.peak) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg=reg)

		if self.capture =='two':
			reg[0] = pars[1]-pars[2]
			reg[1] =pars[4]+pars[5]            
			return FitFuncs.gaussian(x,pars[0:3]) + FitFuncs.gaussian(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg=reg)

		if self.capture =='one':
			reg[0] = pars[1]- pars[2]
			reg[1] =pars[1]+ pars[2]
			return FitFuncs.gaussian(x,pars[0:3]) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg=reg)

		if self.capture =='zero':
			return 0.0 + self.get_bckgrd(x,pars[-3:])

	def get_Xfit(self,x,*pars):
		if self.xray=='five':
			return FitFuncs.threshold(x,pars[0:3]) + FitFuncs.gaussian(x,pars[3:6]) + FitFuncs.gaussian(x,pars[6:9]) + FitFuncs.double_gaus(x,pars[9:12],amp=self.amp,peak=self.peak) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='four':
			return FitFuncs.threshold(x,pars[0:3]) + FitFuncs.gaussian(x,pars[3:6]) + FitFuncs.double_gaus(x,pars[6:9],amp=self.amp,peak=self.peak) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='three':
			return FitFuncs.threshold(x,pars[0:3]) + FitFuncs.double_gaus(x,pars[3:6],amp=self.amp,peak=self.peak) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='zero':
			return FitFuncs.threshold(x,pars[0:3]) + FitFuncs.gaussian(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:])
