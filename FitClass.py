import numpy as np
from config import conf
from scipy.optimize import curve_fit

#define Sn source class
class SnCalibration:

	def __init__(self):
		#initialize conf to determine if performing xray or CE fits
		self.capture = conf['capture']
		self.xray = conf['xray']

		#define CE intensity and energy location ratio 
		self.CEamp = 1.137/5.60
		self.CEpeak = 390.872/387.461

		#define histogram range and double gaussian range
		self.CE1 = [300.0,1400.0]
		self.CE2 = [622.0,700.0]

		#define intensity and energy ratios for xray fit
		self.Xamp = 16.05/79.8
		self.Xpeak = 27.3523/22.59

		#define xray fit range and double gaussian range
		self.X1 = [0.0,200.0]
		self.X2 = [27.0,50.0]

	#threshold gaussian
	def threshold(self,x,pars):
		return pars[0]*np.exp(-1.0 * (x/pars[2])**2.0)

	#gaussian
	def gaussian(self,x,pars):
		return pars[0]*np.exp(-1.0*((x-pars[1])/pars[2])**2.0)

	#double gaussian, dependent on initialization (i.e. fitting xray peaks if conf['capture']=='OFF' and fitting CE peaks if conf['xray']=='OFF')
	def double_gaus(self,x,pars):

		if self.xray=='OFF':
			amp = self.CEamp
			peak = self.CEpeak

		if self.capture=='OFF':
			amp = self.Xamp
			peak = self.Xpeak
			
		return pars[0]*np.exp(-1.0 * ((x-pars[1])/np.abs(pars[2]))**2.0) + (pars[0]*amp)*np.exp(-1.0 * ((x-(pars[1]*peak))/np.abs(pars[2]))**2.0)

	def line1(self,x,pars):
		return pars[1]*x + pars[0]

	def line2(self,x,pars):
		return -pars[0]*x + pars[1]

	def poly(self,x,pars):
		return pars[0] + pars[1]*x +pars[2]*x**2

	#define background funtion, dependent on CE or xray fit
	def get_bckgrd(self,x,pars,reg=None):

		y = np.zeros(x.shape)

		if self.capture == 'zero':
			return self.line1(x,pars[:2])

		if self.capture !='OFF':
			extrap_reg = np.logical_and(x>reg[0],x<reg[1])
			y[extrap_reg] = (self.line2(x[extrap_reg],pars[1:])-self.line1(x[extrap_reg],pars[:2]))/len(extrap_reg)
			y[x<reg[0]] = self.line1(x[x<reg[0]],pars[:2])
			y[x>reg[1]] = self.line2(x[x>reg[1]],pars[1:])
			return y

		elif self.capture == 'OFF':
			return self.poly(x,pars)

	#define multiple gaussian function
	def get_fit(self,x,*pars):
		#define various kinds of CE gaussian fits, based on number of gaussians used in the fit, this also determines the background fit range/location
		reg = {}
		if self.capture =='three':
			reg[0] = pars[1]-pars[2]
			reg[1] = (pars[4]*self.CEpeak)+pars[5]
			return self.gaussian(x,pars[0:3]) + self.double_gaus(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg)

		if self.capture =='two':
			reg[0] = pars[1]-pars[2]
			reg[1] =pars[4]+pars[5]
			return self.gaussian(x,pars[0:3]) + self.gaussian(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg)

		if self.capture =='one':
			reg[0] = pars[1]-pars[2]
			reg[1] =pars[1]+pars[2]
			return self.gaussian(x,pars[0:3]) + pars[-4] + self.get_bckgrd(x,pars[-3:],reg)

		if self.capture =='zero':
			return 0.0 + self.get_bckgrd(x,pars[-3:])


		#define xray gaussian fits, based on height of amplitudes of the various peaks seen
		if self.xray=='five':
			return self.threshold(x,pars[0:3]) + self.gaussian(x,pars[3:6]) + self.gaussian(x,pars[6:9]) + self.double_gaus(x,pars[9:12]) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='four':
			return self.threshold(x,pars[0:3]) + self.gaussian(x,pars[3:6]) + self.double_gaus(x,pars[6:9]) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='three':
			return self.threshold(x,pars[0:3]) + self.double_gaus(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:])

		if self.xray=='zero':
			return self.threshold(x,pars[0:3]) + self.gaussian(x,pars[3:6]) + pars[-4] + self.get_bckgrd(x,pars[-3:])


	#do fit
	def fitter(self,results,bins,pars):
		histogram, bin_edges = np.histogram(results.data()['energy'], bins = bins)

		width = bin_edges[1]-bin_edges[0]

		parameters, errors = curve_fit(self.get_fit, bin_edges[:-1]+width/2, histogram, p0 = pars)

		chi2 = sum((histogram-self.get_fit(bin_edges[:-1]+width/2,*parameters))**2)/(len(bin_edges[:-1]+width/2)-len(parameters))

		return histogram, parameters, chi2, errors


