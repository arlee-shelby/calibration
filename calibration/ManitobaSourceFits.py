#import libraries and imports, need Funcs, new_spectra_peak, and conf to be in same location as this file
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import sys
import Funcs
import csv
import pandas as pd
from FitClass import SnCalibration
from FitClass import CdCalibration
from config import conf

np.set_printoptions(threshold=sys.maxsize)

warnings.simplefilter('ignore')

#define functionality that can change during implementation

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--data', default='/home', help='Path to data directory')
parser.add_argument('-r', '--run', default=1374,type=int, help='Run number')
parser.add_argument('-p', '--pixels', default='None', help='pixels to analyze')
parser.add_argument('-o', '--output', default='output', help='Output file name')
parser.add_argument('-c','--config',default='None',help='Input configuration file path')

args = vars(parser.parse_args())

print('working')

#get configuration of the data from the --config file, should have the pyNab file location, information on the source, pixels, and trap filter parameters

with open(args['config'],'r') as file:
	for line in file:
		if 'nabPy' in line:
			nab_path = line.split()[1]
		if 'slowData' in line:
			slowPath= line.split()[1]
		if 'Cd109' in line:
			CdBool = line.split()[1]
		if 'Sn113' in line:
			SnBool = line.split()[1]
		if 'TrapPars' in line:
			rise, length, decay = line.split()[1:]
		if 'pixels' in line and args['pixels']=='None':
			pixel_list = np.int_(line.split()[1].split(','))
			print(pixel_list)
		if args['pixels']!='None':
			pixel_list = args['pixels']


	print(nab_path, slowPath, CdBool,SnBool,rise,length,decay,pixel_list,pixel_list[0])
	file.close()

run_number = args['run']
directory = args['data']
out_put = args['output']

#import nabPy
sys.path.append(nab_path)
import nabPy as Nab

print(run_number)

#get run
run = Nab.DataRun(directory, run_number, ignoreEventFile = True,subRunMin=0, subRunMax=2)
print('got run')

#set up datafame functionality for outputcsv file
df = {}

run_list = []

ecap_list = []
xray_list = []

chi2_ecap = []
chi2_xray = []

CEspectrum = []
Xspectrum = []

trap_rise = []
trap_length = []
trap_decay = []

CEfit = []
Xfit = []

for i in pixel_list:
	#record run number and trap filter parameters
	run_list.append(run_number)
	trap_rise.append(int(rise))
	trap_length.append(int(length))
	trap_decay.append(int(decay))

	#get results
	results = Funcs.results(run, i, int(rise), int(length), int(decay))
	print('got results')

	if SnBool=='True':
		#generic initialization of fit class
		conf['capture'] = ''
		conf['xray'] = ''
		Sn = SnCalibration()

		#determine counts in CE and X-ray regions
		CEcounts = Funcs.get_counts(results,Sn.CE1[0], Sn.CE1[1])
		Xcounts = Funcs.get_counts(results,Sn.X1[0], Sn.X1[1])s

		#determine general peak amplitudes and locations for guesses for the fit
		##1 guesses 363keV CE peak, #2 guesses 387keV CE peak
		CEpeak1, CEcenter1 = Funcs.get_peak(results,Sn.CE1[0],Sn.CE1[1])
		CEpeak2, CEcenter2 = Funcs.get_peak(results,Sn.CE2[0],Sn.CE2[1])

		#threshold for low energy spectrum is cut off, need to determine location of start of histogram such that the threshold centered around zero is fit properly
		thresh_peak, thresh_start = Funcs.get_peak(results,Sn.X1[0],Sn.X1[1])

		#determine general peak amplitudes and locations for two peaks seen between threshold and xray peaks
		peak1,center1 = Funcs.get_peak(results,thresh_start+5,thresh_start+7)
		peak2, center2 = Funcs.get_peak(results,20,28)

		#determine xray peak amplidude and location guess
		Xpeak, Xcenter = Funcs.get_peak(results,Sn.X2[0],Sn.X2[1])

		#define function to evaluate fit, need to initialize the SnCalibration class with conf['capture'] and conf['xray']
		#also define the parameters used for the CE peak fits based on counts in the histograms
		#raise flag if it fails
		def CEevaluate(SN, conf, CEpeak1,CEcenter1,CEpeak2, CEcenter2,i):
			conf['xray'] = 'OFF'
			bins = np.arange(SN.CE1[0],SN.CE1[1])

			if CEpeak2>=7:
				conf['capture'] = 'three'
				CEfit.append(3)
				pars = [CEpeak1, CEcenter1, 5, CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			elif CEpeak2>=6:
				conf['capture'] = 'two'
				CEfit.append(2)
				print('CE two, pixel:%d'%i)
				pars = [CEpeak1, CEcenter1, 5, CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			elif CEpeak2<=5 and CEpeak1>CEpeak2:
				conf['capture'] = 'one'
				CEfit.append(1)
				print('CE one, pixel:%d'%i)
				pars = [CEpeak1, CEcenter1, 5, 1, 2, 1e-8, 1]

			else:
				conf['capture'] = 'zero'
				CEfit.append(0)
				print('CE zero, not enought counts pixel:%d'%i)
				pars = [2, 1e-8, 1]
				

			Sn  = SnCalibration()
			try:
				return Sn.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed upper fits')
				return 0,0,0,0


		histogram, parameters, chi2, errors = CEevaluate(Sn, conf, CEpeak1, CEcenter1, CEpeak2, CEcenter2, i)
		print(chi2)

		#record fit results and histogram
		CEspectrum.append(histogram)
		ecap_list.append(parameters)
		chi2_ecap.append(chi2)

		#Do same thing for the xray peaks, initialization based on amplitudes of peaks rather than counts
		def Xevaluate(SN, conf, thresh_start, thresh_peak, peak1, center1, peak2, center2, Xpeak, Xcenter, i):
			conf['capture'] = 'OFF'
			bins = np.arange(thresh_start,SN.X1[1])

			if peak2>10 and Xpeak>10:
				conf['xray'] = 'five'
				Xfit.append(5)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]

			elif peak2<10 and Xpeak>10:
				conf['xray'] = 'four'
				Xfit.append(4)
				print('X four, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start+1, peak1, center1, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]

			elif peak2>=10 and Xpeak<=10:
				conf['xray'] = 'three'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, 10, 1, 3, 5]

			elif peak2<=10 and Xpeak<=10:
				conf['xray'] = 'zero'
				Xfit.append(0)
				print('X zero, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, 10, 1, 3, 5]

			Sn  = SnCalibration()
			try:
				return Sn.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed lower fits')
				return 0,0,0,0


		histogram, parameters, chi2, errors = Xevaluate(Sn, conf, thresh_start, thresh_peak, peak1, center1, peak2, center2, Xpeak, Xcenter, i)
		print(chi2)

		Xspectrum.append(histogram)
		xray_list.append(parameters)
		chi2_xray.append(chi2)

	if CdBool=='True':
		conf['capture1'] = ''
		conf['capture2'] = ''
		conf['xray'] = ''
		Cd = CdCalibration()

		#determine general peak amplitudes and locations for guesses for the fit
		CEpeak1, CEcenter1 = Funcs.get_peak(results,Cd.CE1[0],Cd.CE1[1])
		CEpeak2, CEcenter2 = Funcs.get_peak(results,Cd.CE2[0],Cd.CE2[1])

		#threshold for low energy spectrum is cut off, need to determine location of start of histogram such that the threshold centered around zero is fit properly
		thresh_peak, thresh_start = Funcs.get_peak(results,Cd.X1[0],Cd.X1[1])

		# #determine general peak amplitudes and locations for two peaks seen between threshold and xray peaks
		# peak1,center1 = Funcs.get_peak(results,thresh_start+5,thresh_start+7)
		# peak2, center2 = Funcs.get_peak(results,20,28)

		#determine xray peak amplidude and location guess
		Xpeak, Xcenter = Funcs.get_peak(results,Cd.X2[0],Cd.X2[1])

		CEfit_list = []

		def CEevaluate1(CD, conf, CEpeak1,CEcenter1,i):
			conf['xray'] = 'OFF'
			conf['capture2']= 'OFF'

			bins = np.arange(CD.CE1[0],CD.CE1[1])

			if CEpeak1>5:
				conf['capture1'] = 'ON'
				CEfit_list.append(1)
				print('CE peak1 good, pixel:%d'%i)
				pars = [CEpeak1, CEcenter1, 5, 1, 2, 1e-8, 1]
			else:
				conf['capture1'] = 'zero'
				CEfit_list.append(0)
				print('CE peak1 bad, pixel:%d'%i)
				pars = [2, 1e-8, 1]

			Cd = CdCalibration()
			try:
				return Cd.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed CE peak 1')
				return 0,0,0,0

		histogram, parameters, chi2, errors = CEevaluate1(Cd, conf, CEpeak1, CEcenter1, i)
		print(chi2)

		#record fit results and histogram

		CEspectrum.append(histogram)
		ecap_list.append(parameters)
		chi2_ecap.append(chi2)

		def CEevaluate2(CD, conf, CEpeak2, CEcenter2,i):
			conf['xray'] = 'OFF'
			conf['capture1']= 'OFF'
			bins = np.arange(CD.CE2[0],CD.CE2[1])

			if CEpeak2>=7:
				conf['capture2'] = 'two'
				CEfit_list.append(2)
				print('CE peak2 two, pixel:%d'%i)
				pars = [CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			elif CEpeak2>=5:
				conf['capture2'] = 'one'
				CEfit_list.append(1)
				print('CE peak2 one, pixel:%d'%i)
				pars = [CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			else:
				conf['capture2'] = 'zero'
				CEfit_list.append(0)
				print('CE peak2 zero, pixel:%d'%i)
				pars = [2, 1e-8, 1]

			Cd = CdCalibration()
			try:
				return Cd.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed CE peak 2')
				return 0,0,0,0

		histogram, parameters, chi2, errors = CEevaluate2(Cd, conf, CEpeak2, CEcenter2, i)
		print(chi2)

		#record fit results and histogram

		CEspectrum.append(histogram)
		ecap_list.append(parameters)
		chi2_ecap.append(chi2)

		CEfit.append(CEfit_list)

		def Xevaluate(CD, conf, thresh_start, thresh_peak, Xpeak, Xcenter, i):
			conf['capture1'] = 'OFF'
			conf['capture2']= 'OFF'
			bins = np.arange(thresh_start,CD.X1[1])
			if Xpeak>=7:
				conf['xray'] = 'three'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start+1, Xpeak, Xcenter, 4, 10, 1, 3, 5]

			elif Xpeak>=5:
				conf['xray'] = 'two'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start+1, Xpeak, Xcenter, 4, 10, 1, 3, 5]

			else:
				conf['xray'] = 'zero'
				Xfit.append(0)
				print('X zero, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start,10, 1, 3, 5]

			Cd = CdCalibration()
			try:
				return Cd.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed X peak')
				return 0,0,0,0

		histogram, parameters, chi2, errors = Xevaluate(Cd, conf, thresh_start, thresh_peak, Xpeak, Xcenter, i)
		print(chi2)

		Xspectrum.append(histogram)
		xray_list.append(parameters)
		chi2_xray.append(chi2)



#construct pandas DataFrame to store data to csv file
df['run'] = run_list
df['pixel'] = pixel_list

#record slow data
slow_df = pd.read_table(slowPath,delimiter = '|')
b_v = np.array(slow_df[slow_df.columns[5]][slow_df['RunID']==run_number])[0]
df['Bias Voltage'] = [b_v]*len(pixel_list)

for i in range(10,15):
	val = np.array(slow_df[slow_df.columns[i]][slow_df['RunID']==run_number])[0]
	df[slow_df.columns[i]] = [val]*len(pixel_list)

#record trap filter parameters
df['trap rise'] = trap_rise
df['trap length'] = trap_length
df['trap decay'] = trap_decay

#record fit parameters and chi2
df['ecap'] = ecap_list
df['chi2_e'] = chi2_ecap
df['xray'] = xray_list
df['chi2_x'] = chi2_xray

#record number gaussians - should be 3 for CE and 5 for xray for pixels with good counts
df['CE'] = CEfit
df['Xray'] = Xfit

#record histogram data for later analysis
df['CE hist'] = CEspectrum
df['Xray hist'] = Xspectrum

d = pd.DataFrame(df)
d.to_csv('%s%d.csv'%(out_put,run_number),mode = 'w', header = True, index = False)




