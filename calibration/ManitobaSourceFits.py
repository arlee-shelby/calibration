#import libraries and imports, need Funcs, new_spectra_peak, and conf to be in same location as this file
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import sys
import FitFuncs
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
		if 'output' in line:
			out_path = line.split()[1]


	print(nab_path, slowPath, CdBool, SnBool,rise,length,decay,pixel_list,pixel_list[0])
	file.close()

run_number = args['run']
directory = args['data']
out_put = args['output']

#import nabPy
sys.path.append(nab_path)
import nabPy as Nab

print(run_number)

#get run
run = Nab.DataRun(directory, run_number, ignoreEventFile = True)
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
	results = FitFuncs.results(run, i, int(rise), int(length), int(decay))
	print('got results')

	if SnBool=='True':
		#generic initialization of fit class
		conf['capture'] = ''
		conf['xray'] = ''
		Sn = SnCalibration()
		print('Sn')
		#determine counts in CE and X-ray regions
		CEcounts = FitFuncs.get_counts(results,Sn.CE1[0], Sn.CE1[1])
		Xcounts = FitFuncs.get_counts(results,Sn.X1[0], Sn.X1[1])

		#determine general peak amplitudes and locations for guesses for the fit
		##1 guesses 363keV CE peak, #2 guesses 387keV CE peak
		CEpeak1, CEcenter1 = FitFuncs.get_peak(results,Sn.CE1[0],Sn.CE1[1])
		CEpeak2, CEcenter2 = FitFuncs.get_peak(results,Sn.CE2[0],Sn.CE2[1])

		#threshold for low energy spectrum is cut off, need to determine location of start of histogram such that the threshold centered around zero is fit properly
		thresh_peak, thresh_start = FitFuncs.get_peak(results,Sn.X1[0],Sn.X1[1])

		#determine general peak amplitudes and locations for two peaks seen between threshold and xray peaks
		peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
		peak2, center2 = FitFuncs.get_peak(results,20,22)

		#determine xray peak amplidude and location guess
		Xpeak, Xcenter = FitFuncs.get_peak(results,Sn.X2[0],Sn.X2[1])

		#define function to evaluate fit, need to initialize the SnCalibration class with conf['capture'] and conf['xray']
		#also define the parameters used for the CE peak fits based on counts in the histograms
		#raise flag if it fails
		def CEevaluate(SN, conf, CEpeak1,CEcenter1,CEpeak2, CEcenter2,i):
			conf['xray'] = 'OFF'
			bins = np.arange(SN.CE1[0],SN.CE1[1])
			print(CEpeak1,CEpeak2)
			if CEpeak2>=7:
				conf['capture'] = 'three'
				CEfit.append(3)
				pars = [CEpeak1, CEcenter1, 5, CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			elif CEpeak2>=6:
				conf['capture'] = 'two'
				CEfit.append(2)
				print('CE two, pixel:%d'%i)
				pars = [CEpeak1, CEcenter1, 5, CEpeak2, CEcenter2, 5, 1, 2, 1e-8, 1]

			elif CEpeak2<=5 and CEpeak1>CEpeak2 and CEpeak1>2:
				conf['capture'] = 'one'
				CEfit.append(1)
				print('CE one, pixel:%d'%i)
				pars = [CEpeak1, CEcenter1, 6, 1, 2, 1e-8, 1]

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
				histogram, bin_edges = np.histogram(results.data()['energy'], bins = bins) 
				return histogram,np.zeros(len(pars)),0,0


		histogram, parameters, chi2, errors = CEevaluate(Sn, conf, CEpeak1, CEcenter1, CEpeak2, CEcenter2, i)
		print(chi2)

		#record fit results and histogram
		CEspectrum.append(histogram)
		ecap_list.append(parameters)
		chi2_ecap.append(chi2)

		#Do same thing for the xray peaks, initialization based on amplitudes of peaks rather than counts
		def Xevaluate(SN, conf, thresh_start, thresh_peak,peak2, center2, peak1,center1,Xpeak, Xcenter, i):
			conf['capture'] = 'OFF'
# 			bins = np.arange(thresh_start+2,SN.X1[1])
			print(peak2,Xpeak,Xcenter,thresh_start,peak1)
            
			if peak1>306 and Xpeak>48:
				conf['xray'] = 'five'
				Xfit.append(5)
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
				bins = np.arange(thresh_start,SN.X1[1])
				print('X five, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]
        
			elif peak2<15 and Xpeak>18 or peak2/Xpeak<0.14:
				conf['xray'] = 'three'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				bins = np.arange(thresh_start+4,SN.X1[1])
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
				pars = [thresh_peak+400, 0, thresh_start, Xpeak, Xcenter, 6, 1, 1, 3, 5]
    
# 			elif peak2<15 and Xpeak<=20:
# 				conf['xray'] = 'zero'
# 				Xfit.append(2)
# 				print('X zero, pixel:%d'%i)
# # 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
# 				pars = [thresh_peak+400, 0, thresh_start+1, Xpeak, Xcenter, 3, 1, 1, 3, 5]
    
			elif peak2<=26 and Xpeak>=19 or peak2/Xpeak<=0.43 or thresh_start>10.0:
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+3,thresh_start+5)
				bins = np.arange(thresh_start,SN.X1[1])
				conf['xray'] = 'four'
				Xfit.append(4)
				print('X four, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start+1, peak1, center1, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]
                
			elif peak2>=10 and Xpeak>11 and peak2/Xpeak>0.60:
				conf['xray'] = 'five'
				bins = np.arange(thresh_start,SN.X1[1])
				Xfit.append(5)
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
				print('X five, pixel:%d'%i)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]

			elif peak1>10 and Xpeak>11:
				bins = np.arange(thresh_start+4,SN.X1[1])
				conf['xray'] = 'three'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
				pars = [thresh_peak+400, 0, thresh_start, Xpeak, Xcenter, 4, 10, 1, 3, 5]
                
			else:
				bins = np.arange(thresh_start,SN.X1[1])
				conf['xray'] = 'zero'
				Xfit.append(0)
				print('X zero, pixel:%d'%i)
# 				peak1,center1 = FitFuncs.get_peak(results,thresh_start+5,thresh_start+7)
				pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, 10, 1, 3, 5]

			Sn  = SnCalibration()
			try:
				return Sn.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed lower fits')
				histogram, bin_edges = np.histogram(results.data()['energy'], bins = bins)              
				return histogram,np.zeros(len(pars)),0,0


		histogram, parameters, chi2, errors = Xevaluate(Sn, conf, thresh_start, thresh_peak, peak2, center2, peak1,center1,Xpeak, Xcenter, i)
		print(chi2)

		Xspectrum.append(histogram)
		xray_list.append(parameters)
		chi2_xray.append(chi2)

	elif CdBool=='True':
		conf['capture1'] = ''
		conf['capture2'] = ''
		conf['xray'] = ''
		Cd = CdCalibration()
		print('Cd')

		#determine general peak amplitudes and locations for guesses for the fit
		CEpeak1, CEcenter1 = FitFuncs.get_peak(results,Cd.CE1[0],Cd.CE1[1])
		CEpeak2, CEcenter2 = FitFuncs.get_peak(results,Cd.CE2[0],Cd.CE2[1])

		#threshold for low energy spectrum is cut off, need to determine location of start of histogram such that the threshold centered around zero is fit properly
		thresh_peak, thresh_start = FitFuncs.get_peak(results,Cd.X1[0],Cd.X1[1])

		# #determine general peak amplitudes and locations for two peaks seen between threshold and xray peaks
		# peak1,center1 = Funcs.get_peak(results,thresh_start+5,thresh_start+7)
		# peak2, center2 = Funcs.get_peak(results,20,28)

		#determine xray peak amplidude and location guess
		Xpeak, Xcenter = FitFuncs.get_peak(results,Cd.X2[0],Cd.X2[1])

		CEfit_list = []
		histogram_list = []
		parameters_list = []
		chi2e = []

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
				return np.zeros(len(bins)-1),np.zeros(len(pars)),0,0

		histogram1, parameters1, chi21, errors1 = CEevaluate1(Cd, conf, CEpeak1, CEcenter1, i)
		print(chi21)

		#record fit results and histogram

		histogram_list.append(histogram1.tolist())
		parameters_list.append(parameters1.tolist())
		chi2e.append(float(chi21))

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
				return np.zeros(len(bins)-1),np.zeros(len(pars)),0,0

		histogram2, parameters2, chi22, errors2 = CEevaluate2(Cd, conf, CEpeak2, CEcenter2, i)
		print(chi22)

		#record fit results and histogram

		histogram_list.append(histogram2.tolist())
		parameters_list.append(parameters2.tolist())
		chi2e.append(float(chi22))

		CEfit.append(CEfit_list)
		CEspectrum.append(histogram_list)
		ecap_list.append(parameters_list)
		chi2_ecap.append(chi2e)

		def Xevaluate(CD, conf, thresh_start, thresh_peak, Xpeak, Xcenter, i):
			conf['capture1'] = 'OFF'
			conf['capture2']= 'OFF'
			bins = np.arange(thresh_start,CD.X1[1])
			if Xpeak>=50:
				conf['xray'] = 'three'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				pars = [thresh_peak+200, 0, thresh_start+1, Xpeak, Xcenter, 5, 1, 1, 3, 5]

			elif Xpeak>=20:
				conf['xray'] = 'two'
				Xfit.append(3)
				print('X three, pixel:%d'%i)
				pars = [thresh_peak+200, 0, thresh_start+1, Xpeak, Xcenter, 5, 1, 1, 3, 5]

			else:
				conf['xray'] = 'zero'
				Xfit.append(0)
				print('X zero, pixel:%d'%i)
				pars = [thresh_peak+200, 0, thresh_start, 1, 1, 3, 5]

			Cd = CdCalibration()
			try:
				return Cd.fitter(results,bins,pars)

			except Exception as e:
				print(e,i,'failed X peak')
				return np.zeros(len(bins)-1),np.zeros(len(pars)),0,0

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

print(len(CEspectrum),len(Xspectrum),len(CEfit),len(Xfit),len(ecap_list),len(xray_list))

d = pd.DataFrame(df)
d.to_csv('%s%s%d.csv'%(out_path,out_put,run_number),mode = 'w', header = True, index = False)




