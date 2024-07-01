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
from config import conf

warnings.simplefilter('ignore')

#define functionality that can change during implementation

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--data', default='/home', help='Path to data directory')
parser.add_argument('-r', '--run', default=1374,type=int, help='Run number')
parser.add_argument('-p', '--pixels', default='None', help='pixels to analyze')
parser.add_argument('-o', '--output', default='output', help='Output file name')
parser.add_argument('-c','--config',default='None',help='Input configuration file')

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

for i in pixel_list:
	#record run number and trap filter parameters
	run_list.append(run_number)
	trap_rise.append(int(rise))
	trap_length.append(int(length))
	trap_decay.append(int(decay))

	#get results
	results = Funcs.results(run, i, int(rise), int(length), int(decay))

	#generic initialization of fit class
	conf['capture'] = ''
	conf['xray'] = ''
	Sn = SnCalibration()

	#determine counts in CE and X-ray regions
	CEcounts = Funcs.get_counts(results,Sn.CE1[0], Sn.CE1[1])
	Xcounts = Funcs.get_counts(results,Sn.X1[0], Sn.X1[1])

	#determine general peak amplitudes and locations for guesses for the fit
	##1 guesses 363keV CE peak, #2 guesses 387keV CE peak
	CEpeak1, CEcenter1 = Funcs.get_peak(results,Sn.CE1[0],Sn.CE1[1])
	CEpeak2, CEcenter2 = Funcs.get_peak(results,Sn.CE2[0],Sn.CE2[1])

	#threshold for low energy spectrum is cut off, need to determine location of start of histogram such that the threshold centered around zero is fit properly
	thresh_peak, thresh_start = Funcs.get_peak(results,Sn.X1[0],Sn.X1[1])

	#determine general peak amplitudes and locations for two peaks seen between threshold and xray peaks
	peak1,center1 = Funcs.get_peak(results,thresh_start+3,thresh_start+6)
	peak2, center2 = Funcs.get_peak(results,20,28)

	#determine xray peak amplidude and location guess
	Xpeak, Xcenter = Funcs.get_peak(results,Sn.X2[0],Sn.X2[1])

	#define function to evaluate fit, need to initialize the SnCalibration class with conf['capture'] and conf['xray']
	#also define the parameters used for the CE peak fits based on counts in the histograms
	#raise flag if it fails
	def CEevaluate(SN, conf, CEpeak1,CEcenter1,CEpeak2, CEcenter2,i):
		conf['xray'] = 'OFF'
		bins = np.arange(SN.CE1[0],SN.CE2[1])

		if CEcounts>100:
			conf['capture'] = 'three'
			pars = [CEpeak1, CEcenter1, 3, CEpeak2, CEcenter2, 3, 1, 2, 1e-8, 1]

		if CEcounts<25:
			conf['capture'] = 'two'
			pars = [CEpeak1, CEcenter1, 3, CEpeak2, CEcenter2, 3, 1, 2, 1e-8, 1]

		if CEcounts<10:
			conf['capture'] = 'one'
			pars = [CEpeak1, CEcenter1, 3, 1, 2, 1e-8, 1]

		if CEcounts<4:
			conf['capture'] = 'zero'
			pars = [2, 1e-8, 1]

		Sn  = SnCalibration()
		try:
			return Sn.fitter(results,pars,bins)

		except Exception as e:
			print(e,i,'failed upper fits')
			return 0,0,0,0


	histogram, parameters, chi2, errors = CEevaluate(Sn, conf, CEpeak1, CEcenter1, CEpeak2, CEcenter2, i)

	#record fit results and histogram
	CEspectrum.append(histogram)
	ecap_list.append(parameters.tolist())
	chi2_ecap.append(chi2)

	#Do same thing for the xray peaks, initialization based on amplitudes of peaks rather than counts
	def Xevaluate(SN, conf, thresh_start, thresh_peak, peak1, center1, peak2, center2, Xpeak, Xcenter, i):
		conf['capture'] = 'OFF'
		bins = np.arange(thresh_start,SN.X2[1])

		if peak2>60 and Xpeak>60:
			conf['xray'] = 'five'
			pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, Xpeak, Xcenter, 5, 10, 1, 3, 5]

		if peak2<60 and Xpeak>60:
			conf['xray'] = 'four'
			pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, Xpeak, Xcenter, 5, 10, 1, 3, 5]

		if peak2>60 and Xpeak<60:
			conf['xray'] = 'three'
			pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, peak2, center2, 4, 10, 1, 3, 5]

		if peak2<60 and Xpeak<60:
			conf['xray'] = 'zero'
			pars = [thresh_peak+400, 0, thresh_start, peak1, center1, 3, 10, 1, 3, 5]

		Sn  = SnCalibration()
		try:
			return Sn.fitter(results,pars,bins)

		except Exception as e:
			print(e,i,'failed lower fits')
			return 0,0,0,0


	histogram, parameters, chi2, errors = Xevaluate(Sn, conf, thresh_start, thresh_peak, peak1, center1, peak2, center2, Xpeak, Xcenter, i)

	Xspectrum.append(histogram)
	xray_list.append(parameters.tolist())
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

df['trap rise'] = trap_rise
df['trap length'] = trap_length
df['trap decay'] = trap_decay

df['ecap'] = ecap_list
df['chi2_e'] = chi_2_ecap
df['xray'] = xray_list
df['chi2_x'] = chi_2_xray
df['CE hist'] = CEspectrum
df['Xray hist'] = Xspectrum

d = pd.DataFrame(df)
d.to_csv('%s%d.csv'%(out_put,run_number),mode = 'w', header = True, index = False)




