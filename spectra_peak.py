import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import warnings
import sys
import FitFuncs
import ecap_fits
import xray_fits
import csv
import pandas as pd

warnings.simplefilter('ignore')

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-n', '--nabPy', default='/home', help='Path to nabPy')
parser.add_argument('-p', '--path', default='/home', help='Path to data directory')
parser.add_argument('-r', '--run', default=1374,type=int, help='Run number')
parser.add_argument('-pix', '--pixel', default='None', help='Pixel list')
parser.add_argument('-o', '--output', default='output.csv', help='Output file name')
parser.add_argument('-s','--slow',default='None',help='Slow data csv file path')


print('working')
args = vars(parser.parse_args())

nab_path = args['nabPy']
run_number = args['run']
directory = args['path']
out_put = args['output']

sys.path.append('nab_path')
import nabPy as Nab

print(args['pixel'])

run = Nab.DataRun(directory, run_number, ignoreEventFile = True)
print('got run')

if args['pixel']==None:
    pixel_list = get_pixels(run)
else:
    file = open(args['pixel'],'r')
    pixel_list = np.int_(file.read().split(',')).tolist()
    print(pixel_list)


# with open(out_put, 'w') as file:
# headers_ecap = ['amp','center','sigma','slope','left intercept','right intercept','offset']
# headers_xray = ['1 amp','center','sigma','2',' ',' ','3 ',' ',' ','4',' ',' ','5 amp','sigma','quadradic','linear','constant','offset']

# file.write()
    # writer = csv.writer(file,delimiter = ' ')
df = {}
run_list = []
ecap_list = []
xray_list = []
chi_2_ecap = []
chi_2_xray = []

for i in pixel_list:

    run.singleWaves().resetCuts()
    run.singleWaves().defineCut('pixel', '=', i)

    results = run.singleWaves().determineEnergyTiming(method='trap', params=[1250, 50, 1250])

    # print(results.data()['pixel'])

    histogram, bin_edges = np.histogram(results.data()['energy'], bins = np.arange(300, 1400))
    width = bin_edges[1]-bin_edges[0]

    peak, center = FitFuncs.get_peak(results,620,700)

    thresh_peak, thresh_start = FitFuncs.get_peak(results,0,150)
    xray_peak, xray_center = FitFuncs.get_peak(results,27,50)

    histogram_lower, bin_edges_lower = np.histogram(results.data()['energy'], bins = np.arange(thresh_start, 200))
    width_lower = bin_edges_lower[1]-bin_edges_lower[0]

    run_list.append(run_number)

    try:
        parameters, errors = curve_fit(ecap_fits.gaus, bin_edges[:-1]+width/2, histogram,p0 = [60,598,3,peak,center,3,2,1e-8,1,1])
        print('got upper fit',i)
        chi2 = sum((histogram-ecap_fits.gaus(bin_edges[:-1]+width/2,*parameters))**2)/(len(bin_edges[:-1]+width/2)-len(parameters))
        print(chi2)
        print(parameters)

        # file.write('%.3f\t'%chi2)

        # for j in parameters[:-1]:
        #     file.write('%.3f\t'%j)

        # file.write('%.3f\n'%parameters[-1])

        ecap_list.append(parameters.tolist())
        chi_2_ecap.append(chi2)

        # file.write('%.3f\n'%chi2_lower)

        
        # for k in parameters_lower[:-1]:
        #     file.write('%.3f\t'%k)

        # file.write('%.3f\n'%parameters_lower[-1])
        
        # print(parameters_lower)

    except Exception as e:
        print(e,i,'failed upper fits')
        ecap_list.append(0)
        chi_2_ecap.append(0)

        pass

    try:

        parameters_lower, errors_lower = curve_fit(xray_fits.gaus_poly, bin_edges_lower[:-1]+width_lower/2, histogram_lower, p0 = [100000,0,thresh_start+4,xray_peak,xray_center,5,600,thresh_start+4,4.41,300,22,4,300,5,3,3])
        print('got lower fit')

        chi2_lower = sum((histogram_lower-xray_fits.gaus_poly(bin_edges_lower[:-1]+width_lower/2,*parameters_lower))**2)/(len(bin_edges_lower[:-1]+width_lower/2)-len(parameters_lower))
        print(chi2_lower)
        print(parameters_lower)

        xray_list.append(parameters_lower.tolist())
        chi_2_xray.append(chi2_lower)

    except Exception as e:
        print(e,i,'failed lower fit, 1')
        try:
            parameters_lower, errors_lower = curve_fit(xray_fits.gaus_poly, bin_edges_lower[:-1]+width_lower/2, histogram_lower, p0 = [100000,0,thresh_start+4,xray_peak,xray_center,5,600,thresh_start+4,3.41,300,22,5,300,5,3,3])
            print('got lower fit, try 2')
            xray_list.append(parameters_lower.tolist())
            chi_2_xray.append(chi2_lower)

        except Exception as e:
            print(e,i,'failed lower fit, 2')
            xray_list.append(0)
            chi_2_xray.append(0)
            pass

        pass



df['run'] = run_list
df['pixel'] = pixel_list

if args['slow']== None:
    pass
else:
    slow_df = pd.read_table(args['slow'],delimiter = '|')
    b_v = np.array(slow_df['Detector Bias Voltage [V]'][slow_df['RunID']==1374])[0]
    df['Bias Voltage'] = [b_v]*len(pixel_list)

    for i in range(10,15):
        val = np.array(slow_df[slow_df.columns[i]][slow_df['RunID']==1374])[0]
        df[slow_df.columns[i]] = [val]*len(pixel_list)

df['ecap'] = ecap_list
df['chi2_e'] = chi_2_ecap
df['xray'] = xray_list
df['chi2_x'] = chi_2_xray

d = pd.DataFrame(df)
d.to_csv(out_put,mode = 'w', header = True, index = False)

# file.close()

        #for txt file write

            # parameters, errors = curve_fit(ecap_fits.gaus, bin_edges[:-1]+width/2, histogram,p0 = [60,598,3,peak,center,3,2,1e-8,1,1])
            # print('got upper fit',i)
            # chi2 = sum((histogram-ecap_fits.gaus(bin_edges[:-1]+width/2,*parameters))**2)/(len(bin_edges[:-1]+width/2)-len(parameters))
            # print(chi2)

            # # file.write('%.3f\t'%chi2)

            # # for j in parameters[:-1]:
            # #     file.write('%.3f\t'%j)

            # # file.write('%.3f\n'%parameters[-1])

            # print(parameters)
            # parameters_lower, errors_lower = curve_fit(xray_fits.gaus_poly, bin_edges_lower[:-1]+width_lower/2, histogram_lower, p0 = [100000,0,thresh_start+4,xray_peak,xray_center,5,600,thresh_start+4,4.41,300,22,4,300,5,3,3])

            # chi2_lower = sum((histogram_lower-xray_fits.gaus_poly(bin_edges_lower[:-1]+width_lower/2,*parameters_lower))**2)/(len(bin_edges_lower[:-1]+width_lower/2)-len(parameters_lower))
            # file.write('%.3f\n'%chi2_lower)

            
            # # for k in parameters_lower[:-1]:
            # #     file.write('%.3f\t'%k)

            # # file.write('%.3f\n'%parameters_lower[-1])
            
            # # print(parameters_lower)



