This repository is intended to obtain fits to energy spectra for the Manitoba data. It is currently outfitted for the Sn113 calibration data runs. The main script is ManitobaSourceFits.py which will output a csv file with the fit results, histogramed data, and run information. The Sncalibration_analysis.ipynb notebook shows an example of how to read the data, access the relvant fit parameters, verify the fit quality, calculate a calibration curve, and resolution analysis. 

# Set Up
There are a two things you need to specify/change before you can run the script:
- Modify the input_config.txt:
  - Change the first line to your corresponding pyNab src file location 
  - Change the second line to the location of the manitobametadata_filters.csv file (which you can obtain from this repo)
# Run
There are three options you need to provide to be able to run the script correctly:
  - The -d option specifies your data location
  - The -r option specifies the run number and should be changed for your purposes (note: the default is set to 1374)
  - The -c option specifies the location of the input_config.txt file
