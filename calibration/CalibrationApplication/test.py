import sys
nab_path = '/storage/home/hcoda1/4/ashelby8/Manitoba/pyNab/src'
sys.path.append(nab_path)
import nabPy as Nab

directory = '/storage/home/hcoda1/4/ashelby8/scratch/ManitobaData/'
run_number = 1389
pixel_number = 76
trap_rise = 1250
trap_length = 50
trap_decay = 1250

nab_run = Nab.DataRun(directory, run_number, ignoreEventFile = True)

nab_run.singleWaves().resetCuts()
nab_run.singleWaves().defineCut('pixel', '=', pixel_number)

nab_result = run.singleWaves().determineEnergyTiming(method='trap', params=[trap_rise, trap_length, trap_decay])
