#%%
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
# from sqlalchemy import MetaData
import json
import pandas as pd
import numpy as np
import re
from datetime import date
from models.base import Session, engine, Base
from models.run import Run
from models.pixel import Pixel
from controllers.run_controller import RunController
from controllers.pixel_controller import PixelController

# Base.metadata.drop_all(engine)
# Base.metadata.create_all(engine)

Base.metadata.create_all(engine)
session = Session(bind=engine)

# run1374 = Run(1374, -300, 'True', 'False', 121.5, 15, date(2021,1,26), 'False', 'False',0)

# run1389 = Run(1389, -320, 'False', 'False', 124.589, 1, date(2022,3,19), 'False', 'False',30)

# e = [4,6,7,8,9,4,5,6,7,78,89,34.5,6.7,33,6]
# json_obj = json.dumps(e)

# pixel76 = Pixel(76, 'good', 2.34, 1.19, 465.99, json_obj,run1374)

# pixel87 = Pixel(87, 'bad', 2.34, 1.19, 465.99, json_obj,run1374)

# session.add(run1374)
# session.add(run1389)
# session.add(pixel76)
# session.add(pixel87)

# session.commit()
# session.close()

# runs = RunController.get_runs()
# print(runs)
# for run in runs:
#     print(run.id,run.run_number)


cwd = os.getcwd()

file_name = "manitobametadata_filters.csv"

file_path = os.path.join(cwd, file_name)

df = pd.read_table(file_path,delimiter = '|')
match = re.match(r"(\d{4})-(\d{2})-(\d{2})",df[df['RunID']==1072]['Date Time [UTC]'].iloc[0])

# print(df)

# print(type(df[df['RunID']==700]['RunID'].iloc[0]))

# for i in range(len(df['RunID'].unique())):
run_number = 1389
bias_voltage = int(df[df['RunID']==run_number]['Detector Bias Voltage [V]'].iloc[0])
sn_source = str(df[df['RunID']==run_number]['Sn113'].iloc[0])
cd_source = str(df[df['RunID']==run_number]['Cd109'].iloc[0])
average_temp = np.mean(df[df['RunID']==run_number]['Detector Armor Temperature [K]'])
number_subruns = len(df[df['RunID']==run_number])

match = re.match(r"(\d{4})-(\d{2})-(\d{2})",df[df['RunID']==run_number]['Date Time [UTC]'].iloc[0])
year, month, day = match.groups()
run_date = date(int(year),int(month),int(day))

pulser = str(df[df['RunID']==run_number]['Pulser'].iloc[0])
proton = str(df[df['RunID']==run_number]['Proton'].iloc[0])
proton_energy = int(df[df['RunID']==run_number]['Proton Energy'].iloc[0])

directory = '/storage/home/hcoda1/4/ashelby8/scratch/ManitobaData/'

run = Run(run_number, bias_voltage, sn_source, cd_source, average_temp, number_subruns, run_date, pulser, proton, proton_energy, directory, True)

# print(run.nab_run.singleWaves().resetCuts().defineCut('pixel', '=', pixel_number))

pixel = Pixel(run, 76, 1250, 50, 1250)

session.add(run)
session.add(pixel)

session.commit()
session.close()

#%%
# import matplotlib.pyplot as plt
# import numpy as np
# from controllers.pixel_controller import PixelController
# from models.base import Session, engine, Base
# # Base.metadata.create_all(engine)

# hist, bins,width =  PixelController.energy_hist(PixelController.get_pixel_by_number(76),np.arange(0,1000))
# PixelController.get_pixel_by_number(76)

# plt.plot(hist)




# %%
