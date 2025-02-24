import sys
nab_path = '/storage/home/hcoda1/4/ashelby8/Manitoba/pyNab/src'
sys.path.append(nab_path)
import nabPy as Nab
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sqlalchemy import Column, String, Integer, Date, FLOAT, Table
from sqlalchemy.orm import relationship
from models.base import Base
from models.sn import Sn


class Run(Base):
    __tablename__ = 'runs'

    id = Column(Integer, primary_key=True)
    run_number = Column(Integer)
    bias_voltage = Column(Integer)
    sn_source = Column(String)
    cd_source = Column(String)
    average_temp = Column(FLOAT)
    number_subruns = Column(Integer)
    pixels_calibrated = Column(Integer)
    date = Column(Date)
    pulser = Column(String)
    proton = Column(String)
    proton_energy = Column(Integer)
    pixels = relationship('Pixel', back_populates='run')

    def __init__(self, run_number, bias_voltage, sn_source, cd_source, average_temp, number_subruns, date, pulser, proton, proton_energy, directory, event_bool):
        self.run_number = run_number
        self.bias_voltage = bias_voltage
        self.sn_source = sn_source
        self.cd_source = cd_source
        self.average_temp = average_temp
        self.number_subruns = number_subruns
        self.date = date
        self.pulser = pulser
        self.proton = proton
        self.proton_energy = proton_energy
        self.nab_run = Nab.DataRun(directory, run_number, ignoreEventFile = event_bool)
        if sn_source:
            self.source = Sn()
        # if cd_source:
        #     self.source = Cd()

