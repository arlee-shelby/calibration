import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sqlalchemy import Column, Integer, Boolean, String, JSON, FLOAT, ForeignKey
from sqlalchemy.orm import relationship
from models.base import Base

class Pixel(Base):
    __tablename__ = 'pixels'

    id = Column(Integer, primary_key=True)
    run_number = Column(Integer, ForeignKey('runs.run_number'))
    run = relationship('Run', back_populates='pixels')
    pixel_number = Column(Integer)
    calibration_flag = Column(String)
    ecap_fit_chi2 = Column(FLOAT)
    xray_fit_chi2 = Column(FLOAT)
    calibration_chi2 = Column(FLOAT)
    energy = Column(JSON)
    

    def __init__(self, run, pixel_number, trap_rise, trap_length, trap_decay):
        self.pixel_number = pixel_number
        self.run = run
        run.nab_run.singleWaves().resetCuts()
        run.nab_run.singleWaves().defineCut('pixel', '=', pixel_number)
        trap_filter = run.nab_run.singleWaves().determineEnergyTiming(method='trap', params=[trap_rise, trap_length, trap_decay])
        self.trap_filter = trap_filter
        self.energy = trap_filter.data()['energy'].tolist()