from sqlalchemy import Column, Integer, Boolean, String, ARRAY

from base import Base

from sqlalchemy.orm import relationship

class Pixel(Base):
    __tablename__ = 'pixels'

    id = Column(Integer, primary_key=True)
    pixel_number = Column(Integer)
    calibration_flag = Column(String)
    ecap_fit_chi2 = Column(Numeric)
    xray_fit_chi2 = Column(Numeric)
    calibration_chi2 = Column(Numeric)
    energy = Column(ARRAY)
    

    def __init__(self, pixel_number, calibration_flag, calibration_chi2, ecap_fit_chi2, xray_fit_chi2, energy):
        self.pixel_number = pixel_number
        self.calibration_flag = calibration_flag
        self.ecap_fit_chi2 = ecap_fit_chi2
        self.xray_fit_chi2 = xray_fit_chi2
        self.calibration_chi2 = calibration_chi2
        self.energy = energy
