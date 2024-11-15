from sqlalchemy import Column, String, Integer, DateTime, Numeric, Boolean

from base import Base

from sqlalchemy.orm import relationship

class Run(Base):
    __tablename__ = 'runs'
    id = Column(Integer, primary_key=True)
    run_number = Column(Integer)
    bias_voltage = Column(Integer)
    sn_source = Column(Boolean)
    cd_source = Column(Boolean)
    average_temp = Column(Numeric)
    number_subruns = Column(Integer)
    pixels_on = Column(Integer)
    pixels_calibrated = Column(Integer)
    date = Column(DateTime)
    pulser = Column(Boolean)
    proton = Column(Boolean)

    def __init__(self, run_number, bias_voltage, sn_source, cd_source, average_temp, number_subruns, date, pulser, proton):
        self.run_number = run_number
        self.bias_voltage = bias_voltage
        self.sn_source = sn_source
        self.cd_source = cd_source
        self.average_temp = average_temp
        self.number_subruns = number_subruns
        self.date = date
        self.pulser = pulser
        self.proton = proton