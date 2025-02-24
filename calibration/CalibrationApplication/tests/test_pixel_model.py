import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pytest
from models.pixel import Pixel

def test_pixel_valid(session, valid_pixel):
    session.add(valid_pixel)
    session.commit()
    pixel76 = session.query(Pixel).filter_by(pixel_number=76).first()
    assert pixel76.calibration_chi2 == pytest.approx(2.345)