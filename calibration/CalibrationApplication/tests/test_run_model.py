import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pytest
from models.run import Run

def test_run_valid(session, valid_run):
    session.add(valid_run)
    session.commit()
    run1374 = session.query(Run).filter_by(run_number=1374).first()
    assert run1374.bias_voltage == -300




