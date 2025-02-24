import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from models.base import Session
from models.run import Run

session = Session()

class RunController:
    def get_runs():
        return session.query(Run).all()

    def get_run_by_number(num):
        return session.query(Run).filter_by(run_number = num).all()[0]
