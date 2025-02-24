import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pytest
import json
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from datetime import date
from models.base import Base
from models.run import Run
from models.pixel import Pixel

@pytest.fixture(scope="module")
def engine():
    engine = create_engine('sqlite:///:memory:')
    Base.metadata.create_all(engine)
    yield engine
    Base.metadata.drop_all(engine)

@pytest.fixture(scope="module")
def session(engine):
    Session = sessionmaker(bind=engine)
    session = Session()
    yield session
    session.rollback()

@pytest.fixture(scope="module")
def valid_run():
    valid_run = Run(1374, -300, True, False, 121.565, 15, date(2021,1,26), False, False)
    return valid_run

@pytest.fixture(scope="module")
def valid_pixel():
    energy_list = [1,2,5,7,8,8,6,9.0,5.7,3]
    energy_json = json.dumps(energy_list)
    valid_pixel = Pixel(76, 'good', 2.345, 1.199, 465.992, energy_json)
    return valid_pixel