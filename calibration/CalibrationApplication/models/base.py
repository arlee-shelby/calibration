from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os

if os.path.exists('test1.db'):
   os.remove('test1.db')

engine = create_engine('sqlite:///test1.db')
Session = sessionmaker(bind=engine)

Base = declarative_base()
