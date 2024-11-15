from sqlalchemy import create_engine, Column, Integer, String 
from sqlalchemy.orm import sessionmaker 
from sqlalchemy.ext.declarative import declarative_base 

# create a SQLite database engine 
engine = create_engine('sqlite:///example.db') 

# create a session factory 
Session = sessionmaker(bind=engine) 

# create a declarative base 
Base = declarative_base() 

# define a model class 


class User(Base): 
	__tablename__ = 'users'
	id = Column(Integer, primary_key=True) 
	name = Column(String) 
	age = Column(Integer) 



# create the database tables 
Base.metadata.create_all(engine) 

# insert some data 
session = Session() 
session.add_all([ 
	User(name='Alice', age=30), 
	User(name='Bob', age=35), 
	User(name='Charlie', age=40), 
]) 
session.commit() 

# retrieve the data 
users = session.query(User).all() 
for user in users: 
	print(user) 

# close the session 
session.close() 
