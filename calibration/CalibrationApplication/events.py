import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from sqlalchemy import event
from models.pixel import Pixel
from models.run import Run

# @event.listens_for(Child, 'after_insert')
# def update_parent_column(mapper, connection, target):
#     parent = target.parent  # Assuming a relationship between Child and Parent
#     parent.column_to_update = "New value"  # Update the desired column in the parent
#     connection.execute(
#         Parent.__table__.update().
#         where(Parent.id == parent.id).
#         values(column_to_update=parent.column_to_update)
#     )