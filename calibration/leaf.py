import json
import sqlite3
from component import Component

class Leaf(Component):
    def __init__(self, key, value):
        self.key = key
        self.value = value

    def to_json(self):
        return {self.key: self.value}

    def save_to_db(self, cursor):
        cursor.execute(f''' CREATE TABLE IF NOT EXISTS your_table6 (id INTEGER PRIMARY KEY AUTOINCREMENT, {self.key} INTEGER)''')
        cursor.execute(f"INSERT INTO your_table6 ({self.key}) VALUES (?)", (self.value,))