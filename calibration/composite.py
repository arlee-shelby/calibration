import json
from component import Component

class Composite(Component):
    def __init__(self, key):
        self.key = key
        self.children = []

    def add(self, component):
        self.children.append(component)

    def remove(self, component):
        self.children.remove(component)

    def get_children(self):
        return self.children

    def to_json(self):
        result = {}
        for child in self.children:
            result.update(child.to_json())
        return {self.key: result}

    def save_to_db(self, cursor):
        cursor.execute(f''' CREATE TABLE IF NOT EXISTS your_table6 (id INTEGER PRIMARY KEY AUTOINCREMENT, {self.key} INTEGER)''')
        
        cursor.execute(f"INSERT INTO your_table6 ({self.key}) VALUES (?)", (self.value,))