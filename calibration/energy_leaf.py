from leaf import Leaf

class Leaf(Component):
    def __init__(self, key, value):
        self.key = key
        self.value = value

    def to_json(self):
        return {self.key: self.value}

    def save_to_db(self, cursor):
        sql = "INSERT INTO your_table (name, value) VALUES (%s, %s)"
        val = (self.key, json.dumps(self.value))
        cursor.execute(sql, val)