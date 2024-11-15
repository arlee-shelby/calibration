from leaf import Leaf
from composite import Composite
import numpy as np

x = np.arange(4,8)
composite = Composite('pixels')
composite.add(Leaf('pixel_id', 76))
composite.add(Leaf('run_id', 1374))
leaf1 = Leaf('pixel_id', 76)
leaf2 = Leaf('run_id', 1374)

leaf3 = Leaf('pixel_id',97)

json_data = composite.to_json()
print(json_data)

import sqlite3

conn = sqlite3.connect('mydatabase2.db')
cursor = conn.cursor()

cursor.execute('''
    CREATE TABLE IF NOT EXISTS mytable2 (
        id INTEGER PRIMARY KEY AUTOINCREMENT, 
        pixel_id INTEGER, 
        run_id INTEGER
    )
''')

cursor.execute("INSERT INTO mytable2 (pixel_id,run_id) VALUES (?,?)",(leaf1.value,leaf2.value))
conn.commit()
conn.close()

print(type(leaf1.value))
conn = sqlite3.connect('mydatabase3.db')
cursor = conn.cursor()
leaf3.save_to_db(cursor)
leaf1.save_to_db(cursor)
conn.commit()
conn.close()