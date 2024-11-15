from collections import defaultdict
import numpy as np
import json

# class Room():
#     def __init__(self, num):
#         self.room_number = num

#     def key(self):
#         return str(self.room_number)

# class Guest():
#     def __init__(self, name):
#         self.name = name

#     def key(self):
#         return self.name

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))
# room_booking = nested_dict(2, str)

# class Room_Booking():
#     def __init__(self, date):
#         self.date = date


# room1 = Room(1)
# guest1 = Guest("Joe")

# room_booking[room1.key()][guest1.key()] = Room_Booking("some date")

# print(room_booking[room1.key()][guest1.key()])
x = np.arange(0,5).tolist()
y = np.arange(9,15).tolist()
d = defaultdict(list)
for i in x:
    d[0].append(i)
print(d)

p = (76,98,4)
d2 = defaultdict(dict)
for i in p:
    d2[i] = d

data = {
    "run_id":1374,
    "pixel_id":76,
    "energy": x,
    "parameters": y
}

data3 = {
    "run_id":1374,
    "pixel_id":87
}
data2 = {
    "run_id":1375,
    "pixel_id":78,
    "energy": x,
    "parameters": y
}
entreis = {}
print(data)
print(json.dumps(data))
with open('data.json', 'w') as f:
    json.dump(data, f,sort_keys = True, indent = 4,ensure_ascii = False)
    json.dump(data2, f,sort_keys = True, indent = 4,ensure_ascii = False)
    json.dump(data3, f,sort_keys = True, indent = 4,ensure_ascii = False)

# with open('data.json') as data_file:
#     data_loaded = json.load(data_file)

# print(data_loaded)


import sqlite3

conn = sqlite3.connect('mydatabase.db')
cursor = conn.cursor()

# Create table if not exists
cursor.execute('''
    CREATE TABLE IF NOT EXISTS mytable (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        json_col, json_col
    )
''')
print(data3.keys())
dd = json.dumps(data3)
print(dd)
# Insert JSON data
[cursor.execute("INSERT INTO mytable VALUES (?, ?)", (json.dumps(key),json.dumps(data3[key]))) for key in data3.keys()]
conn.commit()
conn.close()

print(data3.keys())
