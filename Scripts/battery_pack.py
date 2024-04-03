# Does nothing because I don't know what to do with it
# Maybe work on this later

import numpy as np
import pandas as pd
class BatteryPack:
    def __init__(self, name, mass):
        self.name = name
        self.mass = mass

def import_battery_csv(csv_file):
    batteries_list = []
    battery_sheet = pd.read_csv(csv_file,index_col=0)

    for i in range(len(battery_sheet.index)):
        batteries_list.append(
            BatteryPack(battery_sheet.index[i]),
            BatteryPack(battery_sheet['Mass'][i])
        )

    return batteries_list

#%%
()
#%%
