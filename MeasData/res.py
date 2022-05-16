import pandas as pd
import numpy as np


main_dict = {"Freq":  290,
"Eccentricity": 0.5021,
"Maj. Axis": 10.4221,
"Min. Axis": 9.0131,
"Deg": 20,
"Freq":  290,
"Eccentricity": 0.3792,
"Maj. Axis": 9.8189,
"Min. Axis": 9.0856,
"Deg": 30,
"Freq":  290,
"Eccentricity": 0.67279,
"Maj. Axis": 13.0607,
"Min. Axis": 9.6627,

"Freq":  290,
"Eccentricity": 0.8099,
"Maj. Axis": 17.5577,
"Min. Axis": 10.2988,
"Deg": 40,
"Freq":  290,
"Eccentricity": 0.83247,
"Maj. Axis": 18.3197,
"Min. Axis": 10.1505,
"Deg": 45,}

eccen = [0.5021, 0.3792, 0.6727, 0.8099, 0.83247]

print(np.mean(eccen))
