import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp
import findOrbitalElements as foe
import orbitalTransforms2 as ot
import compute_params as cp

offset = 2430000.0

print(fsv.date_to_JD(0,23,6,2033) - offset)

print(fsv.find_orbital_elements(fsv.date_to_JD(0,23,3,2032),'Starman'))