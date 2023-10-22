import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
from scipy.optimize import minimize
import FindStateVector as fsv
import PorkchopPlot as pp

### Analyse launch dates over 2 years
## Porkchop Plot for Earth to Starman
earth_to_starman = pp.PorkchopPlot('Earth',   'Starman',  0,1,1,2031, 0,1,10,2033,      120,                 15,24*30,                    60,               "pro",         120)
earth_to_starman.get_plot()
earth_to_starman.get_min_dv()

## Porkchop Plot for Starman to Mars
starman_to_Mars = pp.PorkchopPlot('Starman',   'Mars',  0,1,10,2033, 0,1,8,2036,      120,                 15,24*30,                    60,               "pro",         120)
starman_to_Mars.get_plot()
starman_to_Mars.get_min_dv()

## Porkchop Plot for Mars to Earth
mars_to_earth = pp.PorkchopPlot('Mars',   'Earth',  0,10,10,2034, 0,10,11,2034,      120,                 15,24*30,                    60,               "pro",         120)
mars_to_earth.get_plot()
mars_to_earth.get_min_dv()

### ATTEMPT 1 (depart 31-33) ###
## 'Earth',   'Starman',  0,1,1,2031, 0,1,10,2033
# Minimum delta-v of 3.900596044888159 km/s departing at 2032-03-23 and arriving at 2033-06-23
## 'Starman',   'Mars',  0,1,10,2033, 0,1,8,2036
# Minimum delta-v of 5.148557244135301 km/s departing at 2033-11-22 and arriving at 2034-07-10
## 'Mars',   'Earth',  0,10,10,2034, 0,10,11,2034
# Minimum delta-v of 8.315854140818633 km/s departing at 2034-11-08 and arriving at 2035-10-24

### ATTEMPT 2 (depart 29-31) ###
## 'Earth',   'Starman',  0,1,1,2029, 0,1,10,2031
# Minimum delta-v of 3.894707511253301 km/s departing at 2029-03-18 and arriving at 2030-06-18
## 'Starman',   'Mars',  0,1,8,2030, 0,1,6,2032
# Minimum delta-v of 5.986142496164587 km/s departing at 2031-07-10 and arriving at 2031-11-22
## 'Mars',   'Earth',  0,22,2,2032, 0,22,3,2032
# Minimum delta-v of 9.4369729371202 km/s departing at 2032-03-22 and arriving at 2033-07-16

### ATTEMPT 2 (depart 33-35) ###
## 'Earth',   'Starman',  0,1,1,2033, 0,1,10,2035
# Minimum delta-v of 3.8893148337551584 km/s departing at 2035-03-21 and arriving at 2036-07-02
## 'Starman',   'Mars',  0,1,6,2035, 0,1,4,2037
# Minimum delta-v of 5.1053535317502785 km/s departing at 2037-03-26 and arriving at 2038-10-24
## 'Mars',   'Earth',  0,24,1,2039, 0,24,2,2039
# Minimum delta-v of 9.919899599431957 km/s departing at 2039-02-24 and arriving at 2040-02-21

## Don't worry about this was just testing extra function
# print(earth_to_starman.r_A)
# print(np.shape(earth_to_starman.r_A))

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(earth_to_starman.r_A[0]/(1000.0),earth_to_starman.r_A[1]/(1000.0),earth_to_starman.r_A[2]/(1000.0))
# ax.scatter(earth_to_starman.r_B[0]/(1000.0),earth_to_starman.r_B[1]/(1000.0),earth_to_starman.r_B[2]/(1000.0))
# plt.savefig(f'orbits.png')
# plt.show()