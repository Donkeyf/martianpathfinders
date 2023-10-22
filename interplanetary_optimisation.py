import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp


### Analyse launch dates over 2 years
## Porkchop Plot for Earth to Starman
earth_to_starman = pp.PorkchopPlot('Earth',   'Starman',  0,1,1,2031, 0,1,10,2033,      120,                 15,24*30,                    60,               "pro",         120)
earth_to_starman.get_plot()
earth_to_starman.get_min_dv()

## Porkchop Plot for Starman to Mars
starman_to_Mars = pp.PorkchopPlot('Starman',   'Mars',  0,1,1,2031, 0,1,10,2033,      120,                 15,24*30,                    60,               "pro",         120)
starman_to_Mars.get_plot()
starman_to_Mars.get_min_dv()

## Porkchop Plot for Mars to Earth
mars_to_earth = pp.PorkchopPlot('Mars',   'Earth',  0,1,1,2032, 0,1,10,2034,      120,                 15,24*30,                    60,               "pro",         120)
mars_to_earth.get_plot()
mars_to_earth.get_min_dv()



## Don't worry about this was just testing extra function
# print(earth_to_starman.r_A)
# print(np.shape(earth_to_starman.r_A))

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(earth_to_starman.r_A[0]/(1000.0),earth_to_starman.r_A[1]/(1000.0),earth_to_starman.r_A[2]/(1000.0))
# ax.scatter(earth_to_starman.r_B[0]/(1000.0),earth_to_starman.r_B[1]/(1000.0),earth_to_starman.r_B[2]/(1000.0))
# plt.savefig(f'orbits.png')
# plt.show()