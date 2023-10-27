import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv
import PorkchopPlot as pp

# Use new_jd1 and new_jd2 (outliers removed), where jd1 is departure, and jd2 is arrival
# Use classes like the following: earth_to_mars.new_jd1s

## STARMAN TO MARS (jd1=SD, jd2=MA)

# Define time period (2033 to 2035)
starman_to_mars = pp.PorkchopPlot('Starman',   'Mars',  0,1,10,2033, 0,1,8,2035,      120,                 15,24*30,                    60,               "pro",         120)

# Sorting function that will sort the dv's into an array of smallest to largest for Starman to Mars within a departure range, and sort the jd1 and jd2 arrays (new) in the same order
for i in range(len(starman_to_mars.new_dvs)):
    for j in range(i+1, len(starman_to_mars.new_dvs)):
        if starman_to_mars.new_dvs[i] > starman_to_mars.new_dvs[j]: # swap until in ascending order of values

            sm_dv_temp = starman_to_mars.new_dvs[j]
            starman_to_mars.new_dvs[j] = starman_to_mars.new_dvs[i]
            starman_to_mars.new_dvs[i] = sm_dv_temp

            sm_jd1_temp = starman_to_mars.new_jd_1s[j]
            starman_to_mars.new_jd_1s[j] = starman_to_mars.new_jd_1s[i]
            starman_to_mars.new_jd_1s[i] = sm_jd1_temp

            sm_jd2_temp = starman_to_mars.new_jd_2s[j]            
            starman_to_mars.new_jd_2s[j] = starman_to_mars.new_jd_2s[i]
            starman_to_mars.new_jd_2s[i] = sm_jd2_temp

# print(starman_to_mars.new_dvs)
# print(starman_to_mars.new_jd_1s)
# print(starman_to_mars.new_jd_2s)

## MARS TO EARTH (jd1=MD, jd2=EA)

# Define time period (2033 to 2035)
mars_to_earth = pp.PorkchopPlot('Starman',   'Mars',  0,1,10,2033, 0,1,8,2035,      120,                 15,24*30,                    60,               "pro",         120)

# Sorting function that will sort the dv's into an array of smallest to largest for Mars to Earth within a departure range, and align jd1 and jd2 arrays
for i in range(len(mars_to_earth.new_dvs)):
    for j in range(i+1, len(mars_to_earth.new_dvs)):
        if mars_to_earth.new_dvs[i] > mars_to_earth.new_dvs[j]: # swap until in ascending order of values

            me_dv_temp = mars_to_earth.new_dvs[j]
            mars_to_earth.new_dvs[j] = mars_to_earth.new_dvs[i]
            mars_to_earth.new_dvs[i] = me_dv_temp

            me_jd1_temp = mars_to_earth.new_jd_1s[j]
            mars_to_earth.new_jd_1s[j] = mars_to_earth.new_jd_1s[i]
            mars_to_earth.new_jd_1s[i] = me_jd1_temp

            me_jd2_temp = mars_to_earth.new_jd_2s[j]            
            mars_to_earth.new_jd_2s[j] = mars_to_earth.new_jd_2s[i]
            mars_to_earth.new_jd_2s[i] = me_jd2_temp

# print(mars_to_earth.new_dvs)
# print(mars_to_earth.new_jd_1s)
# print(mars_to_earth.new_jd_2s)

## Now iterate through the dv's and find the first dv where jd2(MA) and jd1(MD) are 3-4 months apart, this will be the min possible dv

# found = False
# i = 0
# j = 0

for i in range(len(mars_to_earth.new_jd_1s)):
    for j in range(len(starman_to_mars.new_jd_2s)):
        difference = abs(mars_to_earth.new_jd_1s[i] - starman_to_mars.new_jd_2s[j])
        if 90 <= difference <= 120:
            found = True
            break
    if found:
        break

# for i in range(len(starman_to_mars.new_dvs)):
#     se_dvs = starman_to_mars.new_dvs[i] + starman_to_mars.new_dvs[i]
#     min_se_dv = min(starman_to_mars.new_dvs[i] + starman_to_mars.new_dvs[i])

if found:
    print("Suitable trajectories found:")
    print("mars_to_earth.new_jd_1s[{}]: {}".format(i, mars_to_earth.new_jd_1s[i]))
    print("starman_to_mars.new_jd_2s[{}]: {}".format(j, starman_to_mars.new_jd_2s[j]))
    print("mars_to_earth.new_dvs[{}]: {}".format(i, mars_to_earth.new_dvs[i]))
    print("starman_to_mars.new_dvs[{}]: {}".format(j, starman_to_mars.new_dvs[j]))
else:
    print("No suitable trajectories found for a 3-4 month stay time.")
