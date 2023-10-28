from PorkchopPlot import PorkchopPlot

##### CHECKING OF OPTIMAL PERIODS TO TRAVEL BETWEEN STARMAN AND MARS #####
## The following dvs were computed between 2027 and 2039, analysing launch dates over 2 year periods
# 2027-10-01 to 2029-08-01: 14.713829889323348
# 2027-29 and 2028-30       14.02744525467152
# 2027-29 and 2029-31       14.030569104630981
# 2028-10-01 to 2030-08-01: 12.335580021467862
# 2028-30 and 2029-31       14.030766890636343
# 2028-30 and 2030-32       14.56184014629024
# 2029-10-01 to 2031-08-01: 14.50881380164731
# 2029-31 and 2031-33       14.31947980658829
# 2030-10-01 to 2032-08-01: 15.191399869653615
# 2030-32 and 2032-34       13.56250889110717
# 2031-10-01 to 2033-08-01: 13.558691335856007     
# 2031-33 and 2032-34       13.558992884033692
# 2031-33 and 2033-35       16.60733755789925
# 2032-10-01 to 2034-08-01: 14.803139543654778
# 2032-34 and 2033-35       12.280214276005227
# 2032-34 and 2034-36       12.278960758338563 ## low
# 2032-34 and 2035-37       21.15204233289733
# 2033-10-01 to 2035-08-01: 12.280344327945683
# 2033-35 and 2034-36       12.279090810279019
# 2033-35 and 2035-37       12.900557420917092
# 2034-10-01 to 2036-08-01: 15.663516183586946
# 2034-36 and 2036-38       12.862498420280257
# 2035-10-01 to 2037-08-01: 14.34086717969144 
# 2035-37 and 2037-39       13.316006085461083
# 2036-10-01 to 2038-08-01: 18.478926383075535
# 2036-38 and 2037-39       13.398959944205519
# 2036-38 and 2038-40       13.35335180261298
# 2037-10-01 to 2039-08-01: 12.99762594071197
# 2037-39 and 2038-40       12.903300442031878
# 2037-39 and 2039-41       15.367273060043996
# 2038-40 and 2040-42       13.836940870075415
# 2039-41 and 2041-43       15.425103202634784
# 2040-42 and 2042-44       11.369656208635345 ## lowest
# 1.10.2041 - 1.8.2043 and 1.10.2042 - 1.8.2044 11.36502777143581 ## lowest
# 2041-43 and 2043-45       14.103415346551383
# 2042-44 and 2044-46       12.857453442471437
# 2043-45 and 2045-47       14.13564726740125
# 2044-46 and 2046-48       14.7382081467543

##### Time period with lowest dv was 2040-42 for departing Starman, and 2042-44 for departing Mars #####
# 2042-10-01 to 2044-08-01

## Porkchop Plot for Earth to Starman
starman_to_mars = PorkchopPlot('Starman',   'Mars',  0,1,8,2040, 0,1,10,2042,      120,                 15,24*30,                    60,               "pro",         120)
# starman_to_mars.get_plot()

## Porkchop Plot for Earth to Starman
mars_to_earth = PorkchopPlot('Mars',   'Earth',  0,1,10,2042, 0,1,8,2044,      120,                 15,24*30,                    60,               "pro",         120)
# mars_to_earth.get_plot()

def min_dvs_about_Mars():
    # Find the smallest combination of dvs for Starman to Mars and Mars to Earth
    # Note: new_jd_2s are arrival dates and new_jd_1s are departure dates (for Mars in the following case)

    min_combined_dv = None # initialise
    min_i = None
    min_j = None

    for i in range(len(mars_to_earth.new_jd_1s)):
        for j in range(len(starman_to_mars.new_jd_2s)):
            date_diff = mars_to_earth.new_jd_1s[i] - starman_to_mars.new_jd_2s[j]

            if 90 <= date_diff <= 120: # only combine the dvs if Mars stay is between 3-4 months
                combined_dv = mars_to_earth.new_dvs[i] + starman_to_mars.new_dvs[j]

                if min_combined_dv is None or combined_dv < min_combined_dv: # if the new combined dv is less than the current minimum, replace
                    min_combined_dv = combined_dv
                    min_i = i # store index of min dv M to E
                    min_j = j # store index of min dv S to M
    
    if min_combined_dv is not None:
        print('STARMAN TO MARS and MARS TO EARTH:')

        print(f'Minimum delta-v of {min_combined_dv} km/s for Starman to Mars and Mars to Earth combined.')
        print(f'Corresponding to a minimum delta-v of {starman_to_mars.new_dvs[min_j]} km/s for departing Starman on {PorkchopPlot.julian_day_number_to_gregorian(int(starman_to_mars.new_jd_1s[min_j]))} and arriving at Mars on {PorkchopPlot.julian_day_number_to_gregorian(int(starman_to_mars.new_jd_2s[min_j]))}')
        print(f'and a minimum delta-v of {mars_to_earth.new_dvs[min_i]} km/s for departing Mars on {PorkchopPlot.julian_day_number_to_gregorian(int(mars_to_earth.new_jd_1s[min_i]))} and arriving at Earth on {PorkchopPlot.julian_day_number_to_gregorian(int(mars_to_earth.new_jd_2s[min_i]))}.')

    else:
        print("No suitable trajectories found from Starman to Mars and Mars to Earth for a 3-4 month stay time at Mars.")

    return min_combined_dv

# se_dv = min_dvs_about_Mars()
# print(se_dv)