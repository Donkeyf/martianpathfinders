from PorkchopPlot import PorkchopPlot

##### CHECKING OF OPTIMAL PERIODS TO TRAVEL BETWEEN STARMAN AND MARS #####
## The following dvs were computed between 2027 and 2039, analysing launch dates over 2 year periods
# 2027-10-01 to 2029-08-01: 10.989313318077789
# 2028-10-01 to 2030-08-01: 12.476533319329398
# 2029-10-01 to 2031-08-01: 11.975404776628531
# 2030-10-01 to 2032-08-01: 11.975404776628531
# 2031-10-01 to 2033-08-01: 13.715518635707486
# 2032-10-01 to 2034-08-01: 12.012644071838682
# 2033-10-01 to 2035-08-01: 10.7315481598824
# 2034-10-01 to 2036-08-01: 10.91158653967717
# 2035-10-01 to 2037-08-01: 11.732436322855172
# 2036-10-01 to 2038-08-01: 13.869721285045102
# 2037-10-01 to 2039-08-01: 13.009896089429802

##### Time period with lowest dv was 2033-2035 so this was selected and run #####

## Porkchop Plot for Earth to Starman
starman_to_mars = PorkchopPlot('Starman',   'Mars',  0,1,10,2033, 0,1,8,2035,      120,                 15,24*30,                    60,               "pro",         120)
starman_to_mars.get_plot()

## Porkchop Plot for Earth to Starman
mars_to_earth = PorkchopPlot('Starman',   'Mars',  0,1,10,2033, 0,1,8,2035,      120,                 15,24*30,                    60,               "pro",         120)
mars_to_earth.get_plot()

def min_dvs_about_Mars():
    # Find the smallest combination of dvs for Starman to Mars and Mars to Earth
    # Note: new_jd_2s are arrival dates and new_jd_1s are departure dates (for Mars in the following case)

    min_combined_dv = None # initialise
    min_i = None
    min_j = None

    for i in range(len(mars_to_earth.new_jd_1s)):
        for j in range(len(starman_to_mars.new_jd_2s)):
            date_diff = abs(mars_to_earth.new_jd_1s[i] - starman_to_mars.new_jd_2s[j])

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
