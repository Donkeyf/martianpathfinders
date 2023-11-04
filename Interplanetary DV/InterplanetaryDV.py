from PorkchopPlot import PorkchopPlot
import CombinedMinDV as combdv

## Porkchop Plot for Earth to Starman (analysing launch dates over 2 years)
earth_to_starman = PorkchopPlot('Earth',   'Starman',  0,1,7,2039, 0,1,2,2041,      120,                 15,24*30,                    60,               "pro",         120)
# earth_to_starman.get_plot()
# earth_to_starman.get_min_dv()

# earth_to_starman.min_dv_toStarman()

def total_dv():

    # Get total Delta-v of all segments (Earth to Starman to Mars to Earth)
    se_dv, jd_starman_dep= combdv.min_dvs_about_Mars() # 2-3. Find the minimum combined dv for Starman to Earth (via Mars)
    # es_dv = earth_to_starman.min_dv
    print(f'date {PorkchopPlot.julian_day_number_to_gregorian(int(jd_starman_dep))}')
    es_dv = None # initialise
    min_i = None

    for i in range(len(earth_to_starman.new_jd_1s)):
        if earth_to_starman.new_jd_2s[i] < jd_starman_dep: # only take if arrival date is earlier than Starman departure date
            if es_dv is None or earth_to_starman.new_dvs[i] < es_dv: # if the new combined dv is less than the current minimum, replace
                es_dv = earth_to_starman.new_dvs[i]
                min_i = i # store index of min dv E to S

    esme_dv = float(es_dv) + float(se_dv)

    print('EARTH TO STARMAN:')
    print(f'Depart Earth at {PorkchopPlot.julian_day_number_to_gregorian(int(earth_to_starman.new_jd_1s[min_i]))} arriving Starman at {PorkchopPlot.julian_day_number_to_gregorian(int(earth_to_starman.new_jd_2s[min_i]))} with {es_dv} km/s')

    print(f'TOTAL INTERPLANETARY DV RESULT: Optimised delta-v of {esme_dv} km/s for all interplanetary segments.')
    return esme_dv

total_dv()
