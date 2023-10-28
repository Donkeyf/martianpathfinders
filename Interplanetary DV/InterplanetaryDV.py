from PorkchopPlot import PorkchopPlot
import CombinedMinDV as combdv

## Porkchop Plot for Earth to Starman (analysing launch dates over 2 years)
earth_to_starman = PorkchopPlot('Earth',   'Starman',  0,1,7,2039, 0,1,3,2041,      120,                 15,24*30,                    60,               "pro",         120)
# earth_to_starman.get_plot()

earth_to_starman.min_dv_toStarman()

def total_dv():        
    # Get total Delta-v of all segments (Earth to Starman to Mars to Earth)
    es_dv = earth_to_starman.min_dv # 1. Find the minimum dv for Earth to Starman
    se_dv = combdv.min_dvs_about_Mars() # 2-3. Find the minimum combined dv for Starman to Earth (via Mars)
    esme_dv = float(es_dv) + float(se_dv)

    print(f'TOTAL INTERPLANETARY DV RESULT: Optimised delta-v of {esme_dv} km/s for all interplanetary segments.')
    return esme_dv

total_dv()