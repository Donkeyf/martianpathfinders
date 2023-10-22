import numpy as np

def fuel_consumption(engine, g, m, delta_v):

    #m = current mass of spacecraft
    #g = specific gravity? 
    #delta_v of the specific maneuver that is being undertaken

    if engine == "roadster":
        Isp = 900

        fuel_used = m*(1 - np.e**(-delta_v/(Isp*g)) )
        new_mass = m - fuel_used

    else:
        Isp = 455 #Isp for liquid hydrogen
        fuel_used = fuel_used = m*(1 - np.e**(-delta_v/(Isp*g)))
        new_mass = m - fuel_used


    return fuel_used, new_mass
