import numpy as np

# Time (Julian date): Time at which to find elements
# Planet: Capitalised string for planet you want to find position of
# eg. Mars Earth
# Returns: 
def find_orbital_elements(JD, planet):
    T0 = (JD - 2451545)/36525 # Number of Julian centuries between J2000 and given time
    a = None # Semi-major axis (Km)
    da_dt = None # (Km/Century)
    e = None # Eccentricity
    de_dt = None # (1/Century)
    i = None # Inclination (Rads)
    di_dt = None # (Rads/Century)
    raan = None # Right ascension of the ascending node (Rads)
    draan_dt = None # (Rads/Century)
    longp = None # Longitude of perihelion (Rads)
    dlongp_dt = None # (Rads/Century)

    if (planet == "Earth"):
        a = 149598261.15044
        da_dt = 840.740033334
        e = 0.01671123
        de_dt = -0.00004392
        i = -0.0000002672
        di_dt = -0.0002259622
        raan = 0
        draan_dt = 0
        longp = 1.796601474
        dlongp_dt = 0.0056421894

    else if (planet == "Mars"):
        a = 227943822.42757
        da_dt = 27630.72671829
        e =  0.09339410
        de_dt = 0.00007882
        i = 0.0322832054
        di_dt = -0.0001419181
        raan = 0.8649771297
        draan_dt = -0.0051063697
        longp = -0.4178951712
        dlongp_dt = 0.0077564331

    parameters = np.array([a, e, i, raan, longp])
    d_parameters = np.array([da_dt, de_dt, di_dt, draan_dt, dlongp_dt])

    def propagate_parameter(x, dx_dt, dt):
        return x + dx_dt*dt

    parameters = propagate_parameter(parameters, d_parameters, T0)

    return parameters

    def main():
        parameters = find_orbital_elements(2451545, "Earth")
        print(parameters)

if __name__ == "__main__":
    main()

    
