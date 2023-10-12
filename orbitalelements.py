import numpy as np

MU_SUN = 1.327 * 10**11 # (km^3/s^2)

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
    L = None # Mean Longitude (Rads)
    dL_dt = None # (Rads/Century)

    # Check angles are in 0-2pi

    def wrap_angle(angle):
        while (angle > 2*np.pi):
            angle -= 2*np.pi
        while (angle < 0):
            angle += 2*np.pi

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
        L =  1.7534375571
        dL_dt = 628.3075779009

    elif (planet == "Mars"):
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
        L = -0.0794723815
        dL_dt = 334.0613016814

    parameters = np.array([a, e, i, raan, longp, L])
    d_parameters = np.array([da_dt, de_dt, di_dt, draan_dt, dlongp_dt, dL_dt])

    def propagate_parameter(x, dx_dt, dt):
        return x + dx_dt*dt

    parameters = propagate_parameter(parameters, d_parameters, T0)
    a = parameters[0]
    e = parameters[1]
    i = parameters[2]
    raan = parameters[3]
    longp = parameters[4]
    L = parameters[5]

    i = wrap_angle(i)
    raan = wrap_angle(raan)
    longp = wrap_angle(longp)
    L = wrap_angle(L)

    h = np.sqrt(MU_SUN * a * (1-e**2)) # Angular Momentum
    argp = longp - raan # Argument of perihelion
    M = L - longp # Mean anomaly
    E = None # Eccentric anomaly

    def kepler_equation(E):
        M = E - e*np.sin(E)
        return M

    init_guess = M
    difference = 100
    tolerance = 10e-12
    E = init_guess
    E1 = None

    # Newton method used to find eccentric anomaly
    while (difference > tolerance):
        F = kepler_equation(E) - M
        F_prime = 1 - e*np.cos(E)
        E1 = E - (F/F_prime)
        difference = np.abs(E1 - E)
        E = E1

    ta = 2 * np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2)) # True anomaly
    if (ta < 0):
        ta = 2*np.pi - ta
    



    
