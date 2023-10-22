import numpy as np

MU_SUN = 1.327 * 10**11 # (km^3/s^2)

# Time (Julian date): Time at which to find elements
# Planet: Capitalised string for planet you want to find the state vector for (or Starman)
# eg. "Mars" "Earth" "Starman"
# Returns: Position vector, Velocity vector (heliocentric ecliptic frame (HEE))
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
        return angle

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
        
    elif (planet == "Starman"):
        # Values are taken from Jan-01-2023
        T0 *= 36525
        T0 += 2451545
        T0 = (JD - date_to_JD(0, 1, 1, 2023))/36525
        a = 1.982407920279165E+08
        da_dt = 80.0555081427*36525
        e = 2.559421363580004E-01
        de_dt = -3.00785672e-7*36525
        i = 1.075055082468543E+00 * np.pi/180
        di_dt = -3.53068568e-7 * np.pi/180*36525
        raan = 3.169096418834785E+02 * np.pi/180
        draan_dt = 0.00004244059 * np.pi/180*36525
        longp = 1.777580453292394E+02 * np.pi/180
        dlongp_dt = 0.000050803 * np.pi/180*36525
        L = 7.789083980704730E+01 * np.pi/180
        dL_dt = -0.64618117421 * np.pi/180*36525

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
    print("Right Ascension: " + str(raan*180/np.pi))
    print("Declination: " + str(i*180/np.pi))

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
        ta += 2*np.pi
        
    perifocal_vectors = elements_to_perifocal(ta, a, e, MU_SUN, h)
    hee_vectors = perifocal_to_hee(perifocal_vectors, i, raan, argp)
    position = hee_vectors[:3]
    velocity = hee_vectors[3:]
    return np.array([position, velocity])
    
def transformation_matrix(i, raan, argp):
    transformer = np.zeros((3, 3))

    transformer[0][0] = (-np.sin(raan) * np.cos(i) * np.sin(argp)
    + np.cos(raan) * np.cos(argp))
    transformer[0][1] = (-np.sin(raan) * np.cos(i) * np.cos(argp)
    - np.cos(raan) * np.sin(argp))
    transformer[0][2] = np.sin(raan) * np.sin(i)

    transformer[1][0] = (np.cos(raan) * np.cos(i) * np.sin(argp)
    + np.sin(raan) * np.cos(argp))
    transformer[1][1] = (np.cos(raan) * np.cos(i) * np.cos(argp)
    - np.sin(raan) * np.sin(argp))
    transformer[1][2] = -np.cos(raan) * np.sin(i)

    transformer[2][0] = np.sin(i) * np.sin(argp)
    transformer[2][1] = np.sin(i) * np.cos(argp)
    transformer[2][2] = np.cos(i)
    
    return transformer

def elements_to_perifocal(ta, a, e, mu, h):
    r = (h**2/mu)*(1/(1+e*np.cos(ta)))

    # Perifocal frame automatically results in w = 0 and dw = 0
    p = r * np.cos(ta)
    q = r * np.sin(ta)
    w = 0

    dp = -mu/h * np.sin(ta)
    dq = mu/h * (e + np.cos(ta))
    dw = 0

    return np.array([p, q, w, dp, dq, dw])


def perifocal_to_hee(perifocal_vectors, i, raan, argp):
    x = np.zeros_like(perifocal_vectors[0])
    y = np.zeros_like(perifocal_vectors[0])
    z = np.zeros_like(perifocal_vectors[0])

    dx = np.zeros_like(perifocal_vectors[0])
    dy = np.zeros_like(perifocal_vectors[0])
    dz = np.zeros_like(perifocal_vectors[0])

    transformation = transformation_matrix(i, raan, argp)

    x, y, z = np.matmul(transformation, [perifocal_vectors[0], perifocal_vectors[1], perifocal_vectors[2]])
    dx, dy, dz = np.matmul(transformation, [perifocal_vectors[3], perifocal_vectors[4], perifocal_vectors[5]])

    return np.array([x, y, z, dx, dy, dz])


def date_to_JD(UT, day, month, year):
    """
    Converts a given date and time into a Julian Date
    
    Args:
        UT (float): Universal time, in hours. (0 <= UT < 24)
        day (int): Day of the month. (1 <= day <= 31)
        month (int): Month of the year. (1 <= month <= 12)
        year (int): Year. (Unbounded)
        
    Returns:
        float: Julian Date equivalent to the given date and time.
    """
    JD = (float) (367*year - int((7)*(year+int((month+9)/12))/4) + int(275*month/9) + day + 1721013.5 + UT/24)
    return JD

if __name__ == "__main__":
    print(find_orbital_elements(date_to_JD(0, 23, 10, 2023), "Starman"))