import numpy as np

def compute_anomalies(t_vec, ma, e, n):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = t_vec * n + ma

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return E - e * np.sin(E)

    # Calculate eccentric anomaly --- Newton's Method
    ea = np.zeros(len(ma_vec))      # initialise eccentric anomaly vector

    for j in range(len(ma_vec)):
        E0 = ma_vec[j]              # set initial guess for E
        tol = 1                     # initialise tolerance value
   
        while tol > 1e-8:           # iterate until tolerance condition (<1e-8) is satisfied
            f = kepler_equation(E0) - ma_vec[j]
            f_prime = 1 - e * np.cos(E0)
            E1 = E0 - f/f_prime     # Newton's method

            tol = abs(E1-E0)        # calculate tolerance
            E0 = E1                 # update guess value for following iteration

        ea[j] = E1                  # update eccentric anomaly vector

    # Calculate true anomaly
    ta = 2 * np.arctan( np.tan( ea/2 ) / np.sqrt( (1-e)/(1+e) ) )

    # Wrap anomalies to -pi:pi
    ta = [i-2*np.pi if i>np.pi else i for i in ta]
    ea = [i-2*np.pi if i>np.pi else i for i in ea]
    ma_vec = [i-2*np.pi if i>np.pi else i for i in ma_vec]

    return ta, ea, ma_vec


def compute_orbital_velocity(h, e, ta, mu):

    v_r = mu/h *  e * np.sin(ta)            # Radial velocity
    v_n = mu/h * (1 + e*np.cos(ta))         # Normal/Tangential velocity

    return v_r, v_n


def elements_to_perifocal(ta, e, mu, h):
    # Calculate perifocal distance in m
    r = h**2/mu * 1/(1+e*np.cos(ta))

    # Compute perifocal coordinates
    p = r * np.cos(ta)
    q = r * np.sin(ta)
    w = np.zeros_like(p)

    # Compute perifocal velocities
    dp = -(mu/h) * np.sin(ta)
    dq = (mu/h) * (e+np.cos(ta))
    dw = np.zeros_like(p)

    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):

    # Transform coordinates to ECI frame
    x = np.zeros_like(p)
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    # Compute coordinates
    for j in range(len(p)):
        perifocal_position = np.matrix([[p[j]],[q[j]],[w[j]]])
        perifocal_velocity = np.matrix([[dp[j]],[dq[j]],[dw[j]]])
        eci_position = np.matmul(perifocal_to_eci_matrix(i, raan, argp), perifocal_position)
        eci_velocity = np.matmul(perifocal_to_eci_matrix(i, raan, argp), perifocal_velocity)
        x[j] = eci_position[0]
        y[j] = eci_position[1]
        z[j] = eci_position[2]
        
        dx[j] = eci_velocity[0]
        dy[j] = eci_velocity[1]
        dz[j] = eci_velocity[2]

    return x, y, z, dx, dy, dz


def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    
    p_to_e = np.matrix([[-np.sin(raan)*np.cos(i)*np.sin(argp)+np.cos(raan)*np.cos(argp), -np.sin(raan)*np.cos(i)*np.cos(argp)-np.cos(raan)*np.sin(argp), np.sin(raan)*np.sin(i)],
              [np.cos(raan)*np.cos(i)*np.sin(argp)+np.sin(raan)*np.cos(argp), np.cos(raan)*np.cos(i)*np.cos(argp)-np.sin(raan)*np.sin(argp), -np.cos(raan)*np.sin(i)], 
              [np.sin(i)*np.sin(argp), np.sin(i)*np.cos(argp), np.cos(i)]])

    return p_to_e


def julian_time_0(y,m,d):

    # Calculate Julian day number at 0h UT
    j_0 = 367*y - int(7 * (y+int((m+9)/12)) / 4) + int(275*m/9) + d + 1721013.5

    # Calculate Julian century
    T_0 = (j_0-2451545) / 36525

    return j_0, T_0


def theta_G(T_0, hour):

    # Calculate Greenwich sidereal time angle 
    
    # Greenwich time (angle) at 0hr UT
    theta_G0 = 100.4606184 + 36000.770044*T_0 + 0.000387933*T_0**2 - 2.583e-8*T_0**3
    if theta_G0 > 360:                  # wrap to 0 < theta < 360
        theta_G0 = theta_G0 - int(theta_G0/360) * 360

    # Greenwich sidereal time (angle) 
    theta_G = theta_G0 + 360.98564724 * hour
    if theta_G > 360:                   # wrap to 0 < theta < 360
        theta_G = theta_G - int(theta_G/360) * 360

    # Convert to radians
    theta_G = theta_G * np.pi/180
    
    return theta_G


def eci_to_ecef_matrix(theta):

    # Calculate transformation matrix from perifocal to ECI frame

    eci_to_ecef = [[np.cos(theta), np.sin(theta), 0],
                   [-np.sin(theta), np.cos(theta), 0],
                   [0, 0, 1]]
    
    return eci_to_ecef


def eci_to_ecef(x, y, z, theta):
    
    # Transform coordinates to ECEF frame
    x_dash = np.zeros_like(x)
    y_dash = np.zeros_like(x)
    z_dash = np.zeros_like(x)

    # Compute coordinates
    for j in range(len(x)):
        eci_position = np.matrix([[x[j]],[y[j]],[z[j]]])
        ecef_position = np.matmul(eci_to_ecef_matrix(theta[j]), eci_position)
        x_dash[j] = ecef_position[0]
        y_dash[j] = ecef_position[1]
        z_dash[j] = ecef_position[2]

    return x_dash, y_dash, z_dash


def ephemeris(x,y,z):

    # Calculate declination and right ascension angles
    dec = np.zeros_like(x)
    ra = np.zeros_like(x)

    for j in range(len(x)):
        r = np.sqrt(x[j]**2 + y[j]**2 + z[j]**2)
        l = x[j]/r
        m = y[j]/r
        n = z[j]/r

        dec[j] = np.arcsin(n)          # declination angle

        # right ascension angle, accounting for quadrant ambiguity
        if m > 0:
            ra[j] = np.arccos(l / np.cos(dec[j]))
        else:
            ra[j] = 2*np.pi - np.arccos(l / np.cos(dec[j]))
        
        # Convert to degrees
        dec[j] = dec[j] * 180/np.pi
        ra[j] = ra[j] * 180/np.pi 

    return dec, ra


def find_orbital_elements(R, V, mu):
    # find the orbital elements when given the state vectors
    # use R1, V1 where applicable for the anomalies at the starting point

    # Distance and Speed
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    ## radial velocity
    v_r = np.dot(R,V) / r

    # Specific Angular Momentum
    H = np.cross(R,V)
    h = np.linalg.norm(H)

    # Inclination
    i = np.arccos(H[2]/h)

    # Node Line
    N = np.cross([0,0,1], H)
    n = np.linalg.norm(N)
    
    # Right Ascension of the Ascending Node
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    # Eccentricity
    E = 1/mu * ((v**2-mu/r)*R - r*v_r*V)
    e = np.linalg.norm(E)

    # Argument of Perigee
    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )
    
    # Semi-major Axis
    a = h**2/mu * 1/(1-e**2)

    # Period
    T = 2*np.pi/np.sqrt(mu) * a**(3/2)

    # Mean Motion
    mm = 2*np.pi/T

    # Anomalies
    ## True anomaly
    if v_r >= 0:
        ta = np.arccos( np.dot(E/e, R/r) )
    elif v_r < 0:
        ta = 2*np.pi - np.arccos( np.dot(E/e, R/r) )

    ## Eccentric Anomaly
    ea = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(ta/2))

    ## Mean Anomaly
    ma = ea - e*np.sin(ea)

    return h, e, a, T, mm, i, raan, argp, ta, ma


def helio_to_eci(Vec_Helio, planet="Earth"):
    # convert vectors between Heliocentric frame to ECI frame
    # only applies for velocities, not positions 

    # obliquity (note negative angle for rotation)
    if planet == "Earth":
        o = -23.44 * np.pi/180
    elif planet == "Mars":
        o = -25.19 * np.pi/180

    # matrix rotation about the x-axis
    R_mat = np.matrix([[1, 0, 0], [0, np.cos(o), np.sin(o)], [0, -np.sin(o), np.cos(o)]])
    Vec_ECI = np.array(np.matmul(R_mat, Vec_Helio))[0]

    return Vec_ECI


# def hyperbolic_orbital_elements(r_p, vinf_vec, mu, direction): 
    # find orbital elements of hyperbolic trajectory given excess velocity and parking orbit periapse

    # calculate excess speed 
    v_inf = np.linalg.norm(vinf_vec)
    
    # calculate beta, the angle between apse line and excess velocity magnitude
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    # rotate excess velocity by beta to get vector in the apse line
    z = np.array([0, 0, 1])
    if direction == 0:
        apse_vec = vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta))
    else:
        apse_vec = vinf_vec * np.cos(beta) + np.cross(-z, vinf_vec) * np.sin(beta) + -z * (np.dot(-z, vinf_vec)) * (1 - np.cos(beta))

    # use unit vector in apse line direction to find periapse position vector
    rp_vec = apse_vec/np.linalg.norm(apse_vec) * r_p

    # calculate specific angular momentum, finding its direction and magnitude in sucession
    h_dir = np.cross(rp_vec, vinf_vec)
    h_unit = h_dir/np.linalg.norm(h_dir)
    h = r_p * np.sqrt(v_inf**2 + (2 * mu)/r_p)
    H = h_unit * h

    # Node Line
    N = np.cross([0,0,1], H)
    n = np.linalg.norm(N)

    # Right Ascension of the Ascending Node
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    # Inclination
    i = np.arccos(H[2]/h)

    # Eccentricity Vector
    E = rp_vec/np.linalg.norm(r_p) * e
    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )

    # True Anomaly of the asymptote
    ta_inf = np.arccos(-1/e)


    return h, e, a, i, raan, argp, ta_inf


def hyperbolic_orbital_elements(r_p, vinf_vec, mu, direction):
    v_inf = np.linalg.norm(vinf_vec)
    
    a = mu/(v_inf ** 2)
    e = r_p/a + 1
    beta = np.arccos(1/e)

    z = np.array([0, 0, 1])

    if direction == 0:                  #Arrival: direction = 0
        apse_vec = (vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta)))
    elif direction == 1:                #Departure: direction = 1
        apse_vec = -(vinf_vec * np.cos(beta) + np.cross(z, vinf_vec) * np.sin(beta) + z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta)))
    elif direction == 2:                #Departure other side: direction = 2
        apse_vec = -(vinf_vec * np.cos(beta) + np.cross(-z, vinf_vec) * np.sin(beta) + -z * (np.dot(z, vinf_vec)) * (1 - np.cos(beta)))

    #periapse position vector
    rp_actual = apse_vec/np.linalg.norm(apse_vec) * r_p

    #specific angular momentum
    h_dir = np.cross(rp_actual, vinf_vec)
    h_unit = h_dir/np.linalg.norm(h_dir)
    h = r_p * np.sqrt(v_inf**2 + (2 * mu)/r_p)
    H = h_unit * h

    #node line
    N = np.cross([0,0,1], H)

    # Right Ascension of the Ascending Node
    n = np.linalg.norm(N)
    if N[1] >= 0:
        raan = np.arccos(N[0]/n)
    elif N[1] < 0:
        raan = 2*np.pi - np.arccos(N[0]/n)

    #inclination
    i = np.arccos(H[2]/h)

    #ta at soi
    ta_inf = np.arccos(-1/e)

    #eccentricity
    E = rp_actual/np.linalg.norm(rp_actual) * e

    #argument of perigee
    if E[2] >= 0:
        argp = np.arccos( np.dot(N,E)/n/e )
    elif E[2] < 0:
        argp = 2*np.pi - np.arccos( np.dot(N,E)/n/e )

    return h, e, a, i, raan, argp, ta_inf


def compute_hyperbolic_anomalies(t_vec, ma, e, h, mu):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly
    ma_vec = t_vec * mu**2/h**3 * (e**2-1)**(3/2) + ma

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(F):
        return e * np.sinh(F) - F

    # Calculate eccentric anomaly --- Newton's Method
    ea = np.zeros(len(ma_vec))      # initialise eccentric anomaly vector

    for j in range(len(ma_vec)):
        F0 = ma_vec[j]              # set initial guess for E
        tol = 1                     # initialise tolerance value
   
        while tol > 1e-8:           # iterate until tolerance condition (<1e-8) is satisfied
            f = kepler_equation(F0) - ma_vec[j]
            f_prime = e * np.cosh(F0) - 1
            F1 = F0 - f/f_prime     # Newton's method

            tol = abs(F1-F0)        # calculate tolerance
            F0 = F1                 # update guess value for following iteration

        ea[j] = F1                  # update eccentric anomaly vector

    # Calculate true anomaly
    ta = 2 * np.arctan( np.sqrt((e+1)/(e-1)) * np.tanh(ea/2) )

    # Wrap anomalies to -pi:pi
    ta = [i-2*np.pi if i>np.pi else i for i in ta]
    ea = [i-2*np.pi if i>np.pi else i for i in ea]
    ma_vec = [i-2*np.pi if i>np.pi else i for i in ma_vec]

    return ta, ea, ma_vec


def circular_parking_orbit(alt, ma, i, raan, argp, mu):
    # simulating a circular parking orbit
    # given its altitude and hyperbolic trajectory orbital elements
     
    if mu == 3.986e5:
        r = 6378 + alt
    elif mu == 4.2828e4:
        r = 3390 + alt

    v = np.sqrt(mu/r)                       # velocity (m/s)

    # orbital parameters
    e = 0
    h = v * r
    T = 2*np.pi*r**(3/2) / np.sqrt(mu)
    n = 2*np.pi/T

    # simulate the orbit
    t_vec = np.linspace(0, T, 10000)
    ta, ea, ma_vec = compute_anomalies(t_vec, ma, e, n)
    p, q, w, dp, dq, dw = elements_to_perifocal(ta, e, mu, h)
    x, y, z, dx, dy, dz = perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp)

    v_vec = [dx, dy, dz]

    return v, x, y, z, v_vec
        


