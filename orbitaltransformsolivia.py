import numpy as np

def compute_anomalies(t_vec, ma, e, n=0):
    # Compute true, eccentric and mean anomalies from mean anomaly and eccentricity

    # Propagate mean anomaly using mean motion (vector over time)
    ma_vec = ma + n * t_vec

    # Define Kepler's equation as a local function (for optimisation)
    def kepler_equation(E):
        return E - e * np.sin(E)

    ea = np.zeros(len(ma_vec)) #initialise eccentric anomaly

    for j in range(len(ma_vec)):
        E0 = ma_vec[j] #this is the initial guess
        tol = 1

        while tol >= 1e-8:
            f = kepler_equation(E0) - ma_vec[j]
            f_prime = 1 - e * np.cos(E0)
            E1 = E0 - f/f_prime
            tol = abs(E1 - E0)
            E0 = E1
        # Calc eccentric anomaly
        ea[j] = E1

    # Calc true anomaly
    ta = 2 * np.arctan(np.tan(ea/2) * np.sqrt((1 + e)/(1 - e)))

    # Wrap anomalies to -pi:pi
    ta = (ta + np.pi) % (2 * np.pi) - np.pi
    ea = (ea + np.pi) % (2 * np.pi) - np.pi
    ma_vec = (ma_vec + np.pi) % (2 * np.pi) - np.pi
    return ta, ea, ma_vec

def compute_orbital_velocity(h, e, ta, mu):
    r = h**2/mu * 1/(1 + e * np.cos(ta))

    v_r = mu / h * e * np.sin(ta) # Radial velocity
    v_n = h/r # Normal/Tangential velocity

    return v_r, v_n

def elements_to_perifocal(ta, a, e, mu, h):
    # Calc perifocal distance in m
    r = h**2/mu * 1/(1 + e * np.cos(ta))

    # Compute perifocal coordinates
    p = r * np.cos(ta) #x-coord
    q = r * np.sin(ta) #y-coord
    w = r * 0 #z-coord is 0 because it is the orbital plane - r * 0 so it remains array of equal length

    # Compute perifocal velocities
    dp = - mu/h * np.sin(ta)
    dq = mu/h * (e + np.cos(ta))
    dw = w

    return p, q, w, dp, dq, dw


def perifocal_to_eci(p, q, w, dp, dq, dw, i, raan, argp):
    
    # Transform coordinates to ECI frame
    x = np.zeros_like(p) #initialising
    y = np.zeros_like(p)
    z = np.zeros_like(p)

    dx = np.zeros_like(p)
    dy = np.zeros_like(p)
    dz = np.zeros_like(p)

    p_to_e = perifocal_to_eci_matrix(i, raan, argp) #calling the function defining the transformation matrix

    # Define vectors of perifocal coordinates/velocities
    perifocal_coord_vec = np.array([p, q, w])
    perifocal_vel_vec = np.array([dp, dq, dw])

    # Compute vectors of ECI coordinates/velocities
    eci_coord_vec = np.dot(p_to_e, perifocal_coord_vec) #dot product of perifocal coords with transformation matrix
    eci_vel_vec = np.dot(p_to_e, perifocal_vel_vec) #perifocal velocities to ECI velocities (via matrix multiplication)

    # Compute coordinates
    x, y, z = eci_coord_vec #ECI coords assigned to array
    dx, dy, dz = eci_vel_vec #ECI velocities assigned to array

    return x, y, z, dx, dy, dz


def perifocal_to_eci_matrix(i, raan, argp):

    # Calculate transformation matrix from perifocal to ECI frame
    c_i = np.cos(i) #rotation angles/components - capture rotation and orientation aspects of transformation
    s_i = np.sin(i) #these are computed to save space in the matrix
    c_raan = np.cos(raan)
    s_raan = np.sin(raan)
    c_argp = np.cos(argp)
    s_argp = np.sin(argp)

    p_to_e = np.array([ #3x3 matrix
        [-s_raan * c_i * s_argp + c_raan * c_argp, -s_raan * c_i * c_argp - c_raan * s_argp, s_raan * s_i],
        [c_raan * c_i * s_argp + s_raan * c_argp, c_raan * c_i * c_argp - s_raan * s_argp, -c_raan * s_i],
        [s_i * s_argp, s_i * c_argp, c_i]
    ])

    return p_to_e

def eci_to_ecef(x, y, z, dx, dy, dz, theta):
    
    # Transform coordinates to ECEF frame
    x_prime = np.zeros_like(x) #initialising
    y_prime = np.zeros_like(x)
    z_prime = np.zeros_like(x)

    dx_prime = np.zeros_like(x)
    dy_prime = np.zeros_like(x)
    dz_prime = np.zeros_like(x)

    # Loop through each time step
    for j in range(len(x)):
        i_to_ef = eci_to_ecef_matrix(theta[j]) #calling the function defining the transformation matrix

        # Define vectors of ECI coordinates/velocities
        eci_coords = np.array([x[j], y[j], z[j]])
        eci_vels = np.array([dx[j], dy[j], dz[j]])

        # Compute vectors of ECEF coordinates/velocities
        ecef_coord_vec = np.dot(i_to_ef, eci_coords) #dot product of ECI coords with transformation matrix
        ecef_vel_vec = np.dot(i_to_ef, eci_vels) #ECI velocities to ECEF velocities (via matrix multiplication)

        # Compute coordinates
        x_prime[j], y_prime[j], z_prime[j] = ecef_coord_vec #ECEF coords assigned to array
        dx_prime[j], dy_prime[j], dz_prime[j] = ecef_vel_vec #ECEF velocities assigned to array

    return x_prime, y_prime, z_prime, dx_prime, dy_prime, dz_prime


def eci_to_ecef_matrix(theta):

    # Calculate rotation matrix from ECI frame to ECEF frame

    i_to_ef = np.array([
        [np.cos(theta), np.sin(theta), 0],
        [-np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])

    return i_to_ef
