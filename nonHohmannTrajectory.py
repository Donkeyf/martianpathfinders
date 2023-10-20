import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt
import FindStateVector as fsv

########################################    CONSTANTS DEFINITIONS   ########################################
mu_e = 3.986004418 * 10**14          # Standard gravitational parameter (Earth) (m^3/s^2)
mu_s = 1.327 * 10**20               # Standard gravitational parameter (Sun) (m^3/s^2)
mu_m = 4.282837 * 10**13            # Standard gravitational parameter (Mars) (m^3/s^2)

r_E = 6378000                      # Radius of Earth (m)
m_E = 5.974 * 10**24               # Mass of Earth (kg)
R_E = 149.6 * 10**9                 # Earth-Sun distance (m)

r_S = 696340000                    # Radius of the Sun (m)
m_S = 1.989 * 10**30               # Mass of the Sun (kg)

r_M = 3397000                      # Radius of Mars (m)
m_M = 5.974 * 10**24                # Mass of Mars (kg)
R_M = 227.9 * 10**9                 # Mars-Sun distance (m)


j2_E = 1.08262668 * 10**(-3)         # J2 Perturbation Constant for Earth
j2_M = 1.9555 * 10**(-3)               # J2 Perturbation Constant for Mars

########################################    Solve Lambert's problem     ########################################
def lamberts_problem(R1, R2, t, mu, case="pro", z=0):
    # z: initial iteration value of solution for z


    # Stumpff functions
    def S(z):
        if z > 0:
            s = (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z)**3
        elif z < 0:
            s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z)**3
        else:
            s = 1/6
        return s

    def C(z):
        if z > 0:
            c = (1 - np.cos(np.sqrt(z)))/z
        elif z < 0:
            c = (np.cosh(np.sqrt(-z)) - 1)/(-z)
        else:
            c = 1/2
        return c    

    # subfunctions
    def y(z):
        return r1 + r2 + A*(z*S(z) - 1)/np.sqrt(C(z))
    
    def F(z,t):
        F = (y(z)/C(z))**(1.5)*S(z) + A*np.sqrt(y(z)) - np.sqrt(mu)*t
        return F
    
    def dF(z):
        if z == 0:
            dF = np.sqrt(2)/40*y(0)**1.5 + A/8*( np.sqrt(y(0)) + A*np.sqrt(1/(2*y(0))) )
        else:
            dF = (y(z)/C(z))**1.5 * ( 1/2/z * (C(z)-3/2*S(z)/C(z)) + 3/4*S(z)**2/C(z) ) + A/8*( 3*S(z)/C(z)*np.sqrt(y(z)) + A*np.sqrt(C(z)/y(z)) )
        return dF


    # (1) calculate position magnitudes
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    # (2) calculate change in true anomaly
    theta = np.arccos(np.dot(R1, R2) / r1 / r2)

    # z-component of positions cross product
    c_z = np.cross(R1, R2)[2]

    # check and update theta's quadrant as necessary
    if case == "pro":
        if c_z < 0:
            theta = 2*np.pi - theta
    elif case == "retro":
        if c_z >= 0:
            theta = 2*np.pi - theta

    # (3) 'A' constant
    A = np.sin(theta) * np.sqrt(r1*r2/(1-np.cos(theta)))

    # (4) solving for z with Newton's Method
    ## initialise z value
    while F(z,t) < 0:
        z += 0.1
        z = round(z, 5)

    ## Newton's method
    tol = 1                             # initialise tolerance value  
    while tol > 1e-8:                   # iterate until tolerance condition (<1e-8) is satisfied
        z_np1 = z - F(z,t)/dF(z)
        tol = abs(z_np1-z)              # update tolerance
        z = z_np1                       # update z for next iteration or for final answer

    # (5) Lagrange functions
    f = 1 - y(z)/r1
    g = A*np.sqrt(y(z)/mu)
    g_dot = 1 - y(z)/r2

    # (6) endpoint velocities of the transfer orbit
    V1 = 1/g * (R2 - f*R1)
    V2 = 1/g * (g_dot*R2 - R1)

    return V1, V2